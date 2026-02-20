
setwd("C:/Users/LORD LAPTOP/Desktop/lncrna/updated")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(biomaRt)
  library(pheatmap)
  library(survival)
  library(survminer)
})

set.seed(123)
Sys.setenv(TAR = "internal")
options(timeout = 120)

# -------------------------
# Helper functions
# -------------------------

# Remove Ensembl version suffix (e.g., ENSG...10.12 -> ENSG...10)
strip_ensembl_version <- function(x) gsub("\\..*", "", x)

# TCGA patient ID from barcode (first 12 chars)
tcga_patient <- function(barcode) substr(barcode, 1, 12)

# For safety: keep one sample per patient (first occurrence)
keep_one_sample_per_patient <- function(se_obj) {
  barcodes <- colnames(se_obj)
  patients <- tcga_patient(barcodes)
  keep_idx <- !duplicated(patients)
  se_obj[, keep_idx]
}

# Fetch annotation from Ensembl via biomaRt
get_ensembl_annot <- function(ensembl_ids) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  annot <- getBM(
    attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = unique(ensembl_ids),
    mart = mart
  )
  annot[!duplicated(annot$ensembl_gene_id), ]
}

# =========================
# 1) DOWNLOAD: pick 50 tumor + 50 normal
# =========================
query_meta <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

files <- getResults(query_meta)

tumor_cases <- files %>%
  filter(sample_type == "Primary Tumor") %>%
  slice_sample(n = 50) %>%
  pull(cases)

normal_cases <- files %>%
  filter(sample_type == "Solid Tissue Normal") %>%
  slice_sample(n = 50) %>%
  pull(cases)

selected_cases <- c(tumor_cases, normal_cases)

query_small <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = selected_cases
)

GDCdownload(query_small, method = "api", files.per.chunk = 1)
se <- GDCprepare(query_small)

# =========================
# 2) DESeq2 Differential Expression (Tumor vs Normal)
# =========================
counts <- assay(se, "unstranded")

condition <- ifelse(se$shortLetterCode == "TP", "Tumor", "Normal")
condition <- factor(condition, levels = c("Normal", "Tumor"))

coldata <- data.frame(row.names = colnames(counts),
                      condition = condition)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Filter: gene expressed (>=10 counts) in at least 10 samples
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# LFC shrink for stable ranking
resLFC <- lfcShrink(
  dds,
  coef = "condition_Tumor_vs_Normal",
  type = "apeglm"
)

# Clean Ensembl IDs
rownames(resLFC) <- strip_ensembl_version(rownames(resLFC))

# =========================
# 3) Annotation + lncRNA filtering
# =========================
annot <- get_ensembl_annot(rownames(resLFC))

res_df <- as.data.frame(resLFC)
res_df$ensembl_gene_id <- rownames(resLFC)

res_annot <- merge(res_df, annot,
                   by = "ensembl_gene_id",
                   all.x = TRUE)

lnc_types <- c(
  "lncRNA", "lincRNA", "antisense",
  "processed_transcript", "sense_intronic", "sense_overlapping"
)

sig_lnc <- res_annot %>%
  filter(gene_biotype %in% lnc_types,
         !is.na(padj),
         padj < 0.05,
         abs(log2FoldChange) > 1)

cat("Significant lncRNAs:", nrow(sig_lnc), "\n")

# =========================
# 4) Volcano plot (lncRNA highlighted)
# =========================
res_annot$significant <- "Not Sig"
res_annot$significant[
  !is.na(res_annot$padj) &
    res_annot$padj < 0.05 &
    abs(res_annot$log2FoldChange) > 1 &
    res_annot$gene_biotype %in% lnc_types
] <- "Significant"

ggplot(res_annot, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  ggtitle("TCGA-BRCA: lncRNA Differential Expression")

# =========================
# 5) PCA + Heatmap
# =========================
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

# Heatmap top 50 sig lncRNAs
top50 <- sig_lnc %>% arrange(padj) %>% slice(1:50)
top50_ids <- top50$ensembl_gene_id

rownames(vsd) <- strip_ensembl_version(rownames(vsd))
common_ids <- intersect(top50_ids, rownames(vsd))

mat <- assay(vsd)[common_ids, , drop = FALSE]
mat_scaled <- t(scale(t(mat)))

annotation_col <- as.data.frame(colData(dds)[, "condition", drop = FALSE])

pheatmap(mat_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Top 50 Differentially Expressed lncRNAs")

# =========================
# 6) SURVIVAL ANALYSIS (Tumor only + one sample per patient)
# =========================

# 6.1) Clinical data
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

clinical$time <- ifelse(!is.na(clinical$days_to_death),
                        clinical$days_to_death,
                        clinical$days_to_last_follow_up)

clinical$event <- ifelse(clinical$vital_status == "Dead", 1, 0)

clinical <- clinical %>%
  mutate(time = as.numeric(time),
         event = as.numeric(event)) %>%
  filter(!is.na(submitter_id),
         !is.na(time), !is.na(event),
         time > 0)

# 6.2) Build VSD for tumor samples only, one sample per patient
tumor_cols <- which(condition == "Tumor")
vsd_tumor <- vsd[, tumor_cols]

# Enforce one tumor sample per patient (important!)
vsd_tumor <- keep_one_sample_per_patient(vsd_tumor)

# 6.3) Map HOTAIR symbol -> Ensembl
# (use annotation table already fetched)
hotair_id <- annot$ensembl_gene_id[annot$external_gene_name == "HOTAIR"][1]

if (is.na(hotair_id) || !hotair_id %in% rownames(vsd_tumor)) {
  stop("HOTAIR Ensembl ID not found in your VSD tumor matrix.")
}

hotair_expr <- assay(vsd_tumor)[hotair_id, ]

expr_df <- data.frame(
  submitter_id = tcga_patient(colnames(vsd_tumor)),
  HOTAIR = as.numeric(hotair_expr)
)

surv_data <- inner_join(expr_df, clinical, by = "submitter_id") %>%
  filter(!is.na(HOTAIR))

cat("Survival N:", nrow(surv_data),
    "Events:", sum(surv_data$event == 1), "\n")

# 6.4) Kaplan–Meier: Median split (robust baseline)
median_expr <- median(surv_data$HOTAIR, na.rm = TRUE)
surv_data$group_median <- ifelse(surv_data$HOTAIR > median_expr, "High", "Low")

fit_km <- survfit(Surv(time, event) ~ group_median, data = surv_data)

ggsurvplot(fit_km,
           data = surv_data,
           pval = TRUE,
           risk.table = TRUE,
           conf.int = TRUE,
           title = "KM Survival – HOTAIR (Tumor-only, 1 sample/patient)",
           legend.title = "HOTAIR",
           legend.labs = c("High", "Low"))

# 6.5) Cox: continuous expression
cox_hotair <- coxph(Surv(time, event) ~ HOTAIR, data = surv_data)
print(summary(cox_hotair))

# 6.6) Optional: optimal cutpoint (may overfit)
cut <- surv_cutpoint(surv_data,
                     time = "time",
                     event = "event",
                     variables = "HOTAIR",
                     minprop = 0.2)

cat_dat <- surv_categorize(cut)
surv_data$group_optimal <- cat_dat$HOTAIR

fit_opt <- survfit(Surv(time, event) ~ group_optimal, data = surv_data)

ggsurvplot(fit_opt,
           data = surv_data,
           pval = TRUE,
           risk.table = TRUE,
           conf.int = TRUE,
           title = "KM (Optimal Cutpoint) – HOTAIR",
           legend.title = "HOTAIR",
           legend.labs = c("Low", "High"))

# =========================
# 7) Survival screening for ALL significant lncRNAs (Tumor-only)
# =========================

# Prepare tumor-only patient IDs once
tumor_patients <- tcga_patient(colnames(vsd_tumor))

# Make a quick lookup: gene symbol -> Ensembl id (first match)
sym2ens <- annot %>%
  filter(!is.na(external_gene_name), !is.na(ensembl_gene_id)) %>%
  group_by(external_gene_name) %>%
  summarise(ensembl_gene_id = first(ensembl_gene_id), .groups = "drop")

sig_symbols <- unique(na.omit(sig_lnc$external_gene_name))

results_tbl <- lapply(sig_symbols, function(gene){
  
  gene_id <- sym2ens$ensembl_gene_id[sym2ens$external_gene_name == gene][1]
  if (is.na(gene_id)) return(NULL)
  if (!gene_id %in% rownames(vsd_tumor)) return(NULL)
  
  expr <- assay(vsd_tumor)[gene_id, ]
  
  tmp <- data.frame(
    submitter_id = tumor_patients,
    expr = as.numeric(expr)
  )
  
  tmp2 <- inner_join(tmp, clinical, by = "submitter_id") %>%
    filter(!is.na(expr), !is.na(time), !is.na(event))
  
  # safeguards
  if (nrow(tmp2) < 30) return(NULL)
  if (sum(tmp2$event == 1) < 10) return(NULL)
  
  fit <- coxph(Surv(time, event) ~ expr, data = tmp2)
  s <- summary(fit)
  
  data.frame(
    gene = gene,
    ensembl_gene_id = gene_id,
    HR = s$conf.int[1, "exp(coef)"],
    HR_low95 = s$conf.int[1, "lower .95"],
    HR_high95 = s$conf.int[1, "upper .95"],
    pvalue = s$coefficients[1, "Pr(>|z|)"]
  )
})

results_tbl <- bind_rows(results_tbl) %>%
  mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
  arrange(FDR)

head(results_tbl, 20)
cat("FDR < 0.05:", sum(results_tbl$FDR < 0.05, na.rm = TRUE), "\n")

# =========================
# 8) Multivariate Cox (age + stage) for HOTAIR
# =========================
age_col   <- intersect(c("age_at_diagnosis", "age", "age_at_initial_pathologic_diagnosis"), names(clinical))[1]
stage_col <- intersect(c("ajcc_pathologic_stage", "pathologic_stage", "tumor_stage"), names(clinical))[1]

clinical2 <- clinical %>%
  transmute(
    submitter_id,
    time, event,
    age = if (!is.na(age_col)) as.numeric(.data[[age_col]]) else NA_real_,
    stage = if (!is.na(stage_col)) as.character(.data[[stage_col]]) else NA_character_
  ) %>%
  filter(!is.na(time), !is.na(event), time > 0) %>%
  mutate(
    stage = gsub("Stage ", "", stage),
    stage = gsub("[A-C]$", "", stage),
    stage = factor(stage)
  )

dat_hotair <- inner_join(expr_df, clinical2, by = "submitter_id") %>%
  filter(!is.na(age), !is.na(stage))

cox_hotair_mva <- coxph(Surv(time, event) ~ HOTAIR + age + stage, data = dat_hotair)
print(summary(cox_hotair_mva))

ggadjustedcurves(cox_hotair_mva, data = dat_hotair, variable = "HOTAIR")
