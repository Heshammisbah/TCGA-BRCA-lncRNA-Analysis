TCGA-BRCA-lncRNA-Analysis/
│
├── README.md
├── LICENSE
├── data/
│   └── (not raw TCGA data — only metadata or small example files)
│
├── scripts/
│   ├── 01_download_TCGA.R
│   ├── 02_differential_expression.R
│   ├── 03_annotation_lncRNA.R
│   ├── 04_survival_analysis.R
│   ├── 05_visualization.R
│
├── results/
│   ├── figures/
│   │   ├── volcano_plot.png
│   │   ├── PCA_plot.png
│   │   ├── heatmap_top50.png
│   │   ├── KM_median.png
│   │   ├── KM_optimal.png
│   │
│   └── tables/
│       ├── DE_lncRNAs.csv
│       ├── Cox_results.csv
│
└── manuscript/
    └── TCGA_BRCA_lncRNA_Manuscript.pdf
