# SNH_AIH_spatial_transcriptome

This repository contains the analysis code and workflows associated with the manuscript:

“Spatial transcriptomics links hepatocyte-macrophage interactions to viral signatures in seronegative hepatitis”
Submitted to Nature Communications, 2025

The repository includes workflows for CosMx single-cell spatial transcriptomics, Visium whole-transcriptome profiling, and metatranscriptomic analysis of explanted liver tissue, along with scripts for data processing, integration, and figure generation.

---

## Repository Structure

```
SNH_spatial_transcriptome/
├── CosMx/
├── Visium/
├── Metatranscriptomics/
├── .gitignore/
├── LICENSE/
└── README.md/
```

---

# CosMx Pre-processing the Data and Figure Generation

This repository provides an **easy-to-use pipeline for pre-processing CosMx data and generating figures**. It is designed for liver disease research (e.g., autoimmune hepatitis (AIH), seronegative (SN) liver disease, and donor (D) samples) and can produce UMAPs, barplots, boxplots, heatmaps, dotplots, and CellChat-based cell-cell interaction visualizations.  

The scripts are structured for **reproducibility**, with figures saved automatically to a dedicated folder. This repository is ideal for immunology researchers looking to visualize cell types, gene expression patterns, and intercellular communication in spatial transcriptomics data.

---

## Repository Structure

The CosMx repository files are structured by the below:

```
├── 2.1.combine_all.r          # Script to combine all samples
├── 2.2.cell_label.r           # Script to label cell types
├── 2.3.singleR.r              # Script to use singleR to identify cell types
├── 2.4.immune.r               # Script to split cells into immune and non-immune to help with cell labelling
├── 2.5.polygon.r              # Script to 
├── 2.6.de.r                   # Script to do differentials
└── combined_plots.r           # Script to generate all figures
```

---

## Data structure 

```
├── data/                          
│   ├── seurat_obj_processed.rds   # Preprocessed Seurat object
├── scripts/                       # Script for workflow
│   └── 2.1.combine_all.r          # Script to combine all samples
│   └── 2.2.cell_label.r           # Script to label cell types
│   └── 2.3.singleR.r              # Script to use singleR to identify cell types
│   └── 2.4.immune.r               # Script to split cells into immune and non-immune to help with cell labelling
│   └── 2.5.polygon.r              # Script to 
│   └── 2.6.de.r                   # Script to do differentials
│   └── combined_plots.r           # Script to generate all figures
└── figures/                       # Output folder for all generated plots
```

- **data/**: Store your preprocessed Seurat object here. Must contain a `final_cell_types` column in `meta.data` for annotated cell types.
- **figures/**: All output figures (UMAP, barplot, boxplot, heatmap, dotplot, CellChat networks) will be automatically saved here.
- **scripts/**: Contains the combined plotting script. Designed to run as a single R script for all figure types.

---

## Requirements

- **R ≥ 4.2**
- Required packages:

```r
install.packages(c("dplyr", "ggplot2", "pheatmap", "patchwork", "ggrepel"))
BiocManager::install(c("SingleR", "CellChat", "Seurat"))
```

•	Input data: a Seurat object with:
	•	Annotated cell types in meta.data$final_cell_types
	•	Expression layers including data and counts

---

## Usage
	1.	Clone the repository:

```
git clone https://github.com/AmberBozward/SNH_AIH_spatial_transcriptomics.git
cd spatial-figures
```

	2.	Place your preprocessed Seurat/CosMx object in data/:

```r
data/seurat_obj_processed.rds
```

	3.	Run the combined plotting script in R:

```r
source("scripts/combined_plots.R")
```

	4.	Check the figures/ folder for all generated outputs:
	•	UMAPs
	•	Barplots
	•	Boxplots
	•	Heatmaps
	•	Dotplots
	•	CellChat network plots

---

## Customisation
	•	Modify genes for plots: Edit the genes_to_plot vector in the script for boxplots and dotplots.
	•	Change cell types: Ensure final_cell_types contains the desired annotations for coloring and subsetting.

---

## Neighbourhood enrichment analysis and Ripley's spatial statistics

All relevant Python code for neighborhood enrichment analysis and Ripley's spatial statistics to understand cellular spatial organisation patterns can be found in the the following project repository:

https://github.com/kyliesavoye/hepatitis-spatial-transcriptomics/

---

# Visium Pre-processing the Data and Figure Generation

This repository provides an **easy-to-use pipeline for pre-processing Visium data and generating figures**. It is designed for liver disease research (e.g., autoimmune hepatitis (AIH), seronegative (SN) liver disease, and donor (D) samples) and can produce UMAPs, barplots, boxplots, heatmaps and dotplots.

---

## 📂 Repository Structure

The repository files are structured by the below:

```
├── .gitignore                 
├── 1.QC/r                                      # Script for inital QC
├── 2a.integration.light.r                      # Script for light integration
├── 2b.integration.strict.r                     # Script for strict integration
├── 3a.cell_types_public.r                      # Script to label cell types
├── 3b.spatial_deconvolution_strict.r           # Script to strictly deconvolute cell types
├── 3c.spatial_deconvolution_light.r            # Script to lightly deconvolute cell types
├── 3d.pseudobulk_linneages_deconvolution.r     # Script to pseudobulk deconvoluted cell types
├── 3f_spatial_deconvolution_lineagges_light.r  # Script to pseudobulk deconvoluted cell types
├── 4a_region_barcodes.r                        # Script to generate barcodes from regions of interest (parenchyma vs non-parenchyma)
├── 4b_region_barcodes_pseudo.r                 # Script to pseudobulk chosen regions
├── 5a.plots.r                                  # Script to generate all figures
├── LICENSE                                     # License information for this repository
└── README.md/                                  # This document
```

---

## ⚙️ Requirements

### R Version
- R >= 4.2.0

### Required R Packages
Make sure the following packages are installed before running the scripts:

- `Seurat` (>= 5.0)
- `ggplot2`
- `dplyr`
- `tidyr`
- `patchwork`
- `Matrix`
- `stringr`
- `cowplot`
- `ggrepel`
- `gridExtra`
- `hdf5r`
- `spatstat` (if using spatial statistics)
- `SeuratDisk` (if working with `.h5Seurat` or converting from AnnData)

You can install missing packages with:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", 
                   "Matrix", "stringr", "cowplot", "ggrepel", "gridExtra"))
```

And from Bioconductor if needed:

```r
install.packages("Seurat")
install.packages("SeuratDisk")
install.packages("spatstat.geom")
install.packages("spatstat.core")
```

### Input Data

The scripts assume 10x Genomics Space Ranger outputs in the following format:

data/
  sample1/
    filtered_feature_bc_matrix.h5
    spatial/
      tissue_hires_image.png
      scalefactors_json.json
  sample2/
    ...

---
## Contributors

Dr Amber Bozward — CosMx and Visium workflows, integration, figure generation

Dr Mahboobeh Behruznia - Metatranscriptomics sequencing analysis

John Cole - CosMx and Visium workflows, integration, figure generation

Chiranjit Das - Figure generation

Kylie Savoye - CosMx neighbourhood analysis, figure generation

Professor Ye Oo - Supervisory role
  
---

## License

This repository is licensed under the MIT License (see LICENSE file).

---

## Citation

If you use this code, please cite:

Bozward et al., “Spatial transcriptomics links hepatocyte-macrophage interactions to viral signatures in seronegative hepatitis”, Nature Communications (2025).
