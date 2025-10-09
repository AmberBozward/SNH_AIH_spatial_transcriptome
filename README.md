# SNH_AIH_spatial_transcriptome

This repository contains the analysis code and workflows associated with the manuscript:

“Spatial transcriptomics links hepatocyte-macrophage interactions to viral signatures in seronegative hepatitis”
Submitted to Nature Communications, 2025

The repository includes workflows for CosMx single-cell spatial transcriptomics, Visium whole-transcriptome profiling, and metatranscriptomic analysis of explanted liver tissue, along with scripts for data processing, integration, and figure generation.

---

## Table of Contents

- [Overall Repository Structure](#overall-repository-structure)
- [CosMx Processing the Data and Figure Generation](#cosmx-processing-the-data-and-figure-generation)
  - [CosMx Repository Structure](#cosmx-repository-structure)
  - [Data Structure](#data-structure)
  - [Requirements](#requirements)
  - [Usage](#usage)
  - [Customisation](#customisation)
  - [Choosing areas on CosMx spatial plot](#choosing-areas-on-cosmx-spatial-plot)
  - [Neighbourhood Enrichment Analysis and Ripleys Spatial Statistics](#neighbourhood-enrichment-analysis-and-ripleys-spatial-statistics)
- [Visium Processing the Data and Figure Generation](#visium-processing-the-data-and-figure-generation)
  - [Visium Repository Structure](#visium-repository-structure)
  - [Visium Requirements](#visium-requirements)
  - [Input Data](#input-data)
- [Metatranscriptomics](#metatranscriptomics)
  - [Metatrascriptomics Repository Structure ](#metatrascriptomics-repository-structure )
  - [Metatranscriptomics Code Usage](#metatranscriptomics-code-usage)
- [Contributors](#contributors)
- [License](#license)
- [Citation](#citation)

---

## Overall Repository Structure

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

# CosMx Processing the Data and Figure Generation

This repository provides an **easy-to-use pipeline for pre-processing CosMx data and generating figures**. It is designed for liver disease research (e.g., autoimmune hepatitis (AIH), seronegative (SN) liver disease, and donor (D) samples) and can produce UMAPs, barplots, boxplots, heatmaps, dotplots, and CellChat-based cell-cell interaction visualizations.  

The scripts are structured for **reproducibility**, with figures saved automatically to a dedicated folder. This repository is ideal for immunology researchers looking to visualize cell types, gene expression patterns, and intercellular communication in spatial transcriptomics data.

---

## CosMx Repository Structure

The CosMx repository files are structured by the below:

```
├── 2.1.combine_all.r          # Script to combine all samples
├── 2.2.cell_label.r           # Script to label cell types
├── 2.3.singleR.r              # Script to use singleR to identify cell types
├── 2.4.immune.r               # Script to split cells into immune and non-immune to help with cell labelling
├── 2.5.polygon.r              # Script to 
├── 2.6.de.r                   # Script to do differentials
├── Choosing_CosMx_areas.py    # Script to chosse areas of interest on spatial CosMx plot
└── CosMx_plots.R          	   # Script to generate all figures
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
│   └── CosMx_plots.r              # Script to generate all figures
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

## Choosing areas on CosMx spatial plot

'Choosing_CosMx_areas.py' is a Python script for interactively selecting and extracting specific spatial regions of interest (ROIs) from CosMx spatial transcriptomics data.

The script allows users to:

	•	Load the CosMx segmentation image and associated spatial coordinates.
	•	Visually select regions within the tissue using polygon or rectangular selection tools.
	•	Extract cell IDs or molecular data corresponding to the selected region.

Dependencies:

	•	Python ≥ 3.8
	•	matplotlib
	•	numpy
	•	pandas
	•	opencv-python (for image handling)

Purpose:
This tool was developed to facilitate manual spatial subsetting of CosMx datasets, enabling region-specific analyses (e.g., comparing parenchyma vs non-parenchyma areas in liver tissue).

---

## Neighbourhood Enrichment Analysis and Ripleys Spatial Statistics

All relevant Python code for neighborhood enrichment analysis and Ripley's spatial statistics to understand cellular spatial organisation patterns can be found in the the following project repository:

https://github.com/kyliesavoye/hepatitis-spatial-transcriptomics/

---

# Visium Processing the Data and Figure Generation

This repository provides an **easy-to-use pipeline for pre-processing Visium data and generating figures**. It is designed for liver disease research (e.g., autoimmune hepatitis (AIH), seronegative (SN) liver disease, and donor (D) samples) and can produce UMAPs, barplots, boxplots, heatmaps and dotplots.

---

## Visium Repository Structure

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

## Visium Requirements

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

# Metatranscriptomics

This repository provides an **easy-to-use pipeline for pre-processing metatranscriptomics data and generating figures**. It is designed for liver disease research (e.g., autoimmune hepatitis (AIH), seronegative (SN) liver disease, and donor (D) samples).

---

## Metatrascriptomics Repository Structure 

└── Metatranscriptomics.R          	   # Script for all analysis and figures

---

## Metatranscriptomics Code Usage 

Barcode trimming and QC with fastp:
fastp -i input_R1.fastq.gz -I input_R2.fastq.gz \
-o output_R1_trimmed.fastq.gz -O output_R2_trimmed.fastq.gz \
--detect_adapter_for_pe  --trim_poly_g  --html fastp_report.html \
--length_required 40 --qualified_quality_phred 20 --thread 10

Creating bowtie2 index and aligning reads to it
bowtie2-build /path/to/GRCh38.fasta GRCh38_index
bowtie2 -x GRCh38_index -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz  --very-sensitive-local -k 100 --score-min L, 0, 1.6 -S output.sam
samtools view -bS output.sam | samtools sort -o output_sorted.bam
samtools index output_sorted.bam 

Running Telescope
telescope assign /path/to/output_sorted.bam/  /path/to/GTF_file/transcripts.gtf --ncpu 12
Note: obtain HERV and L1 annotation file (transcripts.gtf) following the instructions available in: https://github.com/mlbendall/telescope_annotation_db/tree/af3c359/builds/retro.hg38.v1

Host gene quantification with Salmon:
preparing metadata:
grep "^>" <(gunzip -c primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
preparing concatenated transcriptome and genome reference file for index:
zcat gencode.transcripts.fa.gz primary_assembly.genome.fa.gz | gzip -c > gentrome.fa.gz
running Salmon:
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
salmon quant -i salmon_index -l A -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz -p 16 -o salmon_output

Import Salmon into R and run DESeq2 using tximport (see R script for full detail, the main steps are listed below)
use tximport and tx2gene
run DESeq2 pipeline 
visualisation and combining Telescope counts with Salmon gene counts

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
