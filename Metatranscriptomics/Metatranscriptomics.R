BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("biomaRt")
BiocManager::install("EnhancedVolcano")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
# Load libraries
library(tximport)
library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Define paths to quant.sf files
setwd("/Users/m.behruznia@bham.ac.uk/Library/CloudStorage/OneDrive-UniversityofBirmingham/DeSeq/DeSeq")
samples <- read.csv("Book1.csv", stringsAsFactors = TRUE)
quant_files <- file.path(samples$SampleID, "quant.sf")
names(quant_files) <- samples$SampleID




# Download from Ensembl or create manually from Biomart website by downloading TranscriptID and GeneID 
# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve transcript-to-gene mapping
tx2gene <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
                 mart = ensembl)

# Save to CSV
write.csv(tx2gene, "tx2gene.csv", row.names = FALSE)



tx2gene <- read.csv("tx2gene.csv")


# Import data with tximport
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, 
                ignoreTxVersion = TRUE)

# Verify imported data
head(txi$counts)


# DIFERRENTIAL EXPRESSION ANALYSIS ----

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Condition)

# Pre-filtering (optional): Remove low-expression genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2 pipeline
dds <- DESeq(dds)
# After running DESeq I got this message: 1 rows did not converge in beta, labelled 
# in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
# To see what it is:
# non_converged <- rownames(dds)[!mcols(dds)$betaConv]
# print(non_converged) #"ENSG00000164941"


# Extract results for pairwise comparisons
results_AIH_vs_Donor <- results(dds, contrast = c("Condition", "AIH", "Donor"))
results_SN_vs_Donor <- results(dds, contrast = c("Condition", "SN", "Donor"))

# After moving to the next step, I got this error: logical subscript contains NAs
# Inspect the padj column to confirm the presence of NA values 
# sum(is.na(results_AIH_vs_Donor$padj))  # Count the number of NA values
# 51405 I will proceed with filtering genes with very low expression (Line 48)
# After filtering these genes, 12779 genes have NA

# Exclude rows with NA values in the padj

sig_results_AIH <- results_AIH_vs_Donor[!is.na(results_AIH_vs_Donor$padj) & results_AIH_vs_Donor$padj < 0.05, ]
sig_results_SN <- results_SN_vs_Donor[!is.na(results_SN_vs_Donor$padj) & results_SN_vs_Donor$padj < 0.05, ]

# Investigate genes with NA values
# na_genes <- results_AIH_vs_Donor[is.na(results_AIH_vs_Donor$padj), ]
# head(na_genes)  # View the first few rows

# View significant results (PROBLEM WITH NA RESULTS)
# sig_results_AIH <- results_AIH_vs_Donor[results_AIH_vs_Donor$padj < 0.05, ]
# sig_results_SN <- results_SN_vs_Donor[results_SN_vs_Donor$padj < 0.05, ]

# Save results
write.csv(sig_results_AIH, "DE_AIH_vs_Donor.csv")
write.csv(sig_results_SN, "DE_SN_vs_Donor.csv")


# VISUALISATION----


# PCA Plot to visualize sample clustering
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Condition")



# Create volcano plot for AIH vs Donor
library(EnhancedVolcano)
EnhancedVolcano(results_AIH_vs_Donor,
                lab = rownames(results_AIH_vs_Donor),
                x = "log2FoldChange",
                y = "pvalue",
                title = "AIH vs Donor")


# Create volcano plot for SN vs Donor
EnhancedVolcano(results_SN_vs_Donor,
                lab = rownames(results_SN_vs_Donor),
                x = "log2FoldChange",
                y = "pvalue",
                title = "SN vs Donor")




# FUNCTIONAL ANNOTATION ----


sig_genes_SN <- rownames(sig_results_SN)
sig_genes_AIH <- rownames(sig_results_AIH)


enrich_GO_SN <- enrichGO(gene = sig_genes_SN,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENSEMBL",
                      ont = "BP",
                      pvalueCutoff = 0.05)

enrich_GO_AIH <- enrichGO(gene = sig_genes_AIH,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENSEMBL",
                         ont = "BP",
                         pvalueCutoff = 0.05)

head(enrich_GO_SN)
head(enrich_GO_AIH)

write.csv(as.data.frame(enrich_GO_SN), "SN_GO_enrichment_results.csv", row.names = FALSE)
write.csv(as.data.frame(enrich_GO_AIH), "AIH_GO_enrichment_results.csv", row.names = FALSE)

# Barplot of top enriched GO terms
barplot(enrich_GO_SN, showCategory = 10, title = "Top Enriched GO Terms in SN Patients")
barplot(enrich_GO_AIH, showCategory = 10, title = "Top Enriched GO Terms in AIH patients")


# Dotplot visualization
dotplot(enrich_GO_SN, showCategory = 20, title = "Top Enriched GO Terms in SN Patients")


# ggplot visualisation of EnrichGO categories
enrich_GO_SN_2 <- as.data.frame(enrich_GO_SN)

# Filter for the top 20 rows with the smallest p.adjust values
library(dplyr)

# Filter for the top 20 rows with the smallest p.adjust values
top_10 <- enrich_GO_SN_2 %>%
  arrange(p.adjust) %>%    # Sort by p.adjust in ascending order
  slice_head(n = 10)       # Select the top 20 rows

# Pass the filtered data to ggplot
ggplot(data = top_10, aes(x = reorder(Description, p.adjust), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", alpha = .8, width = .4) +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") + # Change color from blue to red
  xlab("") +
  ggtitle("Top 10 Gene Categories in SN Patients") +
  theme_bw()



# Concatenate Salmon output ----

# Load required libraries
library(dplyr)

telescope_samples <- read.csv("Book1.csv", stringsAsFactors=TRUE)
telescope_quant_files <- file.path(samples$SampleID, "telescope_report.tsv")

# Load transcript-to-gene mapping (tx2gene)
# Ensure you have a `tx2gene.csv` file that maps transcript IDs to gene IDs
tx2gene <- read.csv("tx2gene.csv", header = TRUE)  # Adjust path if needed

# Use tximport to import Salmon outputs and summarize at the gene level
# Ensure you're getting the raw counts from Salmon by using countsFromAbundance = "no" to get raw counts
txi <- tximport(files = quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "no")

# Round the counts to integers
gene_counts <- round(txi$counts)

# Extract gene-level counts
# gene_counts <- as.data.frame(txi$counts)

# Save the rounded counts to CSV
write.csv(gene_counts, "salmon_gene_counts.csv", row.names = TRUE)
head(gene_counts)




# Concatenate Telescope output ----

# Ensure telescope_samples dataframe contains SampleID and correct file paths
telescope_samples <- read.csv("Book1.csv", stringsAsFactors = TRUE)
telescope_quant_files <- file.path(telescope_samples$SampleID, "telescope_report.tsv")
names(telescope_quant_files) <- telescope_samples$SampleID

# Load required library
library(dplyr)

# Set directory containing Telescope output files
telescope_output <- "/Users/m.behruznia@bham.ac.uk/Library/CloudStorage/OneDrive-UniversityofBirmingham/DeSeq/DeSeq/telescope_output/"  # Replace with the correct path
file_list <- list.files(telescope_output, full.names = TRUE, pattern = "\\.tsv$")

# Initialize an empty list to store data
telescope_data <- list()

# Loop through each file and extract transcript name and final_count
for (file in file_list) {
  # Extract sample name from the filename before "_telescope_report.tsv"
  sample_name <- gsub("_telescope_report\\.tsv$", "", basename(file))
  
  # Read the file, skipping the first metadata line
  temp_data <- read.csv(file, sep = "\t", header = TRUE, skip = 1)
  
  # Extract relevant columns: transcript and final_count
  temp_data <- temp_data %>%
    select(transcript, final_count) %>%
    rename(!!sample_name := final_count)  # Rename the final_count column to the sample name
  
  # Add the data to the list
  telescope_data[[sample_name]] <- temp_data
}

# Combine all Telescope data into one data frame by transcript names (rows)
merged_telescope_counts <- Reduce(function(x, y) merge(x, y, by = "transcript", all = TRUE),
                                  telescope_data)

# Set row names to transcript names
row.names(merged_telescope_counts) <- merged_telescope_counts$transcript
merged_telescope_counts <- merged_telescope_counts[, -1]  # Remove the transcript column

# Replace NA with 0 for missing values
merged_telescope_counts[is.na(merged_telescope_counts)] <- 0

# Save the merged Telescope counts to a CSV file
write.csv(merged_telescope_counts, "telescope_counts_all_samples.csv", row.names = TRUE)


# DE analysis with both salmon and telescope ----

# Read the Salmon counts and Telescope counts CSV files
salmon_data <- read.csv("salmon_gene_counts.csv", row.names = 1)
telescope_data <- read.csv("telescope_counts_all_samples.csv", row.names = 1)

# Display the first few rows to ensure they are in the correct format
head(salmon_data)
head(telescope_data)

# Create a sample metadata table
sample_metadata <- data.frame(
  sample = colnames(salmon_data),  # Same columns as salmon_data and telescope_data
  condition = c("donor", "SN", "donor", "donor", "donor", "AIH", "SN", "AIH", "AIH", "SN", "AIH", "AIH", "SN", "SN", "SN")  # Example conditions
)
# Check the metadata table
sample_metadata

# Create DESeq2 Dataset for Salmon data
dds_salmon <- DESeqDataSetFromMatrix(countData = salmon_data, colData = sample_metadata, design = ~ condition)

# Create DESeq2 Dataset for Telescope data
dds_telescope <- DESeqDataSetFromMatrix(countData = telescope_data, colData = sample_metadata, design = ~ condition)

# Filter out rows with low counts (e.g., keep genes with at least 10 counts in total across samples)
dds_salmon <- dds_salmon[rowSums(counts(dds_salmon)) >= 10, ]
dds_telescope <- dds_telescope[rowSums(counts(dds_telescope)) >= 10, ]

# Run DESeq2 on Salmon data
# dds_salmon <- DESeq(dds_salmon)

# Run DESeq2 on Telescope data
dds_telescope <- DESeq(dds_telescope)


# Get results for Salmon dataset
res_salmon_AIH_vs_donor <- results(dds_salmon, contrast = c("condition", "AIH", "donor"))
res_salmon_SN_vs_donor <- results(dds_salmon, contrast = c("condition", "SN", "donor"))

# Get results for Telescope dataset
res_telescope_AIH_vs_donor <- results(dds_telescope, contrast = c("condition", "AIH", "donor"))
res_telescope_SN_vs_donor <- results(dds_telescope, contrast = c("condition", "SN", "donor"))

# Filter significant results for AIH vs donor
sig_salmon_AIH_vs_donor <- subset(res_salmon_AIH_vs_donor, padj < 0.05 & abs(log2FoldChange) > 1)

# Filter significant results for SN vs donor
sig_salmon_SN_vs_donor <- subset(res_salmon_SN_vs_donor, padj < 0.05 & abs(log2FoldChange) > 1)

# Filter significant results for AIH vs donor
sig_telescope_AIH_vs_donor <- subset(res_telescope_AIH_vs_donor, padj < 0.05 & abs(log2FoldChange) > 1)

# Filter significant results for SN vs donor
sig_telescope_SN_vs_donor <- subset(res_telescope_SN_vs_donor, padj < 0.05 & abs(log2FoldChange) > 1)

# Save significant results to CSV
# write.csv(as.data.frame(sig_salmon_AIH_vs_donor), "sig_salmon_AIH_vs_donor.csv")
# write.csv(as.data.frame(sig_salmon_SN_vs_donor), "sig_salmon_SN_vs_donor.csv")
write.csv(as.data.frame(sig_telescope_AIH_vs_donor), "sig_telescope_AIH_vs_donor.csv")
write.csv(as.data.frame(sig_telescope_SN_vs_donor), "sig_telescope_SN_vs_donor.csv")


# ggplot volcano ----
create_volcano_plot <- function(results, title) {
  # Add a column to categorize significance
  results$Significance <- "Not Significant"
  results$Significance[results$padj < 0.05 & abs(results$log2FoldChange) > 1] <- "Significant"
  
  ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    theme_minimal() +
    theme(legend.title = element_blank())
}


# Convert DESeq results to data frame for ggplot
res_salmon_AIH_vs_donor_df <- as.data.frame(res_salmon_AIH_vs_donor)
res_salmon_SN_vs_donor_df <- as.data.frame(res_salmon_SN_vs_donor)
res_telescope_AIH_vs_donor_df <- as.data.frame(res_telescope_AIH_vs_donor)
res_telescope_SN_vs_donor_df <- as.data.frame(res_telescope_SN_vs_donor)

# Plot for Salmon dataset
volcano_salmon_AIH <- create_volcano_plot(res_salmon_AIH_vs_donor_df, "Salmon: AIH vs Donor")
volcano_salmon_SN <- create_volcano_plot(res_salmon_SN_vs_donor_df, "Salmon: SN vs Donor")

# Plot for Telescope dataset
volcano_telescope_AIH <- create_volcano_plot(res_telescope_AIH_vs_donor_df, "Telescope: AIH vs Donor")
volcano_telescope_SN <- create_volcano_plot(res_telescope_SN_vs_donor_df, "Telescope: SN vs Donor")

# Display plots
print(volcano_salmon_AIH)
print(volcano_salmon_SN)
print(volcano_telescope_AIH)
print(volcano_telescope_SN)



vsd <- vst(dds_telescope, blind = FALSE)
plotPCA(vsd, intgroup = "condition")



# Create volcano plot for AIH vs Donor
library(EnhancedVolcano)

# Create volcano plot for SN vs Donor
EnhancedVolcano(res_telescope_SN_vs_donor_df,
                lab = rownames(res_telescope_SN_vs_donor_df),
                x = "log2FoldChange",
                y = "pvalue",
                title = "SN vs Donor")


EnhancedVolcano(res_telescope_AIH_vs_donor_df,
                lab = rownames(res_telescope_AIH_vs_donor_df),
                x = "log2FoldChange",
                y = "pvalue",
                title = "AIH vs Donor")



# Combined volcano plot ----
# Add source column
res_salmon_AIH_vs_donor_df$Source <- "Salmon"
res_telescope_AIH_vs_donor_df$Source <- "Telescope"

# Combine the two datasets
combined_AIH_vs_donor <- rbind(res_salmon_AIH_vs_donor_df, res_telescope_AIH_vs_donor_df)


# Repeat the same process for SN
res_salmon_SN_vs_donor_df$Source <- "Salmon"
res_telescope_SN_vs_donor_df$Source <- "Telescope"

combined_SN_vs_donor <- rbind(res_salmon_SN_vs_donor_df, res_telescope_SN_vs_donor_df)


create_combined_volcano_plot <- function(data, title) {
  # Add a column for significance categorization
  data$Significance <- "Not Significant"
  data$Significance[data$padj < 0.05 & abs(data$log2FoldChange) > 1] <- "Significant"
  
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = Source, shape = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Salmon" = "blue", "Telescope" = "green")) +
    scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", color = "Dataset", shape = "Significance") +
    theme_minimal()
}


# Combined Volcano Plot for AIH vs donor
volcano_combined_AIH <- create_combined_volcano_plot(combined_AIH_vs_donor, "Combined Volcano Plot: AIH vs Donor")

# Combined Volcano Plot for SN vs donor
volcano_combined_SN <- create_combined_volcano_plot(combined_SN_vs_donor, "Combined Volcano Plot: SN vs Donor")

# Display plots
print(volcano_combined_AIH)
print(volcano_combined_SN)




# Filter for a single dataset (e.g., AIH vs Donor)
single_dataset <- combined_AIH_vs_donor

# Add significance column and color mapping
single_dataset$Significance <- ifelse(
  single_dataset$padj < 0.05 & abs(single_dataset$log2FoldChange) > 1, 
  "Significant", 
  "Not Significant"
)

single_dataset$Color <- ifelse(
  single_dataset$Significance == "Not Significant", 
  "grey", 
  single_dataset$Source
)

# Create the volcano plot
ggplot(single_dataset, aes(x = log2FoldChange, y = -log10(padj), color = Color, shape = Source)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Salmon" = "green", "Telescope" = "red", "Not Significant" = "grey")
  ) +
  scale_shape_manual(
    values = c("Salmon" = 16, "Telescope" = 17)  # Circle for Salmon, triangle for Telescope
  ) +
  # Add threshold lines for significance and fold change
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + # Fold change thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # P-value threshold
  labs(
    title = "Volcano Plot: AIH vs Donor",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Dataset",
    shape = "Dataset"
  ) +
  theme_minimal()




# Combined plot with enhanced volcano package----

# Convert results to data frames
salmon_AIH_vs_donor_df <- as.data.frame(res_salmon_AIH_vs_donor)
telescope_AIH_vs_donor_df <- as.data.frame(res_telescope_AIH_vs_donor)

# Add a column to indicate the source
salmon_AIH_vs_donor_df$Source <- "Salmon"
telescope_AIH_vs_donor_df$Source <- "Telescope"

# Combine the datasets
combined_AIH_vs_donor <- rbind(salmon_AIH_vs_donor_df, telescope_AIH_vs_donor_df)


# The same for SN and donor
# Convert results to data frames
salmon_SN_vs_donor_df <- as.data.frame(res_salmon_SN_vs_donor)
telescope_SN_vs_donor_df <- as.data.frame(res_telescope_SN_vs_donor)

# Add a column to indicate the source
salmon_SN_vs_donor_df$Source <- "Salmon"
telescope_SN_vs_donor_df$Source <- "Telescope"

# Combine the datasets
combined_SN_vs_donor <- rbind(salmon_SN_vs_donor_df, telescope_SN_vs_donor_df)

# Create a named color vector for the sources
color_mapping <- setNames(
  c("green", "red"), # Colors for Salmon and Telescope
  c("Salmon", "Telescope") # Values in the Source column
)


EnhancedVolcano(
  combined_AIH_vs_donor,
  lab = rownames(combined_AIH_vs_donor), # Gene/transcript names
  x = "log2FoldChange",
  y = "padj",
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value",
  title = "Combined Volcano Plot: AIH vs Donor",
  subtitle = "Salmon and Telescope Data",
  pCutoff = 0.05,
  FCcutoff = 1,
  colCustom = color_mapping[combined_AIH_vs_donor$Source], # Map colors by Source
  colAlpha = 0.8,
  pointSize = 2.0
)



EnhancedVolcano(
  combined_SN_vs_donor,
  lab = rownames(combined_SN_vs_donor), # Gene/transcript names
  x = "log2FoldChange",
  y = "padj",
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value",
  title = "Combined Volcano Plot: SN vs Donor",
  subtitle = "Salmon and Telescope Data",
  pCutoff = 0.05,
  FCcutoff = 1,
  colCustom = color_mapping[combined_SN_vs_donor$Source],, # Custom colors for sources
  colAlpha = 0.8
)


# Add telescope locus name only

# Add a Locus column to the combined dataset
combined_AIH_vs_donor$Locus <- ifelse(
  combined_AIH_vs_donor$Source == "Telescope",
  rownames(combined_AIH_vs_donor), # Replace with actual locus names if stored elsewhere
  NA # No labels for Salmon rows
)

EnhancedVolcano(
  combined_AIH_vs_donor,
  lab = combined_AIH_vs_donor$Locus, # Use Locus column for labels
  x = "log2FoldChange",
  y = "padj",
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value",
  title = "Combined Volcano Plot: AIH vs Donor",
  subtitle = "Salmon and Telescope Data",
  pCutoff = 0.05,
  FCcutoff = 1,
  colCustom = color_mapping[combined_AIH_vs_donor$Source], # Map colors by Source
  colAlpha = 0.8,
  pointSize = 2.0, cutoffLineType = "longdash",
  cutoffLineCol = "black",
  cutoffLineWidth = 0.4,
)


combined_SN_vs_donor$Locus <- ifelse(
  combined_SN_vs_donor$Source == "Telescope",
  rownames(combined_SN_vs_donor), # Replace with actual locus names if stored elsewhere
  NA # No labels for Salmon rows
)

# Define colors based on significance----
combined_AIH_vs_donor$color <- ifelse(
  combined_AIH_vs_donor$padj < 0.05 & abs(combined_AIH_vs_donor$log2FoldChange) > 1,
  color_mapping[combined_AIH_vs_donor$Source], # Use custom colors for significant genes
  "grey" # Use grey for non-significant genes
)


color_mapping <- setNames(
  c("green", "red", "grey"), # Colors for Salmon, Telescope, and Non-Significant
  c("Salmon", "Telescope", "Non-Significant") # Corresponding labels
)

# Add color column based on the mapping
combined_AIH_vs_donor$color <- color_mapping[combined_AIH_vs_donor$color_group]

# Plot with EnhancedVolcano
EnhancedVolcano(
  combined_AIH_vs_donor,
  lab = combined_AIH_vs_donor$Locus, # Use Locus column for labels
  x = "log2FoldChange",
  y = "padj",
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value",
  title = "Combined Volcano Plot: AIH vs Donor",
  subtitle = "Salmon and Telescope Data",
  pCutoff = 0.05,
  FCcutoff = 1,
  colCustom = color_mapping[combined_AIH_vs_donor$color_group], # Use the new color column
  colAlpha = 0.8,
  pointSize = 2.0
)




# same process for SN
combined_SN_vs_donor$color <- ifelse(
  combined_SN_vs_donor$padj < 0.05 & abs(combined_SN_vs_donor$log2FoldChange) > 1,
  color_mapping[combined_SN_vs_donor$Source], # Use custom colors for significant genes
  "grey" # Use grey for non-significant genes
)


combined_SN_vs_donor$color_group <- ifelse(
  combined_SN_vs_donor$padj < 0.05 & abs(combined_SN_vs_donor$log2FoldChange) > 1,
  combined_SN_vs_donor$Source, # Use source for significant genes
  "Non-Significant" # Label for non-significant genes
)

color_mapping_2 <- setNames(
  c("green", "red", "grey"), # Colors for Salmon, Telescope, and Non-Significant
  c("Salmon", "Telescope", "Non-Significant") # Corresponding labels
)

# Add color column based on the mapping
combined_SN_vs_donor$color <- color_mapping_2[combined_SN_vs_donor$color_group]

# Plot with EnhancedVolcano
EnhancedVolcano(
  combined_SN_vs_donor,
  lab = combined_SN_vs_donor$Locus, # Use Locus column for labels
  x = "log2FoldChange",
  y = "padj",
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value",
  title = "Combined Volcano Plot: SN vs Donor",
  subtitle = "Salmon and Telescope Data",
  pCutoff = 0.05,
  FCcutoff = 1,
  colCustom = color_mapping_2[combined_SN_vs_donor$color_group], # Use the new color column
  colAlpha = 0.8,
  pointSize = 2.0)


# barplots ----

# load the library
library(forcats)
library(ggplot2)
library(dplyr)
# Reorder following the value of another column:

data_AIH <- read.csv("sig_telescope_AIH_vs_donor.csv")

data_SN <- read.csv("sig_telescope_SN_vs_donor.csv")

data_AIH %>%
  mutate(X = fct_reorder(X, log2FoldChange),
         color = ifelse(log2FoldChange < 0, "red", "green")) %>%
  ggplot(aes(x = X, y = log2FoldChange, fill = color)) +
  geom_bar(stat = "identity", alpha = .6, width = .4) +
  coord_flip() +
  scale_fill_identity() +
  xlab("") +
  ggtitle("Differentially expressed TEs in AIH patients") + # Add your title here
  theme_bw()


data_SN %>%
  mutate(X = fct_reorder(X, log2FoldChange),
         color = ifelse(log2FoldChange < 0, "red", "green")) %>%
  ggplot(aes(x = X, y = log2FoldChange, fill = color)) +
  geom_bar(stat = "identity", alpha = .6, width = .4) +
  coord_flip() +
  scale_fill_identity() +
  xlab("") +
  ggtitle("Differentially expressed TEs in SN patients") + # Add your title here
  theme_bw()


# PCA based on both datasets
# Load required libraries
library(DESeq2)
library(PCAtools)
library(ggplot2)

# Step 1: Variance Stabilization for Both Datasets
# Variance-stabilized data for Salmon
vsd_salmon <- varianceStabilizingTransformation(dds_salmon, blind = FALSE)

# Variance-stabilized data for Telescope
vsd_telescope <- varianceStabilizingTransformation(dds_telescope, blind = FALSE)

# Step 2: Combine Variance Stabilized Matrices
# Extract variance-stabilized assay data
assay_salmon <- assay(vsd_salmon)
assay_telescope <- assay(vsd_telescope)

# Combine the assays row-wise (genes/TEs as rows, samples as columns)
combined_matrix <- rbind(assay_salmon, assay_telescope)

# Combine metadata (assumes both datasets use the same `colData` structure)
combined_metadata <- rbind(
  colData(vsd_salmon) %>% as.data.frame(),
  colData(vsd_telescope) %>% as.data.frame()
)


colnames(colData(vsd_salmon))
colnames(colData(vsd_telescope))

# Step 1: Rename columns to match
colnames(colData(vsd_salmon)) <- c("sample", "condition")
colnames(colData(vsd_telescope)) <- c("sample", "condition", "sizeFactor")

# Step 2: Add missing column `sizeFactor` to `vsd_salmon`
colData(vsd_salmon)$sizeFactor <- NA

dim(combined_matrix)       # Rows: features (genes/TEs), Columns: samples
dim(combined_metadata)    # Rows: samples, Columns: metadata variables


# Step 3: Combine metadata
combined_metadata <- rbind(
  colData(vsd_salmon) %>% as.data.frame(),
  colData(vsd_telescope) %>% as.data.frame()
)

# Verify the combined metadata
head(combined_metadata)
dim(combined_metadata)  # Check the dimensions


# Combine metadata (assumes both datasets use the same `colData` structure)
combined_metadata <- rbind(
  colData(vsd_salmon) %>% as.data.frame(),
  colData(vsd_telescope) %>% as.data.frame()
)


# Step 3: Perform PCA on the Combined Data
# Create a DESeq2 object for combined data (design ~1 because we are not testing a condition here)
dds_combined <- DESeqDataSetFromMatrix(countData = combined_matrix, colData = combined_metadata, design = ~1)
dds_combined <- DESeq(dds_combined, parallel = TRUE)
vsd_combined <- varianceStabilizingTransformation(dds_combined, blind = FALSE)

# Optional: Filter by variance (remove features with low variance)
removeVar <- 0.3  # Retain top 70% most variable features
pca_obj_comb
