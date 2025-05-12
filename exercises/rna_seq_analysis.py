# /bioinformatics-portfolio/exercises/rna_seq_analysis.py
# This script performs RNA-Seq analysis using the GSE176078 dataset from GEO,
# which provides detailed spatial and single-cell transcriptomics in breast cancer.
# It includes bulk RNA sequencing of primary tumor cells from the three main
# breast cancer subtypes: ER+, HER2+, and TNBC (triple-negative).

import os
import pandas as pd
import numpy as np
import statsmodels.api as sm  # For statistical tests
import matplotlib.pyplot as plt
import seaborn as sns  # For data visualization

# Create the "results" directory if it does not exist
if not os.path.exists("results"):
    os.makedirs("results")

# Path to the RNA-Seq raw counts file
file_path = "/Users/Olivermop/Documents/bioinformatics_portfolio/data/GEO/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt"

# Load RNA-Seq data (tab-delimited, first column is gene names)
data = pd.read_csv(file_path, sep="\t", index_col=0)

# Filter out lowly expressed genes (genes with fewer than 1 count in all samples)
filtered_data = data[(data > 1).sum(axis=1) > 1]

# Display the first few rows and column names to verify proper loading
print(filtered_data.head())
print(filtered_data.columns)

# Separate data into study groups (e.g., 'control' and 'treatment')
# Adjust group assignments based on the original sample metadata
control_samples = filtered_data[['CID3586', 'CID3921', 'CID3941']]  # Example control group
treatment_samples = filtered_data[['CID3948', 'CID3963', 'CID4066']]  # Example treatment group

# Calculate mean expression per group for downstream comparison
control_mean = control_samples.mean(axis=1)
treatment_mean = treatment_samples.mean(axis=1)

# Compute log2 fold change (log2 FC) with a pseudocount to avoid division by zero
log2_fc = np.log2(treatment_mean + 1) - np.log2(control_mean + 1)

# Perform Student’s t-test for each gene to calculate p-values
p_values = []
for gene in filtered_data.index:
    control_values = control_samples.loc[gene]   # Gene expression in control group
    treatment_values = treatment_samples.loc[gene]  # Gene expression in treatment group
    ttest = sm.stats.ttest_ind(control_values, treatment_values, usevar='pooled')
    p_values.append(ttest[1])  # Extract the p-value

# Combine log2 FC and p-values into a DataFrame for further analysis
results = pd.DataFrame({
    'log2_FC': log2_fc,
    'p_value': p_values
})

# Adjust p-values using Benjamini-Hochberg FDR to control false discoveries
results['adjusted_p_value'] = sm.stats.multipletests(results['p_value'], method='fdr_bh')[1]

# Create a volcano plot for differential expression
# Significant genes are marked in red
results['significant'] = results['adjusted_p_value'] < 0.05

plt.figure(figsize=(10, 6))
sns.scatterplot(x=results['log2_FC'], y=-np.log10(results['p_value']),
                hue=results['significant'], palette={True: 'red', False: 'blue'})
plt.title('Volcano Plot')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 p-value')
plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--')  # P-value cutoff line
plt.axvline(x=1, color='gray', linestyle='--')   # Positive log2 FC cutoff
plt.axvline(x=-1, color='gray', linestyle='--')  # Negative log2 FC cutoff
plt.savefig("results/volcano_plot.png")
plt.show()

# Generate a heatmap of the top 20 most significant genes
significant_genes = results[results['significant']].index
top_genes = filtered_data.loc[significant_genes].head(20)

plt.figure(figsize=(10, 8))
sns.heatmap(top_genes, cmap="RdBu_r", linewidths=0.5)
plt.title('Heatmap of Top Significant Genes')
plt.savefig("results/heatmap.png")
plt.show()

# Save a detailed report of the differential expression analysis
output_file_path = "results/differential_expression_results.txt"
try:
    with open(output_file_path, "w") as output_file:
        output_file.write("Differential Expression Analysis (RNA-Seq) Results:\n")
        output_file.write("This analysis compares gene expression between control and treatment conditions.\n")
        output_file.write(
            "Log2 Fold Change (log2 FC) indicates the change in gene expression between the two groups.\n")
        output_file.write(
            "Positive log2 FC means higher expression in treatment; negative means higher expression in control.\n")
        output_file.write(
            "Adjusted p-values are calculated with the Benjamini-Hochberg method to control FDR.\n\n")
        output_file.write("Summary of top results:\n")
        output_file.write(results.head(20).to_string())
        output_file.write("\n\nInterpretation:\n")
        output_file.write(
            "Genes with adjusted p-value < 0.05 are considered significantly differentially expressed.\n")
        output_file.write(
            "The volcano plot highlights the most over- and under-expressed genes, "
            "while the heatmap displays expression patterns of top significant genes.\n")
    print(f"Report saved successfully to: {output_file_path}")
except Exception as e:
    print(f"Error saving report: {e}")
