# /bioinformatics-portfolio/exercises/genetic_variation_analysis.py
# This script uses a VCF file to identify genetic variants (SNPs) and classify them by functional impact.
# It uses `vcfpy` to read the data and classify SNPs based on information contained in the VCF file.
# The analysis results are saved in the "results" folder along with a classified list of detected variants.

import os
import vcfpy
import pandas as pd

# Create the "results" directory if it does not exist
if not os.path.exists("results"):
    os.makedirs("results")

# Path to the VCF file containing genetic data
vcf_file = '/Users/Olivermop/Documents/bioinformatics_portfolio/data/VCF/example_data.vcf'

# Open the VCF file and read the variant records
reader = vcfpy.Reader.from_path(vcf_file)

# Lists to store analysis results
variants = []
classification = []

# Process each variant in the VCF file
for record in reader:
    # Extract variant information
    variant_id = record.ID if record.ID else "N/A"
    chromosome = record.CHROM
    position = record.POS
    ref_allele = record.REF
    alt_alleles = ','.join(record.ALT)

    # Determine if the variant is a SNP (Single Nucleotide Polymorphism)
    if len(ref_allele) == 1 and all(len(str(alt)) == 1 for alt in record.ALT):
        variant_type = "SNP"
    else:
        variant_type = "Other"

    # Classify by functional impact (may vary depending on INFO annotations)
    impact_annotation = record.INFO.get("ANN", ["N/A"])[0]  # ANN field may contain functional annotations
    impact = impact_annotation.split('|')[2] if len(impact_annotation.split('|')) > 2 else "No information"

    # Store the variant data
    variants.append([variant_id, chromosome, position, ref_allele, alt_alleles, variant_type, impact])

# Convert the results to a pandas DataFrame
df_variants = pd.DataFrame(
    variants,
    columns=["ID", "Chromosome", "Position", "Reference Allele", "Alternative Alleles", "Variant Type", "Impact"]
)

# Write a text report with an explanation
with open("results/variant_analysis_report.txt", "w") as output_file:
    output_file.write("Genetic Variation Analysis (SNPs):\n")
    output_file.write("This analysis identifies genetic variants (SNPs) and classifies them by functional impact based on VCF annotations.\n\n")
    output_file.write("Summary of detected variants:\n\n")
    output_file.write(df_variants.to_string(index=False))
    output_file.write("\n\nInterpretation:\n")
    output_file.write("Variants are classified by functional impact using the 'ANN' field annotations in the VCF file.\n")
    output_file.write("SNPs with higher impact may be associated with significant functional effects.\n")

# Also save the analysis as a CSV file
df_variants.to_csv("results/variant_analysis_results.csv", sep="\t", index=False)

print("Analysis complete. Results saved in the 'results' folder.")
