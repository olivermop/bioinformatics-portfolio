# /bioinformatics-portfolio/exercises/vcf_analysis.py
# This script analyzes genetic variants (SNPs) from a VCF file related to breast cancer.
# It uses the cyvcf2 library to read the file and classify variants by their functional impact.

from cyvcf2 import VCF
import os

# Path to the downloaded VCF file containing structural variants (SVs) from the MCF7 cell line
vcf_file = '/Users/Olivermop/Documents/bioinformatics_portfolio/data/VCF/GSM3336911_MCF7-lumpy-sv.vcf.gz'

# Directory where the analysis report will be stored
results_dir = '/Users/Olivermop/Documents/bioinformatics_portfolio/results'
# Create the directory if it doesn't exist
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Path to the output report file for variant analysis
output_report = os.path.join(results_dir, 'vcf_variants_report.txt')

# Create and open the report file
with open(output_report, 'w') as report:
    # Write the report header
    report.write("Genetic Variant Report (SNPs)\n")
    report.write("Analysis based on a VCF file related to breast cancer.\n\n")

    # Read the VCF file using cyvcf2
    vcf_reader = VCF(vcf_file)

    # Initialize counters to classify variants by functional impact
    high_impact = 0      # High-impact variants
    moderate_impact = 0  # Moderate-impact variants
    low_impact = 0       # Low-impact variants
    modifier_impact = 0  # Modifier-impact or no-impact variants

    # Iterate over each variant record in the VCF
    for record in vcf_reader:
        # Check for functional annotations in the INFO field (ANN)
        if 'ANN' in record.INFO:
            annotations = record.INFO['ANN'].split(',')  # Annotations are comma-separated
            for annotation in annotations:
                # The functional impact is the third field (index 2) in the annotation
                impact = annotation.split('|')[2]

                # Classify the variant by its impact
                if impact == "HIGH":
                    high_impact += 1
                elif impact == "MODERATE":
                    moderate_impact += 1
                elif impact == "LOW":
                    low_impact += 1
                else:
                    modifier_impact += 1

    # Write the analysis results to the report
    report.write(f"High-impact variants: {high_impact}\n")
    report.write(f"Moderate-impact variants: {moderate_impact}\n")
    report.write(f"Low-impact variants: {low_impact}\n")
    report.write(f"Modifier-impact variants: {modifier_impact}\n")

# Print confirmation message
print(f"Analysis complete. See report at {output_report}")
