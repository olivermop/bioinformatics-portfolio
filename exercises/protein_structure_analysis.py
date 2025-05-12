# /bioinformatics-portfolio/exercises/protein_structure_analysis.py
# This script analyzes the secondary structure of the BRCA1 protein using the PDB database.

from Bio import PDB
import matplotlib.pyplot as plt
import os

# Directory where results and the downloaded PDB file will be saved
results_dir = '/Users/Olivermop/Documents/bioinformatics_portfolio/results'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Path to the downloaded PDB file (BRCA1: 1T29)
pdb_file = '/Users/Olivermop/Documents/bioinformatics_portfolio/data/PDB/1T29.pdb'

# Load the PDB structure using BioPython
parser = PDB.PDBParser()
structure = parser.get_structure('BRCA1', pdb_file)

# Initialize counters for secondary structure elements
alpha_helices = 0
beta_sheets = 0
other_structures = 0

# Iterate over models, chains, and residues to analyze secondary structure
for model in structure:
    for chain in model:
        for residue in chain:
            # Only consider standard amino acid residues
            if residue.get_id()[0] == ' ':
                res_id = residue.get_id()
                # Identify secondary structure by DSSP codes in res_id
                if 'H' in res_id:
                    alpha_helices += 1
                elif 'E' in res_id:
                    beta_sheets += 1
                else:
                    other_structures += 1

# Visualize the results in a bar chart
labels = ['Alpha Helices', 'Beta Sheets', 'Other Structures']
values = [alpha_helices, beta_sheets, other_structures]

plt.bar(labels, values)
plt.title('BRCA1 Protein Secondary Structure')
plt.ylabel('Number of Residues')
plt.savefig(os.path.join(results_dir, 'brca1_secondary_structure.png'))
plt.show()

# Save the results to a text report
output_file = os.path.join(results_dir, 'protein_structure_report.txt')
with open(output_file, 'w') as report:
    report.write("BRCA1 Protein Secondary Structure Analysis\n")
    report.write(f"Number of residues in alpha helices: {alpha_helices}\n")
    report.write(f"Number of residues in beta sheets: {beta_sheets}\n")
    report.write(f"Number of residues in other structures: {other_structures}\n")

print(f"Analysis complete. See report at {output_file}")
