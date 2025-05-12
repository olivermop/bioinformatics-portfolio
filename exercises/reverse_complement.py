# /bioinformatics-portfolio/exercises/reverse_complement.py
# This script computes the reverse complement of a given DNA sequence.

# Function to obtain the reverse complement of a DNA sequence
def reverse_complement(dna_sequence):
    # Dictionary defining the complement for each DNA nucleotide
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Clean the DNA sequence by removing line breaks and whitespace
    cleaned_sequence = dna_sequence.replace("\n", "").replace(" ", "")

    # Reverse the cleaned sequence and build the complement for each nucleotide
    # A generator is used to iterate over the reversed sequence and fetch the matching complement
    return ''.join(complement[base] for base in reversed(cleaned_sequence))


# The DNA sequence dataset is provided here
dataset = """ DATASET HERE """

# Call the function to get the reverse complement of the dataset
result = reverse_complement(dataset)

# Print the result, which is the reverse complement of the input DNA sequence
print(result)
