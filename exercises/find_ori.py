# bioinformatics-portfolio/exercises/find_ori.py
import re
import os

def find_origin_in_genome(genome_sequence, dnaA_box_pattern="TTATCCACA"):
    """
    Finds the origin of replication in a genome by locating multiple occurrences of the DnaA box pattern.
    :param genome_sequence: Full genome sequence as a string
    :param dnaA_box_pattern: The DnaA box consensus pattern to search for
    :return: The start and end positions of the predicted origin region, or a message if none found
    """
    # Build a regex that allows a single variation at the last position of the DnaA box
    pattern = re.compile(f"({dnaA_box_pattern[:7]}.{dnaA_box_pattern[8]})")

    # Collect all match start positions
    positions = [match.start() for match in pattern.finditer(genome_sequence)]

    if not positions:
        return "No DnaA box found."

    # Estimate the origin region as the span from the first to the last DnaA box occurrence
    start_origin = min(positions)
    end_origin = max(positions)

    return start_origin, end_origin

def read_genome_file(file_path):
    """
    Reads a genome sequence from a file.
    :param file_path: Path to the genome file
    :return: Genome sequence as a single string
    """
    with open(file_path, 'r') as file:
        return file.read().replace('\n', '')

def save_origin_to_file(origin_info, result_file_path):
    """
    Saves the origin position information to a results file.
    :param origin_info: Origin position tuple or message
    :param result_file_path: Path to the output file
    """
    with open(result_file_path, 'w') as file:
        file.write(f"Estimated origin position: {origin_info}\n")

if __name__ == "__main__":
    # Input genome file and output result paths
    genome_file_path = "/Users/Olivermop/Documents/bioinformatics_portfolio/data/DnaABoxes/vibrio_cholerae.txt"
    result_file_path = "/Users/Olivermop/Documents/bioinformatics_portfolio/results/ori_estimation.txt"

    # Load the genome sequence
    genome_sequence = read_genome_file(genome_file_path)

    # Find the origin region
    origin_position = find_origin_in_genome(genome_sequence)

    # Write the result to a file
    save_origin_to_file(origin_position, result_file_path)

    print(f"Estimated origin position saved in {result_file_path}")
