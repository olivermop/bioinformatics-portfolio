# /bioinformatics-portfolio/exercises/calculate_skew.py
# This code calculates the skew function for a DNA sequence and finds the positions where the skew reaches its maximum value.

# Function to calculate the skew of a DNA sequence. Incrementing for each G and decrementing for each C.
def calculate_skew(sequence):
    skew = [0]  # Initialize the skew array with 0 at position 0
    for i in range(len(sequence)):
        if sequence[i] == 'G':
            skew.append(skew[-1] + 1) # Increase skew for each G
        elif sequence[i] == 'C':
            skew.append(skew[-1] - 1)  # Decrease skew for each C
        else:
            skew.append(skew[-1])  # Do not change the skew value for other letters
    return skew

# Function to find the position where the skew reaches its maximum value
def find_max_skew_position(sequence):
    skew = calculate_skew(sequence)
    max_skew_value = max(skew)  # Find the maximum skew value
    max_positions = [i for i, value in enumerate(skew) if value == max_skew_value]  # Positions where the skew is maximum
    return max_positions, max_skew_value

# Example sequence
sequence = "CATTCCAGTACTTCGATGATGGCGTGAAGA"

# Find the position where the skew reaches the maximum value
max_positions, max_skew_value = find_max_skew_position(sequence)


print(f"Maximum skew value: {max_skew_value}in the positions: {max_positions}")
