import numpy as np
import argparse

##############################################
# This program calculates the percent of a genome in contigs >= 1,000,000 base pairs
# Usage: python percent_of_genome_over_1mil.py <input.fa>
##############################################

# Define a function to get data frame of contig names and lengths
    # based on fastalength command from chromosomer package
def get_contig_lengths(fasta):
    contig_lengths = []
    contig_names = []

    with open(fasta, 'r') as f:
        for line in f:
                if line.startswith('>'):
                    contig_names.append(line.rstrip('\n'))
                else:
                    contig_lengths.append(len(line.rstrip('\n')))
    contig_lengths = np.array(contig_lengths)

    return(contig_lengths)

##############################################
# MAIN PROGRAM
# Use the array returned from contig_lengths() to calculate the percent of the genome in contigs >= 1,000,000 base pairs
##############################################

    # 1. Use argparse to read in the path to the fasta file
parser = argparse.ArgumentParser(description="Calculate the percent of a genome in contigs >= 1,000,000 base pairs")
parser.add_argument("fasta_path", help="path to uncompressed fasta file", type=str)
args = parser.parse_args()

fasta_path = args.fasta_path

    # 2. Calculate assembly length and number of contigs

contig_lengths_array = get_contig_lengths(fasta_path)
num_contigs = int(len(contig_lengths_array))
assembly_length = int(sum(contig_lengths_array))


    # 3. Calculate the % of contigs >= 1,000,000 base pairs
contig_lengths_over_1mil = contig_lengths_array[contig_lengths_array >= 1000000]
assembly_length_over_1mil = int(sum(contig_lengths_over_1mil))
assembly_percent_over_1mil = (assembly_length_over_1mil / assembly_length) * 100

    # 4. Print the results
print("You have " + str(num_contigs) + " contigs in your assembly.")
print("The total length of your assembly is " + str(assembly_length) + " base pairs.")
print("The total length of contigs >= 1,000,000 base pairs is " + str(assembly_length_over_1mil) + " base pairs.") 
print("This is " + str(assembly_percent_over_1mil) + "% of your assembly.")
