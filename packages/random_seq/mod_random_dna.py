# 
# Part I:
# Function that writes random DNA sequences following 3 input parameters:
# 1. how many sequences do you want?
# 2. how long do you want the sequences?
# 3. what GC% should sequences be?
# 
# Part II:
# Function that writes a FASTA file from a list of sequences, and names them following a given an identifier.
#
# Part III:
# Script should be made ready to be used through Bash.
#

# Packages
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np

# Function to generate list of random DNA sequences as SeqRecords objects
def random_dna_generator(n, length, gc, batch_id):
    """
    Generates random DNA sequences following input parameters.

    Input:    
    n (int): Number of sequence to generate.
    length (int): Lenght of each sequence (in bp).
    gc (int): Sequence GC% (between 0 and 100).
    batch_id (str): General identifier for sequences (e.g., experiment)
    
    Output:
    List of SeqRecord objects with random DNA sequences

    """
    # Define DNA nucleotides and their probabilities based on GC%
    nts = ['A', 'T', 'G', 'C']
    nts_prob = [(100-gc)/100/2, (100-gc)/100/2, gc/100/2, gc/100/2]

    # Create 2D numpy array with n rows - each row forming a random DNA sequence of specified length and gc%
    random_dna_npa = np.random.choice(nts, size = (n, length), p = nts_prob)
    
    # Convert the numpy array into a list of random DNA sequences
    random_dna_seq = [''.join(seq) for seq in random_dna_npa]

    # Convert list of strings into list of SeqRecords with random DNA sequences, IDs and descriptions
    records = []
    for i, dna_seq in enumerate(random_dna_seq):
        record_id = f'{batch_id}_seq{i+1}'
        record_name = f'randomdna_{batch_id}_seq{i+1}'
        record_description = f'Random {length} bp DNA sequence with {gc}% GC from {batch_id}'
        record = SeqRecord(Seq(dna_seq), id = record_id, name = record_name, description = record_description)
        records.append(record)
    
    return records

# Function to write FASTA file from list of SeqRecords
def seqrecord_fasta_writer(records, filename):
    """Write a list of SeqRecord objects to a FASTA file."""
    with open(filename, 'w') as outfile:
        SeqIO.write(records, outfile, "fasta")

# Funtion to print sequences from list of SeqRecords
def seqrecord_fasta_printer(records):
    """Print SeqRecords from list to the console."""
    for record in records:
        print(f'>{record.id} {record.description}\n{record.seq}')

# Example of usage within Python
# test = random_dna_generator(n=10, length=10000, gc=50, batch_id='test')
# seqrecord_fasta_writer(test, "random_dna_test.fasta")
# seqrecord_fasta_printer(test)






# 
# Initial version did not use numpy arrays:
# for n in range(n_seq):
#     random_seq = ''.join(np.random.choice(nts, size = len_seq, p = nts_prob))
#     print(random_seq)
# 
# Note:
# could have been written using similar random.choice() function from library `random`,
# but the `numpy` function performs better with larger sequences.
# random_seq = ''.join(random.choice(nts, k = len_seq, weights = nts_prob))
#