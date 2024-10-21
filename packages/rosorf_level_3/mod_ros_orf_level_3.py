# Import packages
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import pandas as pd

# PART 1: PREPARE DNA INPUT

# ros-orf level 2 v2.2
def fasta_parser_lvl2(fasta_file):
    """Parse FASTA file with sequences into dictionary with ID and sequence.
    
    Input: path to FASTA file
    Output: dictionary"""

    # Create empty dictionary
    seq_dict = {}

    # Open and parse FASTA file, then store sequneces in dictionary
    with open(fasta_file, 'r') as file:
        for item in SeqIO.parse(file, 'fasta'):
            seq_id = item.id
            seq = str(item.seq).upper()
            seq_dict[seq_id] = seq

    return seq_dict
    
# ros-orf level 3 tentative improvement:
# Don't store input DNA sequences from FASTA in dictionary, instead use list of SeqRecords
def fasta_parser_lvl3(fasta_file):
    """Parse FASTA file with sequences into list of SeqRecords.
    
    Input: path to FASTA file
    Output: list of SeqRecord objects"""


    # Read and parse FASTA file, then store SeqRecords in list
    seqrecords = list(SeqIO.parse(fasta_file, 'fasta'))

    return seqrecords
    
# PART 2: Find putative ORFs on DNA sequences

# ros-orf level 2 v2.2
def orf_finder_lvl2(dna_dict):
    """Find putative ORFs in DNA sequences provided in dataframe,
    and return dataframe with location data for putative ORFs and unique putative ORF ID.
    Note that start and stop codon location is stored as 0-based python index of the first codon nt.
    
    Input: dataframe with DNA ID and sequences    
    Output: dataframe with putative ORF location data"""

    # Set start and stop codons
    start_codon = ['ATG']
    stop_codon = ['TAA', 'TAG', 'TGA']

    # Create empty dictionary to store locations of putative ORFs
    orf_dict = {
        'Source DNA ID': [],
        'Putative ORF ID': [],
        'Strand': [],
        'Start codon index': [],
        'Stop codon index': []
    }

    # Prepare DNA from dataframe 
    for id, seq in dna_dict.items():
        dna_id = id
        dna_seq_sense = seq
        dna_seq_asense = str(Seq(dna_seq_sense).reverse_complement())

        # Search for putative ORFs and fill-out dictionary with data of found ORFs
        for dir, seq in {'Sense': dna_seq_sense, 'Antisense': dna_seq_asense}.items():
            for i in range(len(seq) - 5): 
                if seq[i:i+3] in start_codon:
                    for j in range(i+3, len(seq)-2, 3):
                        if seq[j:j+3] in stop_codon:
                            orf_id = f'{dna_id.lower()}_orf_{len(orf_dict['Putative ORF ID'])+1}'
                            orf_dict['Source DNA ID'].append(dna_id)
                            orf_dict['Putative ORF ID'].append(orf_id)
                            orf_dict['Strand'].append(dir)
                            orf_dict['Start codon index'].append(i)
                            orf_dict['Stop codon index'].append(j)
                            break
                        
    # Return dictionary with location data of putative ORFs
    return pd.DataFrame(orf_dict)


# ros-orf level 3 tentative mprovements:
# numpy arrays used for efficient string searching when looking for start and stop codons
def orf_finder_lvl3_v1(seqrecords):
    """Find putative ORFs in DNA sequences provided in list of SeqRecords,
    and return dataframe with location data for putative ORFs and unique putative ORF ID.
    
    Input: list of SeqRecord objects  
    Output: dataframe with putative ORF location data"""

    # Set start and stop codons
    start_codon = ['ATG']
    stop_codon = ['TAA', 'TAG', 'TGA']
    
    # Create empty dictionary to store locations of putative ORFs
    orf_dict = {
        'Source DNA ID': [],
        'Putative ORF ID': [],
        'Strand': [],
        #'Frame': [],
        'Start codon index': [],
        'Stop codon index': []
    }
    
    # Prepare DNA 
    for record in seqrecords:
        dna_id = record.id
        dna_seq_sense = str(record.seq)
        dna_seq_asense = str(Seq(dna_seq_sense).reverse_complement())
            
        # Create lists of all codon frames
        codons_sense_1 = [dna_seq_sense[i:i+3] for i in range(0, len(dna_seq_sense) - 2, 3)]
        codons_asense_1 =[dna_seq_asense[i:i+3] for i in range(0, len(dna_seq_asense) - 2, 3)]
        codons_sense_2 = [dna_seq_sense[i:i+3] for i in range(1, len(dna_seq_sense) - 2, 3)]
        codons_asense_2 = [dna_seq_asense[i:i+3] for i in range(1, len(dna_seq_asense) - 2, 3)]
        codons_sense_3 = [dna_seq_sense[i:i+3] for i in range(2, len(dna_seq_sense) - 2, 3)]
        codons_asense_3 = [dna_seq_asense[i:i+3] for i in range(2, len(dna_seq_asense) - 2, 3)]

    # Create empty array with shape 2 strands, 3 frames, n codons
    no_codons = len(codons_sense_1)
    codons = np.empty((2, 3, no_codons), dtype=object)

    # Add lists into the empty array
    codons[0, 0, :] = codons_sense_1
    codons[0, 1, :len(codons_sense_2)] = codons_sense_2
    codons[0, 2, :len(codons_sense_3)] = codons_sense_3
    codons[1, 0, :] = codons_asense_1
    codons[1, 1, :len(codons_asense_2)] = codons_asense_2
    codons[1, 2, :len(codons_asense_3)] = codons_asense_3
    #print(codons)
    # Note: lists have different length and frame 1 is longest possible, some array rows will have dummy element

    # Find all start and stop codons
    # this generates an array with same shape as codons but with True or False values depending on match
    starts = np.isin(codons, start_codon)
    stops = np.isin(codons, stop_codon)

    for strand_i in range(2): # 0: sense, 1: antisense
        for frame_i in range(3): # 0: frame 1, 1: frame 2, 2: frame 3
            start_positions = np.where(starts[strand_i, frame_i])[0]
            stop_positions = np.where(stops[strand_i, frame_i])[0]

            for the_start in start_positions:
                valid_stops = stop_positions[stop_positions > the_start]

                if len(valid_stops) > 0:
                    the_stop = valid_stops[0]
                    
                    orf_id = f'{dna_id.lower()}_orf_{len(orf_dict['Putative ORF ID'])+1}'
                    orf_dict['Source DNA ID'].append(dna_id)
                    orf_dict['Putative ORF ID'].append(orf_id)
                    orf_dict['Strand'].append('Sense' if strand_i == 0 else 'Antisense')
                    #orf_dict['Frame'].append(f'+{frame_i+1}') 
                    orf_dict['Start codon index'].append(the_start)
                    orf_dict['Stop codon index'].append(the_stop)

    return pd.DataFrame(orf_dict)

# like orf_finder_npar but input is dictionary
def orf_finder_lvl3_v2(dna_dict):
    """Find putative ORFs in DNA sequences provided in dictionary,
    and return dataframe with location data for putative ORFs and unique putative ORF ID.
    
    Input: dictionary of sequences  
    Output: dataframe with putative ORF location data"""

    # Set start and stop codons
    start_codon = ['ATG']
    stop_codon = ['TAA', 'TAG', 'TGA']
    
    # Create empty dictionary to store locations of putative ORFs
    orf_dict = {
        'Source DNA ID': [],
        'Putative ORF ID': [],
        'Strand': [],
        #'Frame': [],
        'Start codon index': [],
        'Stop codon index': []
    }
    
    # Prepare DNA 
    for seq_id, seq in dna_dict.items():
        dna_id = seq_id
        dna_seq_sense = str(seq)
        dna_seq_asense = str(Seq(dna_seq_sense).reverse_complement())
            
        # Create lists of all codon frames
        codons_sense_1 = [dna_seq_sense[i:i+3] for i in range(0, len(dna_seq_sense) - 2, 3)]
        codons_asense_1 =[dna_seq_asense[i:i+3] for i in range(0, len(dna_seq_asense) - 2, 3)]
        codons_sense_2 = [dna_seq_sense[i:i+3] for i in range(1, len(dna_seq_sense) - 2, 3)]
        codons_asense_2 = [dna_seq_asense[i:i+3] for i in range(1, len(dna_seq_asense) - 2, 3)]
        codons_sense_3 = [dna_seq_sense[i:i+3] for i in range(2, len(dna_seq_sense) - 2, 3)]
        codons_asense_3 = [dna_seq_asense[i:i+3] for i in range(2, len(dna_seq_asense) - 2, 3)]

    # Create empty array with shape 2 strands, 3 frames, n codons
    no_codons = len(codons_sense_1)
    codons = np.empty((2, 3, no_codons), dtype=object)

    # Add lists into the empty array
    codons[0, 0, :] = codons_sense_1
    codons[0, 1, :len(codons_sense_2)] = codons_sense_2
    codons[0, 2, :len(codons_sense_3)] = codons_sense_3
    codons[1, 0, :] = codons_asense_1
    codons[1, 1, :len(codons_asense_2)] = codons_asense_2
    codons[1, 2, :len(codons_asense_3)] = codons_asense_3
    #print(codons)
    # Note: lists have different length and frame 1 is longest possible, some array rows will have dummy element

    # Find all start and stop codons
    # this generates an array with same shape as codons but with True or False values depending on match
    starts = np.isin(codons, start_codon)
    stops = np.isin(codons, stop_codon)

    for strand_i in range(2): # 0: sense, 1: antisense
        for frame_i in range(3): # 0: frame 1, 1: frame 2, 2: frame 3
            start_positions = np.where(starts[strand_i, frame_i])[0]
            stop_positions = np.where(stops[strand_i, frame_i])[0]

            for the_start in start_positions:
                valid_stops = stop_positions[stop_positions > the_start]

                if len(valid_stops) > 0:
                    the_stop = valid_stops[0]
                    
                    orf_id = f'{dna_id.lower()}_orf_{len(orf_dict['Putative ORF ID'])+1}'
                    orf_dict['Source DNA ID'].append(dna_id)
                    orf_dict['Putative ORF ID'].append(orf_id)
                    orf_dict['Strand'].append('Sense' if strand_i == 0 else 'Antisense')
                    #orf_dict['Frame'].append(f'+{frame_i+1}') 
                    orf_dict['Start codon index'].append(the_start)
                    orf_dict['Stop codon index'].append(the_stop)

    return pd.DataFrame(orf_dict)

# Like orf_finder_npar_v1 but reducing nested loops
def orf_finder_npar_v3(seqrecords):
    """Find putative ORFs in DNA sequences provided in list of SeqRecords,
    and return dataframe with location data for putative ORFs and unique putative ORF ID.
    
    Input: list of SeqRecord objects  
    Output: dataframe with putative ORF location data"""

    # Set start and stop codons
    start_codon = ['ATG']
    stop_codon = ['TAA', 'TAG', 'TGA']
    
    # Create empty dictionary to store locations of putative ORFs
    orf_dict = {
        'Source DNA ID': [],
        'Putative ORF ID': [],
        'Strand': [],
        #'Frame': [],
        'Start codon index': [],
        'Stop codon index': []
    }
    
    # Prepare DNA 
    for record in seqrecords:
        dna_id = record.id
        dna_seq_sense = str(record.seq)
        dna_seq_asense = str(Seq(dna_seq_sense).reverse_complement())
            
        # Create lists of all codon frames
        codons_sense_1 = [dna_seq_sense[i:i+3] for i in range(0, len(dna_seq_sense) - 2, 3)]
        codons_asense_1 =[dna_seq_asense[i:i+3] for i in range(0, len(dna_seq_asense) - 2, 3)]
        codons_sense_2 = [dna_seq_sense[i:i+3] for i in range(1, len(dna_seq_sense) - 2, 3)]
        codons_asense_2 = [dna_seq_asense[i:i+3] for i in range(1, len(dna_seq_asense) - 2, 3)]
        codons_sense_3 = [dna_seq_sense[i:i+3] for i in range(2, len(dna_seq_sense) - 2, 3)]
        codons_asense_3 = [dna_seq_asense[i:i+3] for i in range(2, len(dna_seq_asense) - 2, 3)]

    # Create empty array with shape 2 strands, 3 frames, n codons
    no_codons = len(codons_sense_1)
    codons = np.empty((2, 3, no_codons), dtype=object)

    # Add lists into the empty array
    codons[0, 0, :] = codons_sense_1
    codons[0, 1, :len(codons_sense_2)] = codons_sense_2
    codons[0, 2, :len(codons_sense_3)] = codons_sense_3
    codons[1, 0, :] = codons_asense_1
    codons[1, 1, :len(codons_asense_2)] = codons_asense_2
    codons[1, 2, :len(codons_asense_3)] = codons_asense_3
    #print(codons)
    # Note: lists have different length and frame 1 is longest possible, some array rows will have dummy element

    # Find all start and stop codons
    # this generates an array with same shape as codons but with True or False values depending on match
    starts = np.isin(codons, start_codon)
    stops = np.isin(codons, stop_codon)

    for strand_i in range(2): # 0: sense, 1: antisense
        for frame_i in range(3): # 0: frame 1, 1: frame 2, 2: frame 3
            start_positions = np.where(starts[strand_i, frame_i])[0]
            stop_positions = np.where(stops[strand_i, frame_i])[0]

            for the_start in start_positions:
                valid_stops = stop_positions[stop_positions > the_start]

                if len(valid_stops) > 0:
                    the_stop = valid_stops[0]
                    
                    orf_id = f'{dna_id.lower()}_orf_{len(orf_dict['Putative ORF ID'])+1}'
                    orf_dict['Source DNA ID'].append(dna_id)
                    orf_dict['Putative ORF ID'].append(orf_id)
                    orf_dict['Strand'].append('Sense' if strand_i == 0 else 'Antisense')
                    #orf_dict['Frame'].append(f'+{frame_i+1}') 
                    orf_dict['Start codon index'].append(the_start)
                    orf_dict['Stop codon index'].append(the_stop)

    return pd.DataFrame(orf_dict)


# The new fasta parser is not better but the orf finder is maybe 2x better!
# Actually the level 3 orf finder is only 2x better when dealing with large number of DNA sequences input, 
# but worse at longer single sequences 
# Maybe I can use numpy arrays to solve the overlaps problem
# Maybe for the fasta writing section it is still worth keeping seqrecords for the fasta parser
# Maybe write code that plots top 1% of ORF lengths    