# Import packages
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import sys

sys.path.append('C:\\Users\\joanaq\\Documents\\learning-bioinformatics\\projects\\ros-orf\\packages')
from random_seq import seqrecord_fasta_writer


# PART 1: PREPARE DNA INPUT
    
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

# like orf_finder_lvl2 but takes seqrecords as input
def orf_finder_lvl2_v2(seqrecords):
    """Find putative ORFs in DNA sequences provided in list of SeqRecords,
    and return dataframe with location data for putative ORFs and unique putative ORF ID.
    Note that start and stop codon location is stored as 0-based python index of the first codon nt.
    
    Input: List of DNA SeqRecords    
    Output: dataframe with putative ORF location data"""

    # Set start and stop codons
    start_codon = ['ATG']
    stop_codon = ['TAA', 'TAG', 'UGA']

    # Create empty dictionary to store locations of putative ORFs
    orf_dict = {
        'Source DNA ID': [],
        'Putative ORF ID': [],
        'Strand': [],
        'Start codon index': [],
        'Stop codon index': []
    }

    # Prepare DNA from dataframe 
    for dna in seqrecords:
        dna_id = dna.id
        dna_seq_sense = dna.seq
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
# WARNING: IT IS WRONG! DOES NOT FIND ORFs CORRECTLY!!
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

# like orf_finder_lvl3_v1 but input is dictionary, so it is compatible with fasta_parser_lvl2
# WARNING: IT IS WRONG! DOES NOT FIND ORFs CORRECTLY!!
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

# Like orf_finder_npar_v1 but trying further improvements
def orf_finder_npar_v3(seqrecords):
    """Find putative ORFs in DNA sequences provided in list of SeqRecords,
    and return dataframe with location data for putative ORFs and unique putative ORF ID.
    
    Input: list of SeqRecord objects  
    Output: dataframe with putative ORF location data"""

    return



# PART 3: ASSEMBLE DATAFRAME WITH MORE PUTATIVE ORF DATA

def orf_data_collector(dna_dict, orf_locations_df):
    """Gather more data regarding putative ORFs and store it in more comprehensive dataframe.
    This data can be useful for subsequent analysis and plotting.
    Examples of data to collect:
    Location of first and last bp of putative ORF in the DNA duplex (First bp and Last bp),
    an alternative location metric similar to the aforementioned but retaining orientation of the ORF (Start and Stop),
    length of the mRNA (Length of mRNA), size of the peptide (Peptide size),
    whether putative ORF overlaps with others in the DNA duplex (Overlap in DNA)...

    Input: dataframe with locations (0-based index) of putative ORFs
    Output: comprehensive dataframe with more data on putative ORFs"""
    
    # Create expanded dataframe
    orf_extra_df = orf_locations_df.copy()

    # Initiate list for data variables:
    fromhere_dna, tohere_dna = [], [] # First and last bp of ORF in DNA duplex
    start_dna, stop_dna = [], [] # Location of start and stop codons in DNA duplex
    frame = []
    length, pept_size = [], []

    # Collect that data based on strand and position indexes of putative ORF
    for index, row in orf_locations_df.iterrows():
        start_i = row['Start codon index']
        stop_i = row['Stop codon index']
        strand = row['Strand']
        dna_id = row['Source DNA ID']
        dna_seq = dna_dict[dna_id]
        dna_len = len(dna_seq)

        # Calculate size of putative ORF and peptide
        length.append(stop_i - start_i)
        pept_size.append((stop_i - start_i) // 3)

        # Calculate frame (strand dependent)
        if (start_i + 3) % 3 == 0:
            frame.append(f'{strand} +1')
        elif (start_i + 2) % 3 == 0:
            frame.append(f'{strand} +2')
        elif (start_i + 1) % 3 == 0:
            frame.append(f'{strand} +3')

        # ORF location
        if strand == 'Sense':
            fromhere_dna.append(start_i + 1)
            tohere_dna.append(stop_i +3) 
            start_dna.append(start_i + 1)
            stop_dna.append(stop_i +3)
        elif strand == 'Antisense':
            fromhere_dna.append(dna_len - stop_i - 2)
            tohere_dna.append(dna_len - start_i)
            stop_dna.append(dna_len - stop_i - 2)
            start_dna.append(dna_len - start_i)
           
    # Add variables to dataframe
    orf_extra_df['Frame'] = frame
    orf_extra_df['Length'] = length
    orf_extra_df['Peptide size'] = pept_size 
    orf_extra_df['ORF first bp'] = fromhere_dna
    orf_extra_df['ORF last bp'] = tohere_dna
    orf_extra_df['Start bp'] = start_dna
    orf_extra_df['Stop bp'] = stop_dna

    return orf_extra_df     




# PART 4: GENERATE FASTA FILES WITH PUTATIVE ORFS

def orf_seqrecorder(dna_input, orf_df, seq_type, fasta_output=False):
    """
    Input:
    dna_input (list): List of SeqRecords with input DNA
    orf_df (df): Dataframe with putative ORF IDs and locations
    seq_type (str): Specify sequence type for output (DNA, RNA or PEPTIDE)
    fasta_output (bool): Indicate whether you want to output a FASTA file 


    Output:
    List of putative ORF SeqRecords as DNA, RNA or peptide sequences.
    FASTA files (optional)
    """
    # Start empty dictionary for putative ORFs
    orfs = {}

    # Prepare input DNA sequences
    for dna in dna_input:
        dna_id = dna.id
        dna_sense = dna.seq
        dna_asense = str(Seq(dna_sense).reverse_complement())

        # For each DNA sequence...
        filt_orf_df = orf_df[orf_df['Source DNA ID'] == dna_id]

        # Store putative ORFs in the list
        for _, row in filt_orf_df.iterrows():
            if row['Strand'] == 'Sense':
                orf_seq = dna_sense[row['Start codon index']:(row['Stop codon index']+3)]
            elif row['Strand'] == 'Antisense':
                orf_seq = dna_asense[row['Start codon index']:(row['Stop codon index']+3)]
            orfs[row['Putative ORF ID']] = orf_seq
            
        # Create list of SeqRecords with DNA, RNA or peptide sequences
        orf_seqrecords = []
        for orf_id, orf_seq in orfs.items():
            if seq_type.upper() == 'DNA':
                orf_rec = SeqRecord(seq=Seq(orf_seq),
                                    id=orf_id,
                                    description=f"Putative ORF from DNA {row['Source DNA ID']}, {row['Strand']} strand")
                
            elif seq_type.upper() == 'RNA':
                orf_rec = SeqRecord(seq=str(Seq(orf_seq).transcribe()),
                                    id=orf_id,
                                    description=f"Putative mRNA from DNA {row['Source DNA ID']}")
        
            elif seq_type.upper() == 'PEPTIDE':
                orf_rec = SeqRecord(seq=str(Seq(orf_seq).translate()),
                                    id=orf_id,
                                    description=f"Putative peptide from DNA {row['Source DNA ID']}")
            else:
                print('Something is wrong')
            
            orf_seqrecords.append(orf_rec)
        
        if fasta_output:
            fasta_file = f"putativeorfs-{dna_id}-{seq_type.lower()}.fasta"
            seqrecord_fasta_writer(orf_seqrecords, fasta_file)
        
    return orf_seqrecords

#os.chdir('C:/Users/joanaq/Documents/learning-bioinformatics/projects/ros-orf/level-3')
#dna_seqrec = fasta_parser_lvl3('sample-dna-fasta/rosorf_benchm_n10_l1.fasta')
#orf_df = orf_finder_lvl2_v2(dna_seqrec)
#orf_seqrecorder(dna_input=dna_seqrec, orf_df=orf_df, seq_type='peptide', fasta_output=True)
