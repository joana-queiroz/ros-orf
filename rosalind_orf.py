### ORF: Open Reading Frame ###

# Problem
# https://rosalind.info/problems/orf/

# Given:
# A DNA string s of length at most 1 kbp in FASTA format.

# Return:
# Every distinct candidate protein string that can be translated
# from ORFs of s. Strings can be returned in any order.

## My answer ##

codon_table = {'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L', 'UCU':'S',
               'UCC':'S', 'UCA':'S', 'UCG':'S', 'UAU':'Y', 'UAC':'Y',
               'UAA':'stop', 'UAG':'stop', 'UGU':'C', 'UGC':'C', 'UGA':'stop',
               'UGG':'W', 'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
               'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAU':'H',
               'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGU':'R', 'CGC':'R',
               'CGA':'R', 'CGG':'R', 'AUU':'I', 'AUC':'I', 'AUA':'I',
               'AUG':'M', 'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
               'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGU':'S',
               'AGC':'S', 'AGA':'R', 'AGG':'R', 'GUU':'V', 'GUC':'V',
               'GUA':'V', 'GUG':'V', 'GCU':'A', 'GCC':'A', 'GCA':'A',
               'GCG':'A', 'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
               'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

with open('rosalind_orf.txt','r') as file:
    infile = file.read().split('\n')
    DNA = ''.join(infile[1:])
    RNA_sense = DNA.replace('T','U')
    RNA_asense = RNA_sense.replace('A','u').replace('U','a').replace('C','g').replace('G','c').upper()[::-1]
    RNA= {'RNA_sense' : RNA_sense,
          'RNA_asense' : RNA_asense}

#print(infile)
#print(DNA)
#print(RNA_sense)
#print(RNA_antisense)
#print(RNA)

# Create dictionary with codons of all 6 possible frames:
# sense_1, sense_2, sense_3, asense_1, asense_2, asense_3

frames = {}
for strand, seq in RNA.items():
    for frame_start in range(3):
        key = f'{strand}_{frame_start+1}'
        frames[key] = [seq[i:i+3] for i in range(frame_start,len(seq)-2, 3)]


# Find all start codons and stop codons positions in frames dictionary
# Store positions in new dictionary, print to see

start_codon = 'AUG'
stop_codon = {'UAA','UAG','UGA'}


start_i = {}
stop_i = {}
for frame, codons in frames.items():
    start_i[frame]= [i for i, codon in enumerate(codons) if codon == start_codon]
    stop_i[frame]= [i for i, codon in enumerate(codons) if codon in stop_codon]
#print(start_i)
#print(stop_i)


# Select all ORF from start to stop codons w/o an intermediate stop codon
# Translate the ORF codons into a peptide
# Add peptide to respective frame

ORFs = {}
for frame, codons in frames.items():
    starts = start_i[frame]
    stops = stop_i[frame]

    for start in starts:
        for stop in stops:
            if stop > start:
                mid_stops = [pos for pos in stops if start < pos < stop]
                if not mid_stops:
                    orf_codons = codons[start:stop]
                    peptide = ''.join(codon_table[codon] for codon in orf_codons)
                    if frame not in ORFs:
                        ORFs[frame] = []
                    ORFs[frame].append(peptide)
                    break

# Write file with all the unique possible peptides

unique_orfs = set()
for orf_list in ORFs.values():
    unique_orfs.update(orf_list)

with open('rosalind_orf_solv.txt','w') as file:
    for unique_orf in unique_orfs:
        file.write(unique_orf +'\n')
    

#for orf in ORFs.values():
#    print(*orf, sep='\n')


## Learning with sample ##

# Sample dataset
# >Rosalind_99
# AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG

# Sample output:
# MLLGSFRLIPKETLIQVAGSSPCNLS
# M
# MGMTPRLGLESLLE
# MTPRLGLESLLE



## Me trying things ##
        
#start_codon_i = {}
#for frame, codons in frames.items():
#    for codon_i in range(len(codons)):
#        codon = codons[codon_i]
#        if codon == start_codon:
#            start_codon_i[frame] += ' '.join(codon_i)
#print(start_codon_i)
            
        




#start_codon = 'AUG'
#stop_codon = {'UAA','UAG','UGA'}

#ORFs = []
#for frame, codons in frames.items():
#    i = 0
#    while i < len(codons):
#        if codons[i] == start_codon:
#            AAs = []
#            for codon_i in range(i,len(codons)):
#                codon = codons[codon_i]
#                if codon in stop_codon:
#                    if AAs:
#                        ORFs.append(''.join(AAs))
#                    break
#                elif codon in codon_table:
#                    AAs.append(codon_table[codon])
#            i = codon_i
        #else:
        #    i += 1
#print(AAs)
#print(ORFs)
            

#for frame,codons in frames.items():
#    for codon_i in range(0,len(codons)):
#        codon = codons[codon_i]
#        #print(codon)
#        if codon == start_codon:
#            ORF.append(''.join(codon_table[codon]))
#            start_codon_i = codon_i
#            for next_codon_i in range(start_codon_i,len(codons)):
#                next_codon = codons[next_codon_i]
#                if codon_table[next_codon] == 'stop':
#                    break
#                elif next_codon in codon_table:
#                    ORF.append(''.join(codon_table[next_codon]))
#                    print(ORF)
#                else:
#                    print('Error: no codon match')

# frames = {}
#f_nt_start = 0
#f_nt_end = 3
#codons = []

#for f_nt in range(f_nt_start,f_nt_end): # 0, 1, 2
#    for nt in range(f_nt,len(RNA_sense),3):
#        key = f'sense_{f_nt+1}'
#        if len(DNA[nt:nt+3]) == 3:
#            codons.append(RNA_sense[nt:nt+3])
#            frames[key] = codons

#for f_nt in range(f_nt_start,f_nt_end): # 0, 1, 2
#    for nt in range(f_nt,len(RNA_asense),3):
#        key = f'asense_{f_nt+1}'
#        if len(DNA[nt:nt+3]) == 3:
#            codons.append(RNA_asense[nt:nt+3])
#            frames[key] = codons


#ORF = []

#start_codon = 'AUG'
#start_codon_i = 0

#for frame,codons in frames.items():
#    for codon_i in range(0,len(codons)):
#        codon = codons[codon_i]
#        #print(codon)
#        if codon == start_codon:
#            ORF.append(''.join(codon_table[codon]))
#            start_codon_i = codon_i
#            for next_codon_i in range(start_codon_i,len(codons)):
#                next_codon = codons[next_codon_i]
#                if codon_table[next_codon] == 'stop':
#                    break
#                elif next_codon in codon_table:
#                    ORF.append(''.join(codon_table[next_codon]))
#                    print(ORF)
#                else:
#                    print('Error: no codon match')
            




#for nt_i in range(len(RNA) -6 +1):
#    for codon
#        if codon == 'AUG':
#            AA += codon_table[codon]
#            if codon_table[codon] == 'stop':
#                break
#            elif codon in codon_table:
#                AA += codon_table[codon]
#            else:
#                print('Error: no codon match')
    
    
