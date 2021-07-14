import numpy as np
import math
import fasta_parser
import csv
from Bio import SeqIO
from sys import argv
import os

ORFs_fasta = argv[1]
effectors_f = argv[2] #fasta protein
working_directory = argv[3]
non_effectors_f = 'non_effectors.faa'
os.chdir(working_directory)
locus_dic=fasta_parser.parse_ORFs(ORFs_fasta)
###############################################################################

def GC_content(gene):
    ''' returns a portion of GC in the seq '''
    gene = str(gene).upper()
    G_count=gene.count('G')
    C_count=gene.count('C')
    tot_count=G_count+C_count
    if len(gene) != 0:
        return tot_count/len(gene)

def abundance_(type_of_aa, peptide):
    ''' returns the portion of this aa type in the peptide '''
    amount = 0
    for aa in type_of_aa:
        amount += peptide.count(aa)
    return amount/len(peptide)

def interaction_value(interaction, peptide):
    ''' returns the interaction value of the peptide, e.g., hydrophobicity value'''
    value = 0
    for aa in interaction:
        value += peptide.count(aa)*interaction[aa]
    return value/len(peptide)

def aa_profile(peptide):
    ''' returns a vector with a sum of 1 that represents the aa profile of the peptide '''
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_array = np.zeros(20)
    for i in range(20):
        aa_array[i] = abundance_(amino_acids[i], peptide)
    return aa_array

def distance(v1, v2):
    ''' returns the euclidean distance between two vectors '''
    return math.sqrt(sum([(a-b)**2 for a,b in zip(v1,v2)]))

effectors=SeqIO.parse(effectors_f,'fasta')
effectors_profiles = {}
for effector in effectors:
    locus = effector.id
    seq = effector.seq
    effectors_profiles[locus] = aa_profile(seq)

effectors_aa_profile_array = np.array(list(effectors_profiles.values()))
effectors_average_aa_profile = np.average(effectors_aa_profile_array,axis=0)

non_effectors=SeqIO.parse(non_effectors_f,'fasta')
non_effectors_profiles=[]
for non_effector in non_effectors:
    seq=non_effector.seq
    non_effectors_profiles.append(aa_profile(seq))

non_effectors_profile_array = np.array(non_effectors_profiles)
non_effectors_average_profile=np.average(non_effectors_profile_array,axis=0)
    
def similarity_to_effectors_vs_non_effectors(locus):
    peptide = locus_dic[locus].translate(to_stop=True)
    aa_profile_pep = aa_profile(peptide)
    if locus in effectors_profiles:
        effectors_profiles_copy = effectors_profiles.copy()
        effectors_profiles_copy.pop(locus)
        effectors_aa_profile_array_ = np.array(list(effectors_profiles_copy.values()))
        effectors_average_aa_profile_ = np.average(effectors_aa_profile_array_,axis=0)
        dis_from_effectots=distance(aa_profile_pep,effectors_average_aa_profile_)
    else:
        dis_from_effectots=distance(aa_profile_pep,effectors_average_aa_profile)
    dis_from_non_effectors=distance(aa_profile_pep,non_effectors_average_profile)
    return dis_from_effectots-dis_from_non_effectors

## hydrophilicity/ amphiphilicity/ hydrophobicity indexes, for interaction_value     
BLAS910101 = {'A':0.616, 'L':0.943, 'R':0, 'K':0.283, 'N':0.236, 'M':0.738,\
                  'D':0.028, 'F':1, 'C':0.68, 'P':0.711, 'Q':0.251, 'S':0.359,\
                  'E':0.043, 'T':0.45, 'G':0.501, 'W':0.878, 'H':0.165, 'Y':0.88,\
                  'I':0.943, 'V':0.825}
MITS020101 = {'A':0, 'L':0, 'R':2.45, 'K':3.67, 'N':0, 'M':0, 'D':0, 'F':0,\
              'C':0, 'P':0, 'Q':1.25, 'S':0, 'E':1.27, 'T':0, 'G':0,\
              'W':6.93, 'H':1.45, 'Y':5.06, 'I':0, 'V':0}        
KUHL950101 = {'A':0.78, 'L':0.56, 'R':1.58, 'K':1.1, 'N':1.2, 'M':0.66,\
              'D':1.35, 'F':0.47, 'C':0.55, 'P':0.69, 'Q':1.19, 'S':1,\
              'E':1.45, 'T':1.05, 'G':0.68, 'W':0.7, 'H':0.99, 'Y':1,\
              'I':0.47, 'V':0.51}
CIDH920105 = {'A':0.02, 'L':1.14, 'R':-0.42, 'K':-0.41, 'N':-0.77, 'M':1.00,\
              'D':-1.04, 'F':1.35, 'C':0.77, 'P':-0.09, 'Q':-1.10, 'S':-0.97,\
              'E':-1.14, 'T':-0.77, 'G':-0.80, 'W':1.71, 'H':0.26, 'Y':1.11,\
              'I':1.81, 'V':1.13}


amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

with open(f'physical_features.csv','w',newline='') as file:
    file_writer = csv.writer(file)
    header=['locus','GC_content','protein_length','similarity_to_effectors_vs_non_effectors']
    for aa in amino_acids:
        header.append(f'{aa}_full_protein')
        header.append(f'{aa}_N-terminus')
    header.append('hydrophobicity_N-terminus_BLAS910101')
    header.append('amphiphilicity_N-terminus_MITS020101')
    header.append('hydrophilicity_N-terminus_KUHL950101')
    header.append('hydrophobicity_N-terminus_CIDH920105')
    file_writer.writerow(header)
    for locus in locus_dic:
        l=[locus]
        l.append(GC_content(locus_dic[locus]))
        l.append(len(locus_dic[locus].translate(to_stop=True)))
        l.append(similarity_to_effectors_vs_non_effectors(locus))
        for aa in amino_acids:
            l.append(abundance_(aa, locus_dic[locus].translate(to_stop=True)))
            l.append(abundance_(aa, locus_dic[locus].translate(to_stop=True)[:25]))
        l.append(interaction_value(BLAS910101, locus_dic[locus].translate(to_stop=True)[:25]))
        l.append(interaction_value(MITS020101, locus_dic[locus].translate(to_stop=True)[:25]))
        l.append(interaction_value(KUHL950101, locus_dic[locus].translate(to_stop=True)[:25]))
        l.append(interaction_value(CIDH920105, locus_dic[locus].translate(to_stop=True)[:25]))
        file_writer.writerow(l)
        
endfile = open(f'physical_features.done','w')
endfile.close()