from sys import argv
from Bio import SeqIO
import os
import subprocess
import csv
import sys
import effectidor_CONSTANTS
from protein_mmseq import protein_mmseqs_all_vs_all

all_prots = argv[1]
effectors_prots = argv[2]
data_dir = effectidor_CONSTANTS.EFFECTIDOR_DATA
k12_dataset = f'{data_dir}/e.coli_k-12_substr.MG1655.proteins.faa'
if not os.path.exists('blast_outputs'):
    os.makedirs('blast_outputs')
if not os.path.exists('blast_data/non_effectors'):
    os.makedirs('blast_data/non_effectors')
cmd = f'cp {k12_dataset} ./blast_data/non_effectors'
subprocess.check_output(cmd, shell=True)

k12_dataset = './blast_data/non_effectors/e.coli_k-12_substr.MG1655.proteins.faa'
k12_out_file = r'blast_outputs/e_coli.blast'


def parse_blast_out(blast_out, e_val=0.01, min_coverage=0):
    blast_out_dic = {}
    with open(blast_out, 'r') as in_f:
        for row in in_f:
            row = row.split('\t')
            prot_id = row[0]
            e_value = float(row[-2])
            coverage = float(row[2])
            if e_value <= e_val and coverage >= min_coverage:
                if prot_id not in blast_out_dic:
                    blast_out_dic[prot_id] = [[row[1]], float(row[-1])]
                else:
                    if row[1] not in blast_out_dic[prot_id][0]:
                        blast_out_dic[prot_id][0].append(row[1])
    return blast_out_dic


prots_dict = SeqIO.to_dict(SeqIO.parse(all_prots, 'fasta'))
if not os.path.exists('tmp_mmseq'):
    os.makedirs('tmp_mmseq')
protein_mmseqs_all_vs_all(all_prots, k12_dataset, k12_out_file, 'tmp_mmseq')
k12_dict = parse_blast_out(k12_out_file, e_val=10**(-6))

if os.path.exists('blast_outputs/effectors_hits.csv'):  # if the effectors were found based on mmseqs, we need to
    # verify their hit to effectors was better than their hit to non-effectors
    effectors = {}
    with open('blast_outputs/effectors_hits.csv') as in_f:
        reader = csv.reader(in_f)
        for row in reader:
            effectors[row[0]] = float(row[1])
    effectors_to_keep = []
    for effector in effectors:
        if effector not in k12_dict:
            effectors_to_keep.append(effector)
        elif effectors[effector] > k12_dict[effector][1]:
            effectors_to_keep.append(effector)
    if len(effectors_to_keep) < len(effectors.keys()):
        old_effectors = SeqIO.parse(effectors_prots, 'fasta')
        new_effectors = [effector for effector in old_effectors if effector.id in effectors_to_keep]
        SeqIO.write(new_effectors, effectors_prots, 'fasta')
        
k12_out_list = list(k12_dict.keys())

effectors_dict = SeqIO.to_dict(SeqIO.parse(effectors_prots, 'fasta'))
non_effectors_recs = []
for prot in prots_dict:
    if prot not in effectors_dict:
        if prot in k12_out_list:
            non_effectors_recs.append(prots_dict[prot])

SeqIO.write(non_effectors_recs, 'non_effectors.faa', 'fasta')
