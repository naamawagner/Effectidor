from sys import argv
#import fasta_parser
from Bio import SeqIO
import os
import subprocess

all_prots = argv[1]
effectors_prots = argv[2]
k12_dataset= r'/groups/pupko/naamawagner/TF/Genomes/e.coli_k-12_substr.MG1655/e.coli_k-12_substr.MG1655.proteins.faa'
T3SS_dataset = r'/groups/pupko/naamawagner/T3Es_webserver/T3SS_proteins.faa'
if not os.path.exists('blast_outputs'):
    os.makedirs('blast_outputs')
if not os.path.exists('blast_data/non_effectors'):
    os.makedirs('blast_data/non_effectors')
cmd = f'cp {k12_dataset} ./blast_data/non_effectors'
subprocess.check_output(cmd,shell=True)
cmd = f'cp {T3SS_dataset} ./blast_data/non_effectors'
subprocess.check_output(cmd,shell=True)

k12_dataset = './blast_data/non_effectors/e.coli_k-12_substr.MG1655.proteins.faa'
k12_out_file=r'blast_outputs/e_coli.blast'
T3SS_dataset = './blast_data/non_effectors/T3SS_proteins.faa'
T3SS_out_file = r'blast_outputs/T3SS_prots.blast'


def protein_blast_all_vs_all(query,dataset,out,e_val='0.01'):
    import subprocess
    make_data_cmd=f"makeblastdb -in {dataset} -dbtype prot -out {dataset[:-4]}db"
    make_blast_cmd=f'blastp -db {dataset[:-4]}db -query {query} -outfmt 6 -out {out} -evalue {e_val}'
    subprocess.check_output(make_data_cmd,shell=True)
    subprocess.check_output(make_blast_cmd,shell=True)
    
def parse_blast_out(blast_out,e_val=0.01,min_coverage=0):
    blast_out_dic={}
    with open(blast_out,'r') as in_f:
        for row in in_f:
            row=row.split('\t')
            prot_id=row[0]
            e_value=float(row[-2])
            coverage=float(row[2])
            if e_value<=e_val and coverage>=min_coverage:
                if prot_id not in blast_out_dic:
                    blast_out_dic[prot_id]=[1,float(row[-1])]
                else:
                    blast_out_dic[prot_id][0]+=1
    return blast_out_dic
    

prots_dict=SeqIO.to_dict(SeqIO.parse(all_prots,'fasta'))
effectors_dict = SeqIO.to_dict(SeqIO.parse(effectors_prots,'fasta'))

protein_blast_all_vs_all(all_prots,k12_dataset,k12_out_file)
k12_dict = parse_blast_out(k12_out_file,e_val=10**(-6))
k12_out_list=sorted(k12_dict.keys())

protein_blast_all_vs_all(all_prots,T3SS_dataset,T3SS_out_file)
T3SS_dict = parse_blast_out(T3SS_out_file,e_val=10**(-10))
T3SS_out_list=sorted(T3SS_dict.keys())

non_effectors_recs=[]
T3SS_recs = []
for prot in prots_dict:
    #if (prot in k12_out_list or prot in T3SS_out_list) and prot not in effectors_dict:
    if prot not in effectors_dict:
        if prot in k12_out_list and prot not in T3SS_out_list:
            non_effectors_recs.append(prots_dict[prot])
        elif prot not in k12_out_list and prot in T3SS_out_list:
            T3SS_recs.append(prots_dict[prot])
        elif prot in k12_out_list and prot in T3SS_out_list:
            if k12_dict[prot][1] >= T3SS_dict[prot][1]:
                non_effectors_recs.append(prots_dict[prot])
            else:
                T3SS_recs.append(prots_dict[prot])
    '''            
    if prot in k12_out_list and prot not in effectors_dict:
        non_effectors_recs.append(prots_dict[prot])
    elif prot in T3SS_out_list and prot not in effectors_dict:
        T3SS_recs.append(prots_dict[prot])
    '''
SeqIO.write(non_effectors_recs,'non_effectors.faa','fasta')
SeqIO.write(T3SS_recs,'T3SS_proteins.faa','fasta')