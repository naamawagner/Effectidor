from Bio import SeqIO
from sys import argv
import os

effectors = argv[1]
bacterial_proteome = argv[2]
effectors_prots = argv[3]
log_file = argv[4]

if not os.path.exists('blast_outputs'):
    os.makedirs('blast_outputs')

def protein_blast_all_vs_all(query, dataset, out, e_val='0.0001'):
    import subprocess
    make_data_cmd=f"makeblastdb -in {dataset} -dbtype prot -out {dataset[:-4]}db"
    make_blast_cmd=f'blastp -db {dataset[:-4]}db -query {query} -outfmt 6 -out {out} -evalue {e_val}'
    subprocess.check_output(make_data_cmd,shell=True)
    subprocess.check_output(make_blast_cmd,shell=True)
    
def parse_blast_out(blast_out,e_val=10**(-10),min_identity=70):
    blast_out_dic={}
    with open(blast_out,'r') as in_f:
        for row in in_f:
            row=row.split('\t')
            prot_id=row[0]
            hit = row[1]
            e_value=float(row[-2])
            identity=float(row[2])
            if e_value<=e_val and identity>=min_identity:
                if prot_id not in blast_out_dic:
                    blast_out_dic[prot_id]=hit
    return blast_out_dic

protein_blast_all_vs_all(effectors,bacterial_proteome,'blast_outputs/effectorsDB.blast')
effectors_homologs = []
k=70
while len(effectors_homologs)<5 and k>30:
    effectors_homologs = set(list(parse_blast_out('blast_outputs/effectorsDB.blast',min_identity=k).values()))
    with open(log_file,'a') as log:
        log.write(f'k:{k},n:{len(effectors_homologs)}\n{effectors_homologs}\n\n')
    k -= 10
    

recs = SeqIO.parse(bacterial_proteome,'fasta')
effectors_recs = []
for rec in recs:
    if rec.id in effectors_homologs:
        effectors_recs.append(rec)
SeqIO.write(effectors_recs,effectors_prots,'fasta')