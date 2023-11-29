from Bio import SeqIO
from sys import argv
import os
from protein_mmseq import protein_mmseqs_all_vs_all

effectors = argv[1]
bacterial_proteome = argv[2]
effectors_prots = argv[3]
log_file = argv[4]

if not os.path.exists('blast_outputs'):
    os.makedirs('blast_outputs')

    
def parse_blast_out(blast_out,e_val=10**(-10),min_identity=0.7):
    blast_out_dic={}
    with open(blast_out,'r') as in_f:
        best_hits={}
        for row in in_f:
            row=row.split('\t')
            prot_id=row[0]
            hit = row[1]
            e_value=float(row[-2])
            identity=float(row[2])
            if e_value<=e_val and identity>=min_identity:
                if prot_id not in blast_out_dic:
                    blast_out_dic[prot_id]=set([hit])
                else:
                    blast_out_dic[prot_id].add(hit)
                if hit not in best_hits:
                    best_hits[hit]=float(row[-1])
                else:
                    if float(row[-1]) > best_hits[hit]:
                        best_hits[hit]=float(row[-1])
        with open('blast_outputs/effectors_hits.csv','w') as out_f:
            for hit in best_hits:
                out_f.write(f'{hit},{best_hits[hit]}\n')
    return blast_out_dic

if not os.path.exists('tmp_mmseq'):
    os.makedirs('tmp_mmseq')
protein_mmseqs_all_vs_all(effectors,bacterial_proteome,'blast_outputs/effectorsDB.blast', 'tmp_mmseq')

k=0.6
effectors_homologs = set()
homologs = list(parse_blast_out('blast_outputs/effectorsDB.blast',min_identity=k).values())
for h in homologs:
    effectors_homologs = set.union(effectors_homologs,h)
    

recs = SeqIO.parse(bacterial_proteome,'fasta')
effectors_recs = []
for rec in recs:
    if rec.id in effectors_homologs:
        effectors_recs.append(rec)
SeqIO.write(effectors_recs,effectors_prots,'fasta')