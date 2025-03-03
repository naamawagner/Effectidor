import csv
from Bio import SeqIO
from sys import argv
import os
from protein_mmseq import protein_mmseqs_all_vs_all

bacterial_proteome = argv[1]
running_d = argv[2]
chaperones = argv[3]

if not os.path.exists(os.path.join(running_d, 'blast_outputs')):
    os.makedirs(os.path.join(running_d, 'blast_outputs'))

chaperones_recs = SeqIO.parse(chaperones, 'fasta')
chaperones_lens = {chaperone_rec.id: len(chaperone_rec.seq) for chaperone_rec in chaperones_recs}


def parse_blast_out(blast_out, lens_dic, e_val=10**(-10), min_identity=0.5, minimal_coverage=50.0):
    blast_out_dic = {}
    with open(blast_out, 'r') as in_f:
        best_hits = {}
        for row in in_f:
            row = row.split('\t')
            prot_id = row[0]
            hit = row[1]
            e_value = float(row[-2])
            identity = float(row[2])
            coverage = (float(row[3])/lens_dic[prot_id])*100
            prot_id = prot_id.replace('^', '|')
            if e_value <= e_val and identity >= min_identity and coverage >= minimal_coverage:
                if prot_id not in blast_out_dic:
                    blast_out_dic[prot_id] = [hit]
                else:
                    blast_out_dic[prot_id].append(hit)
                if hit not in best_hits:
                    best_hits[hit] = float(row[-1])
                else:
                    if float(row[-1]) > best_hits[hit]:
                        best_hits[hit] = float(row[-1])
        with open(os.path.join(running_d, 'blast_outputs/effectors_hits.csv'), 'w') as out_f:
            for hit in best_hits:
                out_f.write(f'{hit},{best_hits[hit]}\n')
    return blast_out_dic


protein_mmseqs_all_vs_all(chaperones, bacterial_proteome,
                          os.path.join(os.path.join(running_d, 'blast_outputs/chaperones.blast')),
                          os.path.join(os.path.join(running_d, 'tmp_mmseq')))

chaperones_homologs = parse_blast_out(os.path.join(running_d, 'blast_outputs/chaperones.blast'), chaperones_lens)
chaperons_d = {chaperone: chaperones_homologs[chaperone][0] for chaperone in chaperones_homologs}
with open(os.path.join(running_d, 'chaperones.csv'), 'w', newline='') as out_f:
    writer = csv.writer(out_f)
    header = ['chaperone', 'bacterial_protein']
    writer.writerow(header)
    for chaperone in chaperons_d:
        row = [chaperone, chaperons_d[chaperone]]
        writer.writerow(row)
