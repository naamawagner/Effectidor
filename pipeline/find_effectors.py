from Bio import SeqIO
from sys import argv
import os
from protein_mmseq import protein_mmseqs_all_vs_all

effectors = argv[1]
bacterial_proteome = argv[2]
effectors_prots = argv[3]
min_coverage = float(argv[4])
running_d = argv[5]

if not os.path.exists(os.path.join(running_d, 'blast_outputs')):
    os.makedirs(os.path.join(running_d, 'blast_outputs'))
in_f = open(effectors)
content = in_f.read()
replaced_pipes = content.replace('|', '^')  # MMSeq changes the IDs when pipes are available in the rec.id,
# so I replace it now with a character that is not expected to be found in the ID, to replace it back while parsing
# the MMSeq result.
in_f.close()
with open(effectors, 'w') as effectors_f:
    effectors_f.write(replaced_pipes)
effectors_recs = SeqIO.parse(effectors, 'fasta')
effectors_lens = {effector_rec.id: len(effector_rec.seq) for effector_rec in effectors_recs}
    

def parse_blast_out(blast_out, e_val=10**(-10), min_identity=0.7, minimal_coverage=50.0):
    blast_out_dic = {}
    with open(blast_out, 'r') as in_f:
        best_hits = {}
        for row in in_f:
            row = row.split('\t')
            prot_id = row[0]
            hit = row[1]
            e_value = float(row[-2])
            identity = float(row[2])
            coverage = (float(row[3])/effectors_lens[prot_id])*100
            prot_id = prot_id.replace('^', '|')
            if e_value <= e_val and identity >= min_identity and coverage >= minimal_coverage:
                if prot_id not in blast_out_dic:
                    blast_out_dic[prot_id] = {hit}
                else:
                    blast_out_dic[prot_id].add(hit)
                if hit not in best_hits:
                    best_hits[hit] = float(row[-1])
                else:
                    if float(row[-1]) > best_hits[hit]:
                        best_hits[hit] = float(row[-1])
        with open(os.path.join(running_d, 'blast_outputs/effectors_hits.csv'), 'w') as out_f:
            for hit in best_hits:
                out_f.write(f'{hit},{best_hits[hit]}\n')
    return blast_out_dic


if not os.path.exists(os.path.join(running_d, 'tmp_mmseq')):
    os.makedirs(os.path.join(running_d, 'tmp_mmseq'))
protein_mmseqs_all_vs_all(effectors, bacterial_proteome,
                          os.path.join(os.path.join(running_d, 'blast_outputs/effectorsDB.blast')),
                          os.path.join(os.path.join(running_d, 'tmp_mmseq')))

k = 0.6
effectors_homologs = set()
homologs = list(parse_blast_out(os.path.join(running_d, 'blast_outputs/effectorsDB.blast'), min_identity=k,
                                minimal_coverage=min_coverage).values())
for h in homologs:
    effectors_homologs = set.union(effectors_homologs, h)
    

recs = SeqIO.parse(bacterial_proteome, 'fasta')
effectors_recs = []
for rec in recs:
    if rec.id in effectors_homologs:
        effectors_recs.append(rec)
SeqIO.write(effectors_recs, effectors_prots, 'fasta')
