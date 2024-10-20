from Bio import SeqIO
from sys import argv

ORFs_file = argv[1]
effectors_file = argv[2]
all_prots = argv[3]
effectors_prots = argv[4]

def translate_fasta(ORFs_fasta,prots_fasta):
    prots = []
    recs = SeqIO.parse(ORFs_fasta,'fasta')
    for rec in recs:
        if len(rec.seq.translate(to_stop=True)) > 10:
            header = rec.description
            if 'locus_tag' in header:
                for field in header.split():
                    if 'locus_tag' in field:
                        locus = field.split('=')[1].strip(']')
                        break
                rec.id = locus
                rec.description = locus
            else:
                rec.description = rec.id
            rec.seq = rec.seq.translate(to_stop=True)
            prots.append(rec)
    SeqIO.write(prots, prots_fasta, 'fasta')

translate_fasta(ORFs_file, all_prots)
if effectors_file:
    translate_fasta(effectors_file, effectors_prots)