import fasta_parser
import csv
from sys import argv
import os
from protein_mmseq import protein_mmseqs_all_vs_all
ORFs_f = argv[1]
blast_query = argv[2]
blast_dataset_dir = argv[3]
working_directory = argv[4]
blast_out = os.path.join(working_directory, 'blast_outputs')
locus_dic = fasta_parser.parse_ORFs(ORFs_f)
#%%


def parse_blast_out(blast_out_f, e_val=0.01):
    blast_out_dic = {}
    # homologs_dic={}
    with open(blast_out_f, 'r') as in_f:
        for row in in_f:
            row = row.split('\t')
            prot_id = row[0].replace('^', '|')  # pipes were checnged to ^ before. We need to go back to pipes to
            # match with the original file.
            subset_id = row[1].replace('^', '|')  # pipes were checnged to ^ before. We need to go back to pipes to
            # match with the original file.
            e_value = float(row[-2])
            if subset_id != prot_id:
                if e_value <= e_val:
                    if prot_id not in blast_out_dic:
                        blast_out_dic[prot_id] = [subset_id, float(row[-1])]
                        #homologs_dic[prot_id]=[subset_id]
    
    return blast_out_dic


def blast_features(locus, blast_out_dic):
    #blast_out_dic=parse_blast_out(blast_out)
    if locus in blast_out_dic.keys():
        bit_score = blast_out_dic[locus][1]
    else:
        bit_score = 0
    return bit_score


for f in os.listdir(blast_dataset_dir):
    if f.endswith('.faa'):
        in_f = open(os.path.join(blast_dataset_dir, f))
        content = in_f.read()
        replaced_pipes = content.replace('|', '^')  # MMSeq changes the IDs when pipes are available in the rec.id,
        # so I replace it now with a character that is not expected to be found in the ID, to replace it back while
        # parsing the MMSeq result.
        in_f.close()
        with open(os.path.join(blast_dataset_dir, f), 'w') as out_f:
            out_f.write(replaced_pipes)
in_f = open(blast_query)
content = in_f.read()
replaced_pipes = content.replace('|', '^')  # MMSeq changes the IDs when pipes are available in the rec.id,
# so I replace it now with a character that is not expected to be found in the ID, to replace it back while parsing
# the MMSeq result.
in_f.close()
with open(blast_query, 'w') as out_f:
    out_f.write(replaced_pipes)

blast_out_dics = []
datasets = {dataset[:-4]: f'{blast_dataset_dir}/{dataset}' for dataset in os.listdir(blast_dataset_dir) if dataset.endswith('faa')}
# datasets['effectors'] = effectors_prots
datasets_l = sorted(datasets.keys())
if not os.path.exists(blast_out):
    os.makedirs(blast_out)
tmp_mmseq = os.path.join(blast_dataset_dir, 'tmp_mmseq')
if not os.path.exists(tmp_mmseq):
    os.makedirs(tmp_mmseq)
for dataset in datasets_l:
    protein_mmseqs_all_vs_all(blast_query, datasets[dataset], os.path.join(blast_out, f'{dataset}.blast'), tmp_mmseq)
    blast_out_dics.append(parse_blast_out(f'{blast_out}/{dataset}.blast'))
with open(os.path.join(working_directory, 'homology_features.csv'), 'w', newline='') as out_f:
    writer = csv.writer(out_f)
    header = ['locus']
    for dataset in datasets_l:
        header.append(f'homology_to_{dataset}_(bit_score)')
    writer.writerow(header)
    for locus in locus_dic:
        l = [locus]
        for dataset in blast_out_dics:
            bit_score = blast_features(locus, dataset)
            l.append(bit_score)
        writer.writerow(l)

# for further reports
with open(os.path.join(working_directory, 'closest_effector_homologs.csv'), 'w', newline='') as out:
    writer = csv.writer(out)
    header = ['locus', 'Effector_ID']
    writer.writerow(header)
    T3Es_blast_dic = parse_blast_out(f'{blast_out}/T3Es.blast', 10**(-4))
    for locus in T3Es_blast_dic:
        row = [locus, T3Es_blast_dic[locus][0]]
        writer.writerow(row)


endfile = open(os.path.join(working_directory, 'homology_features.done'), 'w')
endfile.close()
