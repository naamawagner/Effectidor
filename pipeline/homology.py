import fasta_parser
import csv
from sys import argv
import os
ORFs_f = argv[1]
blast_query = argv[2]
blast_dataset_dir = argv[3]
effectors_prots = argv[4]
working_directory = argv[5]
blast_out = 'blast_outputs'
os.chdir(working_directory)
locus_dic=fasta_parser.parse_ORFs(ORFs_f)
#%%

def protein_blast_all_vs_all(query,dataset,out,e_val='0.01'):
    import subprocess
    make_data_cmd=f"makeblastdb -in {dataset} -dbtype prot -out {dataset[:-4]}db"
    make_blast_cmd=f'blastp -db {dataset[:-4]}db -query {query} -outfmt 6 -out {out} -evalue {e_val}'
    subprocess.check_output(make_data_cmd,shell=True)
    subprocess.check_output(make_blast_cmd,shell=True)

def parse_blast_out(blast_out,e_val=0.01):
    blast_out_dic={}
    homologs_dic={}
    with open(blast_out,'r') as in_f:
        for row in in_f:
            row=row.split('\t')
            prot_id=row[0]
            subset_id=row[1]
            e_value=float(row[-2])
            if subset_id!=prot_id:
                if e_value<=e_val:
                    if prot_id not in blast_out_dic:
                        blast_out_dic[prot_id]=[1,float(row[-1])]
                        homologs_dic[prot_id]=[subset_id]
                    else:
                        if subset_id not in homologs_dic[prot_id]: 
                            blast_out_dic[prot_id][0]+=1
                            homologs_dic[prot_id].append(subset_id)    
    return blast_out_dic

def blast_features(locus,blast_out_dic):
    #blast_out_dic=parse_blast_out(blast_out)
    if locus in blast_out_dic.keys():
        bit_score,hits = blast_out_dic[locus][1],blast_out_dic[locus][0]
    else:
        bit_score,hits = 0,0
    return bit_score,hits

blast_out_dics=[]
datasets = {dataset[:-4]:f'{blast_dataset_dir}/{dataset}' for dataset in os.listdir(blast_dataset_dir) if dataset.endswith('faa')}
datasets['effectors'] = effectors_prots
datasets_l=sorted(datasets.keys())
if not os.path.exists(blast_out):
    os.makedirs(blast_out)
for dataset in datasets_l:
    protein_blast_all_vs_all(blast_query,datasets[dataset],f'{blast_out}/{dataset}.blast')
    blast_out_dics.append(parse_blast_out(f'{blast_out}/{dataset}.blast'))
with open(f'homology_features.csv','w',newline='') as out_f:
    writer=csv.writer(out_f)
    header=['locus']
    for dataset in datasets_l:
        header.append(f'homology_to_{dataset}_(bit_score)')
        header.append(f'hits_to_{dataset}')
    writer.writerow(header)
    for locus in locus_dic:
        l = [locus]
        for dataset in blast_out_dics:
            bit_score,hits = blast_features(locus,dataset)
            l += [bit_score,hits]
        writer.writerow(l)


endfile = open(f'homology_features.done','w')
endfile.close()