from Bio import SeqIO
import re
import csv
from sys import argv
import os
import fasta_parser
ORFs_file = argv[1]
working_directory = argv[2]
gff_dir = argv[3]
gff_files = [f'{gff_dir}/{file}' for file in os.listdir(gff_dir) if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{gff_dir}/{file}')]
genome_dir = argv[4]
os.chdir(working_directory)
locus_dic = fasta_parser.parse_ORFs(ORFs_file)
#%%

def parse_gff(gff_f):
    locus_area_d ={}
    with open(gff_f) as in_f:
        for line in in_f:
            if line.startswith('#'):
                continue
            line_l = line.strip().split('\t')
            if len(line_l) > 2:
                if line_l[2]=='gene' or line_l[2]=='pseudogene':
                    area = [(int(line_l[3])-1,int(line_l[4]))]
                    if line_l[6]=='+':
                        area.append(1)
                    elif line_l[6]=='-':
                        area.append(-1)
                    features_l = line_l[-1].split(';')
                    for feature in features_l:
                        if feature.startswith('locus_tag='):
                            locus_tag = feature.split('=')[1]
                            locus_area_d[locus_tag]=area
                            break
    return locus_area_d

def get_promoters(gene_area_dic,genome_f):
    promoters_d = {}
    genome = next(SeqIO.parse(genome_f,'fasta')).seq
    for gene in gene_area_dic:
        if gene_area_dic[gene][1]==1:
            start,end = gene_area_dic[gene][0][0],gene_area_dic[gene][0][1]
            if start >= 350:
                promoter = genome[start-350:start]
            else:
                to_add = 350-start
                promoter = genome[-to_add:]+genome[:start]
        else: #-1
            start,end = gene_area_dic[gene][0][0],gene_area_dic[gene][0][1]
            if len(genome)-end >=350:
                promoter = genome[end:end+350].reverse_complement()
            else:
                to_add = 350 - (len(genome)-end)
                promoter = genome[end:]+genome[:to_add]
                promoter = promoter.reverse_complement()
        promoters_d[gene] = promoter
    return promoters_d
                       
hrp_box = '[GT]GGA[AG]C[CT][ATGC]{15,16}CCAC[ATGC]{2}A'
hrp_one_mismatch = '[ATGC]GGA[AG]C[CT][ATGC]{15,16}CCAC[ATGC]{2}A|[GT][ATGC]GA[AG]C[CT][ATGC]{15,16}CCAC[ATGC]{2}A\
                    |[GT]G[ATGC]A[AG]C[CT][ATGC]{15,16}CCAC[ATGC]{2}A|[GT]GG[ATGC][AG]C[CT][ATGC]{15,16}CCAC[ATGC]{2}A\
                    |[GT]GGA[ATGC]C[CT][ATGC]{15,16}CCAC[ATGC]{2}A|[GT]GGA[AG][ATGC][CT][ATGC]{15,16}CCAC[ATGC]{2}A\
                    |[GT]GGA[AG]C[ATGC]{16,17}CCAC[ATGC]{2}A|[GT]GGA[AG]C[CT][ATGC]{16,17}CAC[ATGC]{2}A\
                    |[GT]GGA[AG]C[CT][ATGC]{15,16}C[ATGC]AC[ATGC]{2}A|[GT]GGA[AG]C[CT][ATGC]{15,16}CC[ATGC]C[ATGC]{2}A\
                    |[GT]GGA[AG]C[CT][ATGC]{15,16}CCA[ATGC]{3}A|[GT]GGA[AG]C[CT][ATGC]{15,16}CCAC[ATGC]{3}'

def create_box_one_mismatch(box_list):
    box_one_mismatch = []
    for i in range(len(box_list)):
        li = list(box_list)
        if li[i]!='[ATGC]':
            li[i] = '[ATGC]'
            motif = ''.join(li)
            box_one_mismatch.append(motif)
    box_one_mismatch = '|'.join(box_one_mismatch)
    return box_one_mismatch

exs_box = 'A{5}[ATGC][AT][ATGC][AC][CT][ATGC]{3}[AC][CT]TG[CT]A{2}[GT]'                     
exs_box_l = ['A']*5+['[ATGC]']+['[AT]','[ATGC]','[AC]','[CT]']+['[ATGC]']*3+['[AC]','[CT]','T','G','[CT]']+['A']*2+['[GT]']
exs_box_one_mismatch = create_box_one_mismatch(exs_box_l)

pip_box = 'TTCG[TCG][ATCG]{15}TTCG[TCG]'
pip_box_l = ['T']*2+['C','G','[TCG]']+['[ATCG]']*15+['T']*2+['C','G','[TCG]']
pip_box_one_mismatch = create_box_one_mismatch(pip_box_l)
 
mxiE_box = 'GTATCGT{7}A[ATGC]AG'
mxiE_box_l = ['G','T','A','T','C','G']+['T']*7+['A','[ATGC]','A','G']
mxiE_box_one_mismatch = create_box_one_mismatch(mxiE_box_l)                 

promoters_dicts = []
for gff_f in gff_files:
    name = gff_f.split('/')[-1].split('.')[0]
    genome_f = f'{genome_dir}/{name}.fasta'
    locus_area_d = parse_gff(gff_f)
    promoters_d = get_promoters(locus_area_d,genome_f)
    promoters_dicts.append(promoters_d)
                 
def existence_upstream_to_AUG(locus,pattern):
    for promoters_d in promoters_dicts:
        if locus in promoters_d:
            promoter = promoters_d[locus]
            break
    if re.search(pattern,str(promoter),re.I):
        return 1
    else:
        return 0
    

with open(f'pip_box_features.csv','w',newline='') as f:
    csv_writer = csv.writer(f)
    header=['locus','PIP_box','PIP_box_mismatch','hrp_box','hrp_box_mismatch','mxiE_box','mxiE_box_mismatch','exs_box','exs_box_mismatch']
    csv_writer.writerow(header)
    for locus in locus_dic:
        l=[locus]
        l.append(existence_upstream_to_AUG(locus,pip_box))
        l.append(existence_upstream_to_AUG(locus,pip_box_one_mismatch))
        l.append(existence_upstream_to_AUG(locus,hrp_box))
        l.append(existence_upstream_to_AUG(locus,hrp_one_mismatch))
        l.append(existence_upstream_to_AUG(locus,mxiE_box))
        l.append(existence_upstream_to_AUG(locus,mxiE_box_one_mismatch))
        l.append(existence_upstream_to_AUG(locus,exs_box))
        l.append(existence_upstream_to_AUG(locus,exs_box_one_mismatch))
        csv_writer.writerow(l)
endfile = open('pip_box_features.done','w')
endfile.close()  