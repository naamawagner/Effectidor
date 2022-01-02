from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import math
import re
import csv
import os
import fasta_parser

#%%

def parse_gff(gff_f):
    CDS_l =[]
    RNA = []
    with open(gff_f) as in_f:
        for line in in_f:
            if line.startswith('#'): #header
                continue
            
            line_l = line.strip().split('\t')
            if len(line_l) > 2:
                if line_l[2]=='CDS':
                    features_l = line_l[-1].split(';')
                    for feature in features_l:
                        if feature.startswith('locus_tag='):
                            locus_tag = feature.split('=')[1]
                            CDS_l.append(locus_tag)
                            break
                elif 'RNA' in line_l[2]:
                    features_l = line_l[-1].split(';')
                    for feature in features_l:
                        if feature.startswith('locus_tag='):
                            locus_tag = feature.split('=')[1]
                            RNA.append(locus_tag)
                            break
    return set(CDS_l),set(RNA)

def parse_gff_to_CDS_loc(gff_f):
    regions = []
    circulars = []
    linears = []
    locus_area_d ={}
    with open(gff_f) as in_f:
        for line in in_f:
            if line.startswith('#'): #header
                line_l = line.split()
                if 'sequence-region' in line_l[0]:
                    regions.append(line_l[1])
                    locus_area_d[line_l[1]] = {}
            
            line_l = line.strip().split('\t')
            if len(line_l) > 2:
                if line_l[2] == 'region':
                    if 'Is_circular=true' in line_l[-1]:
                        circulars.append(line_l[0])
                    else:
                        linears.append(line_l[0])
                if line_l[2]=='CDS':
                    region = line_l[0]
                    area = [(int(line_l[3])-1,int(line_l[4]))]
                    if line_l[6]=='+':
                        area.append(1)
                    elif line_l[6]=='-':
                        area.append(-1)
                    features_l = line_l[-1].split(';')
                    for feature in features_l:
                        if feature.startswith('locus_tag='):
                            locus_tag = feature.split('=')[1]
                            locus_area_d[region][locus_tag]=area
                            break
    return locus_area_d,circulars

def get_promoters(gene_area_dic,circular_contigs,genome_f):
    promoters_d = {}
    genome = SeqIO.to_dict(SeqIO.parse(genome_f,'fasta'))
    for region_dic in gene_area_dic:
        for gene in gene_area_dic[region_dic]:
            if gene_area_dic[region_dic][gene][1]==1: #forward
                start,end = gene_area_dic[region_dic][gene][0][0],gene_area_dic[region_dic][gene][0][1]
                if start >= 350:
                    promoter = genome[region_dic].seq[start-350:start]
                else:
                    if region_dic in circular_contigs:
                        to_add = 350-start
                        promoter = genome[region_dic].seq[-to_add:]+genome[region_dic].seq[:start]
                    else:
                        promoter = genome[region_dic].seq[:start]
            else: #-1 reverse
                start,end = gene_area_dic[region_dic][gene][0][0],gene_area_dic[region_dic][gene][0][1]
                if len(genome[region_dic].seq)-end >=350:
                    promoter = genome[region_dic].seq[end:end+350].reverse_complement()
                else:
                    if region_dic in circular_contigs:
                        to_add = 350 - (len(genome[region_dic].seq)-end)
                        promoter = genome[region_dic].seq[end:]+genome[region_dic].seq[:to_add]
                        promoter = promoter.reverse_complement()
                    else:
                        promoter = genome[region_dic].seq[end:].reverse_complement()
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
pip_box_l = ['T']*2+['C','G','[TCG]']+['[ATGC]']*15+['T']*2+['C','G','[TCG]']
pip_box_one_mismatch = create_box_one_mismatch(pip_box_l)
 
mxiE_box = 'GTATCGT{7}A[ATGC]AG'
mxiE_box_l = ['G','T','A','T','C','G']+['T']*7+['A','[ATGC]','A','G']
mxiE_box_one_mismatch = create_box_one_mismatch(mxiE_box_l)  

tts_box = 'GTCAG[TCG]T[TCAG]{4}G[AT][AC]AG[CGT][TAC][ATCG]{3}[CTG]{2}[ATCG]{4}A'
tts_box_l = ['G','T','C','A','G','[TGC]','T']+['[ATGC]']*4+['G','[AT]','[AC]','A','G','[CGT]','[TAC]']+['[ATGC]']*3+['[CTG]']*2+['[ATGC]']*4+['A']
tts_box_one_mismatch = create_box_one_mismatch(tts_box_l)              

def run_hmmer(hmm_profile_path,query,out):
    cmd = f'module load hmmer/hmmer-3.3.2;hmmsearch {hmm_profile_path} {query} > {out}'
    subprocess.check_output(cmd,shell=True)
    
def parse_hmmer(hmmer_out,locus):
    with open(hmmer_out) as hmmer_f:
        for line in hmmer_f:
            if line.startswith(' ') and locus in line:
                line_l = line.split()
                e_val_log , score = -math.log(float(line_l[3])) , float(line_l[4])
                return e_val_log,score
    return 0,0
                
    
def main(ORFs_file, working_directory, gff_dir, genome_dir, PIP=False, hrp=False, mxiE=False, exs=False, tts=False):
    gff_files = [f'{gff_dir}/{file}' for file in os.listdir(gff_dir) if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{gff_dir}/{file}')]
    os.chdir(working_directory)
    locus_dic = fasta_parser.parse_ORFs(ORFs_file)
    #promoters_dicts = []
    locus_area_d,circulars = {},[]
    for gff_f in gff_files:
        #name = gff_f.split('/')[-1].split('.')[0]
        #genome_f = f'{genome_dir}/{name}.fasta'
        locus_area_d1,circulars1 = parse_gff_to_CDS_loc(gff_f)
        locus_area_d.update(locus_area_d1)
        circulars.extend(circulars1)
    genome_recs = []
    for genome_f in os.listdir(genome_dir):
        if os.path.isfile(f'{genome_dir}/{genome_f}'):
            recs = list(SeqIO.parse(f'{genome_dir}/{genome_f}','fasta'))
            genome_recs.extend(recs)
    SeqIO.write(genome_recs,'genome.fasta','fasta')
    promoters_d = get_promoters(locus_area_d,circulars,'genome.fasta')
        #promoters_dicts.append(promoters_d)
     
    '''
    if tts:
        os.makedirs(f'{working_directory}/promoters')
        os.makedirs(f'{working_directory}/hmmer_out')
        for promoters_d in promoters_dicts:
            for locus in promoters_d:
                rec = SeqRecord(promoters_d[locus],id=locus)
                SeqIO.write(rec,f'{working_directory}/promoters/{locus}.fasta','fasta')   '''
        
    def existence_upstream_to_AUG(locus,pattern):
        promoter = promoters_d[locus]
        if re.search(pattern,str(promoter),re.I):
            return 1
        else:
            return 0
        
    with open(f'pip_box_features.csv','w',newline='') as f:
        csv_writer = csv.writer(f)
        header=['locus']
        if PIP:
            header += ['PIP_box','PIP_box_mismatch']
        if hrp:
            header += ['hrp_box','hrp_box_mismatch']
        if mxiE:
            header += ['mxiE_box','mxiE_box_mismatch']
        if exs:
            header += ['exs_box','exs_box_mismatch']
        if tts:
            header += ['tts_box','tts_box_mismatch']
            #header += ['tts_hmmer_e_val_log','tts_hmmer_score']
        csv_writer.writerow(header)
        for locus in locus_dic:
            l=[locus]
            if PIP:
                l.append(existence_upstream_to_AUG(locus,pip_box))
                l.append(existence_upstream_to_AUG(locus,pip_box_one_mismatch))
            if hrp:
                l.append(existence_upstream_to_AUG(locus,hrp_box))
                l.append(existence_upstream_to_AUG(locus,hrp_one_mismatch))
            if mxiE:
                l.append(existence_upstream_to_AUG(locus,mxiE_box))
                l.append(existence_upstream_to_AUG(locus,mxiE_box_one_mismatch))
            if exs:
                l.append(existence_upstream_to_AUG(locus,exs_box))
                l.append(existence_upstream_to_AUG(locus,exs_box_one_mismatch))
            if tts:
                l.append(existence_upstream_to_AUG(locus,tts_box))
                l.append(existence_upstream_to_AUG(locus,tts_box_one_mismatch))
                #hmm_profile_path = r'/groups/pupko/naamawagner/T3Es_webserver/tts.hmm'
                #query = f'{working_directory}/promoters/{locus}.fasta'
                #out_hmmer = f'{working_directory}/hmmer_out/{locus}.txt'
                #run_hmmer(hmm_profile_path,query,out_hmmer)
                #e_val , score = parse_hmmer(out_hmmer,locus)
                #l.append(e_val)
                #l.append(score)
                
            csv_writer.writerow(l)
    endfile = open('pip_box_features.done','w')
    endfile.close()

if __name__ == '__main__':

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_ORFs_path',
                            help='A path to a DNA ORFs file.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_dir_path',
                            help='A path to a folder in which the output files will be created.',
                            type=lambda path: path.rstrip('/'))
        parser.add_argument('gff_path',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'),
                            help='A path to GFF files archive.')
        parser.add_argument('genome_path',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'),
                            help='A path to fasta files archive with full genome records.')
        parser.add_argument('--PIP', help='look for PIP-box in promoters', action='store_true')
        parser.add_argument('--hrp', help='look for hrp-box in promoters', action='store_true')
        parser.add_argument('--mxiE', help='look for mxiE-box in promoters', action='store_true')
        parser.add_argument('--exs', help='look for exs-box in promoters', action='store_true')
        parser.add_argument('--tts',help='look for tts-box in promoters', action='store_true')
        
        args = parser.parse_args()
        ORFs_file = args.input_ORFs_path
        working_directory = args.output_dir_path
        gff_dir = args.gff_path
        genome_dir = args.genome_path
        PIP_flag = args.PIP
        hrp_flag = args.hrp
        mxiE_flag = args.mxiE
        exs_flag = args.exs
        tts_flag = args.tts
        main(ORFs_file, working_directory, gff_dir, genome_dir, PIP=PIP_flag, hrp=hrp_flag, mxiE=mxiE_flag, exs=exs_flag, tts=tts_flag)