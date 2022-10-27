import sys
sys.path.append('/bioseq/effectidor/auxiliaries')
import os
import logging
#import shutil
import Bio.SeqUtils
import effectidor_CONSTANTS as CONSTS  # from /effectidor/auxiliaries
from time import sleep
from auxiliaries import fail,update_html,append_to_html # from /effectidor/auxiliaries
from Bio import SeqIO
from T3Es_wrapper import effectors_learn
import shutil
import re
import subprocess

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def verify_fasta_format(fasta_path,Type,input_name):
    logger.info(f'Validating FASTA format:{fasta_path}')
    file_size = os.path.getsize(fasta_path)
    if int(float(file_size)) == 0:
        return f'{input_name} is empty!'
        
    flag = False
    if Type == 'DNA':
        legal_chars = set(Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() + Bio.SeqUtils.IUPACData.ambiguous_dna_letters)
    else: #Type == 'protein'
        legal_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ*')
        flag = True
    with open(fasta_path) as f:
        line_number = 0
        try:
            line = f.readline()
            line_number += 1
            if not line.startswith('>'):
                return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. First line in fasta {input_name} starts with "{line[0]}" instead of ">".'
            previous_line_was_header = True
            putative_end_of_file = False
            curated_content = f'>{line[1:]}'.replace("|", "_")
            for line in f:
                line_number += 1
                line = line.strip()
                if not line:
                    if not putative_end_of_file: # ignore trailing empty lines
                        putative_end_of_file = line_number
                    continue
                #if putative_end_of_file:  # non empty line after empty line
                #    return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {putative_end_of_file} in fasta {input_name} is empty.'
                if line.startswith('>'):
                    if previous_line_was_header:
                        return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. fasta {input_name} contains an empty record. Both lines {line_number-1} and {line_number} start with ">".'
                    else:
                        previous_line_was_header = True
                        curated_content += f'>{line[1:]}\n'.replace("|", "_")
                        continue
                else:  # not a header
                    previous_line_was_header = False
                    for c in line:
                        if c not in legal_chars:
                            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in fasta {input_name} contains illegal {Type} character "{c}".'
                    curated_content += f'{line}\n'
            if flag: # protein
                recs = SeqIO.parse(fasta_path,'fasta')
                for rec in recs:
                    seq = str(rec.seq)
                    AGCT_count = seq.count('A')+seq.count('G')+seq.count('T')+seq.count('C')+seq.count('N')
                    if AGCT_count >= 0.95*len(seq):
                        f = AGCT_count*100/len(seq)
                        return f'Protein fasta {input_name} seems to contain DNA records (record {rec.id} contains {float("%.2f" % (f))}% DNA characters). Make sure all samples in the file are protein sequences and re-submit.'
        except UnicodeDecodeError as e:
            logger.info(e.args)
            line_number += 1  # the line that was failed to be read
            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in fasta contains one (or more) non <a href="https://en.wikipedia.org/wiki/ASCII" target="_blank">ascii</a> character(s).'
    # override the old file with the curated content
    with open(fasta_path, 'w') as f:
        f.write(curated_content)

def verify_ORFs(ORFs_path):
    logger.info(f'Validating ORFs:{ORFs_path}')
    ORFs_recs = set([rec.id for rec in SeqIO.parse(ORFs_path,'fasta')])
    if len(ORFs_recs) < 1000:
        return f'The ORFs file contains only {str(len(ORFs_recs))} records. Make sure this file contains all the ORFs (open reading frames) in the genome - Effectidor is designed to analyze full genomes and not a sample of genes. Also, make sure this file contains ORFs and not full genome sequence! The full genome sequence can be uploaded in the advanced options.'
    elif len(ORFs_recs) > 10000:
        return f'The ORFs file contains {str(len(ORFs_recs))} records! Effectidor is designed to analyze only one bacterial genome at a time. Make sure your ORFs file contains all the ORFs (open reading frames) in the genome, and only the ORFs of one genome. This number cannot exceed 10,000 ORFs.'

def verify_effectors_f(effectors_path, ORFs_path):
    logger.info(f'Validating effectors:{effectors_path}')
    effectors_recs = [rec.id for rec in SeqIO.parse(effectors_path,'fasta')]
    effectors_set = set([rec.id for rec in SeqIO.parse(effectors_path,'fasta')])
    ORFs_recs = set([rec.id for rec in SeqIO.parse(ORFs_path,'fasta')])
    if not effectors_set.issubset(ORFs_recs):
        not_in_ORFs = ', '.join([rec for rec in effectors_set.difference(ORFs_recs)])
        return f'Illegal effectors records. The following records IDs are in the effectors file and not in the ORFs file:\n{not_in_ORFs}.'
    if len(effectors_set) != len(effectors_recs):
        more_than_once = ','.join([effector for effector in effectors_set if effectors_recs.count(effector)>1])
        return f'Illegal effectors records. The following records IDs appear more than once in the file:\n{more_than_once}.'
    
def verify_genome_one_contig(genome_path, file_name):
    logger.info(f'Validating one contig:{genome_path}')
    recs = list(SeqIO.parse(genome_path,'fasta'))
    if len(recs) > 1: # it must be one contig
        return f'Illegal genome file of {file_name}. The FASTA file contains more than one record (more than one contig).'
    
def verify_zip(file,name):
    logger.info(f'Validating zip:{file}')
    if not file.endswith('.zip'):
        return f'{name} not in zip format. Make sure to upload a zip archive and resubmit your job.'
    #unzip it to a tmp_dir
    if not os.path.exists(f'{"/".join(file.split("/")[:-1])}/zip_tmp'):
        os.makedirs(f'{"/".join(file.split("/")[:-1])}/zip_tmp')
    shutil.unpack_archive(file,'/'.join(file.split('/')[:-1])+'/zip_tmp')
    flag = False
    for f in os.listdir(f'{"/".join(file.split("/")[:-1])}/zip_tmp'):
        new_name = f.replace(' ','_')
        os.rename(f'{"/".join(file.split("/")[:-1])}/zip_tmp/{f}',f'{"/".join(file.split("/")[:-1])}/zip_tmp/{new_name}')
        f = new_name
        if not f.startswith('_') and not f.startswith('.') and os.path.isfile(f'{"/".join(file.split("/")[:-1])}/zip_tmp/{f}'):
            error_msg = verify_fasta_format(f'{"/".join(file.split("/")[:-1])}/zip_tmp/{f}','protein',f)
            if error_msg:
                return f'In {name}:\n{error_msg}'
            flag = True
            subprocess.check_output(f'rm {"/".join(file.split("/")[:-1])}/zip_tmp/{f}',shell=True)
    if flag == False:
        return f'{name} contains no valid files. Make sure to include protein fasta files in this archive!' 
            
    
def validate_set(file,name):
    logger.info(f'Validating set:{file}')
    # by IDs
    recs = [rec.id for rec in SeqIO.parse(file,'fasta')]
    if len(recs) > len(set(recs)):
        non_unique = []
        for rec_id in set(recs):
            if recs.count(rec_id) > 1:
                non_unique.append(rec_id)
        non_unique_ids = '<br>'.join(non_unique)
        return f'{name} contain non unique records. The following records appear more than once:<br><br>{non_unique_ids}'
    # by locus_tags
    recs = SeqIO.parse(file,'fasta')
    locus_tags = []
    for rec in recs:
        header = rec.description
        header_l = header.split()
        for field in header_l:
            if "locus_tag" in field:
                locus = field.split('=')[1].strip(']')
                locus_tags.append(locus)
    if len(locus_tags) > len(set(locus_tags)):
        non_unique = []
        for rec_id in set(locus_tags):
            if locus_tags.count(rec_id) > 1:
                non_unique.append(rec_id)
        non_unique_ids = '<br>'.join(non_unique)
        return f'{name} contain non unique records. The following records appear more than once:<br><br>{non_unique_ids}'
    
def validate_gff(gff_dir,ORFs_f):
    logger.info(f'Validating GFF coverage')
    import fasta_parser
    import pip_box_features
    locus_dic = fasta_parser.parse_ORFs(ORFs_f)
    CDS_set = set()
    RNA_set = set()
    for f in os.listdir(gff_dir):
        with open(f'{gff_dir}/{f}') as in_f:
            content = in_f.read().replace('|','_')
            with open(f'{gff_dir}/{f}','w') as out_f:
                out_f.write(content)
        CDS,RNA = pip_box_features.parse_gff(f'{gff_dir}/{f}',locus_dic)
        CDS_set.update(CDS)
        RNA_set.update(RNA)
    not_in_gff = [locus for locus in locus_dic if (locus not in CDS_set and locus not in RNA_set)]
    if len(not_in_gff)>0:
        return f'''There are records in your ORFs input that are not available in the GFF file. Please revise your input and submit again.<br>
                These records are:<br>{", ".join(not_in_gff)}'''
    if len(RNA_set) > 0:
        recs = SeqIO.parse(ORFs_f,'fasta')
        cds_recs = []
        for rec in recs:
            if rec.id in RNA_set:
                continue
            elif re.search(r'locus_tag=(\w+)',rec.description):
                if re.search(r'locus_tag=(\w+)',rec.description).group(1) in RNA_set:
                    continue
            cds_recs.append(rec) 
        SeqIO.write(cds_recs,ORFs_f,'fasta')
        
def validate_genome_and_gff(gff_dir,genome_dir,ORFs_f):
    logger.info(f'Validating genome and gff')
    import fasta_parser
    import pip_box_features
    locus_dic = fasta_parser.parse_ORFs(ORFs_f)
    gff_files = [f'{gff_dir}/{file}' for file in os.listdir(gff_dir) if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{gff_dir}/{file}')]    
    regions = []
    for gff_f in gff_files:
        locus_area_d,circulars = pip_box_features.parse_gff_to_CDS_loc(gff_f,locus_dic)
        for region in locus_area_d:
            regions.append(region)
        #with open(gff_f) as in_f:
        #    for line in in_f:
        #        if line.startswith('#'): #header
        #            line_l = line.split()
        #            if 'sequence-region' in line_l[0]:
        #                regions.append(line_l[1])
    genome_recs = []
    for genome_f in os.listdir(genome_dir):
        if os.path.isfile(f'{genome_dir}/{genome_f}') and not genome_f.startswith('_') and not genome_f.startswith('.'):
            recs_ids = [rec.id for rec in SeqIO.parse(f'{genome_dir}/{genome_f}','fasta')]
            genome_recs.extend(recs_ids)
    if len(genome_recs) > len(set(genome_recs)):
        non_unique = []
        for ID in genome_recs:
            if genome_recs.count(ID) > 1:
                non_unique.append(ID)
        non_unique = set(non_unique)
        return f'Several contigs appear more than once within the genomic sequences data. These contigs IDs are:<br>{", ".join(non_unique)}<br>Please make sure to upload the genomic sequences such that each sequence will apear only once. These sequences can be uploaded in different FASTA files or in one FASTA file.'
    #return f'IDs: {str(genome_recs)}'
    not_in_genome = []
    for region in regions:
        if region not in genome_recs:
            not_in_genome.append(region)
    if len(not_in_genome)>0:
        return f'The following regions from the GFF file are not found in the full genome file:<br>{", ".join(not_in_genome)}.<br>Make sure the regions names are matching between the GFF and genome files and re-submit the job, or contact us for more information.'
    

def validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path, gff_path, no_T3SS_path, error_path):
    logger.info('Validating input...')
    os.makedirs(f'{output_dir_path}/contigs_ORFs')
    if ORFs_path.endswith('.zip'):     
        shutil.unpack_archive(f'{output_dir_path}/ORFs.zip',f'{output_dir_path}/contigs_ORFs')
        ORFs=[]
        for file in os.listdir(f'{output_dir_path}/contigs_ORFs'):
            if os.path.isfile(f'{output_dir_path}/contigs_ORFs/{file}'):
                error_msg = verify_fasta_format(f'{output_dir_path}/contigs_ORFs/{file}','DNA',f'{file} in ORFs archive')
                if error_msg:
                    error_msg = f'Illegal fasta files in {file} in ORFs archive: {error_msg}<br>This archive is expected to contain only <b>DNA</b> FASTA files.'
                    fail(error_msg,error_path)
                recs = SeqIO.parse(f'{output_dir_path}/contigs_ORFs/{file}','fasta')
                for rec in recs:
                    ORFs.append(rec)
        SeqIO.write(ORFs,f'{output_dir_path}/ORFs.fasta','fasta')
        error_msg = validate_set(f'{output_dir_path}/ORFs.fasta','ORFs records')
        if error_msg:
            fail(error_msg,error_path)
    else:
        shutil.copy(ORFs_path,f'{output_dir_path}/contigs_ORFs')
        error_msg = verify_fasta_format(ORFs_path,'DNA','input ORFs')
        if error_msg:
            error_msg = f'Illegal fasta file in ORFs input: {error_msg}<br>This input is expected to hold a <b>DNA</b> FASTA file.'
            fail(error_msg,error_path)
        error_msg = validate_set(ORFs_path,'ORFs records')
        if error_msg:
            fail(error_msg,error_path)
    ORFs_f = f'{output_dir_path}/ORFs.fasta'
    error_msg = verify_ORFs(ORFs_f)
    if error_msg:
        fail(error_msg,error_path)
    if effectors_path:
        error_msg = verify_fasta_format(effectors_path,'DNA', 'input effectors')
        if error_msg:
            error_msg = f'Illegal effectors file: {error_msg}'
            fail(error_msg, error_path)
        error_msg = verify_effectors_f(effectors_path,f'{output_dir_path}/ORFs.fasta')
        if error_msg:
            fail(error_msg, error_path)
        error_msg = validate_set(effectors_path,'effectors records')
        if error_msg:
            fail(error_msg,error_path)
    if input_T3Es_path:
        error_msg = verify_fasta_format(input_T3Es_path,'protein', 'effectors for homology search')
        if error_msg:
            fail(error_msg,error_path)
        error_msg = validate_set(input_T3Es_path,'effectors records for homology search')
        if error_msg:
            fail(error_msg,error_path)
    if host_proteome:
        error_msg = verify_zip(host_proteome,'Host data')
        if error_msg:
            fail(error_msg,error_path)
    if genome_path:
        # genome
        os.makedirs(f'{output_dir_path}/full_genome')
        if genome_path.endswith('.zip'):
            shutil.unpack_archive(f'{output_dir_path}/genome_sequence.zip',f'{output_dir_path}/full_genome')
        else:
            shutil.move(f'{output_dir_path}/genome.fasta',f'{output_dir_path}/full_genome')
        #genome_files_names = []
        for file in os.listdir(f'{output_dir_path}/full_genome'):
            if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{output_dir_path}/full_genome/{file}'): # discard system files and directories
                error_msg = verify_fasta_format(f'{output_dir_path}/full_genome/{file}','DNA', f'{file} in full genome archive')
                if error_msg:
                    error_msg = f'Illegal fasta files in {file} in full genome archive: {error_msg}'
                    fail(error_msg,error_path)
                '''
                error_msg = verify_genome_one_contig(f'{output_dir_path}/full_genome/{file}',file)
                if error_msg:
                    fail(error_msg,error_path)'''
                #genome_files_names.append('.'.join(file.split('.')[:-1]))
    if gff_path:
        # gff
        os.makedirs(f'{output_dir_path}/gff')
        if gff_path.endswith('.zip'):
            shutil.unpack_archive(f'{output_dir_path}/genome_features.zip',f'{output_dir_path}/gff')
        else:
            shutil.move(f'{output_dir_path}/genome.gff3',f'{output_dir_path}/gff')
        #gff_files_names = ['.'.join(file.split('.')[:-1]) for file in os.listdir(f'{output_dir_path}/gff') if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{output_dir_path}/gff/{file}')]
        #if set(gff_files_names) != set(genome_files_names):
            #in_gff_not_genome = set(gff_files_names).difference(set(genome_files_names))
            #in_genome_not_in_gff = set(genome_files_names).difference(set(gff_files_names))
            #error_msg = 'The names of the files in the gff archive are not matching the names of the files in the full genome archive!'
            #if len(in_gff_not_genome) > 0:
            #    error_msg += f' The following file names are present in the GFF archive but not in the full genome archive: {in_gff_not_genome}'
            #if len(in_genome_not_in_gff) > 0:
            #    error_msg += f' The following file names are present in the full genome archive but not in the GFF archive: {in_genome_not_in_gff}'
            #fail(error_msg,error_path)
        error_msg = validate_gff(f'{output_dir_path}/gff',f'{output_dir_path}/ORFs.fasta')
        if error_msg:
            fail(error_msg,error_path)
    if genome_path and gff_path:
        error_msg = validate_genome_and_gff(f'{output_dir_path}/gff',f'{output_dir_path}/full_genome',f'{output_dir_path}/ORFs.fasta')
        if error_msg:
            fail(error_msg,error_path)
    if no_T3SS_path:
        error_msg = verify_zip(no_T3SS_path,'Close bacteria without T3SS data')
        if error_msg:
            fail(error_msg,error_path)



def main(ORFs_path, output_dir_path, effectors_path, input_T3Es_path, host_proteome, html_path, queue, genome_path, gff_path, no_T3SS, full_genome=False, PIP=False, hrp=False, mxiE=False, exs=False, tts=False, homology_search=False, signal=False):

    error_path = f'{output_dir_path}/error.txt'
    try:
        if html_path:
            run_number = initialize_html(CONSTS, output_dir_path, html_path)
            #final_zip_path = f'{os.path.split(output_dir_path)[0]}/{CONSTS.WEBSERVER_NAME}_{run_number}'
    
        os.makedirs(output_dir_path, exist_ok=True)
    
        tmp_dir = f'{os.path.split(ORFs_path)[0]}/tmp'  # same folder as the input files
        os.makedirs(tmp_dir, exist_ok=True)
        
        validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path, gff_path, no_T3SS, error_path)
        if genome_path and gff_path:
            CIS_elements = True
        else:
            CIS_elements = False
        if gff_path:
            organization = True
        else:
            organization = False
        predicted_table, positives_table, T3SS_table, low_confidence_flag = effectors_learn(error_path, f'{output_dir_path}/ORFs.fasta', effectors_path, output_dir_path, tmp_dir, queue, organization=organization, CIS_elements=CIS_elements, PIP=PIP, hrp=hrp, mxiE=mxiE, exs=exs, tts=tts, homology_search=homology_search, signal=signal)
        
        #if full_genome:
        #    if genome_path and gff_path:
        #        predicted_table, positives_table, T3SS_table, low_confidence_flag = effectors_learn(error_path, f'{output_dir_path}/ORFs.fasta', effectors_path, output_dir_path, tmp_dir, queue, organization=True, CIS_elements=True, PIP=PIP, hrp=hrp, mxiE=mxiE, exs=exs, tts=tts, homology_search=homology_search)
        #    else:
        #        predicted_table, positives_table, T3SS_table, low_confidence_flag = effectors_learn(error_path, f'{output_dir_path}/ORFs.fasta', effectors_path, output_dir_path, tmp_dir, queue, organization=True, homology_search=homology_search)
        #else:
        #    predicted_table, positives_table, T3SS_table, low_confidence_flag = effectors_learn(error_path, f'{output_dir_path}/ORFs.fasta', effectors_path, output_dir_path, tmp_dir, queue, homology_search=homology_search)
    
        if html_path:
            #shutil.make_archive(final_zip_path, 'zip', output_dir_path)
            finalize_html(html_path, error_path, run_number, predicted_table, positives_table, T3SS_table, low_confidence_flag)

    except Exception as e:
        logger.info(f'SUCCEEDED = False')
        logger.info(e)
        logger.info(f"""ORFs_path: {ORFs_path}\noutput_dir_path: {output_dir_path}\n
                    effectors_path:{effectors_path}\nhost_proteome:{host_proteome}""")
        if html_path:
            error_msg = e.args[-1]
            if os.path.exists(error_path):
                with open(error_path) as f:
                    error_msg = f.read()
            edit_failure_html(CONSTS, error_msg, html_path, run_number)
            add_closing_html_tags(html_path, CONSTS, run_number)


def finalize_html(html_path, error_path, run_number, predicted_table, positives_table, T3SS_table, low_confidence_flag):
    succeeded = not os.path.exists(error_path)
    logger.info(f'SUCCEEDED = {succeeded}')
    if succeeded:
        edit_success_html(CONSTS, html_path, run_number, predicted_table, positives_table, T3SS_table, low_confidence_flag)
    else:
        edit_failure_html(CONSTS, error_path, html_path, run_number)
    add_closing_html_tags(html_path, CONSTS, run_number)


def edit_success_html(CONSTS, html_path, run_number, predicted_table, positives_table, T3SS_table, low_confidence_flag):
    update_html(html_path, 'RUNNING', 'FINISHED')
    if low_confidence_flag:
        append_to_html(html_path, f'''
                       <div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify"><h3>
                       <font color="red">WARNING: Due to small positive set to train the classifier on, the predictions are of low quality.
                       </font></h3><br>
                       </div>
                       ''')
    if positives_table:
        append_to_html(html_path, f'''
                <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='left'>
                <a href='out_learning/concensus_predictions_with_annotation.xlsx' target='_blank'>Download predictions file</a>
                <br>
                <a href='out_learning/pseudogenes.xlsx' target='_blank'>Download pseudogenes file</a>
                <br>
                <a href='out_learning/feature_importance.csv' target='_blank'>Download feature importance file</a>
                <br>
                <a href='features.csv' target='_blank'>Download raw features file</a>
                <br><br>
                <h3><b>Positive samples that used to train the model</b></h3>
                {positives_table}
                <br>
                <h3><b>Top 10 predictions, among unlabeled samples</b></h3>
                {predicted_table}
                <br>
                <h3><b>Type 3 secretion system proteins that were found in the genome:</b></h3>
                {T3SS_table}
                </div>
                <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='center'>
                <h3><b>feature importance
                <br>
                <a href='out_learning/feature_importance.png'><img src='out_learning/feature_importance.png'></a>
                <br><br>
                best features comparison - effectors vs non-effectors:
                <br>
                <a href='out_learning/0.png'><img src='out_learning/0.png'></a>
                <a href='out_learning/1.png'><img src='out_learning/1.png'></a>
                <a href='out_learning/2.png'><img src='out_learning/2.png'></a>
                <a href='out_learning/3.png'><img src='out_learning/3.png'></a>
                <a href='out_learning/4.png'><img src='out_learning/4.png'></a>
                <a href='out_learning/5.png'><img src='out_learning/5.png'></a>
                <a href='out_learning/6.png'><img src='out_learning/6.png'></a>
                <a href='out_learning/7.png'><img src='out_learning/7.png'></a>
                <a href='out_learning/8.png'><img src='out_learning/8.png'></a>
                <a href='out_learning/9.png'><img src='out_learning/9.png'></a>
                </b></h3>
                </div>
                ''')
    else:
        append_to_html(html_path, f'''
                <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='left'>
                Unfortunatelly, we could not train a satisfying classifier due to small positive set.<br>The effectors found based on homology are listed in the table bellow:<br>
                {predicted_table}
                </div>
                ''')


def edit_failure_html(CONSTS, error_msg, html_path, run_number):
    update_html(html_path, 'RUNNING', 'FAILED')
    append_to_html(html_path,
                   f'<div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify"><h3>\n'
                   f'<font color="red">{error_msg}</font></h3><br><br>'
                   f'Please make sure your input is OK and then try to re-run your job or '
                   f'<a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME}%20Run%20Number:%20{run_number}">'
                   f'contact us'
                   f'</a> '
                   f'for further information.<br>'
                   f'</div>\n')


def add_closing_html_tags(html_path, CONSTS, run_number):
    FORMER_MSG = 'Effectidor is now processing your request. This page will be automatically updated every {CONSTS.RELOAD_INTERVAL} seconds (until the job is done). You can also reload it manually. Once the job has finished, the output will appear below.'
    update_html(html_path, FORMER_MSG, '')  # remove "web server is now processing your request" message
    update_html(html_path, 'progress-bar-striped active', 'progress-bar-striped')  # stop_progress_bar

    append_to_html(html_path, f'''
            <hr>
                <h4 class=footer>
                    <p align='center'>Questions and comments are welcome! Please
                        <span class="admin_link"> 
                        <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME.upper()}%20Run%20Number%20{run_number}">contact us</a> 
                        </span>
                    </p>
                </h4>
                <div id="bottom_links" align="center"><span class="bottom_link">
                <a href="{CONSTS.WEBSERVER_URL}" target="_blank">Home</a>  
            </span>
        </div>
        <br><br><br>
    </body>
</html>''')

    # have to be last thing that is done
    sleep(2*CONSTS.RELOAD_INTERVAL)
    update_html(html_path, CONSTS.RELOAD_TAGS, f'<!--{CONSTS.RELOAD_TAGS}-->')  # stop refresh



def initialize_html(CONSTS, output_dir_path, html_path):

    path_tokens = output_dir_path.split('/')
    # e.g., "/bioseq/data/results/sincopa/12345678/outputs"
    run_number = path_tokens[path_tokens.index(CONSTS.WEBSERVER_NAME) + 1]

    update_html(html_path, 'QUEUED', 'RUNNING')
    #update_html(html_path, CONSTS.PROGRESS_BAR_ANCHOR, CONSTS.PROGRESS_BAR_TAG)  # add progress bar

    return run_number


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_ORFs_path',
                            help='A path to a DNA ORFs file.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_dir_path',
                            help='A path to a folder in which the output files will be created.',
                            type=lambda path: path.rstrip('/'))
        parser.add_argument('--input_effectors_path', default='',
                            help='A path to a DNA fasta with positive samples. '
                                 'All samples in this file should be in the input ORFs file as well.')
        parser.add_argument('--input_T3Es_path', default='',
                            help='A path to protein fasta with T3Es records of other bacteria.')
        parser.add_argument('--host_proteome_path', default='',
                            help='A path to a zip archive with protein fasta files of host proteome.')
        parser.add_argument('--no_T3SS', default='',
                            help='A path to a zip archive with protein fasta files with related bacteria non T3SS proteomes.')
        parser.add_argument('--genome_path', default='',
                            help='A path to a fasta file with full genome record.')
        parser.add_argument('--gff_path', default='',
                            help='A path to a GFF file.')
        parser.add_argument('--html_path', default=None,
                            help='A path to an html file that will be updated during the run.')
        
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
        
        parser.add_argument('--homology_search', help='search additional effectors based on homology to internal dataset', action='store_true')
        parser.add_argument('--translocation_signal',help='extract translocation signal feature', action='store_true')

        parser.add_argument('--full_genome', help='to extract genome organization features', action='store_true')
        parser.add_argument('--PIP', help='look for PIP-box in promoters', action='store_true')
        parser.add_argument('--hrp', help='look for hrp-box in promoters', action='store_true')
        parser.add_argument('--mxiE', help='look for mxiE-box in promoters', action='store_true')
        parser.add_argument('--exs', help='look for exs-box in promoters', action='store_true')
        parser.add_argument('--tts', help='look for tts-box in promoters', action='store_true')
        parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to', default='pupkolabr')

        args = parser.parse_args()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        #if args.full_genome:
        PIP_flag = args.PIP
        hrp_flag = args.hrp
        mxiE_flag = args.mxiE
        exs_flag = args.exs
        tts_flag = args.tts
        
        main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS, full_genome=args.full_genome, PIP=PIP_flag, hrp=hrp_flag, mxiE=mxiE_flag, exs=exs_flag, tts=tts_flag, homology_search=args.homology_search, signal=args.translocation_signal)
        
       #     if args.homology_search:
       #         main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS, full_genome=True, PIP=PIP_flag, hrp=hrp_flag, mxiE=mxiE_flag, exs=exs_flag, tts=tts_flag, homology_search=True)
       #     else:
       #         main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS, full_genome=True, PIP=PIP_flag, hrp=hrp_flag, mxiE=mxiE_flag, exs=exs_flag, tts=tts_flag)
       # else:
       #     if args.homology_search:
       #         main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS, homology_search=True)
       #     else:
       #         main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS)