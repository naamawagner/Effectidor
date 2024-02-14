import sys
#sys.path.append('/bioseq/effectidor/auxiliaries')
import os
import logging
import Bio.SeqUtils
import effectidor_CONSTANTS as CONSTS  # from /effectidor/auxiliaries
from time import sleep,time
from auxiliaries import fail,update_html,append_to_html # from /effectidor/auxiliaries
from Bio import SeqIO
from T3Es_wrapper import effectors_learn
import shutil
import re
import subprocess

data_dir = CONSTS.EFFECTIDOR_DATA
blast_datasets_dir = f'{data_dir}/blast_data'
scripts_dir = CONSTS.EFFECTIDOR_EXEC

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
        legal_chars = set(Bio.SeqUtils.IUPACData.extended_protein_letters+'*-')
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
                    curated_content = curated_content.replace('.','_')
            if flag: # protein
                recs = SeqIO.parse(fasta_path,'fasta')
                for rec in recs:
                    seq = str(rec.seq)
                    if len(seq)==0:
                        return f'FASTA file {os.path.basename(fasta_path)} in {input_name} input contains empty records!'
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
    if len(ORFs_recs) < 900:
        return f'The ORFs file contains only {str(len(ORFs_recs))} records. Make sure this file contains all the ORFs (open reading frames) in the genome - Effectidor is designed to analyze full genomes and not a sample of genes. Also, make sure this file contains ORFs and not full genome sequence! The full genome sequence can be uploaded in the advanced options.'
    elif len(ORFs_recs) > 10000:
        return f'The ORFs file contains {str(len(ORFs_recs))} records! Every file should contain data of a single bacterial genome. Make sure this file contains all the ORFs (open reading frames) in the genome, and only the ORFs of one genome. This number cannot exceed 10,000 ORFs per genome. If it contains data of multiple genomes, separate them to different files (compressed together in a ZIP archive) such that every file will contain the ORFs of a single genome'

def verify_effectors_f(effectors_path, ORFs_path):
    logger.info(f'Validating effectors:{effectors_path}')
    effectors_recs = [rec.id for rec in SeqIO.parse(effectors_path,'fasta')]
    effectors_set = set([rec.id for rec in SeqIO.parse(effectors_path,'fasta')])
    ORFs_recs = set([rec.id for rec in SeqIO.parse(ORFs_path,'fasta')])
    if not effectors_set.issubset(ORFs_recs):
        not_in_ORFs = ', '.join([rec for rec in effectors_set.difference(ORFs_recs)])
        return f'Illegal effectors records. The following records IDs are in the effectors file and not in the ORFs file:<br>{not_in_ORFs}.<br><br>If these effectors are from other bacterial genomes, they can be supplied (in <b>protein</b> FASTA file) in the advanced options.\nIf these are effectors from the analyzed bacterial genome, they should be identical to the records available in the input ORFs file.'
    if len(effectors_set) != len(effectors_recs):
        more_than_once = ','.join([effector for effector in effectors_set if effectors_recs.count(effector)>1])
        return f'Illegal effectors records. The following records IDs appear more than once in the file:<br>{more_than_once}.'
    
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
                return f'In {name}:<br>{error_msg}'
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
        return f'{name} contains non unique records. The following records appear more than once:<br><br>{non_unique_ids}'
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
        return f'{name} contains non unique records. The following records appear more than once:<br><br>{non_unique_ids}'

def validate_gff_format(gff_f,genome_name=''):
    with open(gff_f) as in_f:
        first_line= in_f.readline()
        if not (first_line.startswith("##gff-version") or first_line.startswith("##sequence-region")):
            return f'''Illegal GFF format! The GFF file for the genome {genome_name} is not in <a href="http://gmod.org/wiki/GFF3" target="_blank">GFF3</a> format!<br>
        The content in the file does not start with a line of ##gff-version or ##sequence-region as must be'''
            
def validate_gff(gff_f,ORFs_f,genome_name=''):
    logger.info(f'Validating GFF')
    error_msg = validate_gff_format(gff_f)
    if error_msg:
        return (error_msg)
    logger.info(f'Validating GFF coverage')
    import fasta_parser
    import pip_box_features
    locus_dic = fasta_parser.parse_ORFs(ORFs_f)
    CDS_set = set()
    RNA_set = set()
    with open(gff_f) as in_f:
        content = in_f.read().replace('|','_')
        content = content.replace('.','_')
        with open(gff_f,'w') as out_f:
            out_f.write(content)
    CDS,RNA = pip_box_features.parse_gff(gff_f,locus_dic)
    logger.info(f'{locus_dic.keys()}')
    logger.info(f'CDS:{CDS}\nRNA:{RNA}')
    CDS_set.update(CDS)
    RNA_set.update(RNA)
    not_in_gff = [locus for locus in locus_dic if (locus not in CDS_set and locus not in RNA_set)]
    if len(not_in_gff)>0:
        return f'''For the genome {genome_name} there are records in your ORFs input that are not available in the corresponding GFF input. Please revise your input and submit again.<br>
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
        
def validate_genome_and_gff(gff_f,genome_f,ORFs_f):
    logger.info(f'Validating genome and gff')
    import fasta_parser
    import pip_box_features
    locus_dic = fasta_parser.parse_ORFs(ORFs_f)
    regions = []
    locus_area_d,circulars = pip_box_features.parse_gff_to_CDS_loc(gff_f,locus_dic)
    for region in locus_area_d:
        regions.append(region)
    genome_recs = []
    recs_ids = [rec.id for rec in SeqIO.parse(genome_f,'fasta')]
    genome_recs.extend(recs_ids)
    if len(genome_recs) > len(set(genome_recs)):
        non_unique = []
        for ID in genome_recs:
            if genome_recs.count(ID) > 1:
                non_unique.append(ID)
        non_unique = set(non_unique)
        return f'Several contigs appear more than once within the genomic sequences data. These contigs IDs are:<br>{", ".join(non_unique)}<br>Please make sure to upload the genomic sequences such that each sequence will appear only once in the file, with unique IDs, matching to the IDs listed in the GFF file.'
    #return f'IDs: {str(genome_recs)}'
    not_in_genome = []
    for region in regions:
        if region not in genome_recs:
            not_in_genome.append(region)
    if len(not_in_genome)>0:
        return f'The following regions from the GFF file are not found in the full genome file:<br>{", ".join(not_in_genome)}.<br>Make sure the regions names are matching between the GFF and genome files and re-submit the job.'


def validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path, gff_path, no_T3SS_path, error_path):
    logger.info('Validating input...')
    if ORFs_path.endswith('.zip'):
        ORFs_genomes = []
        os.makedirs(f'{output_dir_path}/ORFs_tmp')
        shutil.unpack_archive(ORFs_path,f'{output_dir_path}/ORFs_tmp')
        for file in os.listdir(f'{output_dir_path}/ORFs_tmp'):
            if os.path.isfile(f'{output_dir_path}/ORFs_tmp/{file}') and not file.startswith('_') and not file.startswith('.'):
                error_msg = verify_fasta_format(f'{output_dir_path}/ORFs_tmp/{file}','DNA',f'{file} in ORFs archive')
                if error_msg:
                    error_msg = f'Illegal fasta files in {file} in ORFs archive: {error_msg}<br>This archive is expected to contain only <b>DNA</b> FASTA files.'
                    fail(error_msg,error_path)
                error_msg = validate_set(f'{output_dir_path}/ORFs_tmp/{file}', f'File {file} in your ORFs input')
                if error_msg:
                    fail(error_msg, error_path)
                error_msg = verify_ORFs(f'{output_dir_path}/ORFs_tmp/{file}')
                if error_msg:
                    fail(f'For file {file} in your ORFs input: {error_msg}')
                genome_name = '_'.join(file.split('/')[-1].split('.')[0].split(' '))
                ORFs_genomes.append(genome_name)
                os.makedirs(os.path.join(output_dir_path,'Effectidor_runs',genome_name),exist_ok=True)
                shutil.move(os.path.join(output_dir_path,'ORFs_tmp',file),os.path.join(output_dir_path,'Effectidor_runs',genome_name,'ORFs.fasta'))
        ORFs_set = set(ORFs_genomes)
    else:
        error_msg = verify_fasta_format(ORFs_path,'DNA','input ORFs')
        if error_msg:
            error_msg = f'Illegal fasta file in ORFs input: {error_msg}<br>This input is expected to hold a <b>DNA</b> FASTA file.'
            fail(error_msg,error_path)
        error_msg = validate_set(ORFs_path,'Your ORFs input')
        if error_msg:
            fail(error_msg,error_path)
        error_msg = verify_ORFs(ORFs_path)
        if error_msg:
            fail(error_msg, error_path)

    if gff_path:
        # gff
        if gff_path.endswith('.zip'):
            gff_genomes = []
            os.makedirs(f'{output_dir_path}/gff_tmp')
            shutil.unpack_archive(gff_path,f'{output_dir_path}/gff_tmp')
            for file in os.listdir(f'{output_dir_path}/gff_tmp'):
                if os.path.isfile(f'{output_dir_path}/gff_tmp/{file}') and not file.startswith('_') and not file.startswith('.'):
                    genome_name = '_'.join(file.split('/')[-1].split('.')[0].split(' '))
                    gff_genomes.append(genome_name)
                    if genome_name in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                        shutil.move(os.path.join(output_dir_path, 'gff_tmp', file),os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.gff3'))
                        error_msg = validate_gff(os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.gff3'),
                                                 os.path.join(output_dir_path,'Effectidor_runs',genome_name,'ORFs.fasta'),genome_name)
                        if error_msg:
                            fail(error_msg, error_path)
            gff_set = set(gff_genomes)
            if not gff_set == ORFs_set:
                missing_gff = ORFs_set.difference(gff_set)
                if missing_gff: # If there are ORFs files without a corresponding GFF file, we need to crush it. Why? Because for every genome we must calculate the same features. We cannot calculate the GFF features only for some of the genomes.
                    error_msg = f'There are several genomes with missing data.<br>The following genomes did not have a matching gff file:<br>{", ".join(missing_gff)}.<br><br>Make sure all the files that refer to the same genome have matching names, as instructed.'
                    fail(error_msg, error_path)

        else:
            error_msg = validate_gff(gff_path,ORFs_path)
            if error_msg:
                fail(error_msg,error_path)
        if genome_path:
            # genome
            if genome_path.endswith('.zip'):
                os.makedirs(f'{output_dir_path}/full_genome_tmp')
                shutil.unpack_archive(genome_path,f'{output_dir_path}/full_genome_tmp')
                full_genome_names = []
                for file in os.listdir(f'{output_dir_path}/full_genome_tmp'):
                    print(file)
                    if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{output_dir_path}/full_genome_tmp/{file}'): # discard system files and directories
                        error_msg = verify_fasta_format(f'{output_dir_path}/full_genome_tmp/{file}','DNA', f'{file} in full genome archive')
                        if error_msg:
                            error_msg = f'Illegal fasta files in {file} in full genome archive: {error_msg}'
                            fail(error_msg,error_path)
                        genome_name = '_'.join(file.split('/')[-1].split('.')[0].split(' '))
                        print(f'file:{file}, genome:{genome_name}')
                        full_genome_names.append(genome_name)
                        if genome_name in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                            shutil.move(os.path.join(output_dir_path, 'full_genome_tmp', file), os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.fasta'))
                            error_msg = validate_genome_and_gff(os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.gff3'),
                                                                os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.fasta'),
                                                                os.path.join(output_dir_path,'Effectidor_runs',genome_name,'ORFs.fasta'))
                            if error_msg:
                                fail(f'In genome {genome_name}:<br>{error_msg}')
                full_genome_set = set(full_genome_names)
                if not gff_set==full_genome_set==ORFs_set:
                    missing_full_genome = ORFs_set.difference(full_genome_set)
                    if missing_full_genome: # If there are ORFs files without a corresponding full genome file, we need to crush it. Why? Because for every genome we must calculate the same features. We cannot calculate the full genome dependent features only for some of the genomes.
                        error_msg = f'There are several genomes with missing data.<br>The following genomes did not have a matching full genome fasta file:<br>{", ".join(missing_full_genome)}.<br><br>Make sure all the files that refer to the same genome have matching names, as instructed.'
                        fail(error_msg,error_path)
            else:
                error_msg = verify_fasta_format(genome_path, 'DNA', 'full genome')
                if error_msg:
                    fail(error_msg)
                if gff_path:
                    error_msg = validate_genome_and_gff(gff_path, genome_path,
                                                        ORFs_path)
                    if error_msg:
                        fail(error_msg, error_path)

    if effectors_path:
        if effectors_path.endswith('.zip'):
            os.makedirs(f'{output_dir_path}/effectors_tmp')
            shutil.unpack_archive(effectors_path,f'{output_dir_path}/effectors_tmp')
            for file in os.listdir(f'{output_dir_path}/effectors_tmp'):
                if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{output_dir_path}/effectors_tmp/{file}'):
                    error_msg = verify_fasta_format(f'{output_dir_path}/effectors_tmp/{file}', 'DNA', f'input effectors {file}')
                    if error_msg:
                        error_msg = f'Illegal effectors input: {error_msg}'
                        fail(error_msg, error_path)
                    error_msg = validate_set(f'{output_dir_path}/effectors_tmp/{file}', f'File {file} in your effectors input')
                    if error_msg:
                        fail(error_msg,error_path)
                    genome_name = '_'.join(file.split('/')[-1].split('.')[0].split(' '))
                    if genome_name in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                        shutil.move(os.path.join(output_dir_path, 'effectors_tmp', file), os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'effectors.fasta'))
                        error_msg = verify_effectors_f(os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'effectors.fasta'),
                                                       os.path.join(output_dir_path, 'Effectidor_runs', genome_name,'ORFs.fasta'))
                        if error_msg:
                            error_msg = f'In {genome}: {error_msg}'
                            fail(error_msg, error_path)
        else:
            error_msg = verify_fasta_format(effectors_path,'DNA', 'input effectors')
            if error_msg:
                error_msg = f'Illegal effectors file: {error_msg}'
                fail(error_msg, error_path)
            error_msg = verify_effectors_f(effectors_path,ORFs_path)
            if error_msg:
                fail(error_msg, error_path)
            error_msg = validate_set(effectors_path,'Your effectors file')
            if error_msg:
                fail(error_msg,error_path)
    if input_T3Es_path:
        error_msg = verify_fasta_format(input_T3Es_path,'protein', 'effectors for homology search')
        if error_msg:
            fail(error_msg,error_path)
        error_msg = validate_set(input_T3Es_path,'Your effectors file for homology search')
        if error_msg:
            fail(error_msg,error_path)
        os.makedirs(f'{output_dir_path}/blast_data',exist_ok=True)
        cmd = f'cp {blast_datasets_dir}/*.faa {output_dir_path}/blast_data/'
        subprocess.check_output(cmd, shell=True)
        blast_datasets_dir = f'{output_dir_path}/blast_data'
        eff1_recs = SeqIO.parse(f'{output_dir_path}/blast_data/T3Es.faa', 'fasta')
        eff_l = list(eff1_recs)
        seqs = [rec.seq for rec in eff_l]
        eff2_recs = SeqIO.parse(input_T3Es_path, 'fasta')
        for rec in eff2_recs:
            if rec.seq not in seqs:
                eff_l.append(rec)
                seqs.append(rec.seq)
        SeqIO.write(eff_l, f'{output_dir_path}/blast_data/T3Es.faa', 'fasta')

        if os.path.exists(f'{output_dir_path}/Effectidor_runs'):
            for genome in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                shutil.copy(input_T3Es_path,f'{output_dir_path}/Effectidor_runs/{genome}')
    if host_proteome:
        error_msg = verify_zip(host_proteome,'Host data')
        if error_msg:
            fail(error_msg,error_path)
        os.makedirs(f'{output_dir_path}/blast_data/temp_extract', exist_ok=True)
        shutil.unpack_archive(host_proteome, f'{output_dir_path}/blast_data/temp_extract')
        recs = []
        for f in os.listdir(f'{output_dir_path}/blast_data/temp_extract'):
            if not f.startswith('_') and not f.startswith('.') and os.path.isfile(
                    f'{output_dir_path}/blast_data/temp_extract/{f}'):
                f_recs = SeqIO.parse(f'{output_dir_path}/blast_data/temp_extract/{f}', 'fasta')
                for rec in f_recs:
                    recs.append(rec)
        SeqIO.write(recs, f'{output_dir_path}/blast_data/host.faa', 'fasta')
        shutil.rmtree(f'{output_dir_path}/blast_data/temp_extract')

    if no_T3SS_path:
        error_msg = verify_zip(no_T3SS_path,'Close bacteria without T3SS data')
        if error_msg:
            fail(error_msg,error_path)
        os.makedirs(f'{output_dir_path}/blast_data/temp_extract', exist_ok=True)
        shutil.unpack_archive(no_T3SS_path, f'{output_dir_path}/blast_data/temp_extract')
        for file in os.listdir(f'{output_dir_path}/blast_data/temp_extract'):
            if not file.startswith('_') and not file.startswith('.') and os.path.isfile(
                    f'{output_dir_path}/blast_data/temp_extract/{file}'):
                new_name_l = file.replace(' ', '_').split('.')
                new_name = '.'.join(new_name_l[:-1]) + '.' + 'faa'
                os.rename(f'{output_dir_path}/blast_data/temp_extract/{file}',
                          f'{output_dir_path}/blast_data/temp_extract/{new_name}')
                shutil.move(f'{output_dir_path}/blast_data/temp_extract/{new_name}',
                            f'{output_dir_path}/blast_data')
        shutil.rmtree(f'{output_dir_path}/blast_data/temp_extract')

    if os.path.exists(f'{output_dir_path}/blast_data') and os.path.exists(f'{output_dir_path}/Effectidor_runs'):
        for genome in os.listdir(f'{output_dir_path}/Effectidor_runs'):
            os.makedirs(os.path.join(output_dir_path,'Effectidor_runs',genome,'blast_data'),exist_ok=True)
            for f in os.listdir(f'{output_dir_path}/blast_data'):
                shutil.copy(f'{output_dir_path}/blast_data/{f}',os.path.join(output_dir_path,'Effectidor_runs',genome,'blast_data'))
            #shutil.copytree(f'{output_dir_path}/blast_data',f'{output_dir_path}/Effectidor_runs/{genome}',dirs_exist_ok=True)

def cleanup_is_running(queues=('pupkolab','pupkoweb')):
    for q in queues:
        try:
            if subprocess.check_output(f'qstat {q} | grep cleanup_effec',shell=True):
            # a cleanup is currently running.
                return True
        except:
            pass
    # none of the queues is running cleanup (none of them returned True)
    return False

def cleanup_ran_today(path=r'/bioseq/data/results/effectidor/'):
    for f in os.listdir(path):
        if f.endswith('.ER'):
            f_path = os.path.join(path,f)
            ctime = os.stat(f_path).st_ctime
            if time() - ctime < 60*60*24:
                return True
    # cleanup did not run today
    return False

def main(ORFs_path, output_dir_path, effectors_path, input_T3Es_path, host_proteome, html_path, queue, genome_path, gff_path, no_T3SS, PIP=False, hrp=False, mxiE=False, exs=False, tts=False, homology_search=False, signal=False):
    '''
    try:
        if not cleanup_is_running() and not cleanup_ran_today():
            subprocess.call(f'/opt/pbs/bin/qsub /bioseq/effectidor/auxiliaries/remove_old_files.pbs',shell=True)
    except:
        pass
    '''
    error_path = f'{output_dir_path}/error.txt'
    try:
        if html_path:
            run_number = initialize_html(CONSTS, output_dir_path, html_path)
            #final_zip_path = f'{os.path.split(output_dir_path)[0]}/{CONSTS.WEBSERVER_NAME}_{run_number}'
    
        os.makedirs(output_dir_path, exist_ok=True)
    
        tmp_dir = f'{output_dir_path}/tmp'
        os.makedirs(tmp_dir, exist_ok=True)
        
        validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path, gff_path, no_T3SS, error_path)

        find_OGs_cmd = f'module load python/python-anaconda3.6.5;!@#python {os.path.join(scripts_dir,"find_OGs_in_genomes.py")} {output_dir_path}\tfind_OGs_effectidor\n'
        with open(os.path.join(output_dir_path,'find_OGs.cmds'),'w') as out_f:
            out_f.write(find_OGs_cmd)
        cmd = f'{os.path.join(scripts_dir,"q_submitter.py")} {os.path.join(output_dir_path,"find_OGs.cmds")} {output_dir_path} -q {queue}'
        subprocess.check_output(cmd, shell=True)

        if os.path.exists(f'{output_dir_path}/Effectidor_runs'):
            currently_running = []
            for genome in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                genome_ORFs_path = os.path.join(f'{output_dir_path}/Effectidor_runs',genome,'ORFs.fasta')
                genome_output_path = os.path.join(f'{output_dir_path}/Effectidor_runs',genome)
                parameters = f'{error_path} {genome_ORFs_path} {genome_output_path} --queue {queue}'
                if effectors_path:
                    input_effectors_path = os.path.join(output_dir_path, 'Effectidor_runs', genome, 'effectors.fasta')
                    parameters += f' --input_effectors_path {input_effectors_path}'
                if gff_path:
                    gff_file = os.path.join(output_dir_path, 'Effectidor_runs', genome, 'genome.gff3')
                    parameters += f' --gff_file {gff_file}'
                if genome_path:
                    full_genome_f = os.path.join(output_dir_path, 'Effectidor_runs', genome, 'genome.fasta')
                    parameters += f' --full_genome_f {full_genome_f}'
                if PIP:
                    parameters += ' --PIP'
                if hrp:
                    parameters += ' --hrp'
                if mxiE:
                    parameters += ' --mxiE'
                if exs:
                    parameters += ' --exs'
                if tts:
                    parameters += ' --tts'
                if homology_search:
                    parameters += ' --homology_search'
                if signal:
                    parameters += ' --translocation_signal'

                job_cmd = f'module load python/python-anaconda3.6.5;!@#python {os.path.join(scripts_dir,"T3Es_wrapper.py")} {parameters}\tEffectidor_features_{genome}\n'
                cmds_f = os.path.join(output_dir_path, 'Effectidor_runs', genome, 'features_wrapper.cmds')
                with open(cmds_f,'w') as job_f:
                    job_f.write(job_cmd)
                cmd = f'{os.path.join(scripts_dir,"q_submitter.py")} {cmds_f} {os.path.join(output_dir_path, "Effectidor_runs", genome)} -q {queue}'
                subprocess.check_output(cmd,shell=True)
                currently_running.append(genome)
                while len(currently_running) > 4: # features can be extracted for up to 5 genomes at a time, to avoid overload on the cluster.
                    sleep(60)
                    currently_running = [genome for genome in currently_running if not os.path.exists(os.path.join(output_dir_path, "Effectidor_runs", genome, 'features.csv'))]
            while len(currently_running) > 0: # wait until they all finish before proceeding with the next step.
                sleep(60)
                currently_running = [genome for genome in currently_running if not os.path.exists(
                    os.path.join(output_dir_path, "Effectidor_runs", genome, 'features.csv'))]

        else:
            effectors_learn(error_path, ORFs_path, effectors_path, output_dir_path, tmp_dir, queue, gff_path, genome_path, PIP=PIP, hrp=hrp, mxiE=mxiE, exs=exs, tts=tts, homology_search=homology_search, signal=signal)
        # add a check for failed features jobs...
        subprocess.check_output(['python', os.path.join(scripts_dir,'merge_features_for_OGs.py'), output_dir_path])

        # learning step

        low_quality_flag = False
        subprocess.check_output(['python', os.path.join(scripts_dir,'learning.py'), output_dir_path, 'OGs_features.csv'])
        if os.path.exists(f'{output_dir_path}/out_learning/learning_failed.txt'):
            low_quality_flag = True
            #return create_effectors_html(effectors_prots,ORFs_file,working_directory)
            #error_msg = 'Learning failed. It can be due to a small training set, or other reasons. For further details you can contact us.'
            #fail(error_msg,error_path)
        # making final output files and tables
        '''
        in_f = f'{output_dir_path}/out_learning/concensus_predictions.csv'
        out_f_normal = f'{output_dir_path}/out_learning/concensus_predictions_with_annotation.csv'
        out_f_pseudo = f'{output_dir_path}/out_learning/pseudogenes.csv'
        out_f_T3SS = f'{output_dir_path}/out_learning/T3SS.csv'
        annotations = ORFs_file
        out_for_html_normal = f'{output_dir_path}/out_learning/concensus_predictions_with_annotation_for_html.csv'
        out_for_html_pseudo = f'{output_dir_path}/out_learning/pseudogenes_predictions_with_annotation_for_html.csv'
        out_T3SS_for_html = f'{output_dir_path}/out_learning/T3SS_for_html.csv'
        if organization:
            gff_dir = f'{working_directory}/gff'
            add_annotations_to_predictions(in_f,out_f_normal,out_f_pseudo,annotations,out_f_T3SS,gff_dir)
        else:
            add_annotations_to_predictions(in_f,out_f_normal,out_f_pseudo,annotations,out_f_T3SS)
        from csv_to_colored_xlsx_converter import convert_csv_to_colored_xlsx
        convert_csv_to_colored_xlsx(out_f_normal)
        convert_csv_to_colored_xlsx(out_f_pseudo)
        if organization:
            gff_dir = f'{working_directory}/gff'
            add_annotations_to_predictions(in_f,out_for_html_normal,out_for_html_pseudo,annotations,out_T3SS_for_html,gff_dir,line_end='<br>')
        else:
            add_annotations_to_predictions(in_f,out_for_html_normal,out_for_html_pseudo,annotations,out_T3SS_for_html,line_end='<br>')
        predicted_table, positives_table, T3SS_table = make_html_tables(out_for_html_normal,out_T3SS_for_html)
        return predicted_table, positives_table, T3SS_table, low_quality_flag'''

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
                       <font color="red">WARNING: The predictions might be of low quality due to small positive set to train the classifier, or for other reasons.
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
                <h3><b>Positive samples that were used to train the model</b></h3>
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
                   f'<div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify"><h3><br>'
                   f'<font color="red">{error_msg}</font></h3><br><br>'
                   f'Please make sure your input is OK and then try to re-run your job or '
                   f'<a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME}%20Run%20Number:%20{run_number}">'
                   f'contact us'
                   f'</a> '
                   f'for further information.<br>'
                   f'</div><br>')


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
        
        main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path,args.host_proteome_path,
             args.html_path, args.queue_name, args.genome_path,args.gff_path, args.no_T3SS, PIP=PIP_flag, hrp=hrp_flag, mxiE=mxiE_flag,
             exs=exs_flag, tts=tts_flag, homology_search=args.homology_search, signal=args.translocation_signal)


#     if args.homology_search:
       #         main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS, full_genome=True, PIP=PIP_flag, hrp=hrp_flag, mxiE=mxiE_flag, exs=exs_flag, tts=tts_flag, homology_search=True)
       #     else:
       #         main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS, full_genome=True, PIP=PIP_flag, hrp=hrp_flag, mxiE=mxiE_flag, exs=exs_flag, tts=tts_flag)
       # else:
       #     if args.homology_search:
       #         main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS, homology_search=True)
       #     else:
       #         main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS)