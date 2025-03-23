import sys
import os
import logging
import Bio.SeqUtils
from time import sleep, time
from Bio import SeqIO
from T3Es_wrapper import effectors_learn
import shutil
import re
import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from add_annotations_to_predictions import make_html_tables
from csv_to_colored_xlsx_converter import convert_csv_to_colored_xlsx
import fasta_parser
from merge_features_for_OGs import get_ortho_dict

scripts_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = f'{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/data'
blast_datasets_dir = f'{data_dir}/blast_data'

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

ILLEGAL_CHARS = '\\;:,^`~\'\"'


def fail(error_msg, error_file_path):
    with open(error_file_path, 'w') as error_f:
        error_f.write(error_msg + '\n')
    raise Exception(error_msg)


def has_illegal_chars(s):
    # Check if the first word in the string (which is the record ID) contains illegal characters
    record_id = s.split(' ')[0]
    return any(char in ILLEGAL_CHARS for char in record_id)


def verify_fasta_format(fasta_path, Type, input_name):
    logger.info(f'Validating FASTA format:{fasta_path}')
    file_size = os.path.getsize(fasta_path)
    if int(float(file_size)) == 0:
        return f'{input_name} is empty!'

    flag = False
    if Type == 'DNA':
        legal_chars = set(
            Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() + Bio.SeqUtils.IUPACData.ambiguous_dna_letters)
    else:  # Type == 'protein'
        legal_chars = set(Bio.SeqUtils.IUPACData.extended_protein_letters + '*-')
        flag = True
    with open(fasta_path) as f:
        line_number = 0
        try:
            line = f.readline()
            line_number += 1
            while line == '\n':  # skip empty lines
                line = f.readline()
                line_number += 1
            if not line.startswith('>'):
                return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA ' \
                       f'format</a>. First line in fasta {input_name} starts with {line[0]} instead of ">". '
            if Type == 'DNA':
                if has_illegal_chars(line):
                    return f'Illegal format. First line in fasta {input_name} contains an illegal character in its ' \
                           f'first word (one of: {ILLEGAL_CHARS}).'
            previous_line_was_header = True
            putative_end_of_file = False
            curated_content = f'>{line[1:]}'.replace("|", "_")
            for line in f:
                line_number += 1
                line = line.strip()
                if not line:
                    if not putative_end_of_file:  # ignore trailing empty lines
                        putative_end_of_file = line_number
                    continue
                # if putative_end_of_file:  # non-empty line after empty line return f'Illegal <a
                # href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {
                # putative_end_of_file} in fasta {input_name} is empty.'
                if line.startswith('>'):
                    if previous_line_was_header:
                        return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" ' \
                               f'target="_blank">FASTA format</a>. fasta {input_name} contains an empty record. Both' \
                               f' lines {line_number - 1} and {line_number} start with ">". '
                    else:
                        previous_line_was_header = True
                        curated_content += f'>{line[1:]}\n'.replace("|", "_")
                        continue
                else:  # not a header
                    previous_line_was_header = False
                    for c in line:
                        if c not in legal_chars:
                            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" ' \
                                   f'target="_blank">FASTA format</a>. Line {line_number} in fasta {input_name} ' \
                                   f'contains illegal {Type} character "{c}". '
                    curated_content += f'{line}\n'
                    curated_content = curated_content.replace('.', '_')
            if flag:  # protein
                recs = SeqIO.parse(fasta_path, 'fasta')
                for rec in recs:
                    seq = str(rec.seq)
                    if len(seq) == 0:
                        return f'FASTA file {os.path.basename(fasta_path)} in {input_name} input contains empty records'
                    AGCT_count = seq.count('A') + seq.count('G') + seq.count('T') + seq.count('C') + seq.count('N')
                    if AGCT_count >= 0.95 * len(seq):
                        f = AGCT_count * 100 / len(seq)
                        return f'Protein fasta {input_name} seems to contain DNA records (record {rec.id} contains ' \
                               f'{float("%.2f" % f)}% DNA characters). Make sure all samples in the file are protein ' \
                               f'sequences and re-submit. '
        except UnicodeDecodeError as e:
            logger.info(e.args)
            line_number += 1  # the line that was failed to be read
            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA ' \
                   f'format</a>. Line {line_number} in fasta contains one (or more) non <a ' \
                   f'href="https://en.wikipedia.org/wiki/ASCII" target="_blank">ascii</a> character(s). '
    # override the old file with the curated content
    with open(fasta_path, 'w') as f:
        f.write(curated_content)


def verify_ORFs(ORFs_path):
    logger.info(f'Validating ORFs:{ORFs_path}')
    ORFs_recs = list(SeqIO.parse(ORFs_path, 'fasta'))
    if len(ORFs_recs) < 800:
        return f'The ORFs file contains only {str(len(ORFs_recs))} records. Make sure this file contains all the ' \
               f'ORFs (open reading frames) in the genome - Effectidor is designed to analyze full genomes and not ' \
               f'a sample of genes. Also, make sure this file contains ORFs and not full genome sequence! The full ' \
               f'genome sequence can be uploaded in the advanced options. '
    elif len(ORFs_recs) > 10000:
        return f'The ORFs file contains {str(len(ORFs_recs))} records! Every file should contain data of a single ' \
               f'bacterial genome. Make sure this file contains all the ORFs (open reading frames) in the genome, ' \
               f'and only the ORFs of one genome. This number cannot exceed 10,000 ORFs per genome. If it contains ' \
               f'data of multiple genomes, separate them to different files (compressed together in a ZIP archive) ' \
               f'such that every file will contain the ORFs of a single genome. '
    start_code = ['ATG', 'GTG', 'TTG']
    end_code = ['TAA', 'TAG', 'TGA']
    # correct frame if needed and specified before verifying the input contains coding sequences
    corrected_recs = []
    for rec in ORFs_recs:
        header = rec.description
        if 'frame=2' in header:
            rec.seq = rec.seq[1:]
        elif 'frame=3' in header:
            rec.seq = rec.seq[2:]
        corrected_recs.append(rec)
    ORFs_seqs = [rec.seq for rec in corrected_recs]
    coding_estimate = [seq for seq in ORFs_seqs if len(seq) % 3 == 0 and seq[:3] in start_code and seq[-3:] in end_code]
    est_cod_percent = int(100*len(coding_estimate)/len(ORFs_recs))
    if est_cod_percent < 40:
        return f'The ORFs file contains only {str(est_cod_percent)}% of coding sequences with valid start and end '\
               f'codons, and length divided by 3. Make sure this input contains <b>coding sequences</b>.'

    SeqIO.write(corrected_recs, ORFs_path, 'fasta')


def verify_effectors_f(effectors_path, ORFs_path):
    logger.info(f'Validating effectors:{effectors_path}')
    effectors_recs = [rec.id for rec in SeqIO.parse(effectors_path, 'fasta')]
    effectors_set = set([rec.id for rec in SeqIO.parse(effectors_path, 'fasta')])
    ORFs_recs = set([rec.id for rec in SeqIO.parse(ORFs_path, 'fasta')])
    if not effectors_set.issubset(ORFs_recs):
        not_in_ORFs = ', '.join([rec for rec in effectors_set.difference(ORFs_recs)])
        return f'Illegal effectors records. The following records IDs are in the effectors file and not in the ORFs ' \
               f'file:<br>{not_in_ORFs}.<br><br>If these effectors are from other bacterial genomes, they can be ' \
               f'supplied (in <b>protein</b> FASTA file) in the advanced options.\nIf these are effectors from the ' \
               f'analyzed bacterial genome, they should be identical to the records available in the input ORFs file. '
    if len(effectors_set) != len(effectors_recs):
        more_than_once = ','.join([effector for effector in effectors_set if effectors_recs.count(effector) > 1])
        return f'Illegal effectors records. The following records IDs appear more than once in the file:<br>{more_than_once}. '


def verify_genome_one_contig(genome_path, file_name):
    logger.info(f'Validating one contig:{genome_path}')
    recs = list(SeqIO.parse(genome_path, 'fasta'))
    if len(recs) > 1:  # it must be one contig
        return f'Illegal genome file of {file_name}. The FASTA file contains more than one record (more than one ' \
               f'contig). '


def verify_zip(file, name):
    logger.info(f'Validating zip:{file}')
    if not file.endswith('.zip'):
        return f'{name} not in zip format. Make sure to upload a zip archive and resubmit your job.'
    # unzip it to a tmp_dir
    if not os.path.exists(f'{"/".join(file.split("/")[:-1])}/zip_tmp'):
        os.makedirs(f'{"/".join(file.split("/")[:-1])}/zip_tmp')
    shutil.unpack_archive(file, '/'.join(file.split('/')[:-1]) + '/zip_tmp')
    logger.info('Validating zip: created zip_tmp')
    flag = False
    for f in os.listdir(f'{"/".join(file.split("/")[:-1])}/zip_tmp'):
        new_name = f.replace(' ', '_')
        os.rename(f'{"/".join(file.split("/")[:-1])}/zip_tmp/{f}',
                  f'{"/".join(file.split("/")[:-1])}/zip_tmp/{new_name}')
        f = new_name
        if not f.startswith('_') and not f.startswith('.') and os.path.isfile(
                f'{"/".join(file.split("/")[:-1])}/zip_tmp/{f}'):
            error_msg = verify_fasta_format(f'{"/".join(file.split("/")[:-1])}/zip_tmp/{f}', 'protein', f)
            if error_msg:
                return f'In {name}:<br>{error_msg}'
            flag = True
            # subprocess.check_output(f'rm {"/".join(file.split("/")[:-1])}/zip_tmp/{f}', shell=True)
            os.remove(os.path.join("/".join(file.split("/")[:-1]), 'zip_tmp', f))
            logger.info(f'Validating zip: removed {f}')
    if not flag:
        return f'{name} contains no valid files. Make sure to include protein fasta files in this archive!'


def validate_set(file, name):
    logger.info(f'Validating set:{file}')
    # by IDs
    recs = [rec.id for rec in SeqIO.parse(file, 'fasta')]
    if len(recs) > len(set(recs)):
        non_unique = []
        for rec_id in set(recs):
            if recs.count(rec_id) > 1:
                non_unique.append(rec_id)
        non_unique_ids = '<br>'.join(non_unique)
        return f'{name} contains non unique records. The following records appear more than once:<br><br>{non_unique_ids}'
    # by locus_tags
    recs = SeqIO.parse(file, 'fasta')
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


def validate_gff_format(gff_f, genome_name=''):
    with open(gff_f) as in_f:
        first_line = in_f.readline()
        if not (first_line.startswith("##gff-version") or first_line.startswith("##sequence-region")):
            return f'''Illegal GFF format! The GFF file for the genome {genome_name} is not in <a 
            href="http://gmod.org/wiki/GFF3" target="_blank">GFF3</a> format!<br> The content in the file does not 
            start with a line of ##gff-version or ##sequence-region as must be '''


def validate_gff(gff_f, ORFs_f, genome_name=''):
    logger.info(f'Validating GFF: {gff_f}')
    error_msg = validate_gff_format(gff_f)
    if error_msg:
        return error_msg
    logger.info(f'Validating GFF coverage')
    import fasta_parser
    import pip_box_features
    locus_dic = fasta_parser.parse_ORFs(ORFs_f)
    CDS_set = set()
    RNA_set = set()
    with open(gff_f) as in_f:
        content = in_f.read().replace('|', '_')
        content = content.replace('.', '_')
        with open(gff_f, 'w') as out_f:
            out_f.write(content)
    CDS, RNA = pip_box_features.parse_gff(gff_f, locus_dic)
    #logger.info(f'{locus_dic.keys()}')
    #logger.info(f'CDS:{CDS}\nRNA:{RNA}')
    CDS_set.update(CDS)
    RNA_set.update(RNA)
    not_in_gff = [locus for locus in locus_dic if (locus not in CDS_set and locus not in RNA_set)]
    if len(not_in_gff) > 0:
        return f'''For the genome {genome_name} there are records in your ORFs input that are not available in the 
        corresponding GFF input. Please revise your input and submit again.<br> These records are:<br>
        {", ".join(not_in_gff)} '''
    if len(RNA_set) > 0:
        recs = SeqIO.parse(ORFs_f, 'fasta')
        cds_recs = []
        for rec in recs:
            if rec.id in RNA_set:
                continue
            elif re.search(r'locus_tag=(\w+)', rec.description):
                if re.search(r'locus_tag=(\w+)', rec.description).group(1) in RNA_set:
                    continue
            cds_recs.append(rec)
        SeqIO.write(cds_recs, ORFs_f, 'fasta')


def validate_genome_and_gff(gff_f, genome_f, ORFs_f):
    logger.info(f'Validating genome and gff: {gff_f}, {genome_f}')
    import fasta_parser
    import pip_box_features
    locus_dic = fasta_parser.parse_ORFs(ORFs_f)
    regions = []
    locus_area_d, circulars = pip_box_features.parse_gff_to_CDS_loc(gff_f, locus_dic)
    for region in locus_area_d:
        regions.append(region)
    genome_recs = []
    recs_ids = [rec.id for rec in SeqIO.parse(genome_f, 'fasta')]
    genome_recs.extend(recs_ids)
    if len(genome_recs) > len(set(genome_recs)):
        non_unique = []
        for ID in genome_recs:
            if genome_recs.count(ID) > 1:
                non_unique.append(ID)
        non_unique = set(non_unique)
        return f'''Several contigs appear more than once within the genomic sequences data. These contigs IDs are:<br>
                {", ".join(non_unique)}<br>Please make sure to upload the genomic sequences such that each sequence 
                will appear only once in the file, with unique IDs, matching to the IDs listed in the GFF file.'''
        # return f'IDs: {str(genome_recs)}'
    not_in_genome = []
    for region in regions:
        if region not in genome_recs:
            not_in_genome.append(region)
    if len(not_in_genome) > 0:
        return f'''The following regions from the GFF file are not found in the full genome file:<br>
                {", ".join(not_in_genome)}.<br>Make sure the regions names are matching between the GFF 
                and genome files and re-submit the job. '''


def validate_OGs_table(OGs_table, ORFs_f, output_dir_path):
    logger.info(f'Validating OGs table: {OGs_table}, {ORFs_f}')
    genomes_orthogroup_dict = get_ortho_dict(OGs_table)
    logger.info('1')
    flatten_ortho_dict = {key: genomes_orthogroup_dict[dic_name][key] for dic_name in genomes_orthogroup_dict
                          for key in genomes_orthogroup_dict[dic_name]}
    logger.info('2')
    table = pd.read_csv(OGs_table)
    logger.info(f'Have read table successfully')
    if os.path.exists(os.path.join(output_dir_path, 'Effectidor_runs')):
        logger.info(f'multiple genomes run')
        genomes = os.listdir(os.path.join(output_dir_path, 'Effectidor_runs'))
        header = list(table.columns)
        for i in range(1, len(header)):
            header[i] = header[i].replace(' ', '_')
        table.columns = header
        table.to_csv(OGs_table, index=False)
        OGs_genomes = set(header[1:])
        ORFs_genomes = set(genomes)
        if not OGs_genomes == ORFs_genomes:
            genomes_not_in_OGs = ORFs_genomes.difference(OGs_genomes)
            if genomes_not_in_OGs:
                return f'The following genomes are not found in the OGs input: {", ".join(genomes_not_in_OGs)}'

        for genome in genomes:
            locus_dic = fasta_parser.parse_ORFs(os.path.join(output_dir_path, 'Effectidor_runs', genome, 'ORFs.fasta'))
            for locus in locus_dic:
                if locus not in flatten_ortho_dict:
                    return f'The OGs table does not hold locus {locus} from genome {genome}!'
    else:
        logger.info(f'single genome run')
        header = list(table.columns)
        if len(header) > 2:
            return 'The OGs table seems to hold data of more than one genome, while the analysis is on a single genome'
        header[1] = 'genome_ORFs'
        table.columns = header
        table.to_csv(OGs_table, index=False)
        locus_dic = fasta_parser.parse_ORFs(ORFs_f)
        for locus in locus_dic:
            if locus not in flatten_ortho_dict:
                return f'The OGs table does not hold locus {locus}'


def validate_file_type(f, allowed_types):
    if not any([f.endswith(f'.{t}') for t in allowed_types]):
        return 'input not of the allowed types!'


def validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path, gff_path,
                   no_T3SS_path, error_path, OG_input_f):
    logger.info('Validating input...')
    ORFs_types = ['zip', 'txt', 'fasta', 'fna']
    error_msg = validate_file_type(ORFs_path, ORFs_types)
    if error_msg:
        error_msg = f'ORFs {error_msg} The allowed types for this input are {", ".join(ORFs_types)}'
        fail(error_msg, error_path)
    if effectors_path:
        error_msg = validate_file_type(effectors_path, ORFs_types)
        if error_msg:
            error_msg = f'effectors {error_msg} The allowed types for this input are {", ".join(ORFs_types)}'
            fail(error_msg, error_path)
    if input_T3Es_path:
        allowed_types_T3Es = ['txt', 'fasta', 'faa']
        error_msg = validate_file_type(input_T3Es_path, allowed_types_T3Es)
        if error_msg:
            error_msg = f'T3Es {error_msg} The allowed types for this input are {", ".join(allowed_types_T3Es)}'
            fail(error_msg, error_path)
    zip_input = ['zip']
    if host_proteome:
        error_msg = validate_file_type(host_proteome, zip_input)
        if error_msg:
            error_msg = f'Host {error_msg} The allowed types for this input are {", ".join(zip_input)}'
            fail(error_msg, error_path)
    if genome_path:
        error_msg = validate_file_type(genome_path, ORFs_types)
        if error_msg:
            error_msg = f'genome {error_msg} The allowed types for this input are {", ".join(ORFs_types)}'
            fail(error_msg, error_path)
    if gff_path:
        gff_types = ['gff', 'zip', 'txt']
        error_msg = validate_file_type(gff_path, gff_types)
        if error_msg:
            error_msg = f'GFF {error_msg} The allowed types for this input are {", ".join(gff_types)}'
            fail(error_msg, error_path)
    if no_T3SS_path:
        error_msg = validate_file_type(no_T3SS_path, zip_input)
        if error_msg:
            error_msg = f'no_T3SS {error_msg} The allowed types for this input are {", ".join(zip_input)}'
            fail(error_msg, error_path)

    if ORFs_path.endswith('.zip'):
        logger.info('multiple genomes run')
        if gff_path:
            if not gff_path.endswith('.zip'):
                error_msg = 'ORFs input in ZIP and GFF input not! Both inputs must be of the same type (zip for ' \
                            'pan genome, single file otherwise)'
                fail(error_msg, error_path)
        if genome_path:
            if not genome_path.endswith('.zip'):
                error_msg = 'ORFs input in ZIP and genome input not! Both inputs must be of the same type (zip for ' \
                            'pan genome, single file otherwise)'
                fail(error_msg, error_path)
        ORFs_genomes = []
        os.makedirs(f'{output_dir_path}/ORFs_tmp')
        shutil.unpack_archive(ORFs_path, f'{output_dir_path}/ORFs_tmp')
        number_of_genomes = 0
        for file in os.listdir(f'{output_dir_path}/ORFs_tmp'):
            if os.path.isfile(f'{output_dir_path}/ORFs_tmp/{file}') and not file.startswith(
                    '_') and not file.startswith('.'):
                error_msg = verify_fasta_format(f'{output_dir_path}/ORFs_tmp/{file}', 'DNA', f'{file} in ORFs archive')
                if error_msg:
                    error_msg = f'Illegal fasta files in {file} in ORFs archive: {error_msg}<br>This archive is ' \
                                f'expected to contain only <b>DNA</b> FASTA files. '
                    fail(error_msg, error_path)
                error_msg = validate_set(f'{output_dir_path}/ORFs_tmp/{file}', f'File {file} in your ORFs input')
                if error_msg:
                    fail(error_msg, error_path)
                error_msg = verify_ORFs(f'{output_dir_path}/ORFs_tmp/{file}')
                if error_msg:
                    fail(f'For file {file} in your ORFs input: {error_msg}', error_path)
                genome_name = '_'.join(file.split('/')[-1].split('.')[0].split(' '))
                ORFs_genomes.append(genome_name)
                os.makedirs(os.path.join(output_dir_path, 'Effectidor_runs', genome_name), exist_ok=True)
                shutil.move(os.path.join(output_dir_path, 'ORFs_tmp', file),
                            os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'ORFs.fasta'))
                number_of_genomes += 1
        if number_of_genomes == 0:
            error_msg = 'No files were found in the ORFs input! Make sure the files in the ZIP archive are not inside ' \
                        'directories.'
            fail(error_msg, error_path)
        logger.info('deleting temp ORFs')
        shutil.rmtree(f'{output_dir_path}/ORFs_tmp', ignore_errors=True)
        logger.info('deleted temp ORFs')
        ORFs_set = set(ORFs_genomes)
    else:
        logger.info('single genomes run')
        if gff_path:
            if gff_path.endswith('.zip'):
                error_msg = 'ORFs input is a single file and GFF input is in a zip archive! Both inputs must be of ' \
                            'the same type (zip for pan genome, single files otherwise)'
                fail(error_msg, error_path)
        if genome_path:
            if genome_path.endswith('.zip'):
                error_msg = 'ORFs input is a single file and genome input is in a zip archive! Both inputs must be ' \
                            'of the same type (zip for pan genome, single file otherwise)'
                fail(error_msg, error_path)
        error_msg = verify_fasta_format(ORFs_path, 'DNA', 'input ORFs')
        if error_msg:
            error_msg = f'Illegal fasta file in ORFs input: {error_msg}<br>This input is expected to hold a ' \
                        f'<b>DNA</b> FASTA file. '
            fail(error_msg, error_path)
        error_msg = validate_set(ORFs_path, 'Your ORFs input')
        if error_msg:
            fail(error_msg, error_path)
        error_msg = verify_ORFs(ORFs_path)
        if error_msg:
            fail(error_msg, error_path)
    error_msg = validate_OGs_table(OG_input_f, ORFs_path, output_dir_path)
    if error_msg:
        fail(error_msg, error_path)
    if gff_path:
        # gff
        if gff_path.endswith('.zip'):
            gff_genomes = []
            os.makedirs(f'{output_dir_path}/gff_tmp')
            shutil.unpack_archive(gff_path, f'{output_dir_path}/gff_tmp')
            for file in os.listdir(f'{output_dir_path}/gff_tmp'):
                if os.path.isfile(f'{output_dir_path}/gff_tmp/{file}') and not file.startswith(
                        '_') and not file.startswith('.'):
                    genome_name = '_'.join(file.split('/')[-1].split('.')[0].split(' '))
                    gff_genomes.append(genome_name)
                    if genome_name in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                        shutil.move(os.path.join(output_dir_path, 'gff_tmp', file),
                                    os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.gff3'))
                        error_msg = validate_gff(
                            os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.gff3'),
                            os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'ORFs.fasta'), genome_name)
                        if error_msg:
                            fail(error_msg, error_path)
            if len(gff_genomes) == 0:
                error_msg = 'No files were found in the GFF input! Make sure the files in the ZIP archive are not inside directories.'
                fail(error_msg, error_path)
            shutil.rmtree(f'{output_dir_path}/gff_tmp', ignore_errors=True)
            gff_set = set(gff_genomes)
            if not gff_set == ORFs_set:
                missing_gff = ORFs_set.difference(gff_set)
                if missing_gff:  # If there are ORFs files without a corresponding GFF file, we need to crush it.
                    # Why? Because for every genome we must calculate the same features. We cannot calculate the GFF
                    # features only for some genomes.
                    error_msg = f'There are several genomes with missing data.<br>The following genomes did not have ' \
                                f'a matching gff file:<br>{", ".join(missing_gff)}.<br><br>Make sure all the files ' \
                                f'that refer to the same genome have matching names, as instructed. '
                    fail(error_msg, error_path)

        else:
            error_msg = validate_gff(gff_path, ORFs_path)
            if error_msg:
                fail(error_msg, error_path)
        if genome_path:
            # genome
            if genome_path.endswith('.zip'):
                os.makedirs(f'{output_dir_path}/full_genome_tmp')
                shutil.unpack_archive(genome_path, f'{output_dir_path}/full_genome_tmp')
                full_genome_names = []
                for file in os.listdir(f'{output_dir_path}/full_genome_tmp'):
                    # print(file)
                    if not file.startswith('_') and not file.startswith('.') and os.path.isfile(
                            f'{output_dir_path}/full_genome_tmp/{file}'):  # discard system files and directories
                        error_msg = verify_fasta_format(f'{output_dir_path}/full_genome_tmp/{file}', 'DNA',
                                                        f'{file} in full genome archive')
                        if error_msg:
                            error_msg = f'Illegal fasta files in {file} in full genome archive: {error_msg}'
                            fail(error_msg, error_path)
                        genome_name = '_'.join(file.split('/')[-1].split('.')[0].split(' '))
                        # print(f'file:{file}, genome:{genome_name}')
                        full_genome_names.append(genome_name)
                        if genome_name in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                            shutil.move(os.path.join(output_dir_path, 'full_genome_tmp', file),
                                        os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.fasta'))
                            error_msg = validate_genome_and_gff(
                                os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.gff3'),
                                os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'genome.fasta'),
                                os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'ORFs.fasta'))
                            if error_msg:
                                fail(f'In genome {genome_name}:<br>{error_msg}', error_path)
                if len(full_genome_names) == 0:
                    error_msg = 'No files were found in the full genome input! Make sure the files in the ZIP archive are not inside directories.'
                    fail(error_msg, error_path)
                shutil.rmtree(f'{output_dir_path}/full_genome_tmp', ignore_errors=True)
                full_genome_set = set(full_genome_names)
                if not gff_set == full_genome_set == ORFs_set:
                    missing_full_genome = ORFs_set.difference(full_genome_set)
                    if missing_full_genome:  # If there are ORFs files without a corresponding full genome file,
                        # we need to crush it. Why? Because for every genome we must calculate the same features. We
                        # cannot calculate the full genome dependent features only for some of the genomes.
                        error_msg = f'There are several genomes with missing data.<br>The following genomes did not ' \
                                    f'have a matching full genome fas' \
                                    f'ta file:<br>{", ".join(missing_full_genome)}.<br><br>Make sure all the files ' \
                                    f'that refer to the same genome have matching names, as instructed. '
                        fail(error_msg, error_path)
            else:
                error_msg = verify_fasta_format(genome_path, 'DNA', 'full genome')
                if error_msg:
                    fail(error_msg, error_path)
                if gff_path:
                    error_msg = validate_genome_and_gff(gff_path, genome_path,
                                                        ORFs_path)
                    if error_msg:
                        fail(error_msg, error_path)

    if effectors_path:
        if effectors_path.endswith('.zip'):
            os.makedirs(f'{output_dir_path}/effectors_tmp')
            shutil.unpack_archive(effectors_path, f'{output_dir_path}/effectors_tmp')
            for file in os.listdir(f'{output_dir_path}/effectors_tmp'):
                if not file.startswith('_') and not file.startswith('.') and os.path.isfile(
                        f'{output_dir_path}/effectors_tmp/{file}'):
                    error_msg = verify_fasta_format(f'{output_dir_path}/effectors_tmp/{file}', 'DNA',
                                                    f'input effectors {file}')
                    if error_msg:
                        error_msg = f'Illegal effectors input: {error_msg}'
                        fail(error_msg, error_path)
                    error_msg = validate_set(f'{output_dir_path}/effectors_tmp/{file}',
                                             f'File {file} in your effectors input')
                    if error_msg:
                        fail(error_msg, error_path)
                    genome_name = '_'.join(file.split('/')[-1].split('.')[0].split(' '))
                    if genome_name in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                        shutil.move(os.path.join(output_dir_path, 'effectors_tmp', file),
                                    os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'effectors.fasta'))
                        error_msg = verify_effectors_f(
                            os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'effectors.fasta'),
                            os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'ORFs.fasta'))
                        if error_msg:
                            error_msg = f'In {genome_name}: {error_msg}'
                            fail(error_msg, error_path)
        else:
            error_msg = verify_fasta_format(effectors_path, 'DNA', 'input effectors')
            if error_msg:
                error_msg = f'Illegal effectors file: {error_msg}'
                fail(error_msg, error_path)
            error_msg = verify_effectors_f(effectors_path, ORFs_path)
            if error_msg:
                fail(error_msg, error_path)
            error_msg = validate_set(effectors_path, 'Your effectors file')
            if error_msg:
                fail(error_msg, error_path)
    if input_T3Es_path:
        error_msg = verify_fasta_format(input_T3Es_path, 'protein', 'effectors for homology search')
        if error_msg:
            fail(error_msg, error_path)
        error_msg = validate_set(input_T3Es_path, 'Your effectors file for homology search')
        if error_msg:
            fail(error_msg, error_path)
        os.makedirs(f'{output_dir_path}/blast_data', exist_ok=True)
        for f in os.listdir(blast_datasets_dir):
            if f.endswith('.faa'):
                shutil.copy(os.path.join(blast_datasets_dir, f), os.path.join(output_dir_path, 'blast_data'))
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
                shutil.copy(input_T3Es_path, f'{output_dir_path}/Effectidor_runs/{genome}')
    if host_proteome:
        error_msg = verify_zip(host_proteome, 'Host data')
        logger.info(f'called verify_zip {host_proteome} error_msg = {error_msg}')
        if error_msg:
            fail(error_msg, error_path)
        os.makedirs(f'{output_dir_path}/blast_data/temp_extract', exist_ok=True)
        logger.info('creating temp_extract')
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
        logger.info('created host.faa and removed temp_extract')

    if no_T3SS_path:
        logger.info(f'calling verify_zip {no_T3SS_path}')
        error_msg = verify_zip(no_T3SS_path, 'Close bacteria without T3SS data')
        if error_msg:
            fail(error_msg, error_path)
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
            os.makedirs(os.path.join(output_dir_path, 'Effectidor_runs', genome, 'blast_data'), exist_ok=True)
            for f in os.listdir(f'{output_dir_path}/blast_data'):
                logger.info(f'copying {f} to Effectidor_runs')
                shutil.copy(f'{output_dir_path}/blast_data/{f}',
                            os.path.join(output_dir_path, 'Effectidor_runs', genome, 'blast_data'))
            # shutil.copytree(f'{output_dir_path}/blast_data',f'{output_dir_path}/Effectidor_runs/{genome}',dirs_exist_ok=True)
    logger.info('finished validate_input')


def main(ORFs_path, output_dir_path, effectors_path, input_T3Es_path, host_proteome, genome_path,
         gff_path, no_T3SS, OG_input_f, PIP=False, hrp=False, mxiE=False, exs=False, tts=False, homology_search=False,
         signal=False, effectors_coverage='50', cpu='1'):

    error_path = f'{output_dir_path}/error.txt'
    try:

        os.makedirs(output_dir_path, exist_ok=True)

        tmp_dir = f'{output_dir_path}/tmp'
        os.makedirs(tmp_dir, exist_ok=True)

        validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path,
                       gff_path, no_T3SS, error_path, OG_input_f)

        if os.path.exists(f'{output_dir_path}/Effectidor_runs'):
            cmds = []
            for genome in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                genome_ORFs_path = os.path.join(f'{output_dir_path}/Effectidor_runs', genome, 'ORFs.fasta')
                genome_output_path = os.path.join(f'{output_dir_path}/Effectidor_runs', genome)
                genome_error = os.path.join(f'{output_dir_path}/Effectidor_runs', genome, 'error.txt')
                parameters = f'{genome_error} {genome_ORFs_path} {genome_output_path}'
                input_effectors_path = os.path.join(output_dir_path, 'Effectidor_runs', genome, 'effectors.fasta')
                if os.path.exists(input_effectors_path):
                    parameters += f' --input_effectors_path {input_effectors_path}'
                if gff_path:
                    gff_file = os.path.join(output_dir_path, 'Effectidor_runs', genome, 'genome.gff3')
                    parameters += f' --gff_file {gff_file}'
                    # parameters += ' --mobile_genetics_elements'
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
                parameters += f' --coverage {effectors_coverage}'

                cmd_params = f'python {os.path.join(scripts_dir, "T3Es_wrapper.py")} {parameters}'
                cmds.append(cmd_params)
            with ProcessPoolExecutor(max_workers=int(cpu)) as executor:
                for i in range(0, len(cmds)):
                    executor.submit(subprocess.check_output, cmds[i], shell=True)

        else:
            effectors_learn(error_path, ORFs_path, effectors_path, output_dir_path, tmp_dir, gff_path, genome_path,
                            PIP=PIP, hrp=hrp, mxiE=mxiE, exs=exs, tts=tts, homology_search=homology_search,
                            signal=signal, coverage=effectors_coverage)

        # add a check for failed features jobs...
        while not os.path.exists(os.path.join(output_dir_path, 'clean_orthologs_table.csv')):
            # make sure the find_OGs job was finished before proceeding
            sleep(60)
        subprocess.check_output(['python', os.path.join(scripts_dir, 'merge_features_for_OGs.py'), output_dir_path,
                                 OG_input_f])

        annotations_df = pd.read_csv(f'{output_dir_path}/OGs_annotations.csv')
        features_data = pd.read_csv(f'{output_dir_path}/OGs_features.csv')
        try:
            positives_size = features_data['is_effector'].value_counts()['effector']
        except KeyError:
            positives_size = 0
        if positives_size == 0:
            error_msg = 'No effectors were found in your data based on homology! If you know your data should contain '\
                        'type III effectors please supply them. '
            fail(error_msg, error_path)
        elif positives_size < 3:
            effectors = features_data[features_data['is_effector'] == 'effector']
            homologs = pd.read_csv(f'{output_dir_path}/OG_effector_homologs.csv')
            positives_table = annotations_df[annotations_df['OG'].isin(effectors['OG'])].merge(homologs, how='left').to_html(
                index=False, justify='left', escape=False)
            predicted_table = ''

        # learning step
        else:
            subprocess.check_output(
                ['python', os.path.join(scripts_dir, 'learning.py'), output_dir_path, 'OGs_features.csv'])
            preds_df = pd.read_csv(f'{output_dir_path}/out_learning/consensus_predictions.csv')
            ortho_df = pd.read_csv(f'{output_dir_path}/clean_orthologs_table_with_pseudo.csv')
            annotated_preds = preds_df.merge(annotations_df)
            annotated_preds.to_csv(f'{output_dir_path}/out_learning/consensus_predictions_with_annotations.csv',
                                   index=False)
            annotated_preds.merge(ortho_df).to_csv(
                f'{output_dir_path}/out_learning/consensus_predictions_with_annotations_and_ortho_table.csv',
                index=False)
            convert_csv_to_colored_xlsx(
                f'{output_dir_path}/out_learning/consensus_predictions_with_annotations_and_ortho_table.csv')

            predicted_table, positives_table = make_html_tables(
                f'{output_dir_path}/out_learning/consensus_predictions_with_annotations.csv',
                f'{output_dir_path}/OG_effector_homologs.csv')
            if os.path.exists(f'{output_dir_path}/Effectidor_runs'):
                num_genomes = len(os.listdir(f'{output_dir_path}/Effectidor_runs'))
                if num_genomes > 1:
                    subprocess.check_output(['python', os.path.join(scripts_dir, 'phyletic_patterns.py'), output_dir_path])
        low_quality_flag = False
        if os.path.exists(f'{output_dir_path}/out_learning/learning_failed.txt'):
            low_quality_flag = True

        T3SS_data = pd.read_csv(os.path.join(output_dir_path, 'T3SS.csv'), dtype=str)
        T3SS_Table = T3SS_data.dropna()
        T3SS_table = T3SS_Table.to_html(index=False, justify='left', escape=False)

        with open(f'{output_dir_path}/output.html', 'w') as out:
            out.write(f'positives:\n{positives_table}\n<br>\npredicted:\n{predicted_table}<br>T3SS and flagella components:<br>{T3SS_table}')

    except Exception as e:
        logger.info(f'SUCCEEDED = False')
        logger.info(e)
        logger.info(f"""ORFs_path: {ORFs_path}\noutput_dir_path: {output_dir_path}\n
                    effectors_path:{effectors_path}\nhost_proteome:{host_proteome}""")


if __name__ == '__main__':
    from sys import argv

    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input_ORFs_path',
                        help='A path to a DNA ORFs file. Can be either a fasta file or zip holding multiple fasta files.',
                        type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('output_dir_path',
                        help='A path to a folder in which the output files will be created.',
                        type=lambda path: path.rstrip('/'))
    parser.add_argument('OGs_table_path', help='A path to the OGs table - csv file.',
                        type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('--input_effectors_path', default='',
                        help='A path to a DNA fasta with positive samples. '
                             'All samples in this file should be in the input ORFs file as well.')
    parser.add_argument('--input_T3Es_path', default='',
                        help='A path to protein fasta with T3Es records of other bacteria.')
    parser.add_argument('--host_proteome_path', default='',
                        help='A path to a zip archive with protein fasta files of host proteome.')
    parser.add_argument('--no_T3SS', default='',
                        help='A path to a zip archive with protein fasta files with related bacteria non T3SS '
                             'proteomes.')
    parser.add_argument('--genome_path', default='',
                        help='A path to a fasta file with full genome record. Can be either a fasta file or zip holding multiple fasta files.')
    parser.add_argument('--gff_path', default='',
                        help='A path to a GFF file. Can be either a GFF file or zip holding multiple GFF files.')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

    parser.add_argument('--homology_search', help='search additional effectors based on homology to internal dataset',
                        action='store_true')
    parser.add_argument('--translocation_signal', help='extract translocation signal feature', action='store_true')
    parser.add_argument('--PIP', help='look for PIP-box in promoters', action='store_true')
    parser.add_argument('--hrp', help='look for hrp-box in promoters', action='store_true')
    parser.add_argument('--mxiE', help='look for mxiE-box in promoters', action='store_true')
    parser.add_argument('--exs', help='look for exs-box in promoters', action='store_true')
    parser.add_argument('--tts', help='look for tts-box in promoters', action='store_true')
    parser.add_argument('--effectors_coverage', help='coverage percentage cutoff for effectors homologs', default='50')
    parser.add_argument('--cpu', help='maximal number of CPUs to use', default='1')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    PIP_flag = args.PIP
    hrp_flag = args.hrp
    mxiE_flag = args.mxiE
    exs_flag = args.exs
    tts_flag = args.tts

    main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path,
         args.host_proteome_path, args.genome_path, args.gff_path, args.no_T3SS, args.OGs_table_path,
         PIP=PIP_flag, hrp=hrp_flag, mxiE=mxiE_flag, exs=exs_flag, tts=tts_flag, homology_search=args.homology_search,
         signal=args.translocation_signal, effectors_coverage=args.effectors_coverage, cpu=args.cpu)
