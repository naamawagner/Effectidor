import sys
import os
import logging
import Bio.SeqUtils
import effectidor_CONSTANTS as CONSTS  # from /effectidor/auxiliaries
from time import sleep, time
sys.path.append('/lsweb/josef_sites/effectidor/auxiliaries')
from auxiliaries import fail, update_html, append_to_html  # from /effectidor/auxiliaries
from Bio import SeqIO
from T3Es_wrapper import effectors_learn
import shutil
import re
import subprocess
import pandas as pd
from add_annotations_to_predictions import make_html_tables
from csv_to_colored_xlsx_converter import convert_csv_to_colored_xlsx

data_dir = CONSTS.EFFECTIDOR_DATA
blast_datasets_dir = f'{data_dir}/blast_data'
scripts_dir = CONSTS.EFFECTIDOR_EXEC

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

ILLEGAL_CHARS = ':;,\'\"'


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
        if int(float(file_size)) >= 11*(10**6):
            return f'{input_name} is too big! It should contain a single bacterial data. For pan-genome analysis, data of separate genomes should be provided in separate files (see instructions)'
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
                    line = re.sub(f'[{ILLEGAL_CHARS}]', '_', line)
                    # return f'Illegal format. First line in fasta {input_name} contains an illegal character in its ' \
                    #        f'first word (one of: {ILLEGAL_CHARS}).'
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
                        if Type == 'DNA':
                            if has_illegal_chars(line):
                                line = re.sub(f'[{ILLEGAL_CHARS}]', '_', line)
                                # return f'Illegal format. Line {line_number} in fasta {input_name} contains an illegal' \
                                #        f' character in its first word (one of: {ILLEGAL_CHARS}).'
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
                        return f'FASTA file {os.path.basename(fasta_path)} in {input_name} input contains empty records!'
                    AGCT_count = seq.count('A') + seq.count('G') + seq.count('T') + seq.count('C') + seq.count('N')
                    if AGCT_count >= 0.95 * len(seq):
                        f = AGCT_count * 100 / len(seq)
                        return f'Protein fasta {input_name} seems to contain DNA records (record {rec.id} contains ' \
                               f'{float("%.2f" % f)}% DNA characters). Make sure all samples in the file are protein ' \
                               f'sequences and re-submit. '
        except UnicodeDecodeError as e:
            logger.info(e.args)
            line_number += 1  # the line that was failed to be read
            return f'In {input_name}: Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA ' \
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
    start_code = ['ATG', 'GTG', 'TTG', 'atg', 'gtg', 'ttg']
    end_code = ['TAA', 'TAG', 'TGA', 'taa', 'tag', 'tga']
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
    if len(effectors_set) > 150:
        return f'Illegal effectors records. There are {str(len(effectors_set))} records in this input. This is more than can be encoded in a single bacterial genome.'
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
            subprocess.check_output(f'rm {"/".join(file.split("/")[:-1])}/zip_tmp/{f}', shell=True)
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


def validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path, gff_path,
                   no_T3SS_path, error_path):
    logger.info('Validating input...')
    if ORFs_path.endswith('.zip'):
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
                genome_name = re.sub(r'\W', '', rf'{genome_name}')
                ORFs_genomes.append(genome_name)
                os.makedirs(os.path.join(output_dir_path, 'Effectidor_runs', genome_name), exist_ok=True)
                shutil.move(os.path.join(output_dir_path, 'ORFs_tmp', file),
                            os.path.join(output_dir_path, 'Effectidor_runs', genome_name, 'ORFs.fasta'))
                number_of_genomes += 1
        if number_of_genomes == 0:
            error_msg = 'No files were found in the ORFs input! Make sure the files in the ZIP archive are not inside ' \
                        'directories. '
            fail(error_msg, error_path)
        elif number_of_genomes > 1000:
            error_msg = f"There are {str(number_of_genomes)} genomes in your input! Effectidor's limit is 1,000 genomes."
            fail(error_msg, error_path)
        shutil.rmtree(f'{output_dir_path}/ORFs_tmp', ignore_errors=True)
        ORFs_set = set(ORFs_genomes)
    else:
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
                    genome_name = re.sub(r'\W', '', rf'{genome_name}')
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
                error_msg = 'No files were found in the GFF input! Make sure the files in the ZIP archive are not ' \
                            'inside directories. '
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
                        genome_name = re.sub(r'\W', '', rf'{genome_name}')
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
                    error_msg = 'No files were found in the full genome input! Make sure the files in the ZIP archive ' \
                                'are not inside directories. '
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
                    genome_name = re.sub(r'\W', '', rf'{genome_name}')
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
        cmd = f'cp {blast_datasets_dir}/*.faa {output_dir_path}/blast_data/'
        subprocess.check_output(cmd, shell=True)
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


def cleanup_is_running(queues=('pupkolab', 'pupkoweb')):
    for q in queues:
        try:
            if subprocess.check_output(f'qstat {q} | grep cleanup_effec', shell=True):
                # a cleanup is currently running.
                return True
        except:
            pass
    # none of the queues is running cleanup (none of them returned True)
    return False


def cleanup_ran_today(path=r'/bioseq/data/results/effectidor/'):
    for f in os.listdir(path):
        if f.endswith('.ER'):
            f_path = os.path.join(path, f)
            ctime = os.stat(f_path).st_ctime
            if time() - ctime < 60 * 60 * 24:
                return True
    # cleanup did not run today
    return False


def main(ORFs_path, output_dir_path, effectors_path, input_T3Es_path, host_proteome, html_path, queue, genome_path,
         gff_path, no_T3SS, identity_cutoff='50', coverage_cutoff='60', PIP=False, hrp=False, mxiE=False, exs=False,
         tts=False, homology_search=False, signal=False, signalp=False, MGE=True, effectors_coverage='50'):
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
            # final_zip_path = f'{os.path.split(output_dir_path)[0]}/{CONSTS.WEBSERVER_NAME}_{run_number}'

        os.makedirs(output_dir_path, exist_ok=True)

        tmp_dir = f'{output_dir_path}/tmp'
        os.makedirs(tmp_dir, exist_ok=True)

        validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path,
                       gff_path, no_T3SS, error_path)

        find_OGs_cmd = f'module load MMseqs2/May2024;!@#python ' \
                       f'{os.path.join(scripts_dir, "find_OGs_in_genomes.py")} ' \
                       f'{output_dir_path} {identity_cutoff} {coverage_cutoff}\tfind_OGs_effectidor\n'
        with open(os.path.join(output_dir_path, 'find_OGs.cmds'), 'w') as out_f:
            out_f.write(find_OGs_cmd)
        cmd = f'python {os.path.join(scripts_dir,"auxiliaries", "q_submitter_power.py")} {os.path.join(output_dir_path, "find_OGs.cmds")} ' \
              f'{output_dir_path} -q {queue} '
        OGs_job_number = subprocess.check_output(cmd, shell=True).decode('ascii').strip()
        # make sure it was submitted successfully before proceeding
        # TO COMPLETE

        if os.path.exists(f'{output_dir_path}/Effectidor_runs'):
            currently_running = []
            for genome in os.listdir(f'{output_dir_path}/Effectidor_runs'):
                genome_ORFs_path = os.path.join(f'{output_dir_path}/Effectidor_runs', genome, 'ORFs.fasta')
                genome_output_path = os.path.join(f'{output_dir_path}/Effectidor_runs', genome)
                genome_error = os.path.join(f'{output_dir_path}/Effectidor_runs', genome, 'error.txt')
                parameters = f'{genome_error} {genome_ORFs_path} {genome_output_path} --queue {queue}'
                input_effectors_path = os.path.join(output_dir_path, 'Effectidor_runs', genome, 'effectors.fasta')
                if os.path.exists(input_effectors_path):
                    parameters += f' --input_effectors_path {input_effectors_path}'
                if gff_path:
                    gff_file = os.path.join(output_dir_path, 'Effectidor_runs', genome, 'genome.gff3')
                    parameters += f' --gff_file {gff_file}'
                    if MGE:
                        parameters += ' --mobile_genetics_elements'
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
                if signalp:
                    parameters += ' --signalp'
                parameters += f' --coverage {effectors_coverage}'

                job_cmd = f'module load MMseqs2/May2024;!@#python {os.path.join(scripts_dir, "T3Es_wrapper.py")} ' \
                          f'{parameters}\tEffectidor_features_{genome}\n '
                cmds_f = os.path.join(output_dir_path, 'Effectidor_runs', genome, 'features_wrapper.cmds')
                with open(cmds_f, 'w') as job_f:
                    job_f.write(job_cmd)
                cmd = f'python {os.path.join(scripts_dir, "auxiliaries", "q_submitter_power.py")} {cmds_f} ' \
                      f'{os.path.join(output_dir_path, "Effectidor_runs", genome)} -q {queue}'
                subprocess.check_output(cmd, shell=True)
                # make sure it was submitted successfully before proceeding
                # TO COMPLETE
                currently_running.append(genome)
                while len(currently_running) > 4:  # features can be extracted for up to 5 genomes at a time,
                    # to avoid overload the cluster.
                    sleep(30)
                    for Genome in currently_running:
                        if os.path.exists(os.path.join(f'{output_dir_path}/Effectidor_runs', Genome, 'error.txt')):
                            with open(os.path.join(f'{output_dir_path}/Effectidor_runs', Genome, 'error.txt')) as error:
                                msg = error.read()
                            error_msg = f'Oops :(\n In genome {Genome} {msg}'
                            fail(error_msg, error_path)
                    if os.path.exists(os.path.join(output_dir_path, 'error_OGs.txt')):
                        with open(os.path.join(output_dir_path, 'error_OGs.txt')) as f:
                            error_msg = f.read()
                        fail(f'Error in OG calculation: {error_msg}', error_path)
                    currently_running = [genome for genome in currently_running if not
                                        os.path.exists(os.path.join(output_dir_path, "Effectidor_runs",
                                        genome, 'features.csv'))]
            while len(currently_running) > 0:  # wait until they all finish before proceeding with the next step.
                sleep(30)
                currently_running = [genome for genome in currently_running if not os.path.exists(
                    os.path.join(output_dir_path, "Effectidor_runs", genome, 'features.csv'))]

        else:
            effectors_learn(error_path, ORFs_path, effectors_path, output_dir_path, tmp_dir, queue, gff_path,
                            genome_path, PIP=PIP, hrp=hrp, mxiE=mxiE, exs=exs, tts=tts,
                            homology_search=homology_search, signal=signal, signalp=signalp, MGE=MGE,
                            coverage=effectors_coverage)
        while not os.path.exists(os.path.join(output_dir_path, 'clean_orthologs_table.csv')):
            # make sure it didn't fail
            if os.path.exists(os.path.join(output_dir_path, 'error_OGs.txt')):
                with open(os.path.join(output_dir_path, 'error_OGs.txt')) as f:
                    error_msg = f.read()
                fail(f'Error in OG calculation: {error_msg}', error_path)
            # make sure the find_OGs job was finished before proceeding
            sleep(60)
        merge_OGs_cmd = f'module load MMseqs2/May2024;!@#python ' \
                       f'{os.path.join(scripts_dir, "merge_features_for_OGs.py")} ' \
                       f'{output_dir_path}\tmerge_OGs_effectidor\n'
        with open(os.path.join(output_dir_path, 'merge_OGs.cmds'), 'w') as out_f:
            out_f.write(merge_OGs_cmd)
        cmd = f'python {os.path.join(scripts_dir, "auxiliaries", "q_submitter_power.py")} {os.path.join(output_dir_path, "merge_OGs.cmds")} ' \
              f'{output_dir_path} -q {queue} --memory 20'
        subprocess.check_output(cmd, shell=True)
        while not os.path.exists(os.path.join(output_dir_path, 'merge_OGs.done')):
            sleep(30)
        # I put it in a job as for large runs it requires more memory
        #subprocess.check_output(['python', os.path.join(scripts_dir, 'merge_features_for_OGs.py'), output_dir_path])

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
                    subprocess.check_output(['python', os.path.join(scripts_dir, 'phyletic_patterns.py'),
                                             output_dir_path])
        low_quality_flag = False
        if os.path.exists(f'{output_dir_path}/out_learning/learning_failed.txt'):
            low_quality_flag = True

        T3SS_data = pd.read_csv(os.path.join(output_dir_path, 'T3SS.csv'), dtype=str)
        T3SS_Table = T3SS_data.dropna()
        T3SS_table = T3SS_Table.to_html(index=False, justify='left', escape=False)
        if html_path:
            # shutil.make_archive(final_zip_path, 'zip', output_dir_path)
            finalize_html(html_path, error_path, run_number, predicted_table, positives_table, T3SS_table,
                          low_quality_flag, output_dir_path)
        else:
            with open(f'{output_dir_path}/output.html', 'w') as out:
                out.write(f'positives:\n{positives_table}\n<br>\npredicted:\n{predicted_table}<br>T3SS and flagella components:<br>{T3SS_table}')

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


def finalize_html(html_path, error_path, run_number, predicted_table, positives_table, T3SS_table, low_confidence_flag,
                  output_dir_path):
    succeeded = not os.path.exists(error_path)
    logger.info(f'SUCCEEDED = {succeeded}')
    f = open(os.path.join(output_dir_path, f'effectidor_{run_number}.END_OK'), "w")
    f.close()
    if succeeded:
        edit_success_html(CONSTS, html_path, predicted_table, positives_table, T3SS_table,
                          low_confidence_flag, output_dir_path)
    else:
        edit_failure_html(CONSTS, error_path, html_path, run_number)
    add_closing_html_tags(html_path, CONSTS, run_number)


def edit_success_html(CONSTS, html_path, predicted_table, positives_table, T3SS_table,
                      low_confidence_flag, output_dir_path):
    update_html(html_path, 'RUNNING', 'FINISHED')
    if low_confidence_flag:
        append_to_html(html_path, f'''
                       <div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify"><h3>
                       <font color="red">WARNING: The predictions might be of low quality.
                       </font></h3><br>
                       </div>
                       ''')
    if predicted_table:
        append_to_html(html_path, f'''
                <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='left'>
                <a href='out_learning/consensus_predictions_with_annotations_and_ortho_table.xlsx' target='_blank'>Download predictions file</a>
                <br>
                <a href='T3SS.csv' target='_blank'>Download T3SS and flagella components' details</a>
                <br>
                <a href='chaperones.csv' target='_blank'>Download chaperones details</a>
                <br>
                <a href='out_learning/feature_importance.csv' target='_blank'>Download feature importance file</a>
                <br>
                <a href='OGs_features.csv' target='_blank'>Download features file</a>
                <br><br>
                <h3><b>Positive samples that were used to train the model</b></h3>
                {positives_table}
                <br>
                <h3><b>Top 10 predictions, among unlabeled samples</b></h3>
                {predicted_table}
                </div>
                <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='center'>
                ''')
        if os.path.exists(f'{output_dir_path}/Effectidor_runs'):
            num_genomes = len(os.listdir(f'{output_dir_path}/Effectidor_runs'))
            if num_genomes > 1:
                append_to_html(html_path, f'''
                    <br>
                    <h3><b>T3Es presence/absence map</b></h3>
                    <br>
                    <a href='T3Es_presence_absence.png'><img src='T3Es_presence_absence.png' style="max-width: 100%;" /></a>
                    <br>
                    <h3><b>T3SS and flagella components presence/absence map</b></h3>
                    <a href='T3SS_presence_absence.png'><img src='T3SS_presence_absence.png' style="max-width: 100%;" /></a>
                    ''')
                if os.path.exists('chaperones_presence_absence.png'):
                    append_to_html(html_path, f'''
                    <br>
                    <h3><b>chaperones presence/absence map</b></h3>
                    <a href='chaperones_presence_absence.png'><img src='chaperones_presence_absence.png' style="max-width: 100%;" /></a>
                    ''')
        else:
            append_to_html(html_path, f'''
            <br>
            <h3><b>T3SS and flagella components</b></h3>
            {T3SS_table}
            ''')
        append_to_html(html_path, f'''
                <br>
                <h3><b>feature importance</b></h3>
                <br>
                <a href='out_learning/feature_importance.png'><img src='out_learning/feature_importance.png' style="max-width: 100%;" ></a>
                <br><br>
                <h3><b>best features comparison - effectors vs non-effectors:</b></h3>
                <br>
                <a href='out_learning/0.png'><img src='out_learning/0.png' style="max-width: 100%;" ></a>
                <a href='out_learning/1.png'><img src='out_learning/1.png' style="max-width: 100%;" ></a>
                <a href='out_learning/2.png'><img src='out_learning/2.png' style="max-width: 100%;" ></a>
                <a href='out_learning/3.png'><img src='out_learning/3.png' style="max-width: 100%;" ></a>
                <a href='out_learning/4.png'><img src='out_learning/4.png' style="max-width: 100%;" ></a>
                <a href='out_learning/5.png'><img src='out_learning/5.png' style="max-width: 100%;" ></a>
                <a href='out_learning/6.png'><img src='out_learning/6.png' style="max-width: 100%;" ></a>
                <a href='out_learning/7.png'><img src='out_learning/7.png' style="max-width: 100%;" ></a>
                <a href='out_learning/8.png'><img src='out_learning/8.png' style="max-width: 100%;" ></a>
                <a href='out_learning/9.png'><img src='out_learning/9.png' style="max-width: 100%;" ></a>
                </div>
                ''')
    else:
        append_to_html(html_path, f'''
                <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='left'> Unfortunately, we could not 
                train a satisfying classifier due to small positive set.<br>The effectors found based on homology are 
                listed in the table bellow:<br> {positives_table} <br>
                <a href='clean_orthologs_table_with_pseudo.csv' target='_blank'>Download the OGs table</a>
                <br>
                <a href='T3SS.csv' target='_blank'>Download T3SS and flagella components' details</a>
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
    FORMER_MSG = f'Effectidor is now processing your request. This page will be automatically updated every ' \
                 f'{CONSTS.RELOAD_INTERVAL} seconds (until the job is done). You can also reload it manually. Once the ' \
                 f'job has finished, the output will appear below. '
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
    sleep(2 * CONSTS.RELOAD_INTERVAL)
    update_html(html_path, CONSTS.RELOAD_TAGS, f'<!--{CONSTS.RELOAD_TAGS}-->')  # stop refresh


def initialize_html(CONSTS, output_dir_path, html_path):
    path_tokens = output_dir_path.split('/')
    # e.g., "/bioseq/data/results/sincopa/12345678/outputs"
    run_number = path_tokens[path_tokens.index(CONSTS.WEBSERVER_NAME) + 1]

    update_html(html_path, 'QUEUED', 'RUNNING')
    # update_html(html_path, CONSTS.PROGRESS_BAR_ANCHOR, CONSTS.PROGRESS_BAR_TAG)  # add progress bar

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
                        help='A path to a zip archive with protein fasta files with related bacteria non T3SS '
                             'proteomes.')
    parser.add_argument('--genome_path', default='',
                        help='A path to a fasta file with full genome record.')
    parser.add_argument('--gff_path', default='',
                        help='A path to a GFF file.')
    parser.add_argument('--html_path', default=None,
                        help='A path to an html file that will be updated during the run.')

    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

    parser.add_argument('--homology_search', help='search additional effectors based on homology to internal dataset',
                        action='store_true')
    parser.add_argument('--translocation_signal', help='extract translocation signal feature', action='store_true')
    parser.add_argument('--signalp', help='extract SignalP6 feature', action='store_true')
    # parser.add_argument('--mobile_genetics_elements', help='extract distance from mobile genetics elements', action='store_true')
    parser.add_argument('--PIP', help='look for PIP-box in promoters', action='store_true')
    parser.add_argument('--hrp', help='look for hrp-box in promoters', action='store_true')
    parser.add_argument('--mxiE', help='look for mxiE-box in promoters', action='store_true')
    parser.add_argument('--exs', help='look for exs-box in promoters', action='store_true')
    parser.add_argument('--tts', help='look for tts-box in promoters', action='store_true')
    parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to',
                        default='power-pupko')
    parser.add_argument('--identity_cutoff', help='identity percentage cutoff for orthologs and paralogs search',
                        default='50')
    parser.add_argument('--coverage_cutoff', help='coverage percentage cutoff for orthologs and paralogs search',
                        default='60')
    parser.add_argument('--effectors_coverage', help='coverage percentage cutoff for effectors homologs', default='50')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    # if args.full_genome:
    PIP_flag = args.PIP
    hrp_flag = args.hrp
    mxiE_flag = args.mxiE
    exs_flag = args.exs
    tts_flag = args.tts

    main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path,
         args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS,
         identity_cutoff=args.identity_cutoff, coverage_cutoff=args.coverage_cutoff, PIP=PIP_flag, hrp=hrp_flag,
         mxiE=mxiE_flag, exs=exs_flag, tts=tts_flag, homology_search=args.homology_search,
         signal=args.translocation_signal, signalp=args.signalp, effectors_coverage=args.effectors_coverage)
