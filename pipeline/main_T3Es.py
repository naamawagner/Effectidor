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

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def verify_fasta_format(fasta_path,Type,input_name):
    logger.info('Validating FASTA format')
    flag = False
    if Type == 'DNA':
        legal_chars = set(Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() + Bio.SeqUtils.IUPACData.ambiguous_dna_letters)
    else: #Type == 'protein'
        legal_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
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
                if putative_end_of_file:  # non empty line after empty line
                    return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {putative_end_of_file} in fasta {input_name} is empty.'
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
                    AGCT_count = seq.count('A')+seq.count('G')+seq.count('T')+seq.count('C')
                    if AGCT_count >= 0.5*len(seq):
                        return f'Protein fasta {input_name} seems to contain DNA records. Make sure all samples in the file are protein sequences and re-submit.'
        except UnicodeDecodeError as e:
            logger.info(e.args)
            line_number += 1  # the line that was failed to be read
            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in fasta contains one (or more) non <a href="https://en.wikipedia.org/wiki/ASCII" target="_blank">ascii</a> character(s).'
    # override the old file with the curated content
    with open(fasta_path, 'w') as f:
        f.write(curated_content)

def verify_effectors_f(effectors_path, ORFs_path):
    effectors_recs = [rec.id for rec in SeqIO.parse(effectors_path,'fasta')]
    effectors_set = set([rec.id for rec in SeqIO.parse(effectors_path,'fasta')])
    ORFs_recs = set([rec.id for rec in SeqIO.parse(ORFs_path,'fasta')])
    if not effectors_set.issubset(ORFs_recs):
        not_in_ORFs = ','.join([rec for rec in effectors_set.difference(ORFs_recs)])
        return f'Illegal effectors records. The following records IDs are in the effectors file and not in the ORFs file:\n{not_in_ORFs}.'
    if len(effectors_set) != len(effectors_recs):
        more_than_once = ','.join([effector for effector in effectors_set if effectors_recs.count(effector)>1])
        return f'Illegal effectors records. The following records IDs appear more than once in the file:\n{more_than_once}.'
    
def verify_genome_one_contig(genome_path, file_name):
    recs = list(SeqIO.parse(genome_path,'fasta'))
    if len(recs) > 1: # it must be one contig
        return f'Illegal genome file of {file_name}. The FASTA file contains more than one record (more than one contig).'
    
def verify_zip(file,name):
    if not file.endswith('.zip'):
        return f'{name} not in zip format. Make sure to upload a zip archive and resubmit your job.'
    

def validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path, gff_path, no_T3SS_path, error_path):
    logger.info('Validating input...')
    os.makedirs(f'{output_dir_path}/contigs_ORFs')
    if ORFs_path.endswith('.zip'):     
        shutil.unpack_archive(f'{output_dir_path}/ORFs.zip',f'{output_dir_path}/contigs_ORFs')
        ORFs=[]
        for file in os.listdir(f'{output_dir_path}/contigs_ORFs'):
            error_msg = verify_fasta_format(f'{output_dir_path}/contigs_ORFs/{file}','DNA',f'{file} in ORFs archive')
            if error_msg:
                error_msg = f'Illegal fasta files in {file} in ORFs archive: {error_msg}'
                fail(error_msg,error_path)
            recs = SeqIO.parse(f'{output_dir_path}/contigs_ORFs/{file}','fasta')
            for rec in recs:
                ORFs.append(rec)
        SeqIO.write(ORFs,f'{output_dir_path}/ORFs.fasta','fasta')
    else:
        shutil.copy(ORFs_path,f'{output_dir_path}/contigs_ORFs')
        error_msg = verify_fasta_format(ORFs_path,'DNA','input ORFs')
        if error_msg:
            error_msg = f'Illegal fasta file in ORFs file: {error_msg}'
            fail(error_msg,error_path)
    if effectors_path:
        error_msg = verify_fasta_format(effectors_path,'DNA', 'input effectors')
        if error_msg:
            error_msg = f'Illegal effectors file: {error_msg}'
            fail(error_msg, error_path)
        error_msg = verify_effectors_f(effectors_path,f'{output_dir_path}/ORFs.fasta')
        if error_msg:
            fail(error_msg, error_path)
    if input_T3Es_path:
        error_msg = verify_fasta_format(input_T3Es_path,'protein', 'effectors for homology search')
    if host_proteome:
        error_msg = verify_zip(host_proteome,'Host')
        if error_msg:
            fail(error_msg,error_path)
    if genome_path and gff_path:
        # genome
        os.makedirs(f'{output_dir_path}/full_genome')
        shutil.unpack_archive(f'{output_dir_path}/genome_sequence.zip',f'{output_dir_path}/full_genome')
        genome_files_names = []
        for file in os.listdir(f'{output_dir_path}/full_genome'):
            if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{output_dir_path}/full_genome/{file}'): # discard system files and directories
                error_msg = verify_fasta_format(f'{output_dir_path}/full_genome/{file}','DNA', f'{file} in full genome archive')
                if error_msg:
                    error_msg = f'Illegal fasta files in {file} in full genome archive: {error_msg}'
                    fail(error_msg,error_path)
                error_msg = verify_genome_one_contig(f'{output_dir_path}/full_genome/{file}',file)
                if error_msg:
                    fail(error_msg,error_path)
                genome_files_names.append(file.split('.')[0])
        # gff
        os.makedirs(f'{output_dir_path}/gff')
        shutil.unpack_archive(f'{output_dir_path}/genome_features.zip',f'{output_dir_path}/gff')
        gff_files_names = [file.split('.')[0] for file in os.listdir(f'{output_dir_path}/gff') if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{output_dir_path}/gff/{file}')]
        if set(gff_files_names) != set(genome_files_names):
            in_gff_not_genome = set(gff_files_names).difference(set(genome_files_names))
            in_genome_not_in_gff = set(genome_files_names).difference(set(gff_files_names))
            error_msg = 'The names of the files in the gff archive are not matching the names of the files in the full genome archive!'
            if len(in_gff_not_genome) > 0:
                error_msg += f' The following file names are present in the GFF archive but not in the full genome archive: {in_gff_not_genome}'
            if len(in_genome_not_in_gff) > 0:
                error_msg += f' The following file names are present in the full genome archive but not in the GFF archive: {in_genome_not_in_gff}'
            fail(error_msg,error_path)
    if no_T3SS_path:
        error_msg = verify_zip(no_T3SS_path,'no_T3SS')
        if error_msg:
            fail(error_msg,error_path)


def main(ORFs_path, output_dir_path, effectors_path, input_T3Es_path, host_proteome, html_path, queue, genome_path, gff_path, no_T3SS, full_genome=False, PIP=False, hrp=False, mxiE=False, exs=False, tts=False):

    error_path = f'{output_dir_path}/error.txt'
    try:
        if html_path:
            run_number = initialize_html(CONSTS, output_dir_path, html_path)
            #final_zip_path = f'{os.path.split(output_dir_path)[0]}/{CONSTS.WEBSERVER_NAME}_{run_number}'
    
        os.makedirs(output_dir_path, exist_ok=True)
    
        tmp_dir = f'{os.path.split(ORFs_path)[0]}/tmp'  # same folder as the input files
        os.makedirs(tmp_dir, exist_ok=True)
        
        validate_input(output_dir_path, ORFs_path, effectors_path, input_T3Es_path, host_proteome, genome_path, gff_path, no_T3SS, error_path)
        if full_genome:
            if genome_path and gff_path:
                predicted_table, positives_table = effectors_learn(error_path, f'{output_dir_path}/ORFs.fasta', effectors_path, output_dir_path, tmp_dir, queue, organization=True, CIS_elements=True, PIP=PIP, hrp=hrp, mxiE=mxiE, exs=exs, tts=tts)
            else:
                predicted_table, positives_table = effectors_learn(error_path, f'{output_dir_path}/ORFs.fasta', effectors_path, output_dir_path, tmp_dir, queue, organization=True)
        else:
            predicted_table, positives_table = effectors_learn(error_path, f'{output_dir_path}/ORFs.fasta', effectors_path, output_dir_path, tmp_dir, queue)
    
        if html_path:
            #shutil.make_archive(final_zip_path, 'zip', output_dir_path)
            finalize_html(html_path, error_path, run_number, predicted_table, positives_table)

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


def finalize_html(html_path, error_path, run_number, predicted_table, positives_table):
    succeeded = not os.path.exists(error_path)
    logger.info(f'SUCCEEDED = {succeeded}')
    if succeeded:
        edit_success_html(CONSTS, html_path, run_number, predicted_table, positives_table)
    else:
        edit_failure_html(CONSTS, error_path, html_path, run_number)
    add_closing_html_tags(html_path, CONSTS, run_number)


def edit_success_html(CONSTS, html_path, run_number, predicted_table, positives_table):
    update_html(html_path, 'RUNNING', 'FINISHED')
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
                <h3><b>Top 10 predictions, umong unlabeled samples</b></h3>
                {predicted_table}
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
        if args.full_genome:
            PIP_flag = args.PIP
            hrp_flag = args.hrp
            mxiE_flag = args.mxiE
            exs_flag = args.exs
            tts_flag = args.tts
            main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS, full_genome=True, PIP=PIP_flag, hrp=hrp_flag, mxiE=mxiE_flag, exs=exs_flag, tts=tts_flag)
        else:
            main(args.input_ORFs_path, args.output_dir_path, args.input_effectors_path, args.input_T3Es_path, args.host_proteome_path, args.html_path, args.queue_name, args.genome_path, args.gff_path, args.no_T3SS)