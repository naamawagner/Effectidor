def effectors_learn(error_path, ORFs_file, effectors_file, working_directory, tmp_dir,queue,organization=False,pip=False):
    import pandas as pd
    import subprocess
    import os
    import time
    from Bio import SeqIO
    import shutil
    from add_annotations_to_predictions import add_annotations_to_predictions,make_html_tables
    import sys
    sys.path.append('/bioseq/effectidor/auxiliaries')
    from auxiliaries import fail
    # vars
    scripts_dir = '/groups/pupko/naamawagner/T3Es_webserver/scripts/'
    os.chdir(working_directory)
    # input files
    #ORFs_file = 'ORFs.fasta'
    #effectors_file = 'effectors.fasta'
    all_prots = 'translated_ORFs.faa'
    effectors_prots = 'translated_effectors.faa'
    blast_datasets_dir = '/groups/pupko/naamawagner/T3Es_webserver/blast_data'
    
    if not os.path.exists('blast_data'):
        os.makedirs('blast_data')
    if os.path.exists('non_T3SS.zip'):
        #unzip it to the blast datasets dir
        os.makedirs(f'{working_directory}/blast_data/temp_extract')
        shutil.unpack_archive(f'{working_directory}/non_T3SS.zip',f'{working_directory}/blast_data/temp_extract')
        for file in os.listdir(f'{working_directory}/blast_data/temp_extract'):
            if not file.startswith('_') and not file.startswith('.') and os.path.isfile(f'{working_directory}/blast_data/temp_extract/{file}'):
                new_name_l = file.replace(' ','_').split('.')
                new_name = '.'.join(new_name_l[:-1])+'.'+'faa'
                os.rename(f'{working_directory}/blast_data/temp_extract/{file}',f'{working_directory}/blast_data/temp_extract/{new_name}')
                shutil.move(f'{working_directory}/blast_data/temp_extract/{new_name}',f'{working_directory}/blast_data')
    if os.path.exists(f'{working_directory}/blast_data/host.zip'):
        if not os.path.exists(f'{working_directory}/blast_data/temp_extract'):
            os.makedirs(f'{working_directory}/blast_data/temp_extract')
        shutil.unpack_archive(f'{working_directory}/blast_data/host.zip',f'{working_directory}/blast_data/temp_extract')
        recs = []
        for f in os.listdir(f'{working_directory}/blast_data/temp_extract'):
            if not f.startswith('_') and not f.startswith('.') and os.path.isfile(f'{working_directory}/blast_data/temp_extract/{f}'):
                f_recs= SeqIO.parse(f'{working_directory}/blast_data/temp_extract/{f}','fasta')
                for rec in f_recs:
                    recs.append(rec)
        SeqIO.write(recs,f'{working_directory}/blast_data/host.faa','fasta')
    cmd = f'cp {blast_datasets_dir}/*.faa ./blast_data/'
    subprocess.check_output(cmd,shell=True)
    blast_datasets_dir = './blast_data/'
    #tmp_dir = '/groups/pupko/naamawagner/T3Es_webserver/tmp_features/'
    
    # feature extraction step
    
    # translate the input fasta files
    subprocess.check_output(['python',f'{scripts_dir}/translate_fasta.py',ORFs_file,effectors_file,all_prots,effectors_prots])
    if not effectors_file:
        subprocess.check_output(['python',f'{scripts_dir}/find_effectors.py',f'{blast_datasets_dir}/T3Es.faa',all_prots,effectors_prots])
        # make sure effectors were found before proceeding!
        eff_recs = list(SeqIO.parse(effectors_prots,'fasta'))
        if len(eff_recs) == 0:
            error_msg = 'No effectors were found in the data! Make sure you run the analysis on a bacterium with an active T3SS and try to run it again with an effectors file containing all the known effectors in the bacterium.'
            fail(error_msg,error_path)
        elif len(eff_recs) < 5:
            error_msg = f'Not enough effectors were found in the data! Only {len(eff_recs)} were found in the initial homology serach. This is not enough to train a classifier.If you know more effectors are available in the bacterium, try to run it again with an effectors file containing all the known effectors in the bacterium.'
            fail(error_msg,error_path)
    # find and create non effectors fasta file
    subprocess.check_output(['python',f'{scripts_dir}/find_non_effectors.py',all_prots,effectors_prots])
    # creating the features extraction commands
    with open(f'{working_directory}/features_jobs.cmds','w') as jobs_f:
        # sequence features
        jobs_f.write(f'module load python/python-anaconda3.6.5; python {scripts_dir}/sequence_features.py {ORFs_file} {effectors_prots} {working_directory}\tsequence_features\n')
        # homology (blast) features
        jobs_f.write(f'module load python/python-anaconda3.6.5; python {scripts_dir}/homology.py {ORFs_file} {all_prots} {blast_datasets_dir} {effectors_prots} {working_directory}\thomology\n')
        if organization:
            # genome organization features
            ORFs_dir = f'{working_directory}/contigs_ORFs'
            jobs_f.write(f'module load python/python-anaconda3.6.5; python {scripts_dir}/genome_organization.py {ORFs_dir} {effectors_prots} {working_directory}\tgenome_organization\n')
            if pip:
                # PIP-box features
                gff_dir = f'{working_directory}/gff'
                genome_dir = f'{working_directory}/full_genome'
                jobs_f.write(f'module load python/python-anaconda3.6.5; python {scripts_dir}/pip_box_features.py {ORFs_file} {working_directory} {gff_dir} {genome_dir}\tgenome_organization\n')
        # signal peptide
    with open(f'{working_directory}/signalp.cmds','w') as sig_f:
        sig_f.write(f'module load python/python-anaconda3.6.5; python {scripts_dir}/signal_peptide_features.py {ORFs_file} {all_prots} {working_directory}\tsignalP\n')
    subprocess.call(f'/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py {working_directory}/features_jobs.cmds {tmp_dir} -q {queue}',shell=True)
    subprocess.call(f'/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py --cpu 3 {working_directory}/signalp.cmds {tmp_dir} -q {queue}',shell=True)
    
    x = sum([item.endswith('done') for item in os.listdir(working_directory)])
    amount_of_expected_results = 3
    if organization:
        amount_of_expected_results += 1
        if pip:
            amount_of_expected_results += 1
            
    while x < amount_of_expected_results:
        print('number of expected outputs is '+str(amount_of_expected_results)+' , number of existing outputs is '+str(x)+', waiting for additional '+str(amount_of_expected_results-x)+' outputs')
        time.sleep(300)
        num_finished_jobs = sum([item.endswith('ER') for item in os.listdir(f'{working_directory}/tmp')])
        x = sum([item.endswith('done') for item in os.listdir(working_directory)])
        if x < num_finished_jobs:
            break
        
    effectors = SeqIO.to_dict(SeqIO.parse(effectors_prots,'fasta')).keys()
    non_effectors = SeqIO.to_dict(SeqIO.parse('non_effectors.faa','fasta')).keys()
    all_locuses = SeqIO.to_dict(SeqIO.parse(all_prots,'fasta')).keys()
    label_dict = {}
    for locus in all_locuses:
        if locus in effectors:
            label_dict[locus] = [locus,'effector']
        elif locus in non_effectors:
            label_dict[locus] = [locus,'no']
        else:
            label_dict[locus] = [locus, '?']
    label_df = pd.DataFrame.from_dict(label_dict,orient='index',columns=['locus','is_effector'])
    
    files_to_merge =['physical_features.csv','homology_features.csv','signal_p_features.csv']
    done_files = ['physical_features.done','homology_features.done','signal_p_features.done']
    if organization:
        files_to_merge.append('genome_organization_features.csv')
        done_files.append('genome_organization_features.done')
        if pip:
            files_to_merge.append('pip_box_features.csv')
            done_files.append('pip_box_features.done')
    failed_jobs = [done_job[:-5] for done_job in done_files if not os.path.exists(done_job)]
    if len(failed_jobs) > 0:
        failed_str = ', '.join(failed_jobs)
        error_msg = f'The following jobs failed:\n{failed_str}'
        fail(error_msg,error_path)
    merged_df = pd.read_csv(files_to_merge[0])
    for f in files_to_merge[1:]:
        df1 = pd.read_csv(f)
        merged_df = merged_df.merge(df1)
    merged_df = merged_df.merge(label_df)
    merged_df.to_csv('features.csv',index=False)
    
    # learning step
    
    subprocess.check_output(['python',f'{scripts_dir}/learning.py'])
    if os.path.exists(r'out_learning/learning_failed.txt'):
        error_msg = 'Learning failed. It can be due to a small training set, or other reasons. For further details you can contact us.'
        fail(error_msg,error_path)
    # making final output files and tables
    
    in_f = r'out_learning/concensus_predictions.csv'
    out_f_normal = r'out_learning/concensus_predictions_with_annotation.csv'
    out_f_pseudo = r'out_learning/pseudogenes.csv'
    annotations = ORFs_file
    out_for_html_normal = r'out_learning/concensus_predictions_with_annotation_for_html.csv'
    out_for_html_pseudo = r'out_learning/pseudogenes_predictions_with_annotation_for_html.csv'
    add_annotations_to_predictions(in_f,out_f_normal,out_f_pseudo,annotations)
    from csv_to_colored_xlsx_converter import convert_csv_to_colored_xlsx
    convert_csv_to_colored_xlsx(out_f_normal)
    convert_csv_to_colored_xlsx(out_f_pseudo)
    add_annotations_to_predictions(in_f,out_for_html_normal,out_for_html_pseudo,annotations,line_end='<br>')
    predicted_table,positives_table = make_html_tables(out_for_html_normal)
    return predicted_table, positives_table

import os
if __name__ == '__main__':
        from sys import argv
        import logging
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_ORFs_path',
                            help='A path to a DNA ORFs file.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('input_effectors_path',
                            help='A path to a DNA fasta with positive samples. '
                                 'All samples in this file should be in the input ORFs file as well.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_dir_path',
                            help='A path to a folder in which the output files will be created.',
                            type=lambda path: path.rstrip('/'))
        parser.add_argument('--html_path', default=None,
                            help='A path to an html file that will be updated during the run.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))

        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        effectors_learn(args.input_ORFs_path, args.input_effectors_path, args.output_dir_path, f'{args.output_dir_path}/tmp')