def create_effectors_html(effectors_file,ORFs_file,out_dir):
    from Bio import SeqIO
    import csv
    import pandas as pd
    recs = SeqIO.parse(ORFs_file,'fasta')
    locus_annotation={}
    for rec in recs:
        is_locus = False
        annotation = ''
        header = rec.description
        if 'pseudo=true' in header:
            annotation = 'pseudogene'
        header_l = header.split('[')
        for a in header_l:
            if 'locus_tag=' in a:
                locus = a.split('=')[1].strip().strip(']')
                is_locus = True
            elif 'protein=' in a:
                if annotation == 'pseudogene':
                    annotation += ' '+a.split('=')[1].strip().strip(']')
                else:
                    annotation = a.split('=')[1].strip().strip(']')
        if is_locus:
            locus_annotation[locus] = annotation
        else:
            locus_annotation[rec.id] = annotation
    locus_protein={}
    effectors_recs = SeqIO.parse(effectors_file,'fasta')
    for effector in effectors_recs:
        protein = str(effector.seq)
        protein_l = [protein[70*i:70*(i+1)] for i in range(len(protein)//70+1)]
        if protein_l[-1]=='':
            protein_l = protein_l[:-1]
        protein_n = '<br>'.join(protein_l)
        locus_protein[effector.id] = protein_n
    with open(f'{out_dir}/effectors_for_html.csv','w',newline='') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(['Locus tag','Annotation','Protein sequence'])
        for locus in locus_protein:
            writer.writerow([locus,locus_annotation[locus],protein_n])
    data = pd.read_csv(f'{out_dir}/effectors_for_html.csv')
    effectors_table = data.to_html(index=False,justify='left',escape=False)
    return effectors_table,None,None, False
            
def write_sh_file(tmp_dir,path_to_dir_with_fastas,path_for_embedding_files,path_for_output_csv,queue):
    content = f'''#!/bin/bash\n\n#PBS -S /bin/bash\n#PBS -r y\n#PBS -q {queue}@power9\n#PBS -l ncpus=10,mem=15gb\n#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n#PBS -N Effectidor_embedding\n#PBS -e {tmp_dir}\n#PBS -o {tmp_dir}\n\nsource /groups/pupko/alburquerque/miniconda3/etc/profile.d/conda.sh\nconda activate MsaEmbedding\ncd /groups/pupko/alburquerque/MsaTransformerEffectidor/\nPYTHONPATH=$(pwd)\npython extract.py --model_location "/groups/pupko/alburquerque/MsaTransformerEffectidor/esm1_t34_670M_UR50S.pt" --input_dir "{path_to_dir_with_fastas}" --output_dir "{path_for_embedding_files}" --repr_layers 34 --include mean\npython create_effectidor_feats.py --path_to_embedding_files "{path_for_embedding_files}" --path_to_fasta_files "{path_to_dir_with_fastas}" --output_file_path "{path_for_output_csv}"'''
    with open(f'{tmp_dir}/Embedding.sh','w') as out:
        out.write(content)

def effectors_learn(error_path, ORFs_file, effectors_file, working_directory, tmp_dir,queue,organization=False,CIS_elements=False,PIP=False,hrp=False,mxiE=False,exs=False,tts=False,homology_search=False,signal=False):
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
    import effectidor_CONSTANTS
    low_quality_flag = False
    # vars
    scripts_dir = effectidor_CONSTANTS.EFFECTIDOR_EXEC
    data_dir = effectidor_CONSTANTS.EFFECTIDOR_DATA
    os.chdir(working_directory)
    log_file = f'{working_directory}/log.txt'
    # input files
    #ORFs_file = 'ORFs.fasta'
    #effectors_file = 'effectors.fasta'
    all_prots = 'translated_ORFs.faa'
    signal_prot_dir = f'{working_directory}/NTerminal_prots'
    signal_embed_dir = f'{working_directory}/NTerminal_embedding'
    effectors_prots = 'translated_effectors.faa'
    blast_datasets_dir = f'{data_dir}/blast_data'
    
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

    if os.path.exists('user_T3Es.fasta'): # add these records to the T3Es dataset 
        eff1_recs=SeqIO.parse(f'{blast_datasets_dir}/T3Es.faa','fasta')
        eff_l = list(eff1_recs)
        seqs = [rec.seq for rec in eff_l]
        eff2_recs=SeqIO.parse('user_T3Es.fasta','fasta')
        for rec in eff2_recs:
            if rec.seq not in seqs:
                eff_l.append(rec)
                seqs.append(rec.seq)
        SeqIO.write(eff_l,f'{blast_datasets_dir}/T3Es.faa','fasta')
            
    # feature extraction step
    
    # translate the input fasta files
    subprocess.check_output(['python',f'{scripts_dir}/translate_fasta.py',ORFs_file,effectors_file,all_prots,effectors_prots])
    #make N-terminal sequences for Signal search
    if signal:
        if not os.path.exists(signal_prot_dir):
            os.makedirs(signal_prot_dir)
        recs = SeqIO.parse(all_prots,'fasta')
        for rec in recs:
            ID = rec.id
            N_seq = rec.seq[:100]
            rec.seq = N_seq
            file_path = f'{signal_prot_dir}/{ID}.faa'
            SeqIO.write(rec,file_path,'fasta')
        if not os.path.exists(signal_embed_dir):
            os.makedirs(signal_embed_dir)
        write_sh_file(tmp_dir,signal_prot_dir,signal_embed_dir,f'{working_directory}/Embedding_pred.csv',queue)
    
    if not effectors_file:
        subprocess.check_output(['python',f'{scripts_dir}/find_effectors.py',f'{blast_datasets_dir}/T3Es.faa',all_prots,effectors_prots,log_file])
    elif homology_search:
        effectors_prots2 = 'homology_found_effectors.faa'
        subprocess.check_output(['python',f'{scripts_dir}/find_effectors.py',f'{blast_datasets_dir}/T3Es.faa',all_prots,effectors_prots2,log_file])
        eff1 = SeqIO.to_dict(SeqIO.parse(effectors_prots,'fasta'))
        eff2 = SeqIO.to_dict(SeqIO.parse(effectors_prots2,'fasta'))
        eff1.update(eff2)
        SeqIO.write(eff1.values(),effectors_prots,'fasta')

    # make sure effectors were found before proceeding!
    eff_recs = list(SeqIO.parse(effectors_prots,'fasta'))
    if len(eff_recs) == 0:
        error_msg = 'No effectors were found in the data! Make sure you run the analysis on a bacterium with an active T3SS and try to run it again with an effectors file containing all the known effectors in the bacterium.'
        fail(error_msg,error_path)
    elif len(eff_recs) < 3:
        return create_effectors_html(effectors_prots,ORFs_file,working_directory)
        #error_msg = f'Not enough effectors were found in the data! Only {len(eff_recs)} were found in the initial homology serach. This is not enough to train a classifier.If you know more effectors are available in the bacterium, try to run it again with an effectors file containing all the known effectors in the bacterium.'
        #fail(error_msg,error_path)
    # find and create non effectors fasta file
    subprocess.check_output(['python',f'{scripts_dir}/find_non_effectors.py',all_prots,effectors_prots])
    # creating the features extraction commands
    with open(f'{working_directory}/features_jobs.cmds','w') as jobs_f:
        # sequence features
        jobs_f.write(f'module load python/python-anaconda3.6.5; python {scripts_dir}/sequence_features.py {ORFs_file} {effectors_prots} {working_directory}\tsequence_features\n')
        # homology (blast) features
        jobs_f.write(f'module load python/python-anaconda3.6.5; python {scripts_dir}/homology.py {ORFs_file} {all_prots} {blast_datasets_dir} {effectors_prots} {working_directory}\thomology\n')
        if organization:
            gff_dir = f'{working_directory}/gff'
            # genome organization features
            ORFs_dir = f'{working_directory}/contigs_ORFs'
            jobs_f.write(f'module load python/python-anaconda3.6.5; python {scripts_dir}/genome_organization.py {ORFs_dir} {effectors_prots} {working_directory} {gff_dir}\tgenome_organization\n')
            if CIS_elements:
                # PIP-box features
                gff_dir = f'{working_directory}/gff'
                genome_dir = f'{working_directory}/full_genome'
                cmds = f'module load python/python-anaconda3.6.5; python {scripts_dir}/pip_box_features.py {ORFs_file} {working_directory} {gff_dir} {genome_dir}'
                if PIP:
                    cmds += ' --PIP'
                if hrp:
                    cmds += ' --hrp'
                if mxiE:
                    cmds += ' --mxiE'
                if exs:
                    cmds += ' --exs'
                if tts:
                    cmds += ' --tts'
                cmds +='\tcis_regulatory_elements\n'
                jobs_f.write(cmds)
        # signal peptide
    with open(f'{working_directory}/signalp.cmds','w') as sig_f:
        sig_f.write(f'module load python/python-anaconda3.6.5; python {scripts_dir}/signal_peptide_features.py {ORFs_file} {all_prots} {working_directory}\tsignalP\n')
    if not os.path.exists(f'{working_directory}/physical_features.done'):
        subprocess.call(f'/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py {working_directory}/features_jobs.cmds {tmp_dir} -q {queue}',shell=True)
    if not os.path.exists(f'{working_directory}/signal_p_features.done'):
        subprocess.call(f'/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py --cpu 3 {working_directory}/signalp.cmds {tmp_dir} -q {queue}',shell=True)
    if signal:
        if not os.path.exists(f'{working_directory}/Embedding_pred.csv.done'):
            subprocess.call(f'/opt/pbs/bin/qsub {tmp_dir}/Embedding.sh',shell=True)
    
    x = sum([item.endswith('done') for item in os.listdir(working_directory)])
    amount_of_expected_results = 3
    if signal:
        amount_of_expected_results += 1
    if organization:
        amount_of_expected_results += 1
        if CIS_elements:
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
    if signal:
        files_to_merge.append('Embedding_pred.csv')
        done_files.append('Embedding_pred.csv.done')
    if organization:
        files_to_merge.append('genome_organization_features.csv')
        done_files.append('genome_organization_features.done')
        if CIS_elements:
            files_to_merge.append('pip_box_features.csv')
            done_files.append('pip_box_features.done')
    failed_jobs = [done_job.split('.')[0] for done_job in done_files if not os.path.exists(done_job)]
    if len(failed_jobs) > 0:
        failed_str = ', '.join(failed_jobs)
        error_msg = f'Oups :(\nThe following jobs failed:\n\n{failed_str}'
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
        low_quality_flag = True
        #return create_effectors_html(effectors_prots,ORFs_file,working_directory)
        #error_msg = 'Learning failed. It can be due to a small training set, or other reasons. For further details you can contact us.'
        #fail(error_msg,error_path)
    # making final output files and tables
    
    in_f = r'out_learning/concensus_predictions.csv'
    out_f_normal = r'out_learning/concensus_predictions_with_annotation.csv'
    out_f_pseudo = r'out_learning/pseudogenes.csv'
    out_f_T3SS = r'out_learning/T3SS.csv'
    annotations = ORFs_file
    out_for_html_normal = r'out_learning/concensus_predictions_with_annotation_for_html.csv'
    out_for_html_pseudo = r'out_learning/pseudogenes_predictions_with_annotation_for_html.csv'
    out_T3SS_for_html = r'out_learning/T3SS_for_html.csv'
    add_annotations_to_predictions(in_f,out_f_normal,out_f_pseudo,annotations,out_f_T3SS)
    from csv_to_colored_xlsx_converter import convert_csv_to_colored_xlsx
    convert_csv_to_colored_xlsx(out_f_normal)
    convert_csv_to_colored_xlsx(out_f_pseudo)
    add_annotations_to_predictions(in_f,out_for_html_normal,out_for_html_pseudo,annotations,out_T3SS_for_html,line_end='<br>')
    predicted_table, positives_table, T3SS_table = make_html_tables(out_for_html_normal,out_T3SS_for_html)
    return predicted_table, positives_table, T3SS_table, low_quality_flag

import os
if __name__ == '__main__':
        from sys import argv
        import logging
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('error_path',
                            help='A path to the error file.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('input_ORFs_path',
                            help='A path to a DNA ORFs file.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('--input_effectors_path',
                            help='A path to a DNA fasta with positive samples. '
                                 'All samples in this file should be in the input ORFs file as well.',
                            default = '')
        parser.add_argument('output_dir_path',
                            help='A path to a folder in which the output files will be created.',
                            type=lambda path: path.rstrip('/'))
        parser.add_argument('--queue',default='pupkoweb')
        parser.add_argument('--organization', action='store_true',help='calculate organization features')
        parser.add_argument('--CIS_elements',action='store_true',help='extract regulatory elements features')
        parser.add_argument('--PIP', help='look for PIP-box in promoters', action='store_true')
        parser.add_argument('--hrp', help='look for hrp-box in promoters', action='store_true')
        parser.add_argument('--mxiE', help='look for mxiE-box in promoters', action='store_true')
        parser.add_argument('--exs', help='look for exs-box in promoters', action='store_true')
        parser.add_argument('--tts', help='look for tts-box in promoters', action='store_true')
        parser.add_argument('--translocation_signal',help='extract translocation signal feature', action='store_true')
        parser.add_argument('--html_path', default=None,
                            help='A path to an html file that will be updated during the run.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))

        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        effectors_learn(args.error_path,args.input_ORFs_path, args.input_effectors_path, args.output_dir_path, f'{args.output_dir_path}/tmp',args.queue,organization=args.organization,CIS_elements=args.CIS_elements,PIP=args.PIP,hrp=args.hrp,mxiE=args.mxiE,exs=args.exs,tts=args.tts,signal=args.translocation_signal)