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
            
def write_sh_file(tmp_dir,path_to_fasta,path_for_output_csv,queue='power-pupko'):
    content = f'''export HOME=/groups/pupko/yairshimony!@#\
    # Print some information about the job!@#\
    echo "Starting my SLURM job"!@#\
    echo "Job ID: $SLURM_JOB_ID"!@#\
    echo "Running on nodes: $SLURM_JOB_NODELIST"!@#\
    echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"!@#\
    echo "Cuda visible devices: $CUDA_VISIBLE_DEVICES"!@#\
    source ~/miniconda3/etc/profile.d/conda.sh!@#\
    conda activate secretion_signal!@#\
    export PATH=$CONDA_PREFIX/bin:$PATH!@#\
    python ~/secretion_signal_prediction/src/inference/predict_secretion_signal.py --input_fasta_file {path_to_fasta} \
    --output_file {path_for_output_csv} --use_large_model!@#\
    touch {path_for_output_csv}.done!@#\
    touch {os.path.join(tmp_dir,"Embedding_pred.csv.done")}\tEffectidor_embedding'''

    with open(f'{tmp_dir}/Embedding.cmds','w') as out:
        out.write(content)

def effectors_learn(error_path, ORFs_file, effectors_file, working_directory, tmp_dir, queue, gff_file, full_genome_f,
                    PIP=False, hrp=False, mxiE=False, exs=False, tts=False, homology_search=False, signal=False,
                    signalp=False, MGE=True):
    import pandas as pd
    import subprocess
    import os
    import time
    from Bio import SeqIO
    import shutil
    from add_annotations_to_predictions import add_annotations_to_predictions,make_html_tables,create_annotations_f
    import sys
    #sys.path.append('/bioseq/effectidor/auxiliaries')
    from auxiliaries import fail
    import effectidor_CONSTANTS
    low_quality_flag = False
    # vars
    scripts_dir = effectidor_CONSTANTS.EFFECTIDOR_EXEC
    data_dir = effectidor_CONSTANTS.EFFECTIDOR_DATA
    os.chdir(working_directory)
    # input files
    #ORFs_file = 'ORFs.fasta'
    #effectors_file = 'effectors.fasta'
    all_prots = 'translated_ORFs.faa'
    effectors_prots = os.path.join(working_directory, 'translated_effectors.faa')
    blast_datasets_dir = f'{data_dir}/blast_data'

    os.makedirs(tmp_dir,exist_ok=True)
    if not os.path.exists(f'{working_directory}/blast_data/T3Es.faa'): #if no input for effectors homology search was supplied
        os.makedirs(f'{working_directory}/blast_data',exist_ok=True) # this directory was already created if another input for homology search was supplied (host proteome or proteomes of closely related bacteria without T3SS)
        cmd = f'cp {blast_datasets_dir}/*.faa {working_directory}/blast_data/'
        subprocess.check_output(cmd,shell=True)
    cmd = f'cp -r {os.path.join(data_dir,"T3SS_data")} {working_directory}'
    subprocess.check_output(cmd, shell=True)

    blast_datasets_dir = f'{working_directory}/blast_data'
            
    # feature extraction step
    
    # translate the input fasta files
    subprocess.check_output(['python',f'{scripts_dir}/translate_fasta.py',ORFs_file,effectors_file,all_prots,effectors_prots])

    if signal: #make N-terminal sequences for Signal search
        #os.makedirs(signal_prot_dir, exist_ok=True)
        recs = SeqIO.parse(all_prots,'fasta')
        N_terminals = []
        for rec in recs:
            ID = rec.id
            N_seq = rec.seq[:100]
            rec.seq = N_seq
            rec.description = ID
            N_terminals.append(rec)
        N_terminal_file_path = os.path.join(working_directory, 'N_terminals.faa')
        SeqIO.write(N_terminals, N_terminal_file_path, 'fasta')
        write_sh_file(tmp_dir, N_terminal_file_path, f'{working_directory}/Embedding_pred.csv')

    
    if not effectors_file:
        subprocess.check_output(['python', f'{scripts_dir}/find_effectors.py', f'{blast_datasets_dir}/T3Es.faa', all_prots, effectors_prots])
    elif homology_search:
        effectors_prots2 = 'homology_found_effectors.faa'
        subprocess.check_output(['python', f'{scripts_dir}/find_effectors.py', f'{blast_datasets_dir}/T3Es.faa', all_prots, effectors_prots2])
        eff1 = SeqIO.to_dict(SeqIO.parse(effectors_prots, 'fasta'))
        eff2 = SeqIO.to_dict(SeqIO.parse(effectors_prots2, 'fasta'))
        eff1.update(eff2)
        SeqIO.write(eff1.values(), effectors_prots, 'fasta')
    # find and create non effectors fasta file
    subprocess.check_output(['python', f'{scripts_dir}/find_non_effectors.py', all_prots, effectors_prots])
    subprocess.check_output(['python', f'{scripts_dir}/find_T3SS_components.py', working_directory, all_prots, 'T3SS_data'])
    # creating the features extraction commands
    with open(f'{working_directory}/features_jobs.cmds','w') as jobs_f:
        # sequence features
        jobs_f.write(f'''module load MMseqs2/May2024; python {scripts_dir}/sequence_features.py {ORFs_file} {effectors_prots} {working_directory}!@#\
                        touch {os.path.join(tmp_dir,"physical_features.done")}\tsequence_features\n''')
        # sequence similarity features
        jobs_f.write(f'''module load MMseqs2/May2024; python {scripts_dir}/homology.py {ORFs_file} {all_prots} {blast_datasets_dir} {working_directory}!@#\
                        touch {os.path.join(tmp_dir,"homology_features.done")}\thomology\n''')
        if gff_file:
            # genome organization features
            jobs_f.write(f'''module load MMseqs2/May2024; python {scripts_dir}/genome_organization.py {ORFs_file} {effectors_prots} {working_directory} {gff_file}!@#\
                        touch {os.path.join(tmp_dir,"genome_organization_features.done")}\tgenome_organization\n''')
            if MGE:
                jobs_f.write(
                    f'''module load MMseqs2/May2024; python {scripts_dir}/mobile_genetic_elements.py {gff_file} {working_directory}/mobile_genetic_elements.csv!@#\
                    touch {os.path.join(tmp_dir, "mobile_genetic_elements.csv.done")}\tmobile_genetic_elements\n''')
            if full_genome_f:
                # regulatory elements features
                cmds = f'module load MMseqs2/May2024; python {scripts_dir}/pip_box_features.py {ORFs_file} {working_directory} {gff_file} {full_genome_f}'
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
                cmds += f'!@#touch {os.path.join(tmp_dir,"pip_box_features.done")}\tcis_regulatory_elements\n'
                jobs_f.write(cmds)
        # signal peptide
    if signalp:
        with open(f'{working_directory}/signalp.cmds','w') as sig_f:
            sig_f.write(f'''module load MMseqs2/May2024; python {scripts_dir}/signal_peptide_features.py {ORFs_file} {all_prots} {working_directory}!@#\
                            touch {os.path.join(tmp_dir,"signal_p_features.done")}\tsignalP\n''')
        if not os.path.exists(f'{working_directory}/signal_p_features.done'):
            subprocess.call(
                f'{os.path.join(scripts_dir, "q_submitter.py")} --cpu 3 {working_directory}/signalp.cmds {tmp_dir} -q {queue}',
                shell=True)
    if not os.path.exists(f'{working_directory}/physical_features.done'):
        subprocess.call(f'{os.path.join(scripts_dir, "q_submitter.py")} {working_directory}/features_jobs.cmds {tmp_dir} -q {queue}', shell=True)
    if signal:
        if not os.path.exists(f'{working_directory}/Embedding_pred.csv.done'):
            cmd = f'{os.path.join(scripts_dir, "q_submitter.py")} {os.path.join(tmp_dir,"Embedding.cmds")} {tmp_dir} -q {queue} --cpu 10 --memory 15'
            subprocess.call(cmd, shell=True)
    # TO COMPLETE: make sure the jobs were submitted before proceeding
    create_annotations_f(ORFs_file, gff_file, 'annotations.csv', 'pseudogenes.txt')
    x = sum([item.endswith('done') for item in os.listdir(working_directory)])
    # amount_of_expected_results = 3
    amount_of_expected_results = 2
    if signal:
        amount_of_expected_results += 1
    if signalp:
        amount_of_expected_results += 1
    if gff_file:
        amount_of_expected_results += 1
        if MGE:
            amount_of_expected_results += 1
        if full_genome_f:
            amount_of_expected_results += 1
            
    while x < amount_of_expected_results:
        print('number of expected outputs is '+str(amount_of_expected_results)+' , number of existing outputs is '+str(x)+', waiting for additional '+str(amount_of_expected_results-x)+' outputs')
        time.sleep(60)
        num_finished_jobs = sum([item.endswith('done') for item in os.listdir(tmp_dir)])
        x = sum([item.endswith('done') for item in os.listdir(working_directory)])
        if x < num_finished_jobs:
            break
        
    effectors = SeqIO.to_dict(SeqIO.parse(effectors_prots, 'fasta')).keys()
    non_effectors = SeqIO.to_dict(SeqIO.parse('non_effectors.faa', 'fasta')).keys()
    all_loci = SeqIO.to_dict(SeqIO.parse(all_prots, 'fasta')).keys()
    label_dict = {}
    for locus in all_loci:
        if locus in effectors:
            label_dict[locus] = [locus, 'effector']
        elif locus in non_effectors:
            label_dict[locus] = [locus, 'no']
        else:
            label_dict[locus] = [locus, '?']
    label_df = pd.DataFrame.from_dict(label_dict, orient='index', columns=['locus', 'is_effector'])

    files_to_merge = ['physical_features.csv', 'homology_features.csv']
    done_files = ['physical_features.done', 'homology_features.done']
    if signalp:
        files_to_merge.append('signal_p_features.csv')
        done_files.append('signal_p_features.done')
    if signal:
        files_to_merge.append('Embedding_pred.csv')
        done_files.append('Embedding_pred.csv.done')
    if gff_file:
        files_to_merge.append('genome_organization_features.csv')
        done_files.append('genome_organization_features.done')
        if MGE:
            files_to_merge.append('mobile_genetic_elements.csv')
            done_files.append(('mobile_genetic_elements.csv.done'))
        if full_genome_f:
            files_to_merge.append('pip_box_features.csv')
            done_files.append('pip_box_features.done')
    failed_jobs = [done_job.split('.')[0] for done_job in done_files if not os.path.exists(done_job) and os.path.exists(
        os.path.join(tmp_dir, done_job))]
    if len(failed_jobs) > 0:
        failed_str = ', '.join(failed_jobs)
        error_msg = f'The following jobs have failed:\n\n{failed_str}'
        fail(error_msg, error_path)
    merged_df = pd.read_csv(files_to_merge[0])
    for f in files_to_merge[1:]:
        df1 = pd.read_csv(f)
        merged_df = merged_df.merge(df1)
    merged_df = merged_df.merge(label_df)
    merged_df.to_csv('features.csv',index=False)


import os
if __name__ == '__main__':
        from sys import argv
        import logging
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('error_path',
                            help='A path to the error file.')
        parser.add_argument('input_ORFs_path',
                            help='A path to a DNA ORFs file.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('--input_effectors_path',
                            help='A path to a DNA fasta with positive samples. '
                                 'All samples in this file should be in the input ORFs file as well.',
                            default='')
        parser.add_argument('output_dir_path',
                            help='A path to a folder in which the output files will be created.',
                            type=lambda path: path.rstrip('/'))
        parser.add_argument('--queue', default='pupkoweb')
        parser.add_argument('--gff_file', help='A path to a gff file, to calculate genome organization features',
                            default='')
        parser.add_argument('--full_genome_f',
                            help='A path to a full genome fasta file to extract regulatory elements features (used together with the gff file)',
                            default='')
        parser.add_argument('--PIP', help='look for PIP-box in promoters', action='store_true')
        parser.add_argument('--hrp', help='look for hrp-box in promoters', action='store_true')
        parser.add_argument('--mxiE', help='look for mxiE-box in promoters', action='store_true')
        parser.add_argument('--exs', help='look for exs-box in promoters', action='store_true')
        parser.add_argument('--tts', help='look for tts-box in promoters', action='store_true')
        parser.add_argument('--translocation_signal', help='extract translocation signal feature', action='store_true')
        parser.add_argument('--signalp', help='extract SignalP6 feature', action='store_true')
        parser.add_argument('--mobile_genetics_elements', help='extract distance from mobile genetics elements', action='store_true')
        parser.add_argument('--homology_search',
                            help='search additional effectors based on homology to internal dataset',
                            action='store_true')
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        effectors_learn(args.error_path, args.input_ORFs_path, args.input_effectors_path, args.output_dir_path,
                        f'{args.output_dir_path}/tmp', args.queue, args.gff_file, args.full_genome_f, PIP=args.PIP,
                        hrp=args.hrp, mxiE=args.mxiE, exs=args.exs, tts=args.tts, homology_search=args.homology_search,
                        signal=args.translocation_signal, signalp=args.signalp)
