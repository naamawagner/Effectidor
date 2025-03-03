def write_sh_file(tmp_dir, path_to_fasta, path_for_output_csv):
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

    with open(f'{tmp_dir}/Embedding.cmds', 'w') as out:
        out.write(content)


def fail(error_msg, error_file_path):
    with open(error_file_path, 'w') as error_f:
        error_f.write(error_msg + '\n')
    raise Exception(error_msg)


def effectors_learn(error_path, ORFs_file, effectors_file, working_directory, tmp_dir, gff_file, full_genome_f,
                    signal=False, PIP=False, hrp=False, mxiE=False, exs=False, tts=False, homology_search=False,
                    coverage=50):
    import pandas as pd
    import subprocess
    import os
    import time
    from Bio import SeqIO
    import shutil
    from add_annotations_to_predictions import add_annotations_to_predictions, make_html_tables, create_annotations_f
    import sys
    # vars
    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = f'{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/data'
    # os.chdir(working_directory)
    all_prots = os.path.join(working_directory, 'translated_ORFs.faa')
    effectors_prots = os.path.join(working_directory, 'translated_effectors.faa')
    blast_datasets_dir = os.path.join(data_dir, 'blast_data')
    chaperones = os.path.join(data_dir, 'chaperones.faa')

    os.makedirs(tmp_dir, exist_ok=True)
    if not os.path.exists(f'{working_directory}/blast_data/T3Es.faa'):  # if no input for effectors homology search
        # was supplied
        os.makedirs(f'{working_directory}/blast_data', exist_ok=True)  # this directory was already created if
        # another input for homology search was supplied (host proteome or proteomes of closely related bacteria
        # without T3SS)
        # cmd = f'cp {blast_datasets_dir}/*.faa {working_directory}/blast_data/'
        # subprocess.check_output(cmd, shell=True)
        for f in os.listdir(blast_datasets_dir):
            if f.endswith('.faa'):
                shutil.copy(os.path.join(blast_datasets_dir, f), os.path.join(working_directory, 'blast_data'))
    # cmd = f'cp -r {os.path.join(data_dir,"T3SS_data")} {working_directory}'
    # subprocess.check_output(cmd, shell=True)
    shutil.copytree(os.path.join(data_dir, "T3SS_data"), os.path.join(working_directory, "T3SS_data"))

    blast_datasets_dir = f'{working_directory}/blast_data'
            
    # feature extraction step
    
    # translate the input fasta files
    subprocess.check_output(
        ['python', f'{scripts_dir}/translate_fasta.py', ORFs_file, effectors_file, all_prots, effectors_prots])

    if signal: #make N-terminal sequences for Signal search
        #os.makedirs(signal_prot_dir, exist_ok=True)
        recs = SeqIO.parse(all_prots, 'fasta')
        N_terminals = []
        for rec in recs:
            ID = rec.id
            N_seq = rec.seq[:100]
            rec.seq = N_seq
            rec.description = ID
            N_terminals.append(rec)
        N_terminal_file_path = os.path.join(working_directory, 'N_terminals.faa')
        SeqIO.write(N_terminals, N_terminal_file_path, 'fasta')
        signal_script = os.path.join(scripts_dir, r'T3Es_secretion_signal_prediction/src/inference',
                                     'predict_secretion_signal.py')
        signal_out = f'{working_directory}/Embedding_pred.csv'
        cmd = f'python {signal_script} --input_fasta_file {N_terminal_file_path} --output_file {signal_out}'
        if True:  # ADD parameter for large model
            cmd += ' --use_large_model'
        subprocess.check_output(cmd, shell=True)

    if not effectors_file:
        subprocess.check_output(
            ['python', f'{scripts_dir}/find_effectors.py', f'{blast_datasets_dir}/T3Es.faa', all_prots,
             effectors_prots, coverage, working_directory])
    elif homology_search:
        effectors_prots2 = 'homology_found_effectors.faa'
        subprocess.check_output(
            ['python', f'{scripts_dir}/find_effectors.py', f'{blast_datasets_dir}/T3Es.faa', all_prots,
             effectors_prots2, coverage])
        eff1 = SeqIO.to_dict(SeqIO.parse(effectors_prots, 'fasta'))
        eff2 = SeqIO.to_dict(SeqIO.parse(effectors_prots2, 'fasta'))
        eff1.update(eff2)
        SeqIO.write(eff1.values(), effectors_prots, 'fasta')
    # find and create non effectors fasta file
    subprocess.check_output(['python', f'{scripts_dir}/find_non_effectors.py', all_prots, effectors_prots,
                             working_directory])
    subprocess.check_output(
        ['python', f'{scripts_dir}/find_T3SS_components.py', working_directory, all_prots, 'T3SS_data'])
    subprocess.check_output(['python', f'{scripts_dir}/find_chaperones.py', all_prots, working_directory, chaperones])
    # creating the features extraction commands
    Features_cmds = [f'python {scripts_dir}/sequence_features.py {ORFs_file} {working_directory}',
                     f'python {scripts_dir}/homology.py {ORFs_file} {all_prots} {blast_datasets_dir} {working_directory}']
    if gff_file:
        # genome organization features
        Features_cmds.append(f'python {scripts_dir}/genome_organization.py {ORFs_file} {effectors_prots} {working_directory} {gff_file}')
        # mobile genetic elements features
        Features_cmds.append(f'python {scripts_dir}/mobile_genetic_elements.py {gff_file} {working_directory}/mobile_genetic_elements.csv')
        if full_genome_f:
            # regulatory elements features
            cmd = f'python {scripts_dir}/pip_box_features.py {ORFs_file} {working_directory} {gff_file} {full_genome_f}'
            if PIP:
                cmd += ' --PIP'
            if hrp:
                cmd += ' --hrp'
            if mxiE:
                cmd += ' --mxiE'
            if exs:
                cmd += ' --exs'
            if tts:
                cmd += ' --tts'
            Features_cmds.append(cmd)
    if not os.path.exists(f'{working_directory}/physical_features.done'):
        for cmd in Features_cmds:
            subprocess.call(cmd, shell=True)

    create_annotations_f(ORFs_file, gff_file, os.path.join(working_directory, 'annotations.csv'),
                         os.path.join(working_directory, 'pseudogenes.txt'))
    x = sum([item.endswith('done') for item in os.listdir(working_directory)])
    # amount_of_expected_results = 3
    amount_of_expected_results = 2
    # if signal:
    #     amount_of_expected_results += 1
    if gff_file:
        amount_of_expected_results += 2  # Genome organization (effector proximity and MGE)
        if full_genome_f:
            amount_of_expected_results += 1
            
    while x < amount_of_expected_results:
        print(f'number of expected outputs is {str(amount_of_expected_results)} , number of existing outputs is '
              f'{str(x)}, waiting for additional {str(amount_of_expected_results-x)} outputs')
        time.sleep(60)
        x = sum([item.endswith('done') for item in os.listdir(working_directory)])
        
    effectors = SeqIO.to_dict(SeqIO.parse(effectors_prots, 'fasta')).keys()
    non_effectors = SeqIO.to_dict(SeqIO.parse(os.path.join(working_directory, 'non_effectors.faa'), 'fasta')).keys()
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

    files_to_merge = [os.path.join(working_directory, 'physical_features.csv'),
                      os.path.join(working_directory, 'homology_features.csv')]
    done_files = [os.path.join(working_directory, 'physical_features.done'),
                  os.path.join(working_directory, 'homology_features.done')]
    if signal:
        files_to_merge.append(os.path.join(working_directory, 'Embedding_pred.csv'))
    #     done_files.append('Embedding_pred.csv.done')
    if gff_file:
        files_to_merge.append(os.path.join(working_directory, 'genome_organization_features.csv'))
        done_files.append(os.path.join(working_directory, 'genome_organization_features.done'))
        files_to_merge.append(os.path.join(working_directory, 'mobile_genetic_elements.csv'))
        done_files.append(os.path.join(working_directory, 'mobile_genetic_elements.csv.done'))
        if full_genome_f:
            files_to_merge.append(os.path.join(working_directory, 'pip_box_features.csv'))
            done_files.append(os.path.join(working_directory, 'pip_box_features.done'))
    failed_jobs = [done_job.split('.')[0] for done_job in done_files if not os.path.exists(done_job)]
    if len(failed_jobs) > 0:
        failed_str = ', '.join(failed_jobs)
        error_msg = f'The following jobs have failed:\n\n{failed_str}'
        fail(error_msg, error_path)
    merged_df = pd.read_csv(files_to_merge[0], dtype={'locus': str})
    for f in files_to_merge[1:]:
        df1 = pd.read_csv(f, dtype={'locus': str})
        merged_df = merged_df.merge(df1)
    merged_df = merged_df.merge(label_df)
    merged_df.to_csv(os.path.join(working_directory, 'features.csv'), index=False)


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
        parser.add_argument('--homology_search',
                            help='search additional effectors based on homology to internal dataset',
                            action='store_true')
        parser.add_argument('--coverage', help='minimal coverage with effector homolog', default='50')
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        effectors_learn(args.error_path, args.input_ORFs_path, args.input_effectors_path, args.output_dir_path,
                        f'{args.output_dir_path}/tmp', args.gff_file, args.full_genome_f,
                        signal=args.translocation_signal, PIP=args.PIP, hrp=args.hrp, mxiE=args.mxiE, exs=args.exs,
                        tts=args.tts, homology_search=args.homology_search, coverage=args.coverage)
