import argparse
import os
import subprocess
import pandas as pd
from Bio import SeqIO

RUN_WITH_CONDA = False
E_VALUE_CUT_OFF = 1e-10
QUERY_COVERAGE_PERCENTAGE_CUT_OFF = 0.3
FLAGELLA_SYSTEM = "Flagellar"
MMSEQS_OUTPUT_FORMAT = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,qcov,tstart,tend,evalue,bits'
COLUMN_NAMES = ['T3SS_protein', 'bacterial_protein', "identity_percent", "alignment length", "mismatch", "gapopen",
                "query_start", "query_end", "query_ length", 'query_coverage_percentage', "target_ start", "target_end", 'e_value', 'bit_score']
HEADERS = ["Subsystem_T3SS Protein", "Bacterial Protein ID"]


def run_mmseqs(query_files_directory, output_mmseqs, bacterial_proteome, tmp_mmseqs):
    for query_file in os.listdir(query_files_directory):
        query_file_path = os.path.join(query_files_directory, query_file)
        output_path = os.path.join(
            output_mmseqs, query_file.split(".")[0] + ".m8")
        mmseqs_command = f"mmseqs easy-search {query_file_path} {bacterial_proteome} {output_path} {tmp_mmseqs} --format-output {MMSEQS_OUTPUT_FORMAT}"
        if RUN_WITH_CONDA:
            conda_activate_command = ". ~/miniconda3/etc/profile.d/conda.sh; conda activate test;"
            run_mmseqs = conda_activate_command + mmseqs_command
            subprocess.run(run_mmseqs, shell=True)
        else:
            subprocess.run(mmseqs_command, shell=True)


def get_mmseqs_results_dictionary(mmseqs_results_file):
    mmseq_results_dict = {}  # {T3SS_protein: (bacterial_protein, bit_score)}
    df = pd.read_csv(mmseqs_results_file, sep='\t',
                     names=COLUMN_NAMES, header=None)
    filtered_df = df[(df['e_value'] < E_VALUE_CUT_OFF) & (
        df['query_coverage_percentage'] > QUERY_COVERAGE_PERCENTAGE_CUT_OFF)]

    for index, row in filtered_df.iterrows():
        T3SS_protein = row['T3SS_protein']
        bacterial_protein = row['bacterial_protein']
        bit_score = row['bit_score']

        if T3SS_protein not in mmseq_results_dict:
            mmseq_results_dict[T3SS_protein] = (bacterial_protein, bit_score)
        elif bit_score > mmseq_results_dict[T3SS_protein][1]:
            mmseq_results_dict[T3SS_protein] = (bacterial_protein, bit_score)

    return mmseq_results_dict


def get_all_subsystems_dict(output_mmseqs):

    # {subsystem_name: {T3SS_protein: (bacterial_protein, bit_score)}}
    all_subsystems_dict = {}

    for output_mmseqs_file in os.listdir(output_mmseqs):
        subsystem_name = output_mmseqs_file.split(".")[0]
        mmseqs_results_file = os.path.join(output_mmseqs, output_mmseqs_file)

        mmseqs_results_dictionary = get_mmseqs_results_dictionary(
            mmseqs_results_file)

        all_subsystems_dict[subsystem_name] = mmseqs_results_dictionary

    return all_subsystems_dict


def get_T3SS_homologous_bacterial_genes_list(all_subsystems_dict):
    T3SS_homologous_bacterial_genes = []
    for subsystem_dict in all_subsystems_dict.values():
        for bacterial_gene, bit_score in subsystem_dict.values():
            if bacterial_gene not in T3SS_homologous_bacterial_genes:
                T3SS_homologous_bacterial_genes.append(bacterial_gene)
    return T3SS_homologous_bacterial_genes


def get_best_bacterial_T3SS_match_dict(all_subsystems_dict):

    T3SS_homologous_bacterial_genes = get_T3SS_homologous_bacterial_genes_list(
        all_subsystems_dict)
    # {bacterial_gene: (system_gene, subsystem)}
    best_bacterial_T3SS_match = {}

    for homologous_bacterial_gene in T3SS_homologous_bacterial_genes:
        max_bit_score = 0
        best_system_gene = None
        best_subsystem = None
        for subsystem_name, subsystem_dict in all_subsystems_dict.items():
            for subsystem_gene, (bacterial_gene, bit_score) in subsystem_dict.items():
                if homologous_bacterial_gene == bacterial_gene:
                    if bit_score > max_bit_score:
                        max_bit_score = bit_score
                        best_system_gene = subsystem_gene
                        best_subsystem = subsystem_name
        if best_system_gene and best_subsystem:
            best_bacterial_T3SS_match[homologous_bacterial_gene] = (
                best_system_gene, best_subsystem)

    return best_bacterial_T3SS_match


def get_full_bacterial_T3SS_dict(T3SS_data, best_bacterial_T3SS_match_dict):
    # If a homologous bacterial gene is found - {bacterial_gene: (system_gene, subsystem)}. If no homologous bacterial gene is found - {integer: (system_gene, subsystem)}
    full_bacterial_T3SS_dict = {}
    i = 1
    for T3SS_data_file in os.listdir(T3SS_data):
        subsystem_name = T3SS_data_file.split(".")[0]
        if any(subsystem_name in value for value in best_bacterial_T3SS_match_dict.values()):
            path_to_T3SS_data_file = os.path.join(T3SS_data, T3SS_data_file)
            T3SS_proteins_names = set(
                [rec.id for rec in SeqIO.parse(path_to_T3SS_data_file, "fasta")])
            for T3SS_protein in T3SS_proteins_names:
                if (T3SS_protein, subsystem_name) in best_bacterial_T3SS_match_dict.values():
                    for key, value in best_bacterial_T3SS_match_dict.items():
                        if (T3SS_protein, subsystem_name) == value:
                            full_bacterial_T3SS_dict[key] = value
                else:
                    full_bacterial_T3SS_dict[int(i)] = (
                        T3SS_protein, subsystem_name)
                    i = i + 1
    return full_bacterial_T3SS_dict


def write_dict_to_output_file(full_bacterial_T3SS_dict, output_file):
    data = []
    flagella_data = []

    for bacterial_gene, (T3SS_protein, subsystem) in full_bacterial_T3SS_dict.items():
        if isinstance(bacterial_gene, int):
            entry = [f"{subsystem}_{T3SS_protein}", None]
        else:
            entry = [f"{subsystem}_{T3SS_protein}", bacterial_gene]
        if subsystem == FLAGELLA_SYSTEM:
            flagella_data.append(entry)
        else:
            data.append(entry)

    data.extend(flagella_data)
    df = pd.DataFrame(data, columns=HEADERS)
    df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser(
        description='Run mmseqs versus the different T3SS subtypes and process the results files to produce the final report of the T3SS components in the genome')
    parser.add_argument('working_directory', type=str,
                        help='Path to the working directory')
    parser.add_argument('bacterial_proteome', type=str,
                        help='Name of the bacterial proteome file')
    parser.add_argument('T3SS_data', type=str,
                        help='Name of T3SS datasets directory')
    args = parser.parse_args()

    bacterial_proteome = os.path.join(
        args.working_directory, args.bacterial_proteome)
    T3SS_data = os.path.join(args.working_directory, args.T3SS_data)
    tmp_mmseqs = os.path.join(args.working_directory, "tmp_mmseqs")
    output_mmseqs = os.path.join(args.working_directory, "output_mmseqs")
    output_file = os.path.join(args.working_directory, "T3SS.csv")

    os.makedirs(tmp_mmseqs, exist_ok=True)
    os.makedirs(output_mmseqs, exist_ok=True)

    run_mmseqs(T3SS_data, output_mmseqs, bacterial_proteome, tmp_mmseqs)
    all_subsystems_dict = get_all_subsystems_dict(output_mmseqs)
    best_bacterial_T3SS_match_dict = get_best_bacterial_T3SS_match_dict(
        all_subsystems_dict)
    full_bacterial_T3SS_dict = get_full_bacterial_T3SS_dict(
        T3SS_data, best_bacterial_T3SS_match_dict)
    write_dict_to_output_file(full_bacterial_T3SS_dict, output_file)


if __name__ == "__main__":
    main()
