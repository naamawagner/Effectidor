
def protein_mmseqs_all_vs_all(query, dataset, out, mmseq_tmp):
    import subprocess
    command = f'mmseqs easy-search {query} {dataset} {out} {mmseq_tmp} --threads 1'
    subprocess.check_output(command, shell=True)
