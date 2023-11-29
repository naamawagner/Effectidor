
def protein_mmseqs_all_vs_all(query, dataset, out, mmseq_tmp):
    import subprocess
    command = f'mmseqs easy-search {query} {dataset} {out} {mmseq_tmp}'
    subprocess.check_output(command,shell=True)

