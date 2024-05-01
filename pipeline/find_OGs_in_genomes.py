import fasta_parser
from sys import argv
import os
from Bio import SeqIO,SeqRecord
import subprocess
import csv
import shutil
import time

working_directory = argv[1]
os.chdir(working_directory)
if not os.path.exists('genomes_for_Microbializer'):
    os.makedirs('genomes_for_Microbializer')
#Effectidor_features_d = os.path.join(working_directory,'Effectidor_runs')
if os.path.exists('Effectidor_runs'):
    for genome in os.listdir('Effectidor_runs'):
        ORFs_f = os.path.join('Effectidor_runs',genome,'ORFs.fasta')
        locus_dic = fasta_parser.parse_ORFs(ORFs_f)
        recs = [SeqRecord.SeqRecord(locus_dic[locus],id=locus) for locus in locus_dic]
        SeqIO.write(recs,f'genomes_for_Microbializer/{genome}.fasta','fasta')
else:
    ORFs_f = 'ORFs.fasta'
    locus_dic = fasta_parser.parse_ORFs(ORFs_f)
    recs = [SeqRecord.SeqRecord(locus_dic[locus], id=locus) for locus in locus_dic]
    SeqIO.write(recs, f'genomes_for_Microbializer/genome_ORFs.fasta', 'fasta')


# Run Microbializer block
queue = 'power-pupko'
sh_file_content = f'''#!/bin/bash -x\n#PBS -S /bin/bash\n#PBS -q {queue}\n#PBS -o {working_directory}\n#PBS -e {working_directory}\n#PBS -N Microbializer_for_Effectidor\n
#PBS -r y\nhostname\necho job_name: Microbializer_for_Effectidor\n\nsource /powerapps/share/miniconda3-4.7.12/etc/profile.d/conda.sh\n
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer\nexport PATH=$CONDA_PREFIX/bin:$PATH\n\n
python /groups/pupko/yairshimony/microbializer_prod/pipeline/main.py --contigs_dir {working_directory}/genomes_for_Microbializer/ --output_dir output_OGs --add_orphan_genes_to_ogs --only_calc_ogs --inputs_are_annotated_genomes --bypass_number_of_genomes_limit\n'''
with open('search_OGs.pbs','w') as pbs_f:
    pbs_f.write(sh_file_content)
cmd = 'qsub search_OGs.pbs'
run_number = subprocess.check_output(cmd,shell=True).decode('ascii').strip()

while not os.path.exists(f'{run_number}.ER'): #while the job hasn't finished
    time.sleep(300)

# move Microbializer output to the working directory for future use
Microbializer_output_f = 'M1CR0B1AL1Z3R_output_OGs/07a_final_table/final_orthologs_table.csv'
shutil.copy(Microbializer_output_f,working_directory)

#subprocess.call("rm -r M1CR0B1AL1Z3R_output_OGs",shell=True)
#subprocess.call("rm -r output_OGs",shell=True)

with open('final_orthologs_table.csv') as ortho_table:
    ortho_reader = csv.reader(ortho_table)
    header = next(ortho_reader)
    with open(os.path.join(working_directory,'clean_orthologs_table.csv'),'w',newline='') as out_F:
        writer = csv.writer(out_F)
        writer.writerow(header)
        for row in ortho_reader:
            for i in range(1,len(header)):
                row[i]=row[i].replace(f'{header[i]}:','')
            writer.writerow(row)
