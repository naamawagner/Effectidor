import fasta_parser
from sys import argv
import os
from Bio import SeqIO, SeqRecord
import subprocess
import csv
import shutil
import time
import sys
sys.path.append('/lsweb/josef_sites/effectidor/auxiliaries')
from auxiliaries import fail

working_directory = argv[1]
identity_cutoff = argv[2]  # default: '50'
coverage_cutoff = argv[3] # default: '60'
error_path = os.path.join(working_directory, 'error_OGs.txt')
#scripts_dir = '/groups/pupko/naamawagner/T3Es_webserver/scripts/OGs/scripts/'
scripts_dir = '/lsweb/josef_sites/effectidor/auxiliaries/'
os.chdir(working_directory)
if not os.path.exists('genomes_for_Microbializer'):
    os.makedirs('genomes_for_Microbializer')

if os.path.exists('Effectidor_runs'):
    for genome in os.listdir('Effectidor_runs'):
        ORFs_f = os.path.join('Effectidor_runs', genome, 'ORFs.fasta')
        locus_dic = fasta_parser.parse_ORFs(ORFs_f)
        recs = [SeqRecord.SeqRecord(locus_dic[locus], id=locus) for locus in locus_dic]
        SeqIO.write(recs, f'genomes_for_Microbializer/{genome}.fasta', 'fasta')
else:
    ORFs_f = 'ORFs.fasta'
    locus_dic = fasta_parser.parse_ORFs(ORFs_f)
    recs = [SeqRecord.SeqRecord(locus_dic[locus], id=locus) for locus in locus_dic]
    SeqIO.write(recs, f'genomes_for_Microbializer/genome_ORFs.fasta', 'fasta')


# Run Microbializer block
cmds_file_content = f'''source /lsweb/pupko/microbializer/miniconda3/etc/profile.d/conda.sh!@#\
conda activate /lsweb/pupko/microbializer/miniconda3/envs/microbializer!@#\
export PATH=$CONDA_PREFIX/bin:$PATH!@#!@#\
python /lsweb/pupko/microbializer/pipeline/main.py --contigs_dir \
{working_directory}/genomes_for_Microbializer/ --output_dir output_OGs --add_orphan_genes_to_ogs true --only_calc_ogs true \
--inputs_fasta_type orfs --bypass_number_of_genomes_limit true --identity_cutoff {identity_cutoff} --coverage_cutoff \
{coverage_cutoff} --account_name pupkoweb-users -q pupkoweb\tMicrobializer_for_Effectidor'''
with open('search_OGs.cmds', 'w') as cmds_f:
    cmds_f.write(cmds_file_content)
cmd = f'python {os.path.join(scripts_dir, "q_submitter_power.py")} search_OGs.cmds {working_directory} --memory 4'
subprocess.call(cmd, shell=True)

Microbializer_output_f = 'M1CR0B1AL1Z3R_output_OGs/05a_orthogroups/orthogroups.csv'
while not os.path.exists(Microbializer_output_f):
    if os.path.exists('M1CR0B1AL1Z3R_output_OGs/error.txt'):
        with open('M1CR0B1AL1Z3R_output_OGs/error.txt') as error:
            error_msg = error.read()
        fail(error_msg, error_path)
    # while the job hasn't finished
    time.sleep(60)

# move Microbializer output to the working directory for future use
shutil.copy(Microbializer_output_f, working_directory)

# subprocess.call("rm -r M1CR0B1AL1Z3R_output_OGs",shell=True)
# subprocess.call("rm -r output_OGs",shell=True)

with open('orthogroups.csv') as ortho_table:
    ortho_reader = csv.reader(ortho_table)
    header = next(ortho_reader)
    with open(os.path.join(working_directory, 'clean_orthologs_table.csv'), 'w', newline='') as out_F:
        writer = csv.writer(out_F)
        writer.writerow(header)
        for row in ortho_reader:
            for i in range(1, len(header)):
                row[i] = row[i].replace(f'{header[i]}:', '')
            writer.writerow(row)
