import fasta_parser
from sys import argv
import os
from Bio import SeqIO,SeqRecord
import subprocess
import csv
import shutil

working_directory = argv[1]
os.chdir(working_directory)
if not os.path.exists('genomes_for_Microbializer'):
    os.makedirs('genomes_for_Microbializer')
Effectidor_features_d = os.path.join(working_directory,'Effectidor_runs')

for genome in os.listdir(Effectidor_features_d):
    ORFs_f = os.path.join(Effectidor_features_d,genome,'ORFs.fasta')
    locus_dic = fasta_parser.parse_ORFs(ORFs_f)
    recs = [SeqRecord.SeqRecord(locus_dic[locus],id=locus) for locus in locus_dic]
    SeqIO.write(recs,f'genomes_for_Microbializer/{genome}.fasta','fasta')
 
# Run Microbializer block
'''TO COMPLETE'''

# move Microbializer output to the working directory for future use
Microbializer_output_f = 'M1CR0B1AL1Z3R_output_OGs/06e_final_table/final_orthologs_table.csv'
shutil.copy(Microbializer_output_f,working_directory)

#subprocess.call("rm -r M1CR0B1AL1Z3R_output_OGs",shell=True)
#subprocess.call("rm -r output_OGs",shell=True)

with open(os.path.join(working_directory,'final_orthologs_table.csv')) as ortho_table:
    ortho_reader = csv.reader(ortho_table)
    header = next(ortho_reader)
    with open(os.path.join(working_directory,'clean_orthologs_table.csv'),'w',newline='') as out_F:
        writer = csv.writer(out_F)
        writer.writerow(header)
        for row in ortho_reader:
            for i in range(1,len(header)):
                row[i]=row[i].replace(f'{header[i]}:','')
            writer.writerow(row)
