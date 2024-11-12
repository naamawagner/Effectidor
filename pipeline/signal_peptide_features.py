import os
import fasta_parser
import csv
from sys import argv
import subprocess

ORFs_fasta = argv[1]
proteins_seq = argv[2]
working_directory = argv[3]
os.chdir(working_directory)
locus_dic = fasta_parser.parse_ORFs(ORFs_fasta)
###############################################################################

signal_p_out = f"signalp_results"

def make_signal():
    cmd = f"module load signalp6-6.0; signalp6 --fastafile {proteins_seq} --output_dir {signal_p_out} --format none"
    subprocess.check_output(cmd, shell=True)
    
def parse_signalp_out(out_f):
    sigp_dic = {}
    with open(out_f, 'r') as result_f:
        for line in result_f:
            if not line.startswith('#'):
                existence = 0
                result = line.split('\t')
                locus = result[0]
                if result[1] == "SP":
                    existence = 1
                sigp_dic[locus] = existence
    return sigp_dic
                        
make_signal()
signalp_dic = parse_signalp_out(os.path.join(signal_p_out, 'prediction_results.txt'))

with open(f'signal_p_features.csv', 'w', newline='') as sigp_f:
    f_writer = csv.writer(sigp_f)
    header = ['locus', 'SecYEG_existence_(SignalP_6)']
    f_writer.writerow(header)
    for locus in locus_dic:
        # while locus not in signalp_dic:  # in case the process stopped before finishing
        #     make_signal()
        #     signalp_dic = parse_signalp_out(signal_p_out)
        l = [locus]
        score = signalp_dic[locus]
        l.append(score)
        f_writer.writerow(l)
    
endfile = open('signal_p_features.done', 'w')
endfile.close()
    