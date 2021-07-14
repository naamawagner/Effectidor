import os
import fasta_parser
import csv
from sys import argv
import subprocess

ORFs_fasta = argv[1]
proteins_seq = argv[2]
working_directory = argv[3]
os.chdir(working_directory)
locus_dic=fasta_parser.parse_ORFs(ORFs_fasta)
###############################################################################

if not os.path.exists("signalp_results"):
    os.makedirs("signalp_results")
signal_p_out = f"signalp_results/results_summary.txt"

def make_signal():
    cmd=f"module load signalp/4.1; signalp -f summary -t gram- {proteins_seq} > {signal_p_out}"
    subprocess.check_output(cmd,shell=True)
    
def parse_signalp_out(out_f):
    sigp_dic = {}
    with open(out_f,'r') as result_f:
        for line in result_f:
            if line.startswith('Name'):
                existence,D=0,0
                result = line.split()
                locus = result[0].split('=')[1]
                if result[1][3:] == "'YES'":
                    existence = 1
                for segment in result:
                    if segment.startswith('D='):
                        D = float(segment[2:])
                sigp_dic[locus]=(existence,D)
    return sigp_dic
                        
make_signal()
signalp_dic = parse_signalp_out(signal_p_out)

with open(f'signal_p_features.csv','w',newline='') as sigp_f:
    f_writer = csv.writer(sigp_f)
    header = ['locus','SecYEG_existence_(SignalP)','SecYEG_score_(SignalP)']
    f_writer.writerow(header)
    for locus in locus_dic:
        while locus not in signalp_dic: # in case the process stopped before finishing
            make_signal()
            signalp_dic = parse_signalp_out(signal_p_out)
        l = [locus]
        score,D = signalp_dic[locus]
        l.append(score)
        l.append(D)
        f_writer.writerow(l)
    
endfile = open('signal_p_features.done','w')
endfile.close()
    