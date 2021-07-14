import fasta_parser
import csv
from sys import argv
import os
ORFs_dir = argv[1]
effectors_file = argv[2]
working_directory = argv[3]
os.chdir(working_directory)
locus_dics=[fasta_parser.parse_ORFs(f'{ORFs_dir}/{ORFs_file}') for ORFs_file in os.listdir(ORFs_dir)]
effectors_dict = fasta_parser.parse_ORFs(effectors_file)
#%%

all_locuses_in_genome=[list(locus_dic.keys()) for locus_dic in locus_dics]
binary_lsts=[] # effectors locations list
for contig in all_locuses_in_genome:
    binary_l = []
    for locus in contig:
        if locus in effectors_dict:
            binary_l.append(1)
        else:
            binary_l.append(0)
    binary_lsts.append(binary_l)

def closest_effector(locus):
    for i in range(len(all_locuses_in_genome)):
        if locus in all_locuses_in_genome[i]:
            contig = all_locuses_in_genome[i]
            binary_l = binary_lsts[i]
            gene_place = contig.index(locus)
            upstream_closest = 0
            for i in range(1,len(contig)-gene_place):
                if binary_l[gene_place+i]==1:
                    upstream_closest = i
                    break
            if upstream_closest==0:
                for i in range(gene_place):
                    if binary_l[i]==1:
                        upstream_closest=i+len(contig)-gene_place
                        break
            downstream_closest = 0
            for i in range(1,len(contig)):
                if binary_l[gene_place-i] == 1:
                    downstream_closest = i
                    break
            if upstream_closest != 0 and downstream_closest != 0: # there is a neighbor effector
                return min(upstream_closest,downstream_closest)
            elif binary_l[gene_place] == 1: # there is no neighbor effector but this is an effector
                return len(contig)
            else: # this is not an effector nor it has a neighbor effector
                return 5000


def effectors_in_neighbors(locus,k_neighbors):
    for i in range(len(all_locuses_in_genome)):
        if locus in all_locuses_in_genome[i]:
            contig = all_locuses_in_genome[i]
            binary_l = binary_lsts[i]
            gene_place = contig.index(locus)
            if len(contig)-1 <= k_neighbors:
                return sum(binary_l[:gene_place])+sum(binary_l[gene_place:])
            
            if gene_place >= k_neighbors: # enough genes upstream
                down = sum(binary_l[gene_place-k_neighbors:gene_place])
            else:
                more_to_count = k_neighbors-gene_place # to add from the end which is further upstream to the start (circular genome)
                down = sum(binary_l[:gene_place])+sum(binary_l[-more_to_count:])
            if len(binary_l)-gene_place-1 >= k_neighbors: # enough genes downstream
                up = sum(binary_l[gene_place+1:gene_place+1+k_neighbors])
            else:
                more_to_count = k_neighbors-(len(binary_l)-gene_place-1) # to add from the start which is further downstream to the end
                up = sum(binary_l[gene_place+1:])+sum(binary_l[:more_to_count])
            return down+up

with open(f'genome_organization_features.csv','w',newline='') as f:
    csv_writer = csv.writer(f)
    header=['locus','distance_from_closest_effector']
    for k in [5,10,15,20,25,30]:
        header.append(f'effectors_in_closest_{str(k)}_ORFs')
    csv_writer.writerow(header)
    for locus_dic in locus_dics:
        for locus in locus_dic:
            l=[locus]
            l.append(closest_effector(locus))
            for k in [5,10,15,20,25,30]:
                l.append(effectors_in_neighbors(locus,k))
            csv_writer.writerow(l)
endfile = open('genome_organization_features.done','w')
endfile.close()