import fasta_parser
import csv
from sys import argv
import os
from pip_box_features import parse_gff_to_CDS_loc
ORFs_f = argv[1]
effectors_file = argv[2]
working_directory = argv[3]
gff_f = argv[4]
os.chdir(working_directory)
locus_dic=fasta_parser.parse_ORFs(ORFs_f)
effectors_dict = fasta_parser.parse_ORFs(effectors_file,DNA=False)
#%%

locus_area_d,circulars = parse_gff_to_CDS_loc(gff_f,locus_dic)

    
linear_contigs = []
circular_contigs = []
for region in locus_area_d:
    genes_in_order = sorted(locus_area_d[region].keys(),key=lambda k:locus_area_d[region][k][0])
    if region in circulars:
        circular_contigs.append(genes_in_order)
    else:
        linear_contigs.append(genes_in_order)
        
#all_locuses_in_genome=[list(locus_dic.keys()) for locus_dic in locus_dics]
binary_linear_lsts=[] # effectors locations list
for contig in linear_contigs:
    binary_l = []
    for locus in contig:
        if locus in effectors_dict:
            binary_l.append(1)
        else:
            binary_l.append(0)
    binary_linear_lsts.append(binary_l)
    
binary_circular_lsts=[] # effectors locations list
for contig in circular_contigs:
    binary_l = []
    for locus in contig:
        if locus in effectors_dict:
            binary_l.append(1)
        else:
            binary_l.append(0)
    binary_circular_lsts.append(binary_l)

def closest_effector(locus):
    for i in range(len(circular_contigs)):
        if locus in circular_contigs[i]: #circular
            contig = circular_contigs[i]
            binary_l = binary_circular_lsts[i]
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
            else: # it has no neighbor effector
                return None #it will be filled with the median of this feature later (before learning)
    for i in range(len(linear_contigs)):#linear
        if locus in linear_contigs[i]:
            contig = linear_contigs[i]
            binary_l = binary_linear_lsts[i]
            gene_place = contig.index(locus)
            upstream_closest = 0
            for i in range(1,len(contig)-gene_place):
                if binary_l[gene_place+i]==1:
                    upstream_closest = i
                    break
            downstream_closest = 0
            for i in range(1,gene_place):
                if binary_l[gene_place-i] == 1:
                    downstream_closest = i
                    break
            if upstream_closest != 0 and downstream_closest != 0: # there are neighbor effectors both upstream and downstream
                return min(upstream_closest,downstream_closest)
            elif upstream_closest != 0: # there is a neighbor effector upstream
                return upstream_closest
            elif downstream_closest != 0: # there is a neighbor effector downstream
                return downstream_closest
            else: # it has no neighbor effector
                return None #it will be filled with the median of this feature later (before learning)


def effectors_in_neighbors(locus,k_neighbors):
    for i in range(len(circular_contigs)):
        if locus in circular_contigs[i]: #circular
            contig = circular_contigs[i]
            binary_l = binary_circular_lsts[i]
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
        
    for i in range(len(linear_contigs)):#linear
        if locus in linear_contigs[i]:
            contig = linear_contigs[i]
            binary_l = binary_linear_lsts[i]
            gene_place = contig.index(locus)
            if len(contig)-1 <= k_neighbors:
                return sum(binary_l[:gene_place])+sum(binary_l[gene_place:])
            if gene_place >= k_neighbors: # enough genes upstream
                down = sum(binary_l[gene_place-k_neighbors:gene_place])
            else:
                down = sum(binary_l[:gene_place])
            if len(binary_l)-gene_place-1 >= k_neighbors: # enough genes downstream
                up = sum(binary_l[gene_place+1:gene_place+1+k_neighbors])
            else:
                up = sum(binary_l[gene_place+1:])
            return down+up
                

with open(f'genome_organization_features.csv','w',newline='') as f:
    csv_writer = csv.writer(f)
    header=['locus','distance_from_closest_effector']
    for k in [5,10,15,20,25,30]:
        header.append(f'effectors_in_closest_{str(k)}_ORFs')
    csv_writer.writerow(header)
    for locus in locus_dic:
        l=[locus]
        l.append(closest_effector(locus))
        for k in [5,10,15,20,25,30]:
            l.append(effectors_in_neighbors(locus,k))
        csv_writer.writerow(l)
endfile = open('genome_organization_features.done','w')
endfile.close()