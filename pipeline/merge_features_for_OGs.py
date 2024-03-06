from sys import argv
import os
import csv
import pandas as pd
import math

working_directory = argv[1]
os.chdir(working_directory)

Effectidor_features_d = os.path.join(working_directory,'Effectidor_runs')

def get_ortho_dict(wd,ortho_f='clean_orthologs_table.csv'):
    '''from the output file of Microbializer, present in the working directory,
    this function creates a dictionary that maps in every genome the locuses to OGs'''
    genomes_orthogroup_dict = {}
    with open(os.path.join(wd,ortho_f)) as ortho_table:
        ortho_reader = csv.reader(ortho_table)
        header = next(ortho_reader)
        for row in ortho_reader:
            OG = row[0]
            for i in range(1,len(header)):
                genome = header[i]
                if genome not in genomes_orthogroup_dict:
                    genomes_orthogroup_dict[genome]={}
                for locus in row[i].split(';'):
                    genomes_orthogroup_dict[genome][locus]=OG
    return genomes_orthogroup_dict
   
genomes_orthogroup_dict = get_ortho_dict(working_directory)

# prepare the orthologs table for final merge:
with open(os.path.join(working_directory, 'clean_orthologs_table.csv')) as in_f:
    full_content = in_f.read()
if os.path.exists(Effectidor_features_d):
    genomes = os.listdir(Effectidor_features_d)
    for genome in genomes:
        pseudogenes_path = os.path.join(Effectidor_features_d,genome,'pseudogenes.txt')
        with open(pseudogenes_path) as pseudo_f:
            pseudogenes = pseudo_f.read().split('\n')
            for pseudo in pseudogenes:
                full_content = full_content.replace(f'{pseudo},',f'{pseudo}(pseudogene),')
                full_content = full_content.replace(f'{pseudo};', f'{pseudo}(pseudogene);')
else:
    with open('pseudogenes.txt'):
        pseudogenes = pseudo_f.read().split('\n')
        for pseudo in pseudogenes:
            full_content = full_content.replace(f'{pseudo},', f'{pseudo}(pseudogene),')
            full_content = full_content.replace(f'{pseudo};', f'{pseudo}(pseudogene);')
full_content = full_content.replace('OG_name,','OG,')
with open('clean_orthologs_table_with_pseudo.csv','w') as out_f:
    out_f.write(full_content)

def combine_all_genomes_data(out_f_path,in_f_name):
    with open(out_f_path,'w',newline='') as out_path:
        writer = csv.writer(out_path)
        if os.path.exists(Effectidor_features_d):
            genomes = os.listdir(Effectidor_features_d)
            first_genome = genomes[0]
            with open(os.path.join(Effectidor_features_d, first_genome, f'{in_f_name}.csv')) as in_f:
                reader = csv.reader(in_f)
                header = ['OG'] + next(reader)
                writer.writerow(header)
                for row in reader:
                    if row[0] in genomes_orthogroup_dict[first_genome]:
                        row = [genomes_orthogroup_dict[first_genome][row[0]]] + row
                        writer.writerow(row)
            for genome in genomes[1:]:
                with open(os.path.join(Effectidor_features_d, genome, f'{in_f_name}.csv')) as in_f:
                    reader = csv.reader(in_f)
                    next(reader)
                    for row in reader:
                        if row[0] in genomes_orthogroup_dict[genome]:
                            row = [genomes_orthogroup_dict[genome][row[0]]] + row
                            writer.writerow(row)
        else:
            with open(f'{in_f_name}.csv') as in_f:
                reader = csv.reader(in_f)
                header = ['OG'] + next(reader)
                writer.writerow(header)
                for row in reader:
                    if row[0] in genomes_orthogroup_dict['genome_ORFs']:
                        row = [genomes_orthogroup_dict['genome_ORFs'][row[0]]] + row
                        writer.writerow(row)

combine_all_genomes_data('full_data.csv','features')
combine_all_genomes_data('full_OGs_annotations.csv','annotations')

annot_df = pd.read_csv('full_OGs_annotations.csv')
grouped_annot = annot_df.groupby(['OG'])['annotation'].agg(pd.Series.mode)
grouped_annot.to_csv('OGs_annotations.csv')

def label(iterable_arg):
    '''define a label of an OG, based on the labels of its members'''
    if any(iterable_arg=='effector'):
        return 'effector'
    elif any(iterable_arg=='no'):
        return 'no'
    else:
        return '?'   
        

df = pd.read_csv('full_data.csv')
features = list(df.columns[2:-1])
#defining manipulation per feature in the transformation to OGs
minimum = [feature for feature in features if 'distance_from_closest_effector' in feature]
maximum = [feature for feature in features if ('T3_signal' in feature or '_box' in feature)]
for feature in maximum:
    features.remove(feature)
for feature in minimum:
    features.remove(feature)

#%% transformation to OGs
grouped = df.groupby('OG').agg({**{feature:['mean'] for feature in features},\
                    **{feature:['mean','min','median'] for feature in minimum},\
                    **{feature:['max'] for feature in maximum},\
                    **{'is_effector':label}})

#%% adding the "similarity_to_effectors_vs_non_effectors" features for the OGs
aa_freqs = [feature for feature in features if 'full_protein' in feature]
effectors = grouped[grouped.is_effector.label=='effector'][aa_freqs]
non_effectors = grouped[grouped.is_effector.label=='no'][aa_freqs]
effectors_mean = effectors.mean()
non_effectors_mean = non_effectors.mean()
    
def aa_profile(OG):
    OG_aa_freqs = grouped.loc[OG][aa_freqs]
    if OG in effectors.index:
        effector_1 = effectors.drop(index=OG)
        effector_1_mean = effector_1.mean()
        dis_to_eff = math.sqrt(sum([(a-b)**2 for a,b in zip(OG_aa_freqs.values,effector_1_mean.values)]))
    else:
        dis_to_eff = math.sqrt(sum([(a-b)**2 for a,b in zip(OG_aa_freqs.values,effectors_mean.values)]))
    if OG in non_effectors.index:
        non_effector_1 = non_effectors.drop(index=OG)
        non_effector_1_mean = non_effector_1.mean()
        dis_to_non_eff = math.sqrt(sum([(a-b)**2 for a,b in zip(OG_aa_freqs.values,non_effector_1_mean.values)]))
    else:
        dis_to_non_eff = math.sqrt(sum([(a-b)**2 for a,b in zip(OG_aa_freqs.values,non_effectors_mean.values)]))
    
    return dis_to_eff - dis_to_non_eff

OGs = list(grouped.index)
aa_profiles_dis =[aa_profile(OG) for OG in OGs]
grouped.insert(len(grouped.columns)-1,'similarity_to_effectors_vs_non_effectors',aa_profiles_dis)
 
#%% flatten hirarchical columns and indices to a simple features table
updated_features = grouped.reset_index()
updated_features.columns = ['_'.join(col) for col in updated_features.columns.values]
updated_features.columns = [col.replace('label','_').strip('_') for col in updated_features.columns]
updated_features.to_csv('OGs_features.csv',index=False)

