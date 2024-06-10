import pandas as pd
import os
import numpy as np
import csv
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sys import argv

working_directory = argv[1]


def create_effector_phyletic_patterns(in_file, phyletic_csv, phyletic_text, gene_order, start_col=4, delimiter=',', 
                                      T3Es=False):
    '''
    in_file: path to csv file with the data to transform to phyletic pattern, where the columns are genomes, and the
    rows are genes/OGs
    phyletic_csv: path to a csv file that will be created with the phyletic pattern. it will be used later for the
    presence/absense map
    phyletic_text: a path to the output phyletic pattern  that will be created, in txt format
    gene_order: a path to an output txt file that will hold a list of the gene order in the phyletic pattern
    start_col: position of the first column in the input file that represents a genome
    delimiter: delimiter to use to read the dataframe (, for csv and \t for tsv)
    T3Es: whether the input is T3Es prediction (this input requires dropping columns that do not exist elsewhere)
    '''
    df = pd.read_csv(in_file, sep=delimiter)    
    header = df.columns
    species = header[start_col:]
    df.set_index(header[0], inplace=True)
    df.dropna(how='all', subset=species, inplace=True)
    df.fillna('0', inplace=True)
    df[species] = np.where(df[species] != '0', '1', '0')  # put 1 wherever there's a locus
    if T3Es:
        predicted = df[(df['is_effector'] == '?') & (df['score'] > 0.5)][species]
        positives = df[df.is_effector == 'yes'][species]
        positives, predicted = positives.T, predicted.T
        phyletic_df = positives.merge(predicted, left_index=True, right_index=True)
    else:  # T3SS
        phyletic_df = df.T
    phyletic_df.to_csv(phyletic_csv)
    phyletic_sequence = phyletic_df.apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    with open(phyletic_text, 'w') as out_phyletic:
        for i in range(len(phyletic_sequence)):
            out_phyletic.write(f'>{phyletic_sequence.index[i]}\n')
            out_phyletic.write(f'{phyletic_sequence[i]}\n')
    with open(gene_order, 'w') as out_gene_order:
        out_gene_order.write(','.join(phyletic_df.columns))


# generating clustered presence/absence map
# working_directory = r'C:\Users\TalPNB2\Downloads'


def create_presence_absence_map(phyletic_csv, x_label, out_fig_path, colors=['silver', 'lightseagreen'], 
                                base_font_size=10):
    df = pd.read_csv(phyletic_csv)
    df.set_index(df.columns[0], inplace=True)
    fig_width, fig_height = df.shape[1]//2.5, df.shape[0]//2.5
    font_size = base_font_size * (min(fig_width, fig_height) / 10)  # Adjust font size based on figure size

    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
    clustermap = sns.clustermap(df, cmap=cmap, figsize=(fig_width, fig_height), dendrogram_ratio=0.1,
                                tree_kws={"linewidths": 0.}, cbar_pos=(0.01, 0.8, 0.01, 0.1))
    clustermap.ax_heatmap.xaxis.set_ticks_position('top')
    clustermap.ax_heatmap.xaxis.set_label_position('top')
    plt.setp(clustermap.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha='left', rotation_mode='anchor')
    clustermap.ax_heatmap.yaxis.set_ticks_position('left')
    clustermap.ax_heatmap.yaxis.set_label_position('left')
    colorbar = clustermap.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([0.25, 0.75])
    colorbar.set_ticklabels(['Absent', 'Present'])
    for label in colorbar.ax.yaxis.get_ticklabels():
        label.set_fontsize(font_size/2)
    # clustermap.ax_heatmap.set_title('T3Es presence/absence map across analyzed genomes', fontsize=font_size)
    clustermap.ax_heatmap.set_xlabel(x_label, fontsize=font_size)
    clustermap.ax_heatmap.set_ylabel('Genome', fontsize=font_size)
    clustermap.ax_heatmap.xaxis.labelpad = 25
    clustermap.ax_heatmap.yaxis.labelpad = 25
    # plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    plt.savefig(out_fig_path)


# T3Es:
in_f_T3Es = os.path.join(working_directory, 'out_learning', 'consensus_predictions_with_annotations_and_ortho_table.csv')
phyletic_csv_T3Es = os.path.join(working_directory, 'PresenceAbsence_T3Es.csv')
phyletic_text_T3Es = os.path.join(working_directory, 'Effectors_phyletic_pattern.txt')
gene_order_T3Es = os.path.join(working_directory, 'T3Es_order_in_PhyleticPattern.txt')
create_effector_phyletic_patterns(in_f_T3Es, phyletic_csv_T3Es, phyletic_text_T3Es, gene_order_T3Es, T3Es=True)
create_presence_absence_map(phyletic_csv_T3Es, 'T3E (ortholog group)', os.path.join(working_directory,
                                                                                    'T3Es_presence_absence.png'))
# T3SS:
in_file_T3SS = os.path.join(working_directory, 'T3SS.csv')
phyletic_csv_T3SS = os.path.join(working_directory, 'PresenceAbsence_T3SS.csv')
phyletic_text_T3SS = os.path.join(working_directory, 'T3SS_phyletic_pattern.txt')
gene_order_T3SS = os.path.join(working_directory, 'T3SS_order_in_PhyleticPattern.txt')
create_effector_phyletic_patterns(in_file_T3SS, phyletic_csv_T3SS, phyletic_text_T3SS, gene_order_T3SS, start_col=1)
create_presence_absence_map(phyletic_csv_T3SS, 'T3SS component', os.path.join(working_directory,
                                                                              'T3SS_presence_absence.png'))
