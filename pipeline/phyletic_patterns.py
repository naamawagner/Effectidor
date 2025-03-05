import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sys import argv

working_directory = argv[1]


def create_phyletic_patterns(in_file, phyletic_csv, phyletic_text, gene_order, start_col=4, delimiter=',',
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


def create_presence_absence_map(phyletic_csv, x_label, out_fig_path, colors=None,
                                base_font_size=10, col_cluster=True):
    if colors is None:
        colors = ['silver', 'lightseagreen']
    df = pd.read_csv(phyletic_csv)
    df.set_index(df.columns[0], inplace=True)

    # Determine figure dimensions based on data dimensions
    n_rows, n_cols = df.shape
    fig_width = max(n_cols * 0.5, 10)
    fig_height = max(n_rows * 0.5, 10)

    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
    clustermap = sns.clustermap(df, cmap=cmap, figsize=(fig_width, fig_height),
                                tree_kws={"linewidths": 0.}, col_cluster=col_cluster,
                                vmin=0, vmax=1)
    # clustermap = sns.clustermap(df, cmap=cmap, figsize=(fig_width, fig_height),
    #                            tree_kws={"linewidths": 0.}, cbar_pos=(0.01, 0.8, 0.01, 0.1))
    clustermap.ax_heatmap.xaxis.set_ticks_position('top')
    clustermap.ax_heatmap.xaxis.set_label_position('top')

    plt.setp(clustermap.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha='left', rotation_mode='anchor')
    clustermap.ax_heatmap.yaxis.set_ticks_position('left')
    clustermap.ax_heatmap.yaxis.set_label_position('left')

    # Customize colorbar
    colorbar = clustermap.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([0.25, 0.75])
    colorbar.set_ticklabels(['Absent', 'Present'])
    tick_font_size = clustermap.ax_heatmap.xaxis.get_majorticklabels()[0].get_fontsize()

    clustermap.ax_heatmap.xaxis.labelpad = 20 + tick_font_size * 2  # Add padding to the x-axis label to ensure axis
    # label does not overlap with ticks' labels
    clustermap.ax_heatmap.yaxis.labelpad = 20 + tick_font_size * 2  # Add padding to the y-axis label to ensure axis
    # label does not overlap with ticks' labels
    # Adjust font size based on figure size, ensuring it does not go below the tick_font_size
    font_size = max(base_font_size * min(fig_width, fig_height) / 10, base_font_size, tick_font_size)

    for label in colorbar.ax.yaxis.get_ticklabels():
        label.set_fontsize(font_size)

    clustermap.ax_heatmap.set_xlabel(x_label, fontsize=font_size)
    clustermap.ax_heatmap.set_ylabel('Genome', fontsize=font_size)

    plt.tight_layout()
    plt.savefig(out_fig_path)


# T3Es:

in_f_T3Es = os.path.join(working_directory, 'out_learning',
                         'consensus_predictions_with_annotations_and_ortho_table.csv')
phyletic_csv_T3Es = os.path.join(working_directory, 'PresenceAbsence_T3Es.csv')
phyletic_text_T3Es = os.path.join(working_directory, 'Effectors_phyletic_pattern.txt')
gene_order_T3Es = os.path.join(working_directory, 'T3Es_order_in_PhyleticPattern.txt')
create_phyletic_patterns(in_f_T3Es, phyletic_csv_T3Es, phyletic_text_T3Es, gene_order_T3Es, T3Es=True)
create_presence_absence_map(phyletic_csv_T3Es, 'T3E (ortholog group)', os.path.join(working_directory,
                                                                                    'T3Es_presence_absence.png'))

# T3SS:

in_file_T3SS = os.path.join(working_directory, 'T3SS.csv')
phyletic_csv_T3SS = os.path.join(working_directory, 'PresenceAbsence_T3SS.csv')
phyletic_text_T3SS = os.path.join(working_directory, 'T3SS_phyletic_pattern.txt')
gene_order_T3SS = os.path.join(working_directory, 'T3SS_order_in_PhyleticPattern.txt')
create_phyletic_patterns(in_file_T3SS, phyletic_csv_T3SS, phyletic_text_T3SS, gene_order_T3SS, start_col=1)
create_presence_absence_map(phyletic_csv_T3SS, 'T3SS component', os.path.join(working_directory,
                                                                              'T3SS_presence_absence.png'),
                            col_cluster=False)

# chaperones:

in_f_chaperones = os.path.join(working_directory, 'chaperones.csv')

df = pd.read_csv(in_f_chaperones)
n_rows, n_cols = df.shape
if n_rows > 1 and n_cols > 1:

    phyletic_csv_chaperones = os.path.join(working_directory, 'PresenceAbsence_chaperones.csv')
    phyletic_text_chaperones = os.path.join(working_directory, 'chaperones_phyletic_pattern.txt')
    gene_order_chaperones = os.path.join(working_directory, 'chaperones_order_in_PhyleticPattern.txt')
    create_phyletic_patterns(in_f_chaperones, phyletic_csv_chaperones, phyletic_text_chaperones, gene_order_chaperones,
                             start_col=1)
    create_presence_absence_map(phyletic_csv_chaperones, 'chaperone', os.path.join(working_directory,
                                                                                   'chaperones_presence_absence.png'))
