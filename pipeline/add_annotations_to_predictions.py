from Bio import SeqIO
import csv
import pandas as pd
import os


def create_annotations_f(ORFs_f, gff_f='', out_f='', pseudo_f=''):
    '''
    The annotation is extracted from the gff file if given and otherwise from the ORFs fasta file.
    It is then saved to the out_f (csv format)
    '''
    locus_annotation = {}
    pseudo_genes = []
    if gff_f:
        with open(gff_f) as gff:
            for line in gff:
                if not line.startswith('#'):
                    row = line.strip().split('\t')
                    if len(row) == 9:
                        if row[2] == 'CDS':
                            locus, annot = '', ''
                            features = row[-1].split(';')
                            for fe in features:
                                if 'locus_tag=' in fe:
                                    locus = fe.split('=')[1]
                                elif 'product=' in fe:
                                    annot = fe.split('=')[1]
                            if 'pseudo=true' in row[-1]:
                                pseudo_genes.append(locus)
                                # annot = 'pseudogene ' + annot
                            locus_annotation[locus] = annot
    else:  # get the annotations from the fasta file
        recs = SeqIO.parse(ORFs_f, 'fasta')
        for rec in recs:
            is_locus = False
            annotation = ''
            header = rec.description
            header_l = header.split('[')  # split by '[' and not by spaces as the annotation fields contain spaces
            for a in header_l:
                if 'locus_tag=' in a:
                    locus = a.split('=')[1].strip('] ')
                    is_locus = True
                elif 'protein=' in a:
                    annotation = a.split('=')[1].strip('] ')

            if is_locus:
                locus_annotation[locus] = annotation
                if 'pseudo=true' in header:
                    pseudo_genes.append(locus)
            else:
                locus_annotation[rec.id] = annotation
                if 'pseudo=true' in header:
                    pseudo_genes.append(rec.id)
    with open(out_f, 'w', newline='') as out:
        writer = csv.writer(out)
        header = ['locus', 'annotation']
        writer.writerow(header)
        for locus in locus_annotation:
            writer.writerow([locus, locus_annotation[locus]])
    with open(pseudo_f, 'w') as pseudo:
        pseudo.write('\n'.join(pseudo_genes))


def add_annotations_to_predictions(in_f, out_f_normal, out_f_pseudo, annotations_fasta,
                                   out_f_T3SS, gff_d='', line_end='\n'):
    locus_annotation = {}
    if gff_d:
        for f in os.listdir(gff_d):
            gff_f = os.path.join(gff_d, f)
            with open(gff_f) as gff:
                for line in gff:
                    if not line.startswith('#'):
                        locus, annot = '', ''
                        row = line.strip().split('\t')
                        if len(row) == 9:
                            if row[2] == 'CDS':
                                features = row[-1].split(';')
                                for fe in features:
                                    if 'locus_tag=' in fe:
                                        locus = fe.split('=')[1]
                                    elif 'product=' in fe:
                                        annot = fe.split('=')[1]
                                locus_annotation[locus] = annot
                        
    recs = SeqIO.parse(annotations_fasta, 'fasta')
    T3SS_recs = SeqIO.to_dict(SeqIO.parse('T3SS_proteins.faa', 'fasta'))
    T3_hits_f = 'T3SS_hits.csv'
    T3_d = {}
    with open(T3_hits_f) as f:
        reader = csv.reader(f)
        for row in reader:
            T3_d[row[0]] = row[1].replace('\t', ', ')
    locus_prot = {}
    for rec in recs:
        dna_seq = rec.seq
        protein = str(dna_seq.translate(to_stop=True))
        protein_l = [protein[70*i:70*(i+1)] for i in range(len(protein)//70+1)]
        if protein_l[-1] == '':
            protein_l = protein_l[:-1]
        protein_n = line_end.join(protein_l)
        is_locus = False
        annotation = ''
        header = rec.description
        if 'pseudo=true' in header:
            annotation = 'pseudogene'
        header_l_for_prot = header.split('[')
        header_l_for_locus = header.split()
        for a in header_l_for_locus:
            if 'locus_tag=' in a:
                locus = a.split('=')[1].strip(']')
                is_locus = True
        for a in header_l_for_prot:
            if 'protein=' in a:
                if annotation == 'pseudogene':
                    annotation += ' '+a.split('=')[1].strip().strip(']')
                else:
                    annotation = a.split('=')[1].strip().strip(']')
        if is_locus:
            if locus not in locus_annotation:
                locus_annotation[locus] = annotation
            locus_prot[locus] = protein_n
        else:
            if rec.id not in locus_annotation:
                locus_annotation[rec.id] = annotation
            locus_prot[rec.id] = protein_n
    with open(in_f) as f:
        reader = csv.reader(f)
        with open(out_f_normal, 'w', newline='') as out_normal:
            writer_normal = csv.writer(out_normal)
            with open(out_f_pseudo, 'w', newline='') as out_pseudo:
                writer_pseudo = csv.writer(out_pseudo)
                with open(out_f_T3SS, 'w', newline='') as out_T3SS:
                    writer_T3SS = csv.writer(out_T3SS)
                    header = next(reader)
                    header.append('Annotation')
                    header.append('Protein sequence')
                    writer_normal.writerow(header)
                    writer_pseudo.writerow(header)
                    writer_T3SS.writerow(['Locus tag']+header[-2:]+['Hits to T3SS proteins (see Data)'])
                    for row in reader:
                        row.append(locus_annotation[row[0]])
                        row.append(locus_prot[row[0]])
                        if row[0] in T3SS_recs:
                            writer_T3SS.writerow([row[0]]+row[-2:]+[T3_d[row[0]]])
                        elif 'pseudogene' in row[-2]:
                            writer_pseudo.writerow(row)
                        else:
                            writer_normal.writerow(row)


def make_html_tables(predictions_f, effector_homolog):
    data = pd.read_csv(predictions_f)
    predicted = data[data.is_effector == '?']
    positives = data[data.is_effector == 'yes']
    header = list(data.columns)
    header[1] = 'Score'
    predicted.columns = header
    positives.columns = header
    predicted = predicted.drop(columns=['is_effector'])
    positives = positives.drop(columns=['is_effector'])
    effector_homolog_df = pd.read_csv(effector_homolog)
    predicted_table = predicted.head(10).merge(effector_homolog_df, how='left').to_html(
        index=False, justify='left', escape=False)
    positives_table = positives.merge(effector_homolog_df, how='left').to_html(
        index=False, justify='left', escape=False)
    return predicted_table, positives_table
    