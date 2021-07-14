from Bio import SeqIO
import csv
import pandas as pd

def add_annotations_to_predictions(in_f,out_f_normal,out_f_pseudo,annotations,line_end='\n'):
    recs = SeqIO.parse(annotations,'fasta')
    locus_annotation={}
    locus_prot = {}
    for rec in recs:
        dna_seq = rec.seq
        protein = str(dna_seq.translate(to_stop=True))
        protein_l = [protein[70*i:70*(i+1)] for i in range(len(protein)//70+1)]
        if protein_l[-1]=='':
            protein_l = protein_l[:-1]
        protein_n = line_end.join(protein_l)
        is_locus = False
        annotation = ''
        header = rec.description
        if 'pseudo=true' in header:
            annotation = 'pseudogene'
        header_l = header.split('[')
        for a in header_l:
            if 'locus_tag=' in a:
                locus = a.split('=')[1].strip().strip(']')
                is_locus = True
            elif 'protein=' in a:
                if annotation == 'pseudogene':
                    annotation += ' '+a.split('=')[1].strip().strip(']')
                else:
                    annotation = a.split('=')[1].strip().strip(']')
        if is_locus:
            locus_annotation[locus] = annotation
            locus_prot[locus] = protein_n
        else:
            locus_annotation[rec.id] = annotation
            locus_prot[rec.id] = protein_n
    with open(in_f) as f:
        reader = csv.reader(f)
        with open(out_f_normal,'w',newline='') as out_normal:
            writer_normal = csv.writer(out_normal)
            with open(out_f_pseudo,'w',newline='') as out_pseudo:
                writer_pseudo = csv.writer(out_pseudo)
                header = next(reader)
                header.append('Annotation')
                header.append('Protein sequence')
                writer_normal.writerow(header)
                writer_pseudo.writerow(header)
                for row in reader:
                    row.append(locus_annotation[row[0]])
                    row.append(locus_prot[row[0]])
                    if 'pseudogene' in row[-2]:
                        writer_pseudo.writerow(row)
                    else:
                        writer_normal.writerow(row)


def make_html_tables(predictions_f):
    data = pd.read_csv(predictions_f)
    predicted = data[data.is_effector=='?']
    positives = data[data.is_effector=='yes']
    header = list(data.columns)
    header[1]='Score'
    header[0]='Locus tag'
    predicted.columns = header
    positives.columns = header
    predicted.drop(columns=['is_effector'], inplace=True)
    positives.drop(columns=['is_effector'], inplace=True)
    predicted_table = predicted.head(10).to_html(index=False,justify='left',escape=False)
    positives_table = positives.to_html(index=False,justify='left',escape=False)
    return predicted_table,positives_table
    