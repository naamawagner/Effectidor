from Bio import SeqIO

def parse_ORFs(ORFs_fasta_file,DNA=True):
    ORFs_dic = {}
    recs = SeqIO.parse(ORFs_fasta_file,'fasta')
    if DNA:
        for rec in recs:
            seq = rec.seq
            if len(seq.translate(to_stop=True))>10:
                header = rec.description
                header_l = header.split()
                is_locus = False
                for field in header_l:
                    if "locus_tag" in field:
                        locus = field.split('=')[1].strip(']')
                        ORFs_dic[locus] = seq
                        is_locus = True
                        break
                if not is_locus:
                    ORFs_dic[rec.id] = seq
    else:
        for rec in recs:
            seq = rec.seq
            header = rec.description
            header_l = header.split()
            is_locus = False
            for field in header_l:
                if "locus_tag" in field:
                    locus = field.split('=')[1].strip(']')
                    ORFs_dic[locus] = seq
                    is_locus = True
                    break
            if not is_locus:
                ORFs_dic[rec.id] = seq
    return ORFs_dic