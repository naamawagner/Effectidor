import os
import csv


def parse_gff_for_MGE(gff_path: str) -> tuple[dict, dict, dict]:
    mobile_element = {'integrase', 'recombinase', 'transposase'}
    exclude_element = {'reca', 'recb', 'recc', 'recd', 'rada', 'ruvc', 'recf', 'recbcd',
                       'uvrd'}
    # Talk to tal, pretty much Rec a-z should be excluded,
    # but it's not the official names for most of them.
    mobile_genetic_elements = dict()
    genomic_components = dict()
    orfs = dict()

    with open(gff_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue
            parts: list = line.strip().split('\t')
            if len(parts) > 2:
                attributes_dict: dict = {key_value.split('=')[0]: key_value.split('=')[1] for
                                         key_value in parts[-1].split(';') if '=' in key_value}
                sequence_id: str = parts[0]
                if parts[2] == 'CDS':
                    product: str = attributes_dict.get('product', '').lower()
                    locus_tag: str = attributes_dict.get('locus_tag', 'Unknown')
                    if any(element in product for element in mobile_element) and not any(
                            non_mge in product for non_mge in exclude_element):
                        if sequence_id not in mobile_genetic_elements:
                            mobile_genetic_elements[sequence_id] = {}
                        mobile_genetic_elements[sequence_id][locus_tag] = {
                            'MGE': product,
                            'start': int(parts[3]),
                            'end': int(parts[4])
                        }
                    if sequence_id not in orfs:
                        orfs[sequence_id] = dict()
                    orfs[sequence_id][locus_tag] = {'start': int(parts[3]), 'end': int(parts[4])}

                if sequence_id not in genomic_components:  # for new genomic component we need to capture its length and circular status
                    length = int(parts[4])
                    is_circular = bool(attributes_dict.get('Is_circular', 'false').capitalize())
                    genomic_components[sequence_id] = {'length': length, 'is_circular': is_circular}

    return mobile_genetic_elements, orfs, genomic_components


def find_closest_mge(output_path: str, dict_orfs: dict, genomic_components: dict,
                     sequence_id: str, mobile_genetic_elements: dict):
    full_output_path = ''.join([output_path[:-4], '_full.csv'])
    closest_mge_all = {}
    genome_length: int = genomic_components[sequence_id]['length']
    is_circular: bool = genomic_components[sequence_id]['is_circular']
    if sequence_id not in mobile_genetic_elements:
        for orf_locus in dict_orfs[sequence_id].keys():
            closest_mge = (orf_locus, None, None, None)
            closest_mge_all[orf_locus]: tuple = closest_mge
    else:
        for orf_locus in dict_orfs[sequence_id].keys():
            orf: dict = dict_orfs[sequence_id][orf_locus]
            orf_start: int = orf['start']
            orf_end: int = orf['end']

            min_distance = float('inf')
            closest_mge = (orf_locus, None, None, None)

            for mge_locus_tag, mge in mobile_genetic_elements[sequence_id].items():
                mge_start: int = mge['start']
                mge_end: int = mge['end']
                if mge_locus_tag != orf_locus:
                    distance1: int = abs(orf_end - mge_start)
                    distance2: int = abs(mge_end - orf_start)
                    if is_circular:
                        wrap_around_distance1: int = genome_length - distance1
                        wrap_around_distance2: int = genome_length - distance2
                        current_distance: int = min(distance1, wrap_around_distance1, distance2, wrap_around_distance2)
                    else:
                        current_distance: int = min(distance1, distance2)

                    if current_distance < min_distance:
                        min_distance = current_distance
                        closest_mge = (orf_locus, mge_locus_tag, current_distance, mge['MGE'])
                else:
                    continue

            closest_mge_all[orf_locus]: tuple = closest_mge

    file_exists: bool = os.path.exists(output_path)

    with open(output_path, 'a', newline='') as csvfile:
        if not file_exists:
            fieldnames: list[str] = ['locus', 'distance_to_mobile_genetic_element']
            csv.writer(csvfile).writerow(fieldnames)

        for closest_mge in closest_mge_all.values():
            csv.writer(csvfile).writerow((closest_mge[0], closest_mge[2]))

    full_file_exists: bool = os.path.exists(full_output_path)

    with open(full_output_path, 'a', newline='') as full_csvfile:
        if not full_file_exists:
            all_fieldnames: list[str] = ['locus', 'MGE_locus_tag', 'distance', 'MGE']
            csv.writer(full_csvfile).writerow(all_fieldnames)

        for closest_mge in closest_mge_all.values():
            csv.writer(full_csvfile).writerow(closest_mge)


def main(gff_path: str, output_path: str):
    mobile_elements: dict
    dict_orfs: dict
    genomic_components: dict
    mobile_elements, dict_orfs, genomic_components = parse_gff_for_MGE(gff_path)
    for sequence_id in dict_orfs.keys():
        find_closest_mge(output_path, dict_orfs, genomic_components, sequence_id, mobile_elements)
    endfile = open(f'{output_path}.done', 'w')
    endfile.close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Process GFF files and output ORF-MGE proximity data.")

    # Define arguments
    parser.add_argument('gff_file', help="Path to the GFF file",
                        type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('output_path', help="Output CSV file path for locus_tag and distance")

    args = parser.parse_args()

    main(args.gff_file, args.output_path)
