# Effectidor
## https://effectidor.tau.ac.il/
A pipeline to predict Type III effectors in bacterial genomes, using machine-learning.

This repository contains a stand-alone version for your use.

### Requirements
Clone the repository, create a conda environment using env.yml, and activate it.

In a designated running directory (<b><i>output_dir_path</i></b>) save the required input:
1. <b><i>ORFs</i></b> - ORFs fasta files; fasta files with all the open reading frames (ORFs) in the genome, one file per genome. Multiple files (representing multiple genomes) should be archived in a single zip, not inside directories.
2. <b><i>OG_table</i></b> - Orthologs table; the table should be in csv format, such that every column holds the information of a different genome, and every row of an ortholog group. Paralogs should be separated by a semi column (;). The identifier of the genes in the ortholog groups are the ORFs IDs, or, if downloaded from NCBI, the locus_tag.

Additional input may include:
1. <b><i>input_effectors_path</i></b> - effectors data; a subset of the ORFs input, representing the known effectors from the genome. Up to one file per genome. For multiple genome analysis, these files must be zipped, even if a file is supplied just for one of the genomes. The names of the files must match the corresponding ORFs input.
2. <b><i>gff_path</i></b> - GFF annotation file(s). For computing genome organization features, a GFF annotation file should be supplied for each of the genomes. These files should be zipped (not inside directories) and the names of the files should match the corresponding ORFs files.
3. <b><i>genome_path</i></b> - Full genome file(s). For searching for regulatory elements (as detailed in the usage) in the promoter of the genes, in addition to the GFF input, you need to supply full genome fasta file(s). For multiple genome analysis the files should be zipped (not inside directories) and the names of the files should match the corresponding ORFs files.
4. <b><i>input_T3Es_path</i></b> - Your own proven T3Es data from any bacteria, in protein fasta file. This input will be used in addition to Effectidor's T3Es data for homology searches.
5. <b><i>host_proteome_path</i></b> - Host proteome, in protein fasta format. This input must be zipped (not inside directories) and thus may include multiple files.
6. <b><i>no_T3SS</i></b> - Proteomes of closely related bacteria without the T3SS. A protein fasta file per representative bacterium, all files must be zipped together (not inside directories).
### Usage
After activating the environment, run the following script:

python pipeline/main_T3Es.py <b><i>ORFs</i></b> <b><i>output_dir_path</i></b> <b><i>OG_table</i></b> [--input_effectors_path <b><i>input_effectors_path</i></b>] [--input_T3Es_path <b><i>input_T3Es_path</i></b>] [--host_proteome_path <b><i>host_proteome_path</i></b>] [--no_T3SS <b><i>no_T3SS</i></b>] [--genome_path <b><i>genome_path</i></b>] [--gff_path <b><i>gff_path</i></b>]

Additional optional flags:

--homology_search - if you supplied the <b><i>input_effectors_path</i></b> input, this flag will indicate that this input is partial (per genome) and that additional T3Es should be added based on homology to the T3Es data, if found in the genome.

--translocation_signal - to add a prediction of the type III secretion signal using a pre-trained embedding model. Increases the accuracy with a price of longer running time.

One of the following regulatory elements (that will be searched in the promoters, if you use a GFF+full genome inputs):
1. --PIP (for PIP-box)
2. --hrp (for Hrp-box)
3. --mxiE (for MxiE-box)
4. --exs (for exs-box)
5. --tts (for tts-box)

### Example data
in the <i>Example</i> directory you can find example data. To run a pan-genome analysis, on all the reference genomes of <i>Xanthomonas</i>, run the following command:<br><br>
python pipeline/main_T3Es.py Example/Pan_genome/ORFs.zip Example/Pan_genome/ Example/Pan_genome/OG_table.csv --genome_path Example/Pan_genome/genome.zip --gff_path Example/Pan_genome/GFF.zip --PIP<br>
this command will include genome organization and regulatory elements (PIP-box) features in the analysis. To also include the translocation signal feature, add the flag --translocation_signal.

