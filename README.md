# ExInAtor2

# Authors

Roberta Esposito<sup>\*1,2,3</sup>, Andrés Lanzós<sup>\*1,2,4</sup>, Taisia Polidori<sup>1,2</sup>, Hugo Guillen-Ramirez<sup>5,6</sup>,
Bernard Merlin<sup>1,2</sup>, Lia Mela<sup>1,2</sup>, Eugenio Zoni<sup>2,9</sup>, Isabel Büchi<sup>2,8</sup>, Lusine Hovhannisyan
<sup>2,7</sup>, Finn McCluggage<sup>10,11</sup>, Matúš Medo<sup>2,7</sup>, Giulia Basile<sup>1,2</sup>, Dominik F. Meise<sup>1,2</sup>,
Sunandini Ramnarayanan<sup>5,6</sup>, Sandra Zwyssig<sup>1,2</sup>, Corina Wenger<sup>1,2</sup>, Kyriakos Schwarz
<sup>1,2</sup>, Adrienne Vancura<sup>1,2</sup>, Nuria Bosch-Guiteras<sup>1,2,4</sup>, Marianna Kruithof-de Julio<sup>2,9</sup>, Yitzhak
Zimmer<sup>2,7</sup>, Michaela Medová<sup>2,7</sup>, Deborah Stroka<sup>2,8</sup>, Archa Fox<sup>10,11</sup>, Rory Johnson
<sup>1,2,5,6</sup>

1.Department of Medical Oncology, Inselspital, Bern University Hospital, University of Bern, 3010 Bern, Switzerland.

2.Department for BioMedical Research, University of Bern, 3008 Bern, Switzerland

3.Institute of Genetics and Biophysics "Adriano Buzzati-Traverso", CNR, 80131 Naples, Italy.

4.Graduate School of Cellular and Biomedical Sciences, University of Bern, 3012 Bern, Switzerland.

5.School of Biology and Environmental Science, University College Dublin, Dublin D04 V1W8, Ireland.

6.Conway Institute for Biomolecular and Biomedical Research, University College Dublin, Dublin D04 V1W8,
Ireland.

7.Department of Radiation Oncology, Inselspital, Bern University Hospital and University of Bern, Bern, Switzerland

8.University Clinic of Visceral Surgery and Medicine, Bern University Hospital, Inselspital, Department of Biomedical Research, University of Bern, Bern, Switzerland.

9.Department of Urology, Inselspital, Bern University Hospital, Bern, Switzerland.

10.School of Molecular Sciences, University of Western Australia, Crawley, Western Australia, Australia.

11.School of Human Sciences, University of Western Australia, Crawley, Western Australia, Australia.

\*Equal contribution

# Description 

ExInAtor2 is an improvement of ExinAtor that can be [found here](https://github.com/alanzos/ExInAtor).

ExInAtor2 detects positive selection using the following two signatures:

1. Recurrence (RE) that compares the exonic mutation rate to that of the local background

2. Functional impact (FI), that compares estimated functional impact of mutations to background, both in exonic regions. 

# Requirements 

Python
---
Tested on python version 3.6.13/3.7.2

*Python packages required:* regex (tested on version 2016.06.24), numpy (tested on version 1.19.5), pytabix (required only for functional impact, tested on version 0.1)

R
---
Tested on R version 3.5.1

Bedtools
---
Tested on Bedtools version 2.29.2/2.26.0

EPACTS (for functional impact only)
---
Tested on version 3.4.2

# Installation 

Estimated time for installation will greatly vary on your network speed and on whether software Like R or Python is installed on your system. 

1. Download all files from Github Repository 

2. Download FASTA for hg19 [from here](https://www.dropbox.com/s/a6vthezotm6iaih/Genome_v19.fasta.gz?dl=0). <b> If running on FASTA files from other sources, make sure it is unmasked and follows the ENSEMBL FASTA format (UCSC will throw errors as the chromosomes are named differently). </b>

3. Download comprehensive gene annotation V19 (GRCh37 for the test input to work) [from here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz). (Filtered Long non-coding RNA GTF provided in the input folder, however the comprehensive gene annotation must be downloaded.) 

4. *For running functional impact script,* download CADD scores for hg19 [from here](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz) and its index [from here](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi) 

5. Some Input files are stored as .zip files to reduce size. Make sure to uncompress them before running ExInAtor2 on it (for example, Biliary_AdenoCA.bed.zip, black_list_regions_final.bed.zip, gencode.v19.long_noncoding_RNAs.gtf.zip) 
 
# Inputs

<b> For recurrence: </b>

*Mandatory arguments:*

-i <--input_file>: Input file with mutations (SNVs) in BED format

-o <--output_folder>: Path to folder where all output files will be stored

-f <--fasta_file>: FASTA file of the human genome

-g <--gtf_file>: Path to Gencode GTF containing information only on long non-coding RNAs. 

-t <--true_set>: Path to true set of cancer-driver lncRNAs (in this case from the [CLC](https://academic.oup.com/narcancer/article/3/2/zcab013/6225859))

-k <--kmers_file>: Txt file containing all the possible trinucleotides.

-z <--chr_sizes>: Two column file with name and lengh of the chromosomes. These lengths have to match the ones in the FASTA file indicated previously

-y <--chr_sizes_long>: Three column file with name, start and end of the chromosomes. These lengths have to match the ones in the FASTA file indicated previously

-w <--whole_genome>: Path to Gencode GTF containing information about all genes in the genome. 

*Optional arguments:*

-b <--background_size>: The extension length of the background region that includes all introns

-s <--simulations>: Number of iterations to perform simulations, 10000 recommended

-c <--cores>: Number of cores in the computer to use for the simulations

-e <--exc_regions>: BED file with regions from the genome to ignore (such as those with low mappability, high repetitive sequences, etc)


<b> For functional impact: </b>

*Mandatory arguments:*

-m <--mutations>: Input file with mutations (SNVs) in BED format

-o <--output_folder>: Path to folder where all output files will be stored

-f <--fasta_file>: FASTA file of the human genome

-g <--gtf_file>: GTF file of the Gencode annotation to use when creating exons

-z <--chr_sizes>: Two column file with name and lengh of the chromosomes. These lengths have to match the ones in the FASTA file indicated previously

-s <--scores_file>: TSV file with scores having these five columns: #Chr    Pos     Ref     Alt     RawScore        PHRED

*Optional arguments:*

-i <--iterations>: Number of iterations to perform simulations, 10000 recommended

-c <--cores>: Number of cores in the computer to use for the simulations

-t <--true_set>: This is a TXT file with ENSGs of your true positives genes, that will be used in the R script, from the [CLC](https://academic.oup.com/narcancer/article/3/2/zcab013/6225859). 

-e <--exc_regions>: BED file with regions from the genome to ignore (such as those with low mappability, high repetitive sequences, etc)

# Outputs 

<b> For recurrence: </b>

The output files for recurrence are stored in Output_REC. Several intermediary files will be generated. The output files included in the repository are:

1. qqplot.png: QQplot of the expected and observed pvalues

2. venn.png: Venn Diagram of all genes and genes in the true positive set. 

3. precision.png: Plot with genes ranked by p value on the x axis and percentage of genes in the true set on the y axis

4. "results_Biliary-AdenoCA_fixed2_10000_100_12.txt" - The final output containing the p-vals and q-vals of input genes. 

5. "table_probabilities_R.txt" - R table with probabilies.

<b> For functional impact: </b>

Same as above, only the output will be stored in the folder Output_FI. 

# Example

1. Unzip .bed, .fasta and .gtf files before running command: 

unzip -c ~/Exinator2/Inputs/Biliary-AdenoCA.bed.zip > ~/Exinator2/Inputs/Biliary-AdenoCA.bed

unzip -c ~/Exinator2/Inputs/gencode.v19.long_noncoding_RNAs.gtf.zip > ~/Exinator2/Inputs/gencode.v19.long_noncoding_RNAs.gtf

unzip -c ~/gencode.v19.annotation.gtf.zip > ~/gencode.v19.annotation.gtf

gunzip -c ~/Genome_v19.fasta.gz > ~/Genome_v19.fasta

2. Make sure to download additional files mentioned in the Installation section (Genome fasta file "Genome_v19.fasta", CADD score file for functional impact: whole_genome_SNVs.tsv.gz)

3. Full paths are required when running ExInAtor

<b> For recurrence: </b>

python ~/Exinator2/recurrence_script_main.py -b 10000 -f ~/Genome_v19.fasta -o ~/Exinator2/Output_REC2 -g ~/Exinator2/Inputs/gencode.v19.long_noncoding_RNAs.gtf -i ~/Exinator2/Inputs/Biliary-AdenoCA.bed -k ~/Exinator2/Inputs/3mers.txt -w ~/gencode.v19.annotation.gtf -z ~/Exinator2/Inputs/chromosomes_tab.bed -t ~/Exinator2/Inputs/CLC2_final_extended.txt -y ~/Exinator2/Inputs/chromosomes_long_tab.bed -e ~/Exinator2/Inputs/black_list_regions_final.bed -s 100

<b> For functional impact: </b>

python ~/Exinator2/functionalimpact_script_main.py -c 1 -i 100 -m ~/Exinator2/Inputs/Biliary-AdenoCA.bed -o ~/Exinator2/Output_FI2 -g ~/Exinator2/Inputs/gencode.v19.long_noncoding_RNAs.gtf -f ~/Genome_v19.fasta -z ~/Exinator2/Inputs/chromosomes_tab.bed -s ~/whole_genome_SNVs.tsv.gz -t ~/Exinator2/Inputs/CLC2_final_extended.txt -e ~/Exinator2/Inputs/black_list_regions_final.bed

Estimated runtime for 100 simulations (as stated in the example code on the example input) on a system with 4 cores and 32GB RAM: 

1. For Recurrence: ~10 minutes 
2. For Functional Impact: ~8 minutes
