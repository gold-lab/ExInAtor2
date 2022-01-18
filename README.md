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

# Requirements 

Python 

R

Bedtools 

# Installation 

1. Download all files from Github Repository 

2. Download FASTA for hg19 [from here](https://www.dropbox.com/s/a6vthezotm6iaih/Genome_v19.fasta.gz?dl=0)

3. <b> For running functional impact script,</b> download CADD scores for hg19 [from here](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz) and its index [from here](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi) 

<b> For recurrence: </b>

python /home/Exinator2/recurrence_script_main.py -b 10000 -f /home/Exinator2/Inputs/Genome_v19.fasta -o /home/Exinator2/Output_REC -g /home/Exinator2/Inputs/gencode.v19.long_noncoding_RNAs.gtf -i /home/Exinator2/Inputs/Biliary-AdenoCA.bed -k /home/Exinator2/Inputs/3mers.txt -w /home/Exinator2/Inputs/gencode.v19.annotation.gtf -z /home/Exinator2/Inputs/chromosomes_tab.bed -t /home/Exinator2/Inputs/CLC2_final_extended.txt -y /home/Exinator2/Inputs/chromosomes_long_tab.bed -e /home/Exinator2/Inputs/black_list_regions_final.bed -s 100

<b> For functional impact: </b>

python /home/Exinator2/functionalimpact_script_main.py -c 1 -i 100 -m /home/Exinator2/Inputs/Biliary-AdenoCA.bed -o /home/Exinator2/Output_FI -g /home/Exinator2/Inputs/gencode.v19.long_noncoding_RNAs.gtf -f /home/Exinator2/Inputs/Genome_v19.fasta -z /home/Exinator2/Inputs/chromosomes_tab.bed -s /home/Exinator2/Inputs/whole_genome_SNVs.tsv.gz -t /home/Exinator2/Inputs/CLC2_final_extended.txt -e /home/Exinator2/Inputs/black_list_regions_final.bed

# Inputs

<b> For recurrence: </b>

Mandatory arguments: 

-i <--input_file>: Input file with mutations (SNVs) in BED format

-o <--output_folder>: Path to folder where all output files will be stored

-f <--fasta_file>: FASTA file of the human genome

-g <--gtf_file>: Path to GTF containing lncRNAs

-t <--true_set>: Path to true set of cancer-driver lncRNAs (in this case from the CLC)

-k <--kmers_file>: Txt file containing all the possible trinucleotides.

-z <--chr_sizes>: Two column file with name and lengh of the chromosomes. These lengths have to match the ones in the FASTA file indicated previously

-y <--chr_sizes_long>: Two column file with name, start and end of the chromosomes. These lengths have to match the ones in the FASTA file indicated previously

-w <--whole_genome>: Path to whole genome GTF

Optional arguments: 

-b <--background_size>: the extension length of the background region that includes all introns

-s <--simulations>: Number of iterations to perform simulations, 10000 recommended

-c <--cores>: Number of cores in the computer to use for the simulations

-r <--reg_regions>: 

-e <--exc_regions>: BED file with regions from the genome to ignore (such as those with low mappability, high repetitive sequences, etc)


<b> For functional impact: </b>

Mandatory arguments: 

-m <--mutations>: Input file with mutations (SNVs) in BED format

-o <--output_folder>: Path to folder where all output files will be stored

-f <--fasta_file>: FASTA file of the human genome

-g <--gtf_file>: GTF file of the Gencode annotation to use when creating exons

-z <--chr_sizes>: Two column file with name and lengh of the chromosomes. These lengths have to match the ones in the FASTA file indicated previously

-s <--scores_file>: TSV file with scores having these five columns: #Chr    Pos     Ref     Alt     RawScore        PHRED

Optional arguments:

-i <--iterations>: Number of iterations to perform simulations, 10000 recommended

-c <--cores>: Number of cores in the computer to use for the simulations

-t <--true_set>: This is a TXT file with ENSGs of your true positives genes, that will be used in the R script.

-e <--exc_regions>: BED file with regions from the genome to ignore (such as those with low mappability, high repetitive sequences, etc)

# Outputs 

<b> For recurrence: </b>

The output files for recurrence are stored in Output_REC. Several intermediary files will be generated. The output files included in the repository are:

1. QQplot

2. Venn Daigram 

3. Precision plot 

4. "results_Biliary-AdenoCA_fixed2_10000_100_12.txt" - The final output containing the p-vals and q-vals of input genes. 

5. "table_probabilities_R.txt" - R table with probabilies.

<b> For functional impact: </b>

Same as above, only the output will be stored in the folder Output_FI. 
