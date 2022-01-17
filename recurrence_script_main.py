# coding: utf-8
import sys, getopt
import os
import time
import multiprocessing as mp
import random
import string
import math
import numpy
from math import *
import unittest
from random import randint
from sys import stdout
import os.path
import subprocess
import regex as re
start_time = time.time()
seed_number=12
random.seed(seed_number)
numpy.random.seed(seed_number)

if len(sys.argv)==1 or "-h" in sys.argv: 
	print("Mandatory arguments: -i <--input_file> -o <--output_folder> -f <--fasta_file> -g <--gtf_file> -t <--true_set> -k <--kmers_file> -z <--chr_sizes> -y <--chr_sizes_long> -w <--whole_genome>")
	print("Optional arguments: -b <--background_size> -s <--simulations> -c <--cores> -r <--reg_regions> -e <--exc_regions>")
	sys.exit()
n=1
while n<len(sys.argv):
	arg = sys.argv[n]
	if arg == "-i" or arg == "--input_file":
		input_file = sys.argv[n+1]
	elif arg == "-o" or arg == "--output_folder":
		output_folder = sys.argv[n+1]
	elif arg == "-f" or arg == "--fasta_file":
		fasta_file = sys.argv[n+1]
	elif arg == "-g" or arg == "--gtf_file":
		gtf_file = sys.argv[n+1]
	elif arg == "-k" or arg == "--kmers_file":
		kmers_file = sys.argv[n+1]
	elif arg == "-s" or arg == "--simulations":
		simulations = int(sys.argv[n+1])
	elif arg == "-z" or arg == "--chr_sizes":
		chr_sizes = sys.argv[n+1]
	elif arg == "-y" or arg == "--chr_sizes_long":
		chr_sizes_long = sys.argv[n+1]	
	elif arg == "-c" or arg == "--cores":
		cores = int(sys.argv[n+1])
	elif arg == "-b" or arg == "--background_size":
		background_size = sys.argv[n+1]
	elif arg == "-w" or arg == "--whole_genome":
		whole_genome = sys.argv[n+1]
	elif arg == "-t" or arg == "--true_set":
		true_set = sys.argv[n+1]
	elif arg == "-r" or arg == "--reg_regions":
		reg_regions = sys.argv[n+1]
	elif arg == "-e" or arg == "--exc_regions":
		exc_regions = sys.argv[n+1]				
	n+=2
if 'input_file' not in globals():
	print("Please indicate \"-i\" or \"--input_file\"")
	print("type \"test.py -h\" for help")
	sys.exit()
if 'output_folder' not in globals():
	print("Please indicate \"-o\" or \"--output_folder\"")
	print("type \"test.py -h\" for help")
	sys.exit()
if 'fasta_file' not in globals():
	print("Please indicate \"-f\" or \"--fasta_file\"")
	print("type \"test.py -h\" for help")
	sys.exit()
if 'gtf_file' not in globals():
	print("Please indicate \"-g\" or \"--gtf_file\"")
	print("type \"test.py -h\" for help")
	sys.exit()
if 'kmers_file' not in globals():
	print("Please indicate \"-k\" or \"--kmers_file\"")
	print("type \"test.py -h\" for help")
	sys.exit()
if 'whole_genome' not in globals():
	print("Please indicate \"-w\" or \"--whole_genome\"")
	print("type \"v1.py -h\" for help")
	sys.exit()
if 'chr_sizes' not in globals():
	print("Please indicate \"-z\" or \"--chr_sizes\"")
	print("type \"v1.py -h\" for help")
	sys.exit()
if 'chr_sizes_long' not in globals():
	print("Please indicate \"-y\" or \"--chr_sizes_long\"")
	print("type \"v1.py -h\" for help")
	sys.exit()	
if 'true_set' not in globals():
	print("Please indicate \"-t\" or \"--true_set\"")
	print("type \"v1.py -h\" for help")
	sys.exit()	
if 'simulations' not in globals():
	simulations=2000
if 'cores' not in globals():
	cores=12
if 'background_size' not in globals():
	background_size="10000"
if 'reg_regions' not in globals():
	reg_regions="NO"
if 'exc_regions' not in globals():
	exc_regions="NO"	
print("Arguments correctly inserted: %.0f seconds " % (time.time() - start_time))
del arg


# Initial setup


if not os.path.exists(output_folder):
    os.makedirs(output_folder)
os.chdir(output_folder)
input_folder="/".join(input_file.split("/")[:-1])+"/"
input_file=input_file.split("/")[-1]
print("Starting the analysis of "+input_file.split(".")[0]+": %.0f seconds " % (time.time() - start_time))
if not os.path.exists(output_folder+"/"+input_file.split(".")[0]):
	os.makedirs(output_folder+"/"+input_file.split(".")[0])
filename=input_file.split(".")[0]
os.chdir(output_folder+"/"+input_file.split(".")[0])
input_file=input_folder+input_file


# Create a BED file of the merged exons, and introns, for each lncRNA gene.


if not os.path.isfile("table_kmer_counts.txt"):
	subprocess.call("cat "+whole_genome+' | tr --delete \; | tr --delete \\" | grep -v \\# | grep -v chrM | awk \'$3==\"exon\"\' | awk \'{split($10,array,\".\")}{print $10\"\t\"$4-1\"\t\"$5\"\t\"$1\"\t.\t\"$7}\' | sort -T ~/andres/ -k1,1 -k2,2n > exons_n2.bed', shell=True)
	subprocess.call("mergeBed -i exons_n2.bed -c 4,6 -o distinct | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t.\t\"$5}\'> merged_exons_n2.bed", shell=True)
	subprocess.call("cat "+gtf_file+' | tr --delete \; | tr --delete \\" | grep -v \\# | grep -v chrM | awk \'$3==\"exon\"\' | awk \'{split($10,array,\".\")}{print $10\"\t\"$4-1\"\t\"$5\"\t\"$1\"\t.\t\"$7}\' | sort -T ~/andres/ -k1,1 -k2,2n > exons_n.bed', shell=True)
	subprocess.call("mergeBed -i exons_n.bed -c 4,6 -o distinct | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t.\t\"$5}\'> merged_exons_n.bed", shell=True)
	subprocess.call("cat "+gtf_file+' | tr --delete \; | tr --delete \\" | grep -v \\# | grep -v chrM | awk \'$3==\"gene\"\' | awk \'{split($10,array,\".\")}{print $1\"\t\"$4-1\"\t\"$5\"\t\"$10\"\t.\t\"$7}\' | sort -T ~/andres/ -k1,1 -k2,2n > genes_n.bed', shell=True)
	subprocess.call("cat "+whole_genome+' | tr --delete \; | tr --delete \\" | grep -v \\# | grep -v chrM | awk \'$3==\"gene\"\' | awk \'{split($10,array,\".\")}{print $1\"\t\"$4-1\"\t\"$5\"\t\"$10\"\t.\t\"$7}\' | sort -T ~/andres/ -k1,1 -k2,2n > genes_n2.bed', shell=True)
	
	if exc_regions!="NO":
		#subprocess.call("subtractBed -a genes_n2.bed -b "+exc_regions+" > genes_n2_clean.bed", shell=True)
		subprocess.call("subtractBed -a merged_exons_n2.bed -b "+exc_regions+" > merged_exons_n2_clean.bed", shell=True)
		#subprocess.call("subtractBed -a genes_n.bed -b "+exc_regions+" > genes_n_clean.bed", shell=True)
		subprocess.call("subtractBed -a merged_exons_n.bed -b "+exc_regions+" > merged_exons_n_clean.bed", shell=True)	
		#subprocess.call("mv genes_n2_clean.bed genes_n2.bed", shell=True)
		subprocess.call("mv merged_exons_n2_clean.bed merged_exons_n2.bed", shell=True)
		#subprocess.call("mv genes_n_clean.bed genes_n.bed", shell=True)
		subprocess.call("mv merged_exons_n_clean.bed merged_exons_n.bed", shell=True)
		subprocess.call("subtractBed -a "+chr_sizes_long+" -b "+exc_regions+" > chr_sizes_long_clean.bed", shell=True)
		chr_sizes_long="chr_sizes_long_clean.bed"
	
	if reg_regions=="YES":
	
		# Create promoters
		subprocess.call("bedtools flank -i genes_n2.bed -g "+chr_sizes+" -l 200 -r 0 -s > promoters_2.bed", shell=True)
		subprocess.call("bedtools flank -i genes_n.bed -g "+chr_sizes+" -l 200 -r 0 -s > promoters.bed", shell=True)
		subprocess.call("cat promoters_2.bed >> merged_exons_n2.bed", shell=True)
		subprocess.call("cat promoters.bed >> merged_exons_n.bed", shell=True)
		subprocess.call("awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t.\t\"$6}\' genes_n2.bed | sort -T ~/andres/ -k1,1 -k2,2n > genes_n2_temp.bed", shell=True)
		subprocess.call("awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t.\t\"$6}\' genes_n.bed | sort -T ~/andres/ -k1,1 -k2,2n > genes_n_temp.bed", shell=True)
		subprocess.call("awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t.\t\"$6}\' promoters_2.bed | sort -T ~/andres/ -k1,1 -k2,2n > promoters_2_temp.bed", shell=True)
		subprocess.call("awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t.\t\"$6}\' promoters.bed | sort -T ~/andres/ -k1,1 -k2,2n > promoters_temp.bed", shell=True)
		subprocess.call("cat promoters_2_temp.bed >> genes_n2_temp.bed", shell=True)
		subprocess.call("cat promoters_temp.bed >> genes_n_temp.bed", shell=True)
		subprocess.call("sort -T ~/andres/ -k1,1 -k2,2n -o genes_n2_temp.bed genes_n2_temp.bed", shell=True)
		subprocess.call("sort -T ~/andres/ -k1,1 -k2,2n -o genes_n_temp.bed genes_n_temp.bed", shell=True)
		subprocess.call("mergeBed -i genes_n2_temp.bed -c 4,6 -o distinct | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t.\t\"$5}\' | sort -T ~/andres/ -k1,1 -k2,2n> genes_n2.bed", shell=True)
		subprocess.call("mergeBed -i genes_n_temp.bed -c 4,6 -o distinct | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t.\t\"$5}\' | sort -T ~/andres/ -k1,1 -k2,2n> genes_n.bed", shell=True)
					
		# Create splice sites
		
		subprocess.call("cat "+gtf_file+' | tr --delete \; | tr --delete \\" | grep -v \\# | grep -v chrM | awk \'$3==\"gene\"\' | awk \'{print $10\"\t\"$4-1\"\t\"$5\"\t\"$1\"\t.\t\"$7}\' | sort -T ~/andres/ -k1,1 -k2,2n > genes_n_for_ss.bed', shell=True)
		subprocess.call("cat "+whole_genome+' | tr --delete \; | tr --delete \\" | grep -v \\# | grep -v chrM | awk \'$3==\"gene\"\' | awk \'{print $10\"\t\"$4-1\"\t\"$5\"\t\"$1\"\t.\t\"$7}\' | sort -T ~/andres/ -k1,1 -k2,2n > genes_n2_for_ss.bed', shell=True)
		subprocess.call("mergeBed -i exons_n2.bed -c 4,6 -o distinct | sort -T ~/andres/ -k1,1 -k2,2n > merged_exons_n2_for_ss.bed", shell=True)
		subprocess.call("mergeBed -i exons_n.bed -c 4,6 -o distinct | sort -T ~/andres/ -k1,1 -k2,2n > merged_exons_n_for_ss.bed", shell=True)
		subprocess.call("subtractBed -a genes_n_for_ss.bed -b merged_exons_n_for_ss.bed > introns_n_for_ss.bed", shell=True)
		subprocess.call("subtractBed -a genes_n2_for_ss.bed -b merged_exons_n2_for_ss.bed > introns_n2_for_ss.bed", shell=True)
		subprocess.call("awk \'{print $4\"\t\"$2+20\"\t\"$3-20\"\t\"$1\"\t\"$5\"\t\"$6}\' introns_n_for_ss.bed | awk \'$3>$2\' > introns_n_for_ss_reduced.bed", shell=True)
		subprocess.call("awk \'{print $4\"\t\"$2+20\"\t\"$3-20\"\t\"$1\"\t\"$5\"\t\"$6}\' introns_n2_for_ss.bed | awk \'$3>$2\' > introns_n2_for_ss_reduced.bed", shell=True)
		subprocess.call("sort -T ~/andres/ -k1,1 -k2,2n -o introns_n_for_ss_reduced.bed introns_n_for_ss_reduced.bed", shell=True)
		subprocess.call("sort -T ~/andres/ -k1,1 -k2,2n -o introns_n2_for_ss_reduced.bed introns_n2_for_ss_reduced.bed", shell=True)
		subprocess.call("bedtools slop -i introns_n_for_ss_reduced.bed -g "+chr_sizes+" -l 14 -r 0 -s | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t\"$5\"\t\"$6}\' > introns_n_for_ss_reduced_2.bed", shell=True)
		subprocess.call("bedtools slop -i introns_n2_for_ss_reduced.bed -g "+chr_sizes+" -l 14 -r 0 -s | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t\"$5\"\t\"$6}\' > introns_n2_for_ss_reduced_2.bed", shell=True)
		subprocess.call("subtractBed -a introns_n_for_ss.bed -b introns_n_for_ss_reduced_2.bed | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t\"$5\"\t\"$6}\' > splice_sites.bed", shell=True)
		subprocess.call("subtractBed -a introns_n2_for_ss.bed -b introns_n2_for_ss_reduced_2.bed | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1\"\t\"$5\"\t\"$6}\' > splice_sites_2.bed", shell=True)
		subprocess.call("cat splice_sites.bed >> merged_exons_n.bed", shell=True)
		subprocess.call("cat splice_sites_2.bed >> merged_exons_n2.bed", shell=True)
		subprocess.call("sort -T ~/andres/ -k1,1 -k2,2n -o merged_exons_n.bed merged_exons_n.bed", shell=True)
		subprocess.call("sort -T ~/andres/ -k1,1 -k2,2n -o merged_exons_n2.bed merged_exons_n2.bed", shell=True)


	#background
	
	subprocess.call("bedtools slop -i genes_n.bed -g "+chr_sizes+" -b "+background_size+" > genes_extended_n.bed", shell=True)
	
	# Get non-extended introns of genes
	subprocess.call("subtractBed -a genes_n.bed -b merged_exons_n.bed > introns_n.bed", shell=True)
	subprocess.call("subtractBed -a introns_n.bed -b merged_exons_n2.bed > introns_n2.bed", shell=True)
	if exc_regions!="NO":
		subprocess.call("subtractBed -a introns_n2.bed -b "+exc_regions+" > introns_n2_clean.bed", shell=True)	
		subprocess.call("mv introns_n2_clean.bed introns_n2.bed", shell=True)

	intronic_regions_dir={}
	intronic_regions_dir2={}
	file = open ("genes_n.bed","r")
	for line in file:
		chr,start,end,gene,dot,strand=line.strip().split("\t")
		if gene not in intronic_regions_dir:
			intronic_regions_dir[gene]=0
		if gene not in intronic_regions_dir2:
			intronic_regions_dir2[gene]=str(strand)
	file.close()
	
	intronic_regions_file=open("introns_n2.bed","r")

	for line in intronic_regions_file:
		chr,start,end,gene,dot,strand=line.strip().split("\t")
		start=int(start)
		end=int(end)
		length=end-start
		gene=str(gene)
		if gene not in intronic_regions_dir:
			intronic_regions_dir[gene]=0
		if gene not in intronic_regions_dir2:
			intronic_regions_dir2[gene]=str(strand)
		intronic_regions_dir[gene]+=length
	intronic_regions_file.close()	
	
	subprocess.call("subtractBed -a "+chr_sizes_long+" -b merged_exons_n.bed > chr_sizes_n.bed", shell=True)
	subprocess.call("subtractBed -a chr_sizes_n.bed -b merged_exons_n2.bed > chr_sizes_n2.bed", shell=True)
	if exc_regions!="NO":
		subprocess.call("subtractBed -a chr_sizes_n2.bed -b "+exc_regions+" > chr_sizes_n2_clean.bed", shell=True)	
		subprocess.call("mv chr_sizes_n2_clean.bed chr_sizes_n2.bed", shell=True)
	no_exonic_regions_file=open("chr_sizes_n2.bed","r")
	no_exonic_regions_dir={}
	for line in no_exonic_regions_file:
		chr,start,end=line.strip().split()
		chr=str(chr)
		if chr not in no_exonic_regions_dir: no_exonic_regions_dir[chr]=[]
		no_exonic_regions_dir[chr].append([int(start),int(end)])
	no_exonic_regions_file.close()		
	genes_file=open("genes_n.bed","r")
	background_file=open("introns_n2.bed","a")
	for line in genes_file:
		chr,start,end,gene,dot,strand=line.strip().split()
		chr=str(chr)
		start=int(start)
		end=int(end)
		list=no_exonic_regions_dir[chr]
		list_left=[item for item in list if item[1] <= start]
		list_left.reverse()
		background_size_temp=0
		background_size_max=int(background_size)
		i=0
		j=len(list_left)
		while background_size_temp<=background_size_max and i<j:
			background_size_temp+=list_left[i][1]-list_left[i][0]
			if background_size_temp>background_size_max:
				background_size_extra=background_size_temp-background_size_max
				if (list_left[i][0]+background_size_extra)==list_left[i][1]:
					line_to_print=str(chr)+"\t"+str(list_left[i][0]+background_size_extra-1)+"\t"+str(list_left[i][1])+"\t"+str(gene)+"\t.\t"+intronic_regions_dir2[gene]+"\n"
				else:
					line_to_print=str(chr)+"\t"+str(list_left[i][0]+background_size_extra)+"\t"+str(list_left[i][1])+"\t"+str(gene)+"\t.\t"+intronic_regions_dir2[gene]+"\n"
			else:
				line_to_print=str(chr)+"\t"+str(list_left[i][0])+"\t"+str(list_left[i][1])+"\t"+str(gene)+"\t.\t"+intronic_regions_dir2[gene]+"\n"
			background_file.write(line_to_print)
			i+=1
				
		list_right=[item for item in list if item[0] >= end]
		background_size_temp=0
		background_size_max=int(background_size)
		i=0
		j=len(list_right)
		while background_size_temp<=background_size_max and i<j:
			background_size_temp+=list_right[i][1]-list_right[i][0]
			if background_size_temp>background_size_max:
				background_size_extra=background_size_temp-background_size_max
				if (list_right[i][1]-background_size_extra)==list_right[i][0]:
					line_to_print=str(chr)+"\t"+str(list_right[i][0])+"\t"+str(list_right[i][1]-background_size_extra+1)+"\t"+str(gene)+"\t.\t"+intronic_regions_dir2[gene]+"\n"
				else:
					line_to_print=str(chr)+"\t"+str(list_right[i][0])+"\t"+str(list_right[i][1]-background_size_extra)+"\t"+str(gene)+"\t.\t"+intronic_regions_dir2[gene]+"\n"
			else:
				line_to_print=str(chr)+"\t"+str(list_right[i][0])+"\t"+str(list_right[i][1])+"\t"+str(gene)+"\t.\t"+intronic_regions_dir2[gene]+"\n"
			background_file.write(line_to_print)
			i+=1
	background_file.close()
						
	subprocess.call('sort -T ~/andres/ -k 1,1 -k2,2n introns_n2.bed > introns.bed', shell=True)
	subprocess.call('sort -T ~/andres/ -k 1,1 -k2,2n exons_n.bed > exons.bed', shell=True)
	subprocess.call('sort -T ~/andres/ -k 1,1 -k2,2n merged_exons_n.bed > merged_exons.bed', shell=True)
	subprocess.call('sort -T ~/andres/ -k 1,1 -k2,2n genes_extended_n.bed > genes.bed', shell=True)
	print("\tBED files created for merged exons and introns: %.0f seconds " % (time.time() - start_time))

	# Count mutation in merged exons and introns

	subprocess.call('intersectBed -a merged_exons.bed -b '+input_file+' -sorted -wb | awk \'{print $1"\t"$2-1"\t"$3+1"\t"$4"\t"$5"\t"$6}\' > merged_exons_mutations.bed', shell=True)
	subprocess.call('intersectBed -a introns.bed -b '+input_file+' -sorted -wb | awk \'{print $1"\t"$2-1"\t"$3+1"\t"$4"\t"$5"\t"$6}\' > introns_mutations.bed', shell=True)
	print("\tMutations counted in merged exons mutations and introns: %.0f seconds " % (time.time() - start_time))


	# Get fasta for merged exons mutations and introns mutations


	subprocess.call("sed 's/chr//g' merged_exons_mutations.bed > merged_exons_mutations_for_fasta.bed", shell=True)
	subprocess.call("sed 's/chr//g' introns_mutations.bed > introns_mutations_for_fasta.bed", shell=True)
	subprocess.call(' bedtools getfasta -fi '+fasta_file+' -bed merged_exons_mutations_for_fasta.bed -fo merged_exons_mutations.fa -name -s', shell=True)
	subprocess.call(' bedtools getfasta -fi '+fasta_file+' -bed introns_mutations_for_fasta.bed -fo introns_mutations.fa -name -s', shell=True)
	print("\tFasta created for merged exons mutations and introns mutations: %.0f seconds " % (time.time() - start_time))


	# Get genes and kmers


	genes=[]
	file = open ("genes.bed","r")
	for line in file:
		line=line.rstrip()
		line=line.split("\t")
		genes.append(line[3])
	file.close()
	file = open (kmers_file,"r")
	kmers=[]
	for line in file:
		line=line.rstrip()
		kmers.append(line)
	file.close()
	print("\tGenes and kmers gathered: %.0f seconds " % (time.time() - start_time))


	# Define functions to read, count and print kmers


	def count_kmers(list):
		file = open (list[0],"r")
		kmers=list[1]
		order=list[2]
		exons={}
		header = None
		for line in file:
			line=line.rstrip()
			if line.startswith('>'):
				if header is not None:
					if header not in exons:
						exons[header]={}
					if header in exons:
						for motif in kmers:
							if motif not in exons[header]:
								exons[header][motif]=0
							if motif in exons[header]:
								exons[header][motif]+=len(re.findall(motif, sequence, overlapped=True))
				header = line.split("::")[0][1:]
				sequence=""
			else:
				sequence=sequence+line
		if header is not None:
			if header not in exons:
				exons[header]={}
			if header in exons:
				for motif in kmers:
					if motif not in exons[header]:
						exons[header][motif]=0
					if motif in exons[header]:
						exons[header][motif]+=len(re.findall(motif, sequence, overlapped=True))
		file.close()
		return((order,exons))
	def print_kmers(list):
		file=list[1]
		list=list[0]
		file=open(file,"w")
		for header in list:
			for motif in list[header]:
				line_to_print=header+"\t"+motif+"\t"+str(list[header][motif])+"\n"
				file.write(line_to_print)
		file.close()
		return()
	def read_kmers(file_name):
		dir={}
		file=open(file_name,"r")
		for line in file:
			line=line.rstrip()
			line=line.split("\t")
			header=line[0]
			motif=line[1]
			count=line[2]
			if header not in dir:
				dir[header]={}
			if motif not in dir[header]:
				dir[header][motif]=count
		file.close()
		return(dir)


	# Count and print kmers


	if not os.path.isfile(output_folder+"/exonic_kmers.txt"):
		if not os.path.isfile(output_folder+"/intronic_kmers.txt"):
			subprocess.call("bedtools slop -i merged_exons.bed -g "+chr_sizes+" -b 1 | sed 's/chr//g' > merged_exons_for_fasta2.bed", shell=True)
			subprocess.call("bedtools slop -i introns.bed -g "+chr_sizes+" -b 1 | sed 's/chr//g' > introns_for_fasta2.bed", shell=True)
			subprocess.call('bedtools getfasta -fi '+fasta_file+' -bed merged_exons_for_fasta2.bed -fo merged_exons.fa -name -s', shell=True)
			subprocess.call('bedtools getfasta -fi '+fasta_file+' -bed introns_for_fasta2.bed -fo introns.fa -name -s', shell=True)
			pool = mp.Pool(processes=cores)
			results = pool.map(count_kmers, [["introns.fa",kmers,2],["merged_exons.fa",kmers,1],["merged_exons_mutations.fa",kmers,3],["introns_mutations.fa",kmers,4]])
		else:
			subprocess.call("bedtools slop -i merged_exons.bed -g "+chr_sizes+" -b 1 | sed 's/chr//g' > merged_exons_for_fasta2.bed", shell=True)
			subprocess.call('bedtools getfasta -fi '+fasta_file+' -bed merged_exons_for_fasta2.bed -fo merged_exons.fa -name -s', shell=True)
			pool = mp.Pool(processes=cores)
			results = pool.map(count_kmers, [["merged_exons.fa",kmers,1],["merged_exons_mutations.fa",kmers,3],["introns_mutations.fa",kmers,4]])
			introns=read_kmers(output_folder+"/intronic_kmers.txt")
	else:
		if not os.path.isfile(output_folder+"/intronic_kmers.txt"):
			subprocess.call("bedtools slop -i introns.bed -g "+chr_sizes+" -b 1 | sed 's/chr//g' > introns_for_fasta2.bed", shell=True)
			subprocess.call('bedtools getfasta -fi '+fasta_file+' -bed introns_for_fasta2.bed -fo introns.fa -name -s', shell=True)
			pool = mp.Pool(processes=cores)
			results = pool.map(count_kmers, [["introns.fa",kmers,2],["merged_exons_mutations.fa",kmers,3],["introns_mutations.fa",kmers,4]])
			exons=read_kmers(output_folder+"/exonic_kmers.txt")
		else:
			pool = mp.Pool(processes=cores)
			results = pool.map(count_kmers, [["merged_exons_mutations.fa",kmers,3],["introns_mutations.fa",kmers,4]])
			exons=read_kmers(output_folder+"/exonic_kmers.txt")
			introns=read_kmers(output_folder+"/intronic_kmers.txt")
	for element in results:
		if element[0]==1: exons=element[1]
		if element[0]==2: introns=element[1]
		if element[0]==3: exons_mutations=element[1]
		if element[0]==4: introns_mutations=element[1]
			
	if not os.path.isfile(output_folder+"/exonic_kmers.txt"):
		if not os.path.isfile(output_folder+"/intronic_kmers.txt"):
			pool = mp.Pool(processes=cores)
			variable = pool.map(print_kmers, [[exons,output_folder+"/exonic_kmers.txt"],[introns,output_folder+"/intronic_kmers.txt"]])
		else:
			print_kmers(exons,output_folder+"/exonic_kmers.txt")
	else:
		if not os.path.isfile(output_folder+"/intronic_kmers.txt"):
			print_kmers(introns,output_folder+"/intronic_kmers.txt")
	print("\tKmers counted: %.0f seconds " % (time.time() - start_time))


	# Combine all the information before simulation


	file=open("table_kmer_counts.txt","w")
	for gene in genes:
		for motif in kmers:
			if gene not in exons:
				a="0"
				b="0"
			if gene in exons:
				if motif not in exons[gene]:
					a="0"
					b="0"
				if motif in exons[gene]:
					b=str(exons[gene][motif])
					if gene not in exons_mutations:
						a="0"
					if gene in exons_mutations:
						if motif in exons_mutations[gene]:
							a=str(exons_mutations[gene][motif])
						if motif not in exons_mutations[gene]:
							a="0"
			if gene not in introns:
				c="0"
				d="0"
			if gene in introns:
				if motif not in introns[gene]:
					c="0"
					d="0"
				if motif in introns[gene]:
					d=str(introns[gene][motif])
					if gene not in introns_mutations:
						c="0"
					if gene in introns_mutations:
						if motif in introns_mutations[gene]:
							c=str(introns_mutations[gene][motif])
						if motif not in introns_mutations[gene]:
							c="0"
			line_to_print=gene+"\t"+motif+"\t"+a+"\t"+b+"\t"+c+"\t"+d+"\n"
			file.write(line_to_print)
	file.close()
	del results, exons, introns, exons_mutations, introns_mutations
	subprocess.call('cp table_kmer_counts.txt table_kmer_counts_backup.txt', shell=True)
	subprocess.call('grep -v -P "\.\t" table_kmer_counts_backup.txt > table_kmer_counts.txt', shell=True)
	print("\tAll information combined: %.0f seconds " % (time.time() - start_time))


# Get calculated kmers


#if not os.path.isfile('results_'+filename+"_"+str(background_mode)+"_"+str(background_size)+"_"+str(simulations)+"_"+str(seed_number)+'.txt'):

if not os.path.isfile('table_probabilities_R.txt'):
	file=open("table_kmer_counts.txt","r")
	table_kmer_counts={}
	input_file=[]
	real_exonic_mutations={}
	simulated_exonic_mutations={}
	i=1
	for line in file:
		line=line.rstrip()
		line=line.split("\t")
		gene=line[0]
		motif=line[1]
		ex_mut=int(line[2])
		ex_len=int(line[3])
		in_mut=int(line[4])
		in_len=int(line[5])
		tot_mut=ex_mut+in_mut
		tot_len=ex_len+in_len
		if tot_mut>0:
			if gene in table_kmer_counts:
				table_kmer_counts[gene].append([motif,tot_mut,tot_len,ex_len])
			if gene not in table_kmer_counts:
				table_kmer_counts[gene]=[[motif,tot_mut,tot_len,ex_len]]
			if gene not in real_exonic_mutations:
				real_exonic_mutations[gene]=0
			real_exonic_mutations[gene]+=ex_mut
		i+=1
	file.close()
	for gene in table_kmer_counts:
		input_file.append([gene,table_kmer_counts[gene]])
	print("\tTable_kamer_count.txt parsed: %.0f seconds " % (time.time() - start_time))


	# Define funcitons to simulate


	def f2(list):
		tot_mut=list[0]
		tot_len=list[1]
		ex_len=list[2]
		axx=numpy.random.randint(1,tot_len+1,tot_mut).tolist()
		ex_simulated_count= len([x for x in range(len(axx)) if axx[x] <= ex_len])
		return ex_simulated_count
	def simulate(list2):
		numpy.random.seed(seed_number)
		atime=time.time()
		sets=0
		gene=list2[0]
		list2=list2[1]
		results=[0] * simulations
		for list in list2:
			motif=list[0]
			tot_mut=list[1]
			tot_len=list[2]
			ex_len=list[3]
			result=[]
			list=[[tot_mut,tot_len,ex_len]]*simulations
			result=map(f2,list)
			results=[x + y for x, y in zip(results, result)]
		simulated_exonic_mutations=results
		value1=len([i for i in simulated_exonic_mutations if i >= real_exonic_mutations[gene]])
		del simulated_exonic_mutations
		pval=(value1)/float(simulations)
		if pval<=0.1:
			sets+=1
			results=[0] * simulations
			for list in list2:
				motif=list[0]
				tot_mut=list[1]
				tot_len=list[2]
				ex_len=list[3]
				result=[]
				list=[[tot_mut,tot_len,ex_len]]*simulations
				result=map(f2,list)
				results=[x + y for x, y in zip(results, result)]
			simulated_exonic_mutations=results
			value2=len([i for i in simulated_exonic_mutations if i >= real_exonic_mutations[gene]])+value1
			del simulated_exonic_mutations
			pval=(value2)/float(simulations*2)
		if pval<=0.01:
			sets+=1
			for i in range(10):
				results=[0] * simulations*10
				for list in list2:
					motif=list[0]
					tot_mut=list[1]
					tot_len=list[2]
					ex_len=list[3]
					result=[]
					list=[[tot_mut,tot_len,ex_len]]*simulations*10
					result=map(f2,list)
					results=[x + y for x, y in zip(results, result)]
				simulated_exonic_mutations=results
				value2+=len([i for i in simulated_exonic_mutations if i >= real_exonic_mutations[gene]])
				del simulated_exonic_mutations
			value3=value2
			pval=(value3)/float(simulations*102)
		if pval<=0.001:
			sets+=1
			for i in range(100):
				results=[0] * simulations*100
				for list in list2:
					motif=list[0]
					tot_mut=list[1]
					tot_len=list[2]
					ex_len=list[3]
					result=[]
					list=[[tot_mut,tot_len,ex_len]]*simulations*100
					result=map(f2,list)
					results=[x + y for x, y in zip(results, result)]
				simulated_exonic_mutations=results
				value2+=len([i for i in simulated_exonic_mutations if i >= real_exonic_mutations[gene]])
				del simulated_exonic_mutations
			value3=value2
			pval=(value3)/float(simulations*10102)
		line_to_print=gene+"\t"+str(real_exonic_mutations[gene])+"\t"+str(pval)+"\n"
		#if sets>1: print(gene, sets, time.time()-atime, len(list2), pval)
		return(line_to_print)

	# Start simulation
	def simulated_second_round(a):
		value=0
		for i in range(10):
			results=[0] * 100
			for list in table_kmer_counts[gene]:
				motif=list[0]
				tot_mut=list[1]
				tot_len=list[2]
				ex_len=list[3]
				result=[]
				list=[[tot_mut,tot_len,ex_len]]*100
				result=map(f2,list)
				results=[x + y for x, y in zip(results, result)]
			value+=len([i for i in results if i >= real_exonic_mutations[gene]])
		return(value)

	start_time = time.time()
	file2=open("table_probabilities_R.txt","w")
	pool = mp.Pool(processes=cores)
	results = pool.map(simulate, input_file)
	for element in results:
		file2.write(element)
	print("\tSimulation finished: %.0f seconds " % (time.time() - start_time))                                                                                            
	results=None
	file2.close()
	del real_exonic_mutations, simulated_exonic_mutations
	print("\tPvalues calculated: %.0f seconds " % (time.time() - start_time))



# R statistics
subprocess.call('Rscript /home/Exinator2/recurrence_script_statistics_v2.r '+true_set, shell=True)
subprocess.call('mv table_probabilities_R_final.txt results_'+filename+"_"+"_"+str(background_size)+"_"+str(simulations)+"_"+str(seed_number)+'.txt', shell=True)
subprocess.call('mv qqplot.png '+output_folder+'/qqplot_'+filename+"_"+"_"+str(background_size)+"_"+str(simulations)+"_"+str(seed_number)+'.png', shell=True)
subprocess.call('mv venn.png '+output_folder+'/venn_'+filename+"_"+"_"+str(background_size)+"_"+str(simulations)+"_"+str(seed_number)+'.png', shell=True)
subprocess.call('mv precision.png '+output_folder+'/precision_'+filename+"_"+"_"+str(background_size)+"_"+str(simulations)+"_"+str(seed_number)+'.png', shell=True)
print("\tQvalues calculated: %.0f seconds " % (time.time() - start_time))
