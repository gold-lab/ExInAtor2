# coding: utf-8

# Load modules

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
import tabix
start_time = time.time()
seed_number=12
random.seed(seed_number)
numpy.random.seed(seed_number)

# Load arguments

if len(sys.argv)==1 or "-h" in sys.argv: 
	print("Mandatory arguments: -m <--mutations> -o <--output_folder> -f <--fasta_file> -g <--gtf_file> -z <--chr_sizes> -s <--scores_file>")
	print("Optional arguments: -i <--iterations> -c <--cores> -t <--true_set> -e <--exc_regions>")
	sys.exit()
n=1
while n<len(sys.argv):
	arg = sys.argv[n]
	if arg == "-m" or arg == "--mutations_file": # BED file with mutations.
		mutations_file = sys.argv[n+1]
	elif arg == "-o" or arg == "--output_folder": # Folder where all output files will be saved.
		output_folder = sys.argv[n+1]
	elif arg == "-f" or arg == "--fasta_file": # FASTA file of the human genome.
		fasta_file = sys.argv[n+1]
	elif arg == "-g" or arg == "--gtf_file": # GTF file of the Gencode annotation to use when creating exons.
		gtf_file = sys.argv[n+1]
	elif arg == "-i" or arg == "--iterations": # Number of iterations to perform simulations, 10000 recommended.
		iterations = int(sys.argv[n+1])
	elif arg == "-z" or arg == "--chr_sizes": # Two column file with name and lengh of the chromosomes. These lengths have to match the ones in the FASTA file indicated previously.
		chr_sizes = sys.argv[n+1]
	elif arg == "-c" or arg == "--cores": # Number of cores in the computer to use for the simulations.
		cores = int(sys.argv[n+1])
	elif arg == "-s" or arg == "--scores_file": # TSV file with scores having these five columns: #Chr    Pos     Ref     Alt     RawScore        PHRED
		scores_file = sys.argv[n+1]
	elif arg == "-t" or arg == "--true_set": # This is a TXT file with ENSGs of your true positives genes, that will be used in the R script.
		true_set = sys.argv[n+1]
	elif arg == "-e" or arg == "--exc_regions": # BED file with regions from the genome to ignore (such as those with low mappability, high repetitive sequences, etc).
		exc_regions = sys.argv[n+1]	
	n+=2
if 'mutations_file' not in globals():
	print("Please indicate \"-m\" or \"--mutations_file\"")
	print("type \"exinator3.py -h\" for help")
	sys.exit()
if 'output_folder' not in globals():
	print("Please indicate \"-o\" or \"--output_folder\"")
	print("type \"exinator3.py -h\" for help")
	sys.exit()
if 'fasta_file' not in globals():
	print("Please indicate \"-f\" or \"--fasta_file\"")
	print("type \"exinator3.py -h\" for help")
	sys.exit()
if 'gtf_file' not in globals():
	print("Please indicate \"-g\" or \"--gtf_file\"")
	print("type \"exinator3.py -h\" for help")
	sys.exit()
if 'scores_file' not in globals():
	print("Please indicate \"-s\" or \"--scores_file\"")
	print("type \"exinator3.py -h\" for help")
	sys.exit()
if 'chr_sizes' not in globals():
	print("Please indicate \"-z\" or \"--chr_sizes\"")
	print("type \"exinator3.py -h\" for help")
	sys.exit()
if 'true_set' not in globals():
	print("No true set provided, no venn or precision plots will be generated")
if 'iterations' not in globals():
	iterations=10000
	print("No number of iterations provided,",iterations,"will be run")
if 'cores' not in globals():
	cores=6
	print("No cores provided,",cores,"will be used")
if 'exc_regions' not in globals():
	exc_regions="NO"
	print("No regions in the Genome will be ignored")	
print("Arguments correctly inserted: %.0f seconds " % (time.time() - start_time))
del arg

# Initial setup: creates output folder and sets the path of the input folder to the one where the mutation bed file is present

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
os.chdir(output_folder)
input_folder="/".join(mutations_file.split("/")[:-1])+"/"
mutations_file=mutations_file.split("/")[-1]
print("Starting the analysis of "+mutations_file.split(".")[0]+": %.0f seconds " % (time.time() - start_time))
if not os.path.exists(output_folder+"/"+mutations_file.split(".")[0]):
	os.makedirs(output_folder+"/"+mutations_file.split(".")[0])
filename=mutations_file.split(".")[0]
os.chdir(output_folder+"/"+mutations_file.split(".")[0])
resetdir=output_folder+"/"+mutations_file.split(".")[0] 
mutations_file=input_folder+mutations_file

# Create merged exons files and file that contains regions to obtain CADD scores for

if not os.path.isfile("merged_exons_mutations.bed"):
	
	subprocess.call("cat "+gtf_file+' | tr --delete \; | tr --delete \\" | grep -v \\# | awk \'$3==\"exon\"\' | awk \'{split($10,var,".");print var[1]\"\t\"$4-1\"\t\"$5\"\t\"$1}\' | sortBed -i - > exons_n.bed', shell=True)
	subprocess.call("mergeBed -i exons_n.bed -c 4 -o distinct | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1}\' | sortBed -i - > merged_exons_n.bed", shell=True)
	subprocess.call('cat exons_n.bed | awk \'{print $4\"\t\"$2\"\t\"$3\"\t\"$1}\' | sortBed -i - > exons.bed', shell=True)
	subprocess.call('cat merged_exons_n.bed | sortBed -i - > merged_exons.bed', shell=True)
	
	if exc_regions!="NO":
		subprocess.call("subtractBed -a merged_exons.bed -b "+exc_regions+" | sortBed -i - > merged_exons_clean.bed", shell=True)
		subprocess.call("mv merged_exons_clean.bed merged_exons.bed", shell=True)
	
	subprocess.call('awk \'{print $1":"$2"-"$3"\t"$4}\' merged_exons.bed | sed \'s/chr//g\' > regions_to_get_scores.txt',shell=True)
	print("\tBED files created: %.0f seconds " % (time.time() - start_time))
	

# Count mutations in exonic regions

	subprocess.call('intersectBed -a merged_exons.bed -b '+mutations_file+' -sorted -wb | awk \'{print $1"\t"$2-1"\t"$3+1"\t"$4}\' | awk \'$3-$2==3\' > merged_exons_mutations.bed', shell=True)
	print("\tMutations counted: %.0f seconds " % (time.time() - start_time))

# Get CADD scores for exonic regions 

if not os.path.isfile(output_folder+"/merged_exons_scores.txt"):
	
	regions=open("regions_to_get_scores.txt","r")
	tb = tabix.open(scores_file)
	out_file=open("merged_exons_scores.bed","w")
	base=""
	sum_value=[]
	for line in regions:
		records = tb.querys(line.strip().split("\t")[0])
		gene=line.strip().split("\t")[1]
		for record in records:
			new_chro,new_pos,wild,mut,value,value2=record
			new_base=str(new_chro)+str(new_pos)
			if base=="": 
				base=new_base
				chro=new_chro
				pos=new_pos
			if base==new_base: sum_value.append(float(value))
			else: 
				out_file.write("chr"+str(chro)+"\t"+str(int(pos)-1)+"\t"+str(pos)+"\t"+gene+"\t"+str(sum(sum_value)/len(sum_value))+"\n")
				base=new_base
				chro=new_chro
				pos=new_pos
				sum_value=[float(value)]
	out_file.write("chr"+str(chro)+"\t"+str(int(pos)-1)+"\t"+str(pos)+"\t"+gene+"\t"+str(sum(sum_value)/len(sum_value))+"\n")
	regions.close()
	out_file.close()			

	
	#subprocess.call('xargs -a regions_to_get_scores.txt -I {} tabix '+scores_file+' {} | awk \'{print "chr"$1"\t"$2-1"\t"$2"\t"$5}\'  > merged_exons_scores.bed', shell=True)
	print("\tScores gathered: %.0f seconds " % (time.time() - start_time))

# Get fasta files and extract sequences of exonic regions

	subprocess.call("bedtools slop -i merged_exons_scores.bed -g "+chr_sizes+" -b 1 | sed 's/chr//g' "+'| awk \'{print $1"\t"$2"\t"$3"\t"$1"&&"$2"&&"$3"&&"$4"&&"$5}\''+" > merged_exons_scores_for_fasta.bed", shell=True)
	#subprocess.call('sed \'s/chr//g\' merged_exons_scores_for_fasta.bed | awk \'{print $1"\t"$2"\t"$3"\t"$1"&&"$2"&&"$3"&&"$4"&&"$5}\' > merged_exons_scores_for_fasta2.bed', shell=True)
	subprocess.call('bedtools getfasta -fi '+fasta_file+' -bed merged_exons_scores_for_fasta.bed -fo merged_exons_scores.fa -name', shell=True)
	print("\tFasta created: %.0f seconds " % (time.time() - start_time))

# Prepare for iterations

	input_fasta=open("merged_exons_scores.fa","r")
	output_presim=open("pre_simulation.txt","w")
	for line in input_fasta:
		if ">" in line:
			line=line.rstrip().split("&&")
			chrm=line[0].replace(">", "chr")
			start=int(line[1])
			end=int(line[2])
			score=line[4].split("::")[0]
			gene=line[3]
		else:
			line2=chrm+"&&"+str(start)+"&&"+str(end)+"&&"+gene+"\t"+str(score)+"\t"+line.rstrip()+"\n"
			output_presim.write(line2)
	input_fasta.close()
	output_presim.close()
	subprocess.call('sort pre_simulation.txt > '+output_folder+'/pre_simulation2.txt ', shell=True)
	subprocess.call("sed 's/\&\&/\t/g' "+output_folder+'/pre_simulation2.txt | awk \'{print $4"\t"$5"\t"$6}\' > '+output_folder+'/merged_exons_scores.txt ', shell=True)

if not os.path.isfile("merged_exons_mutations_scores.txt"):	
	subprocess.call("sed 's/\t/\&\&/g' merged_exons_mutations.bed | sort > merged_exons_mutations.txt", shell=True)
	subprocess.call('join merged_exons_mutations.txt '+output_folder+'/pre_simulation2.txt | sed \'s/&&/\t/g\' | awk \'{print $4"\t"$5"\t"$6}\' > merged_exons_mutations_scores.txt ', shell=True)
	
	print("\tPrepared for iterations: %.0f seconds " % (time.time() - start_time))

if not os.path.isfile("results_"+filename+"_"+str(iterations)+"_"+str(seed_number)+'.txt'):
	
	# Load files for iterations

	input_mut=open("merged_exons_mutations_scores.txt","r")
	mutations={}
	genes=[]
	for line in input_mut:
		line=line.rstrip().split("\t")
		if "NA" not in line[1]:
			gene=str(line[0])
			score=float(line[1])
			tri=str(line[2])
			if gene not in genes: genes.append(gene)
			if gene not in mutations: mutations[gene]={}
			if tri not in mutations[gene]: mutations[gene][tri]=[]
			mutations[gene][tri].append(score)
	input_mut.close()

	input_nuc=open(output_folder+"/merged_exons_scores.txt","r")
	nucleotides={}
	for line in input_nuc:
		line=line.rstrip().split("\t")
		if "NA" not in line[1]:
			gene=str(line[0])
			score=float(line[1])
			tri=str(line[2])
			if gene not in nucleotides: nucleotides[gene]={}
			if tri not in nucleotides[gene]: nucleotides[gene][tri]=[]
			nucleotides[gene][tri].append(score)
	input_nuc.close()

	simulation_data=[]
	for gene in genes:
		simulation_data.append([gene,mutations[gene],nucleotides[gene]])
	print("\tFiles loaded for iterations: %.0f seconds " % (time.time() - start_time))

	# Start iterations

	def simulate(list_sim):
		numpy.random.seed(seed_number)
		gene=list_sim[0]
		mutations=list_sim[1]
		nucleotides=list_sim[2]
		real_values=[]
		simulated_means=[]
		real_mean="NA"
		pval="NA"
		start_time2 = time.time()
		
		
		if "gsdgfdsgfdgdfs" not in gene:
			for i in range(iterations):
				simulated_values=[]
				for tri in mutations:
					if i==0: real_values.extend(mutations[tri])
					tri_mut_len=len(mutations[tri])
					tri_nuc_values=numpy.array(nucleotides[tri])
					simulated_values.extend(tri_nuc_values[numpy.random.randint(0,len(nucleotides[tri]),tri_mut_len)].tolist())
				simulated_means.append(sum(simulated_values)/len(simulated_values))
			real_mean=sum(real_values)/float(len(real_values))
			value1=len([i for i in simulated_means if i >= real_mean])
			pval=(value1)/float(len(simulated_means))
			
			if pval<=0.1:
				#print(gene, 0.001)
				for i in range(iterations*10):
					simulated_values=[]
					for tri in mutations:
						tri_mut_len=len(mutations[tri])
						tri_nuc_values=numpy.array(nucleotides[tri])
						simulated_values.extend(tri_nuc_values[numpy.random.randint(0,len(nucleotides[tri]),tri_mut_len)].tolist())
					simulated_means.append(sum(simulated_values)/len(simulated_values))
				value1=len([i for i in simulated_means if i >= real_mean])
				pval=(value1)/float(len(simulated_means))
				
			if pval<=0.01:
				#print(gene, 0.0001)
				for i in range(iterations*100):
					simulated_values=[]
					for tri in mutations:
						tri_mut_len=len(mutations[tri])
						tri_nuc_values=numpy.array(nucleotides[tri])
						simulated_values.extend(tri_nuc_values[numpy.random.randint(0,len(nucleotides[tri]),tri_mut_len)].tolist())
					simulated_means.append(sum(simulated_values)/len(simulated_values))
				value1=len([i for i in simulated_means if i >= real_mean])
				pval=(value1)/float(len(simulated_means))
			
			if pval<=0.001:
				#print(gene, 0.00001)
				for i in range(iterations*1000):
					simulated_values=[]
					for tri in mutations:
						tri_mut_len=len(mutations[tri])
						tri_nuc_values=numpy.array(nucleotides[tri])
						simulated_values.extend(tri_nuc_values[numpy.random.randint(0,len(nucleotides[tri]),tri_mut_len)].tolist())
					simulated_means.append(sum(simulated_values)/len(simulated_values))
				value1=len([i for i in simulated_means if i >= real_mean])
				pval=(value1)/float(len(simulated_means))
				
			if pval<=0.0001:
				#print(gene, 0.00001)
				for i in range(iterations*10000):
					simulated_values=[]
					for tri in mutations:
						tri_mut_len=len(mutations[tri])
						tri_nuc_values=numpy.array(nucleotides[tri])
						simulated_values.extend(tri_nuc_values[numpy.random.randint(0,len(nucleotides[tri]),tri_mut_len)].tolist())
					simulated_means.append(sum(simulated_values)/len(simulated_values))
				value1=len([i for i in simulated_means if i >= real_mean])
				pval=(value1)/float(len(simulated_means))
		
		result=[gene,len(real_values),real_mean,value1,float(len(simulated_means)),pval]
		return(result)
	pool = mp.Pool(processes=cores)
	results = pool.map(simulate, simulation_data)
	file2=open("table_probabilities_R.txt","w")
	for element in results:
		line_to_print=element[0]+"\t"+str(element[1])+"\t"+str(element[2])+"\t"+str(element[3])+"\t"+str(element[4])+"\t"+str(element[5])+"\n"
		file2.write(line_to_print)
	file2.close()
	print("\titerations finished: %.0f seconds " % (time.time() - start_time))       

# Run R script

# For setting the directory for R script 


subprocess.call('Rscript /home/Exinator2/functionalimpact_script_statistics.r '+true_set, shell=True)
subprocess.call('mv table_probabilities_R_final.txt results_'+filename+"_"+str(iterations)+"_"+str(seed_number)+'.txt', shell=True)
subprocess.call('mv qqplot.png '+output_folder+'/qqplot_'+filename+"_"+str(iterations)+"_"+str(seed_number)+'.png', shell=True)
subprocess.call('mv venn.png '+output_folder+'/venn_'+filename+"_"+str(iterations)+"_"+str(seed_number)+'.png', shell=True)
subprocess.call('mv precision.png '+output_folder+'/precision_'+filename+"_"+str(iterations)+"_"+str(seed_number)+'.png', shell=True)

print("\tR script ran: %.0f seconds " % (time.time() - start_time))       
