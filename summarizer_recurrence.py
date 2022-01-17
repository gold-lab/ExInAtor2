#!/usr/bin/python

import sys
import os

# Folder where the results per sample are saved
path = sys.argv[1]
simulations = "2000"

# Gather all the genes per sample with the Qvalues
all_gene_list={}
all_gene_list2={}
ID_list=[]
for sample in os.listdir(path):
	if "." not in sample:
		for file in os.listdir(path+sample+"/"):
			if str(simulations)+"_12" in file:
				all_gene_list[sample]={}
				all_gene_list2[sample]={}
				file = open (path+sample+"/"+file,"r")
				next(file)
				for line in file:
					line=line.rstrip()
					gene=line.split(" ")[0]
					if "::" in gene: gene=gene.split("::")[3]
					if "." in gene: gene=gene.split(".")[0]
					all_gene_list[sample][gene]=float(line.split(" ")[3])
					all_gene_list2[sample][gene]=float(line.split(" ")[2])
					if gene not in ID_list:
						ID_list.append(gene)
			
# Make the summary, joining the results per sample, per gene and also joining all the genes in each sample in one file
gene_list = []
gene_list2 = []
total_counter={}
percentage_counter={}
for gene in ID_list:
	total_counter[gene]={}
	total_counter[gene]["0.1"]=0
	total_counter[gene]["0.01"]=0
	total_counter[gene]["0.001"]=0
	total_counter[gene]["0.0001"]=0
	for sample in all_gene_list:
		if sample  not in percentage_counter:
			percentage_counter[sample]={}
			percentage_counter[sample]["0.1"]=0
			percentage_counter[sample]["0.01"]=0
			percentage_counter[sample]["0.001"]=0
			percentage_counter[sample]["0.0001"]=0
		if gene in all_gene_list[sample]:
				if float(all_gene_list[sample][gene]) <= 0.1:
					total_counter[gene]["0.1"]+=1
					percentage_counter[sample]["0.1"]+=1
				if float(all_gene_list[sample][gene]) <= 0.01:
					total_counter[gene]["0.01"]+=1
					percentage_counter[sample]["0.01"]+=1
				if float(all_gene_list[sample][gene]) <= 0.001:
					total_counter[gene]["0.001"]+=1
					percentage_counter[sample]["0.001"]+=1
				if float(all_gene_list[sample][gene]) <= 0.0001:
					total_counter[gene]["0.0001"]+=1
					percentage_counter[sample]["0.0001"]+=1
				gene_list.append([gene,sample,all_gene_list[sample][gene]])
				gene_list2.append([gene,sample,all_gene_list2[sample][gene]])
		else:
			gene_list.append([gene,sample,"NA"])
			gene_list2.append([gene,sample,"NA"])

# Print the results for all genes in all samples
out_file =  open (path+"all_genes_combined_long_"+simulations+"_12.txt","w")
out_file.write("Ensemble_id\tsample\tqvalue\n")

for gene in gene_list:
	line_to_print= gene[0]+"\t"+gene[1]+"\t"+str(gene[2])+"\n"
	out_file.write(line_to_print)
	
out_file =  open (path+"all_genes_combined_long_"+simulations+"_12_pvalue.txt","w")
out_file.write("Ensemble_id\tsample\tpvalue\n")

for gene in gene_list2:
	line_to_print= gene[0]+"\t"+gene[1]+"\t"+str(gene[2])+"\n"
	out_file.write(line_to_print)

# Print the results per each gene, this is, the number of times each gene is detected in all the samples
out_file1 =  open (path+"summary_per_gene_long_"+simulations+"_12.txt","w")
out_file1.write("Ensemble_id\tsamples_with_Q-value_<=0.1\tsamples_with_Q-value_<=0.01\tsamples_with_Q-value_<=0.001\tsamples_with_Q-value_<=0.0001\n")

for gene in sorted(total_counter, key=lambda x:total_counter[x]['0.1'],reverse=True):
	line_to_print= gene+"\t"+str(total_counter[gene]["0.1"])+"\t"+str(total_counter[gene]["0.01"])+"\t"+str(total_counter[gene]["0.001"])+"\t"+str(total_counter[gene]["0.0001"])+"\n"
	out_file1.write(line_to_print)

# Print the results per each sample, this is, the number of genes that are detected in each sample	
out_file2 =  open (path+"summary_per_sample_long_"+simulations+"_12.txt","w")
out_file2.write("Sample\tQ-value_<=0.1\tQ-value_<=0.01\tQ-value_<=0.001\tQ-value_<=0.0001\n")

for sample in sorted(percentage_counter, key=lambda x:percentage_counter[x]['0.1'],reverse=True):
	line_to_print= sample+"\t"+str(percentage_counter[sample]["0.1"])+"\t"+str(percentage_counter[sample]["0.01"])+"\t"+str(percentage_counter[sample]["0.001"])+"\t"+str(percentage_counter[sample]["0.0001"])+"\n"
	out_file2.write(line_to_print)
