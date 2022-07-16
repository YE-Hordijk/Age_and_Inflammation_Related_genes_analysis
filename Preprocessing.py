#Preprocessing.py

#Automatic updating
import requests
import gzip
from zipfile import ZipFile
import shutil
import os
import re
from datetime import datetime
import psutil #for checking if enough diskspace
import operator #for sorting

#General packages
import GeneralFunctions as GF
from GeneralFunctions import st
from Parameters import P

import numpy as np
from numpy import log as ln
import pandas as pd
import math as mt
import statistics
from itertools import islice #for selecting dictionary items
#import sys



#*******************************************************************************


#df_RNA_seq = pd.DataFrame([])

class PrC:
	dict_samples = {}
	df_subjects = pd.DataFrame([])
	df_RNA_seq = pd.DataFrame([])

#*******************************************************************************
#*****************************Functions*****************************************
#*******************************************************************************

#function that decides if this column contains a sample of interest and thus must be read
def read_this_sample(index):
	if index in PrC.dict_samples or index == "Name" or index == "Description":
		return True
	return False
#*******************************************************************************
#function that only reads and stores the columns of interest in a pandas dataframe
def strip_data_file(filename, lambda_function, header, separator):
	print("\nReading", filename ,"...")
	if (lambda_function == 0):
		data = pd.read_csv(filename, header=header, usecols=["SAMPID","SMTSD"], sep=separator).query('SMTSD == "{}"'.format(P.tissue)) #Change the sample type here
	else:
		data = pd.read_csv(filename, header=header, usecols=(lambda x: lambda_function(x)), sep=separator) #'\t'
	print('Reading complete')
	return data
	


#******************************************************************************* ->preprocessing
"""
#function that only reads and stores the columns of interest in a pandas dataframe
def strip_data_file(filename, lambda_function, header, separator):
	print("\nReading", filename ,"...")
	if (lambda_function == 0):
		data = pd.read_csv(filename, header=header, usecols=["SAMPID","SMTSD"], sep=separator).query('SMTSD == "{}"'.format(P.tissue)) #Change the sample type here
	else:
		data = pd.read_csv(filename, header=header, usecols=(lambda x: lambda_function(x)), sep=separator) #'\t'
	print('Reading complete')
	return data
"""
#*******************************************************************************
#Creating a file with metadata about the samples (linking age to sample) for making a PCA plot in R
def make_metadata(RNAseq, subjects):
	temp_subjects = PrC.df_subjects["Age".upper()].to_dict() #making dictionary where SubjectID is key and Age is value

	metadata = {}
	"""
	for column in RNAseq:
		if ((column != "Name") and (column != "Description")):
			col = column.split("-")[0]+"-"+column.split("-")[1] #extracting the SubjectID from the SampleName
			if temp_subjects[col] <= '40-49':		metadata[column] = "Young (20-49)"	#is actually: 50-59
			elif temp_subjects[col] >= '60-69':	metadata[column] = "Old (60-79)"		#is actually: 50-59
			elif P.use_middle_age:								metadata[column] = "Middle (50-59)"
	return metadata
	"""
	max_young = P.YOUNG[-2:]
	min_old = P.OLD[:2]
	
	for column in RNAseq:
		if ((column != "Name") and (column != "Description")):
			col = column.split("-")[0]+"-"+column.split("-")[1] #extracting the SubjectID from the SampleName
			if temp_subjects[col] <= str(max_young[0]+'0-'+max_young):		metadata[column] = "Young ("+P.YOUNG+")"	#is actually: 50-59
			elif temp_subjects[col] >= str(min_old+'-'+min_old[0]+'9'):	metadata[column] = "Old ("+P.OLD+")"		#is actually: 50-59
			elif P.use_middle_age:								metadata[column] = "Middle ("+P.MIDDLE+")"
	
	return metadata
#*******************************************************************************
#This function adds information about a suject to the dataframe
def make_subject_row(info, df_RNA_seq): #Make the row with subject data connected to each sample
	temp_subjects = PrC.df_subjects[info.upper()].to_dict()
	new_row = {'Name':'Subject '+info, 'Description':info}
	for column in df_RNA_seq:
		if ((column != "Name") and (column != "Description")):
			col = column.split("-")[0]+"-"+column.split("-")[1] #extracting the SubjectID from the SampleName
			#col = column[:10]
			#if(column[9] == "-"):
			#	col = column[:9]
			new_row[column] = str(temp_subjects[col])
	return new_row
#*******************************************************************************
def add_new_row(new_row):
	new_row = pd.DataFrame(new_row, index=[0])
	df_RNA_seq2 = PrC.df_RNA_seq.append(new_row, ignore_index=True)
	df_RNA_seq2 = pd.concat([df_RNA_seq2.reindex([len(df_RNA_seq2)-1]) , df_RNA_seq2.iloc[0:len(df_RNA_seq2)-1]])
	return df_RNA_seq2

################################################################################
#******************************************************************************#
#******************************************************************************#
#*****************THIS IS WHERE THE CODE STARTS********************************#
#******************************************************************************#
#******************************************************************************#
################################################################################

def preprocessing():
	#update_files = P.update_files
	#GENE_SELECTION = P.GENE_SELECTION
	#tissue = P.tissue
	#select_on_genes = P.select_on_genes
	#use_middle_age = P.use_middle_age
	RNACountsFileName = "RNAseqfile_Tissue="+P.tissue+".csv"
	pathwaygenes = {}
	
	PrC.df_RNA_seq = pd.DataFrame([]) #emptying dataframe


	help_name = P.GENE_SELECTION #variable for nameing files and folders
	if P.use_middle_age:	help_name += "_With-MiddleAge_METADATA.txt"
	else: 							help_name += "_No-MiddleAge_METADATA.txt"
	#__________________Updating of downloading source files_________________________
	GF.ensure_file('https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz','GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz', P.update_files)
	GF.ensure_file('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', 'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', P.update_files)
	GF.ensure_file('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', P.update_files)


	#__________________Reading and converting the datafiles_________________________
	#Reading SUBJECTS dataframe
	PrC.df_subjects = GF.readfile("Source/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", '\t', 0, "dataframe")
	PrC.df_subjects.columns = ['SUBJID','SEX','AGE','DTHHRDY'] #MISSCHIEN WEGHALEN, LIJKT OVERBODIG
	PrC.df_subjects = pd.DataFrame(PrC.df_subjects) #MISSCHIEN WEGHALEN, LIJKT OVERBODIG
	PrC.df_subjects = PrC.df_subjects.set_index('SUBJID') #Setting subjectID as index-column
	PrC.df_subjects.sort_values(by=['AGE','SEX'], inplace=True) #sort first by Age and then by Sex

	#Reading SAMPLES dictionary and only selecting samples from a specific tissue
	df_samples = strip_data_file("Source/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 0, 0, '\t')
	lst_samples = df_samples['SAMPID'].tolist() #sample IDs from dataframe to list
	PrC.dict_samples = {lst_samples[i] for i in range(0, len(lst_samples), 1)} #list to dictionary

	#Reading the RNA seq data (already filtering on samples by tissuetype (columns))
	#--Option 1) the countset has previously been filtered on samples, use this saved vversion instead for speed
	if ("Countset_Samplefiltered" in os.listdir('./'+P.experiment_name)) and (RNACountsFileName in os.listdir('./'+P.experiment_name+"/Countset_Samplefiltered/")): #Does this folder already exist
		PrC.df_RNA_seq = strip_data_file(P.experiment_name+"/Countset_Samplefiltered/"+RNACountsFileName, read_this_sample, 0, '\t')
	#--Option 2) This file has never been read and saved
	else: PrC.df_RNA_seq = strip_data_file("Source/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", read_this_sample, 2, '\t')

	print(PrC.df_RNA_seq)
	


	#________________Saving in between for faster runs nexttime_____________________
	if "Countset_Samplefiltered" not in os.listdir('./'+P.experiment_name): #Does this folder already exist
			os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Countset_Samplefiltered for intermediate saving. Next run will be faster."))
			print("Writing ", RNACountsFileName, " to a file")
			PrC.df_RNA_seq.to_csv("./"+P.experiment_name+"/Countset_Samplefiltered/"+RNACountsFileName,index=False, sep='\t')
	elif RNACountsFileName not in os.listdir('./'+P.experiment_name+"/Countset_Samplefiltered/"): #File does not yet exist
		print("Writing ", RNACountsFileName, " to a file for intermediate saving. Next run will be faster.")
		PrC.df_RNA_seq.to_csv("./"+P.experiment_name+"/Countset_Samplefiltered/"+RNACountsFileName,index=False, sep='\t')


	#________________Reading the genelist___________________________________________
	if ("Genelists" not in os.listdir('./'+P.experiment_name)) and (P.GENE_SELECTION+"_Genelist.txt" not in os.listdir('./'+P.experiment_name+"Genelists")): #checking if file exists
		exit("Genelist is required but does not exist. Provide a gene list ",P.experiment_name+"/"+P.GENE_SELECTION+"_Genelist.txt")
	else: #If sure the file exists
		f = open("./"+P.experiment_name+"/Genelists/"+P.GENE_SELECTION+"_Genelist.txt", "r+")
		for i in f: pathwaygenes[i.split("\n")[0]] = None
	
	""" #CODE FOR CHECKING TE DOUBLE GENES IN TEH DATASET
	geneocc = {}
	for i in PrC.df_RNA_seq["Description"]:
		if i in geneocc:
			geneocc[i] += 1
		else: geneocc[i] = 1

	#geneocc2 = {k: v for k, v in sorted(geneocc.items(), key=lambda item: item[1])}
	geneocc2 = dict( sorted(geneocc.items(), key=operator.itemgetter(1),reverse=True))
	#print(geneocc2)
	tel = 0
	for i in geneocc2:
		#print(i, ": ", geneocc2[i])
		tel += 1
		if tel > 198: break
	
	nroccurences = 0
	tell = 0
	print("<__> geneselection: ", P.GENE_SELECTION)
	print("<__> number of gene in pathwyagenes: ", len(pathwaygenes))
	for i in pathwaygenes:
		tell += 1
		if i in geneocc2:
			nroccurences += geneocc2[i]
			if geneocc2[i] > 1:
				input(i+": "+str(geneocc2[i]))
	
	print(nroccurences)
	print("Number of genes in pathwaygenes: ", tell)
	"""
	#________________Filtering on gene list (rows)__________________________________

	if (P.select_on_genes):
		print("Filtering on genes")
		PrC.df_RNA_seq = PrC.df_RNA_seq[PrC.df_RNA_seq.Description.isin(pathwaygenes.keys())]
		print("Filtering on genes COMPLETE")

	print(st.YELLOW,st.BOLD,PrC.df_RNA_seq.shape[0],"genes left in de dataset", st.RST)


	################################################################################
	############ Preparing Files for R (for normalization and visualization) #######

	#_1___________Create and write metafile (contains info about the samples)________
	if "Files_for_R" not in os.listdir('./'+P.experiment_name):
		os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Files_for_R"))

	print("Writing the metadata...")
	meta = make_metadata(PrC.df_RNA_seq, PrC.df_subjects) #Create a dictionary with sample ids as keys and age groups as values
	df_meta = pd.DataFrame.from_dict(meta, orient="index", columns=['age']) #making dataframe with sampids as index and 1 column "age"
	f = open(P.experiment_name+"/Files_for_R/"+help_name, 'w')
	#print(df_meta)
	f.writelines(str("\tage"+"\n"))
	for key,value in meta.items():
		f.writelines(str(key+"\t"+value+'\n'))
	f.close()

	#_2______Filtering out middle age if this age group is not required______________
	if not P.use_middle_age:
		print("Filtering out the middle age group (", P.MIDDLE,")...")
		for col in PrC.df_RNA_seq.iloc[0:,2:]:
			if ((col != "Name") and (col != "Description") and (col not in meta) ): #has this subject already been filtered out from the metadatafile because middle-age
				PrC.df_RNA_seq.drop(col, axis=1, inplace=True)


	#_3______Adding new rows with subject information_______________________________
	new_row = make_subject_row('sex', PrC.df_RNA_seq)
	PrC.df_RNA_seq = add_new_row(new_row)
	new_row = make_subject_row('DTHHRDY', PrC.df_RNA_seq)
	PrC.df_RNA_seq = add_new_row(new_row)
	new_row = make_subject_row('age', PrC.df_RNA_seq)
	PrC.df_RNA_seq = add_new_row(new_row)

	PrC.df_RNA_seq = PrC.df_RNA_seq.set_index('Name')


	#_4___________Writing readcountsfile for normalization__________________________

	if P.select_on_genes: help_name = help_name[:-12] + "unnormalized.txt" #Changing the help_name for naming files and folders
	else: 
		if P.use_middle_age: tempvar = "_With-MiddleAge"
		else: tempvar = "_No-MiddleAge"
		help_name = " NoGeneSelection"+tempvar+"_unnormalized.txt" #no selection on genes
	PrC.df_RNA_seq.drop(['Description'], axis=1, inplace=True)
	print("Writing", help_name, "to a file")
	PrC.df_RNA_seq.iloc[3:, 0:].to_csv(P.experiment_name+"/Files_for_R/"+help_name,index=True, sep='\t')
	




