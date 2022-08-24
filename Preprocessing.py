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
import matplotlib.pyplot as plt


#*******************************************************************************
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
		data = pd.read_csv(filename, header=header, usecols=["SAMPID","SMTSD"], sep=separator).query('SMTSD == "{}"'.format(P.tissue)) #Selecting only samples of specified tissue
	else:
		data = pd.read_csv(filename, header=header, usecols=(lambda x: lambda_function(x)), sep=separator) #'\t'
	print('Reading complete')
	return data
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
	#df_RNA_seq2 = PrC.df_RNA_seq.append(new_row, ignore_index=True)
	#df_RNA_seq2 = pd.concat([df_RNA_seq2.reindex([len(df_RNA_seq2)-1]) , df_RNA_seq2.iloc[0:len(df_RNA_seq2)-1]])
	df_RNA_seq2 = pd.concat([new_row, PrC.df_RNA_seq], ignore_index=True)
	return df_RNA_seq2

################################################################################
#******************************************************************************#
#******************************************************************************#
#*****************THIS IS WHERE THE CODE STARTS********************************#
#******************************************************************************#
#******************************************************************************#
################################################################################

def preprocessing():

	RNACountsFileName = "RNAseqfile_Tissue="+P.tissue+".csv"
	pathwaygenes = {}
	PrC.df_RNA_seq = pd.DataFrame([]) #emptying dataframe


	help_name = P.GENE_SELECTION #variable for nameing files and folders
	if P.use_middle_age:	help_name += "_With-MiddleAge_METADATA.txt"
	else: 							help_name += "_No-MiddleAge_METADATA.txt"

	#__________________Updating of downloading source files_________________________
	GF.ensure_file('https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz','GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz', P.update_files)
	GF.ensure_file('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', 'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', P.update_files)
	if P.Use_AllSammple_Data: GF.ensure_file('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', P.update_files)
	else: GF.ensure_file('https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_'+P.tissue.lower().replace(" ", "_")+'.gct.gz', 'gene_reads_2017-06-05_v8_'+P.tissue.lower().replace(" ", "_")+'.gct.gz', P.update_files)


	#__________________Reading and converting the datafiles_________________________
	#Reading SUBJECTS dataframe
	PrC.df_subjects = GF.readfile("Source/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", '\t', 0, "dataframe")
	PrC.df_subjects.columns = ['SUBJID','SEX','AGE','DTHHRDY'] #MISSCHIEN WEGHALEN, LIJKT OVERBODIG
	PrC.df_subjects = pd.DataFrame(PrC.df_subjects) #MISSCHIEN WEGHALEN, LIJKT OVERBODIG
	PrC.df_subjects = PrC.df_subjects.set_index('SUBJID') #Setting subjectID as index-column
	PrC.df_subjects.sort_values(by=['AGE','SEX'], inplace=True) #sort first by Age and then by Sex

	#Reading SAMPLES dictionary and only selecting samples from a specific tissue
	if P.Use_AllSammple_Data:
		df_samples = strip_data_file("Source/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 0, 0, '\t')
		lst_samples = df_samples['SAMPID'].tolist() #sample IDs from dataframe to list
		PrC.dict_samples = {lst_samples[i] for i in range(0, len(lst_samples), 1)} #list to dictionary

	#Reading the RNA seq data (already filtering on samples by tissuetype (columns))
	if P.Use_AllSammple_Data:
		#--Option 1) the countset has previously been filtered on samples, use this saved vversion instead for speed
		if ("Countset_Samplefiltered" in os.listdir('./'+P.experiment_name)) and (RNACountsFileName in os.listdir('./'+P.experiment_name+"/Countset_Samplefiltered/")): #Does this folder already exist
			PrC.df_RNA_seq = strip_data_file(P.experiment_name+"/Countset_Samplefiltered/"+RNACountsFileName, read_this_sample, 0, '\t')
		#--Option 2) This file has never been read and saved
		else: PrC.df_RNA_seq = strip_data_file("Source/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", read_this_sample, 2, '\t')
	if not P.Use_AllSammple_Data: #this is faster because the dataframe is much smaller and contains only the desired sample
		PrC.df_RNA_seq = GF.readfile('Source/gene_reads_2017-06-05_v8_'+P.tissue.lower().replace(" ", "_")+'.gct', '\t', 2, 'dataframe').iloc[:, 1:]



	#________________Saving in between for faster runs nexttime_____________________
	if P.Use_AllSammple_Data:
		if "Countset_Samplefiltered" not in os.listdir('./'+P.experiment_name): #The intermediate save folder does not exist jet 
				os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Countset_Samplefiltered"))
				print("Writing ", RNACountsFileName, " to a file for intermediate saving. Next run will be faster.")
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
	
	
	print(PrC.df_RNA_seq)

	#_2______Filtering out middle age if this age group is not required______________
	if not P.use_middle_age and not P.visualize_data:
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



	#_5_____________________Visualizing the data__________________________________
	if P.visualize_data:
		# 1 = male, 2 = female

		ages = {'20-29':{'thisage':0, 'DTHHRDY':[0]*6, 'sex':[0,0]}, 
						'30-39':{'thisage':0, 'DTHHRDY':[0]*6, 'sex':[0,0]}, 
						'40-49':{'thisage':0, 'DTHHRDY':[0]*6, 'sex':[0,0]}, 
						'50-59':{'thisage':0, 'DTHHRDY':[0]*6, 'sex':[0,0]}, 
						'60-69':{'thisage':0, 'DTHHRDY':[0]*6, 'sex':[0,0]}, 
						'70-79':{'thisage':0, 'DTHHRDY':[0]*6, 'sex':[0,0]}}


		for i in PrC.df_RNA_seq:
			ages[PrC.df_RNA_seq.loc['Subject age',i]]['thisage'] += 1
			if PrC.df_RNA_seq.loc['Subject DTHHRDY', i] == "nan": PrC.df_RNA_seq.loc['Subject DTHHRDY', i] = 5
			else: PrC.df_RNA_seq.loc['Subject DTHHRDY', i] = int(float(PrC.df_RNA_seq.loc['Subject DTHHRDY', i]))
			ages[PrC.df_RNA_seq.loc['Subject age',i]]['DTHHRDY'][int(PrC.df_RNA_seq.loc['Subject DTHHRDY', i])] += 1
			ages[PrC.df_RNA_seq.loc['Subject age',i]]['sex'][int(PrC.df_RNA_seq.loc['Subject sex', i])-1] += 1

		#NORMAL AGE DISTRIBUTION
		plt.figure(figsize=[10,10])
		plt.bar(ages.keys(), [k for k in [p['thisage'] for p in ages.values()]], edgecolor = "black")
		plt.title('Nr '+P.tissue+' samples per agegroup', fontsize=20)
		for i in range(0, len(ages.keys())):
			plt.text(i, [k for k in [p['thisage'] for p in ages.values()]][i]+1 ,[k for k in [p['thisage'] for p in ages.values()]][i], color='black', fontweight=600, ha="center", fontsize=18)
		plt.xlabel('Agegroup', fontsize=18)
		plt.ylabel('Nr '+P.tissue+' samples', fontsize=18)
		plt.xticks(fontsize=15)
		plt.yticks(fontsize=15)
		#plt.show()
		plt.savefig(P.experiment_name+'/Normal_Age_distribution.png')
		plt.close()


		#SEX DISTRIBUTION OVER AGE
		plt.figure(figsize=[10,10])
		plt.bar(ages.keys(), [k for k in [p['sex'][0] for p in ages.values()]], color='cornflowerblue', edgecolor = "black")
		plt.bar(ages.keys(), [k for k in [p['sex'][1] for p in ages.values()]], bottom=[k for k in [p['sex'][0] for p in ages.values()]], color='hotpink', edgecolor = "black")
		for i in range(0, len(ages.keys())):
			plt.text(i, [k for k in [p['thisage'] for p in ages.values()]][i]+2 ,[k for k in [p['thisage'] for p in ages.values()]][i], color='black', fontweight=600, ha="center", fontsize=18)
		plt.title('Nr '+P.tissue+' samples of each sex per agegroup', fontsize=20)
		plt.xlabel('Agegroup', fontsize=18)
		plt.ylabel('Nr '+P.tissue+' samples per sex', fontsize=18)
		plt.xticks(fontsize=15)
		plt.yticks(fontsize=15)
		#Legend
		colors = {'Men':'cornflowerblue', 'Women':'hotpink'}         
		labels = list(colors.keys())
		handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
		plt.legend(handles, labels, loc='upper left', fontsize='xx-large')

		#plt.show()
		plt.savefig(P.experiment_name+'/Sex_Distribution_over_Age.png', bbox_inches = 'tight')
		plt.close()


		#DTHHRDY DISTRIBUTION OVER AGE
		plt.figure(figsize=[10, 10])
		colors = ['#7A67EE', '#FF8247', '#4EEE94', '#D8BFD8', '#EE0000', '#B7B7B7'] 
		counthelp = [k for k in [p['DTHHRDY'][0] for p in ages.values()]]
		plt.bar(ages.keys(), [k for k in [p['DTHHRDY'][0] for p in ages.values()]], color=colors[0], edgecolor = "black")
		for i in range(1, len(ages.keys())):
			plt.bar(ages.keys(), [k for k in [p['DTHHRDY'][i] for p in ages.values()]], bottom=counthelp, color=colors[i], edgecolor = "black")
			counthelp = np.add(counthelp, [k for k in [p['DTHHRDY'][i] for p in ages.values()]])
		for i in range(0, len(ages.keys())):
			plt.text(i, [k for k in [p['thisage'] for p in ages.values()]][i]+2 ,[k for k in [p['thisage'] for p in ages.values()]][i], color='black', fontweight=600, ha="center", fontsize=18)
		plt.title('Nr '+P.tissue+' samples of each sex per agegroup', fontsize=20)
		plt.xlabel('Agegroup', fontsize=18)
		plt.ylabel('Nr '+P.tissue+' samples per sex', fontsize=18)
		plt.xticks(fontsize=15)
		plt.yticks(fontsize=15)
		#Legend
		colors = {'Ventilator Case': colors[0], 'Violent & fast': colors[1], 'Fast of natural cause':colors[2], 'Intermediate':colors[3], 'Slow':colors[4], 'Unknown':colors[5]}         
		labels = list(colors.keys())
		handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
		plt.legend(handles, labels, loc='upper left', fontsize='xx-large')
		#plt.show()
		plt.savefig(P.experiment_name+'/Death Circumstances_DIstribution_over_Age.png', bbox_inches = 'tight')
		plt.close()


		print(PrC.df_RNA_seq)
		smtsd = ["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Adrenal Gland", "Artery - Aorta", "Artery - Coronary", "Artery - Tibial", "Bladder", "Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate (basal ganglia)", "Brain - Cerebellar Hemisphere", "Brain - Cerebellum", "Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Hippocampus", "Brain - Hypothalamus", "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)", "Brain - Substantia nigra", "Brain - Spinal cord (cervical c-1)", "Brain - Substantia nigra", "Breast - Mammary Tissue", "Cervix - Ectocervix", "Cervix - Endocervix", "Colon - Transverse", "Esophagus - Mucosa", "Esophagus - Muscularis", "Fallopian Tube", "Heart - Atrial Appendage", "Heart - Left Ventricle", "Kidney - Cortex", "Kidney - Medulla", "Liver", "Lung", "Muscle - Skeletal", "Nerve - Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)", "Spleen", "Stomach", "Testis", "Thyroid", "Transformed fibroblasts", "Uterus", "Vagina", "Whole Blood"]
		smts = ["Adipose Tissue", "Adrenal Gland", "Bladder", "Blood", "Blood Vessel", "Bone Marrow", "Brain", "Breast", "Cervix Uteri", "Colon", "Esophagus", "Fallopian Tube", "Heart", "Kidney", "Liver", "Lung", "Muscle", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina"]



		fik = GF.readfile("Source/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", '\t', 0, 																																																						"dataframe")
		fik = fik.loc[:, ['SAMPID', 'SMTS', 'SMTSD', 'SMAFRZE']]
		fik = fik.loc[fik['SMAFRZE'] == 'RNASEQ']
		print(fik)
		
		samples_per_tissue = {}
		for tis in smtsd:
			samples_per_tissue[tis] = fik.loc[fik['SMTSD																																																'] == tis].shape[0]
		print(fik)
		
		samples_per_tissue = {k: v for k, v in sorted(samples_per_tissue.																																items(), key=lambda item: item[1], reverse=True)}
		print(samples_per_tissue)
		#exit()																																																																				
		for i in samples_per_tissue:
			print(i, " - ", samples_per_tissue[i])
		
		# Samples per detailed tissue
		plt.figure(figsize=[20,10])
		plt.bar(samples_per_tissue.keys(),samples_per_tissue.values(), edgecolor = "black")
		plt.title('Nr samples per detailed tissuetype', fontsize=20)
		for i in range(0, len(samples_per_tissue.keys())):
			#plt.text(i, [k for k in samples_per_tissue[i].values()][i]+1 , [k for k in samples_per_tissue[i].values()][i], color='black', fontweight=1000, ha="center")
			plt.text(i, list(samples_per_tissue.values())[i]+10, list(samples_per_tissue.values())[i], color='black', fontweight=1000, ha="center", fontsize=13)
		#plt.xlabel('Detailed tissuetypes', fontsize=22)
		plt.xticks(rotation=-90)
		plt.ylabel('Nr samples', fontsize=22)
		
		
		plt.gcf().canvas.draw() #PLAGIAAT
		tl = plt.gca().get_xticklabels() #PLAGIAAT
		maxsize = max([t.get_window_extent().width for t in tl]) #PLAGIAAT
		s = maxsize/plt.gcf().dpi*150+4*0.1 #PLAGIAAT
		margin =  0.031
		plt.gcf().subplots_adjust(left=margin, right=1-margin) #PLAGIAAT
		plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1]) #PLAGIAAT
		
		plt.margins(x=0)# or ax.margins(x=0)
		plt.xticks(fontsize=19)
		plt.yticks(fontsize=19)
		
		
		#plt.rcParams["figure.figsize"] = [10, 20]
		#plt.rcParams["figure.autolayout"] = True
		
		
		#plt.show()
		plt.savefig(P.experiment_name+'/Detailed_Samples.png', bbox_inches = 'tight')
		plt.close()




		while True:
			choice = input("You have just visualized the data and saved barplots in your folder where the code is. Visualizing meant not excluding the middle group are you sure you want to continue? (yes/no)")
			if choice == "no": exit()
			elif choice == "yes": break
		




