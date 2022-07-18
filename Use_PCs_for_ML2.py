#Use_PCs_for_ML2.py

#General packages
import GeneralFunctions as GF
from GeneralFunctions import st
from Parameters import P

import numpy as np
from numpy import log as ln
import os
import math
import pandas as pd
import math as mt
import statistics
from itertools import islice #for selecting dictionary items
import matplotlib
from matplotlib import pyplot as plt
from collections import OrderedDict
#import sys

#For machinelearning
np.random.seed(0)
#Machine learning packages
from sklearn.model_selection import train_test_split
from sklearn import svm, tree, metrics #Support vector machines # DT #for accuracy calculation
from sklearn.ensemble import RandomForestRegressor, BaggingRegressor	#Random Forrest
from sklearn.tree import DecisionTreeClassifier 
from sklearn.tree import DecisionTreeRegressor

from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
import sklearn



#Setting global vriables
help_name = ""
df_subjects = pd.DataFrame([])
df_RNA_seq = pd.DataFrame([])

############################ FUNCTIONS #########################################
def set_vars():
	global help_name
	global df_subjects
	global df_RNA_seq
	help_name = ""
	if P.select_on_genes: help_name = P.GENE_SELECTION #variable for nameing files and folders
	else: help_name = "NoGeneSelection"
	if P.use_middle_age:	help_name += "_With-MiddleAge"
	else: 								help_name += "_No-MiddleAge"

	df_subjects = pd.DataFrame([])
	df_RNA_seq = pd.DataFrame([])



#*******************************************************************************
#******************************* Functions *************************************
#*******************************************************************************
#This function was copied from https://stackoverflow.com/questions/41592661/determining-the-most-contributing-features-for-svm-classifier-in-sklearn
def f_importances(coef, names):
	imp = coef
	imp,names = zip(*sorted(zip(imp,names)))
	from matplotlib.pyplot import figure
	figure(figsize=(6, 9), dpi=80)
	plt.barh(range(len(names)), imp, align='center')
	plt.yticks(range(len(names)), names, fontsize=8)
	plt.savefig(P.experiment_name+'/Compare_outliers/Outlier_Important_genes/'+P.GENE_SELECTION+"["+P.MODEL+"]"+'_FeatureImportances.pdf')

	plt.close()
	plt.show()
	#exit()

#*******************************************************************************
def NameDescription_linking(df):
	temp = {}
	for key, value in df.iterrows():
		temp[value[0]] = value[1]
	f = open("NameDescription.txt", 'w')
	for key in temp:
		f.writelines(str(key+"\t"+temp[key]+'\n'))
	f.close()
	return temp
#*******************************************************************************
#function that decides if this column contains a sample of interest and thus must be read
def read_this_sample(index):
	if index in dict_samples or index == "Name" or index == "Description":
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
#*******************************************************************************
# This funcion is to link the Description to the Name column in the outputed normalized file from R
def integrate_normalized_data(Normalized_and_selected_file_name): 
	#df_normalized_data = GF.readfile("Genes_from_Papers/With MiddleAge/Normalized_genes_from_papers.txt", "\t", 0, "dataframe")
	df_normalized_data = GF.readfile(Normalized_and_selected_file_name, "\t", 0, "dataframe")
	df_normalized_data.rename(columns={'Unnamed: 0':'Name'}, inplace=True )

	if "NameDescription.txt" not in os.listdir('.'): #NameDescription contains the collumns Name and Description
		df = strip_data_file("Source/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", read_this_sample, 2, '\t')
		tempo = NameDescription_linking(df) #making and saving the linking dictionary 
	else: tempo = GF.readfile("NameDescription.txt", "\t", 0, "dictionary") #make a dictionairy with the "Name"-column as keys, and the "Description"-column as values

	df_normalized_data['Description'] = df_normalized_data['Name'].map(lambda x: tempo[x]) #Adding the "Description"-column

	#Reorganizing the columns
	cols = list(df_normalized_data.columns.values) #making a list of all the collumns
	cols.pop(cols.index('Description')) #removing column "Description" from list
	cols.pop(cols.index('Name')) #removing column "Name" from list
	df_normalized_data = df_normalized_data[['Name', 'Description']+cols] #First "Name" then "Description" then the other columns
	
	return df_normalized_data
#*******************************************************************************
def log_adjustment(data):
	ones = pd.DataFrame(int(1), index=data.index, columns=data.columns) #creating dataframe of right size filled with ones
	#dataa2 = pd.DataFrame(dataa+ones)
	data = pd.DataFrame(ones.add(data, fill_value=0))
	data = np.log2(data)
	return data
#*******************************************************************************
#This function adds information about a suject to the dataframe
def make_subject_row(info): #Make the row with subject data connected to each sample
	temp_subjects = df_subjects[info.upper()].to_dict()
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
	df_RNA_seq2 = df_RNA_seq.append(new_row, ignore_index=True)
	df_RNA_seq2 = pd.concat([df_RNA_seq2.reindex([len(df_RNA_seq2)-1]) , df_RNA_seq2.iloc[0:len(df_RNA_seq2)-1]])
	return df_RNA_seq2
#*******************************************************************************
def write_latex_line(w, bold, f):
	if w[1] == 1: w[1] = "Young"
	elif w[1] == 2: w[1] = "Middle"
	elif w[1] == 3: w[1] = "Old"

	if bold:
		f.write("\\begin{table}[H]\n")
		f.write("\t\centering\n\t\small\n")
		#tab_cols = str("c|"*len(w))[:-1] #without side lines
		tab_cols = str("|c"*len(w))+"|" #with side lines
		f.write(str("\t\\begin{}{}\n").format("{tabular}","{"+tab_cols+"}"))
		f.write("\t\t\hline\n")
	f.write('\t\t')
	for i in range(len(w)):
		if not isinstance(w[i], (str)):
			w[i] = round(w[i], 2)
			w[i] = str(w[i])
		if bold: f.write("\\textbf{}".format("{"+w[i]+"}"))
		else: f.write(w[i])
		if (w[i] == w[-1]) and (w[i] != " "): 
			f.write(" \\\\")
			if bold:	f.write(" \hline\n")
			else:			f.write("\n")
		else: f.write(" & ")

#*******************************************************************************
def barplot_expressiondifference(Genes, ExpressionDifference, Agegroup, Threshold):
	plt.rcParams["figure.figsize"] = [10, 8]
	plt.rcParams["figure.autolayout"] = True
	plt.bar(Genes, ExpressionDifference)
	plt.title('Expression difference between normal '+Agegroup+' and outlier '+Agegroup+'\n (Negative values mean the normal group is downregulated compared to the outlier group)')
	plt.xlabel('Genes in dataset \''+P.GENE_SELECTION+"\'")
	plt.ylabel('Expression difference')
	
	#plt.xticks(rotation=90, fontsize=2)#ticks=None, labels=None)
	plt.tick_params( #removing the xticks (https://stackoverflow.com/questions/12998430/remove-xticks-in-a-matplotlib-plot)
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=False) # labels along the bottom edge are off
	plt.xticks(ticks=None, labels=None)
	plt.axhline(y=Threshold, color='green')
	plt.axhline(y=-Threshold, color='green')
	plt.axhline(y=0, color='k')

	ticklist = list(range(math.floor(min(ExpressionDifference)), math.ceil(max(ExpressionDifference))))
	for i in ticklist: 
		if i > 0 and i == round(Threshold, 1): ticklist.remove(i)
		elif i < 0 and i == round((-1*Threshold), 1): ticklist.remove(i)
	plt.yticks(ticklist + [Threshold, (-1*Threshold)])
	plt.savefig(P.experiment_name+"/Compare_outliers/"+P.GENE_SELECTION+"_"+Agegroup+".png")
	#plt.show()
	plt.clf()
	plt.close()
	
################################################################################
#******************************************************************************#
#******************************************************************************#
#******************************* MAIN CODE ************************************#
#******************************************************************************#
#******************************************************************************#
################################################################################



def use_pcs_for_ml2():
	#________________________Setting variables____________________________________
	set_vars()
	global df_subjects
	global df_RNA_seq
	global help_name
	
	#_______________________Checking directory____________________________________
	if "Files_from_R" not in os.listdir('./'+P.experiment_name):
		exit("No PCfiles generated by R")
	
	print("help_name: ", help_name)
	
	#__________________Reading and converting the datafiles_______________________
	#Reading PCs dataframe
	df = pd.read_csv(P.experiment_name+"/Files_from_R/"+help_name+"_PCs.txt", delimiter="\t")
	df.rename({'Unnamed: 0': 'SubjectID'}, axis=1, inplace=True) #Renaming subjectID column
	df.set_index('SubjectID', inplace=True) #Setting the subjectID column as index
	
	#Reading a metafile to dictionary so that we have info[age, sex, DTHHRDY] about the subjects
	subjectinfo = {}
	with open(P.experiment_name+"/Files_for_R/"+help_name+"_METADATA.txt") as f:
		for line in f:
			(k, v) = line.split("\t")
			if k != "": subjectinfo[k] = v[:-1]

	df.insert(0, 'Agegroup', [subjectinfo[x] for x in subjectinfo.keys()]) #Adding the agegroup column to the dataframe
	df.sort_values('Agegroup', axis='rows', inplace=True, ascending=False)
	df = df.iloc[:, 0:2] #only keeping PC1 to work with

	AgeRightWrong = {}

	for subject, row in df.iterrows():
		AgeRightWrong[subject] = {}
		if df["Agegroup"][subject] == "Young (20-49)": 
			AgeRightWrong[subject]["Agegroup"] = "Young (20-49)"
			if df["PC1"][subject] < 0: AgeRightWrong[subject]["inowngroup"]=True#correctgroup = "right"
			else: AgeRightWrong[subject]["inowngroup"]=False#correctgroup = "wrong"

		if df["Agegroup"][subject] == "Old (60-79)": 
			AgeRightWrong[subject]["Agegroup"]="Old (60-79)"
			if df["PC1"][subject] > 0: AgeRightWrong[subject]["inowngroup"]=True #correctgroup = "right"
			else: AgeRightWrong[subject]["inowngroup"]=False #correctgroup = "wrong"
		#AgeRightWrong[subject] = correctgroup

	#Reading SAMPLES dictionary and only selecting samples from a specific tissue
	df_samples = strip_data_file("Source/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 0, 0, '\t')
	lst_samples = df_samples['SAMPID'].tolist() #sample IDs from dataframe to list
	dict_samples = {lst_samples[i] for i in range(0, len(lst_samples), 1)} #list to dictionary

	#_______Reading the normalized files_________________________
	df_RNA_seq = integrate_normalized_data(P.experiment_name+"/Files_from_R/"+help_name+"_NORMALIZED.txt")

	#_________________________Creating a random baseline____________________________
	if P.random_baseline: #fill the matrix with random numbers for machine learning baseline
		x = df_RNA_seq.iloc[0:, 2:].shape[0] # x = number of rows 
		y = df_RNA_seq.iloc[0:, 2:].shape[1] # y = number of columns
		df_RNA_seq.iloc[0:,2:]  = pd.DataFrame(np.random.randint(1,28,size=(x, y))) #creating random dataframe with same size
		print("Random matrix for baseline: ", df_RNA_seq)


	#________Adding new rows with subject information_______________________________
	df_RNA_seq_copy = df_RNA_seq.copy()


	#______________Finding expressiondifference between normal &outliers__________

	#1. Add row Agegroup using the dictionary
	AgeRightWrong["Name"]={"Agegroup":"Agegroup", "inowngroup":"Inowngroup"}
	df_RNA_seq_copy = df_RNA_seq_copy.drop(columns="Description")
	df_RNA_seq_copy.loc[-0.5] = [AgeRightWrong[x]["Agegroup"] for x in df_RNA_seq_copy.columns]
	
	#2.Add row "inowngroup (True/False)"
	df_RNA_seq_copy.loc[-0.2] = [AgeRightWrong[x]["inowngroup"] for x in df_RNA_seq_copy.columns]
	df_RNA_seq_copy = df_RNA_seq_copy.sort_index().reset_index(drop=True)
	df_RNA_seq_copy.set_index('Name', inplace=True) #Setting the Name column as index

	#3. Split the dataframe in 1: YoungAgeGroup and 2: OldAgeGroup
	YoungAgeGroup = df_RNA_seq_copy.loc[:, df_RNA_seq_copy.loc['Agegroup'] == 'Young (20-49)'].copy()
	OldAgeGroup = df_RNA_seq_copy.loc[:, df_RNA_seq_copy.loc['Agegroup'] == 'Old (60-79)'].copy()

	#making a copy that can later be used for machinelearning
	copyYoungAgeGroup = YoungAgeGroup.copy()
	copyOldAgeGroup = OldAgeGroup.copy()

	#4. Calculating the sums
	#pd.options.mode.chained_assignment = None  # default='warn' #TURNING OFF WANRINGS
	YoungAgeGroup['SumNormal'] = YoungAgeGroup.iloc[2:, :].loc[:, (YoungAgeGroup.loc['Inowngroup'] == True)].sum(axis=1)
	YoungAgeGroup['SumInWrongGroup'] = YoungAgeGroup.iloc[2:, :].loc[:, (YoungAgeGroup.loc['Inowngroup'] == False)].sum(axis=1)
	OldAgeGroup['SumNormal'] = OldAgeGroup.iloc[2:, :].loc[:, (OldAgeGroup.loc['Inowngroup'] == True)].sum(axis=1)
	OldAgeGroup['SumInWrongGroup'] = OldAgeGroup.iloc[2:, :].loc[:, (OldAgeGroup.loc['Inowngroup'] == False)].sum(axis=1)
	print("Number of Young samples that are in the correct agegroup: ",YoungAgeGroup.loc[:, (YoungAgeGroup.loc['Inowngroup'] == True)].shape[1])
	print("Number of Young samples that are in the wrong agegroup: ",YoungAgeGroup.loc[:, (YoungAgeGroup.loc['Inowngroup'] == False)].shape[1])
	print("Number of Old samples that are in the correct agegroup: ",OldAgeGroup.loc[:, (OldAgeGroup.loc['Inowngroup'] == True)].shape[1])
	print("Number of Old samples that are in the wrong agegroup: ",OldAgeGroup.loc[:, (OldAgeGroup.loc['Inowngroup'] == False)].shape[1])
	
	#5. Calculating the mean of the sums
	YoungAgeGroup['MeanNormal'] = YoungAgeGroup["SumNormal"].iloc[2:]/YoungAgeGroup.loc[:, (YoungAgeGroup.loc['Inowngroup'] == True)].shape[1]
	YoungAgeGroup['MeanInWrongGroup'] = YoungAgeGroup["SumInWrongGroup"].iloc[2:]/YoungAgeGroup.loc[:, (YoungAgeGroup.loc['Inowngroup'] == False)].shape[1]
	OldAgeGroup['MeanNormal'] = OldAgeGroup["SumNormal"].iloc[2:]/OldAgeGroup.loc[:, (OldAgeGroup.loc['Inowngroup'] == True)].shape[1]
	OldAgeGroup['MeanInWrongGroup'] = OldAgeGroup["SumInWrongGroup"].iloc[2:]/OldAgeGroup.loc[:, (OldAgeGroup.loc['Inowngroup'] == False)].shape[1]

	#6. Remove the other expression columns
	YoungAgeGroup = YoungAgeGroup.loc[:, ['MeanNormal', 'MeanInWrongGroup']]
	OldAgeGroup = OldAgeGroup.loc[:, ['MeanNormal', 'MeanInWrongGroup']]

	#7. For each dataframe add the column "difference"
	YoungAgeGroup["Difference"] = YoungAgeGroup.loc[:, 'MeanNormal'] - YoungAgeGroup.loc[:, 'MeanInWrongGroup']
	OldAgeGroup["Difference"] = OldAgeGroup.loc[:, 'MeanNormal'] - OldAgeGroup.loc[:, 'MeanInWrongGroup']

	#8. Add the gene description again
	YoungAgeGroup.insert(0, 'Description', ["-", "-"] + list(df_RNA_seq['Description'])) #Adding the agegroup column to the dataframe
	OldAgeGroup.insert(0, 'Description', ["-", "-"] + list(df_RNA_seq['Description'])) #Adding the agegroup column to the dataframe

	#9. Sort on difference
	YoungAgeGroup.iloc[2:, :] = YoungAgeGroup.iloc[2:, :].sort_values('Difference', ascending=False, key=abs)
	OldAgeGroup.iloc[2:, :] = OldAgeGroup.iloc[2:, :].sort_values('Difference', ascending=False, key=abs)

	#10. uUe the threshold line to make a list of genes that are down or upregulated in the two groups
	#print(YoungAgeGroup.iloc[2:(P.NrFoundRelevantGenes+2), [0,3]]) #printing the relevant genes
	YoungThreshold = YoungAgeGroup.iloc[(P.NrFoundRelevantGenes+1), 3]
	OldThreshold = OldAgeGroup.iloc[(P.NrFoundRelevantGenes+1), 3]

	# 11. Writing most important genes to filename
	if "Compare_outliers" not in os.listdir(P.experiment_name): #Making a folder for the machinelearning results
		os.mkdir(os.path.join(os.getcwd(), P.experiment_name+"/Compare_outliers"))
	if "Outlier_determining_genes" not in os.listdir(P.experiment_name+"/Compare_outliers"): #Making a folder for the machinelearning results
		os.mkdir(os.path.join(os.getcwd(), P.experiment_name+"/Compare_outliers/Outlier_determining_genes"))
	YoungAgeGroup.iloc[2:(P.NrFoundRelevantGenes+2), [0,3]].to_csv(P.experiment_name+"/Compare_outliers/Outlier_determining_genes/Young_"+P.GENE_SELECTION+"_Outlier_determining_genes.txt", index=False, header=False, sep="\t")
	OldAgeGroup.iloc[2:(P.NrFoundRelevantGenes+2), [0,3]].to_csv(P.experiment_name+"/Compare_outliers/Outlier_determining_genes/Old_"+P.GENE_SELECTION+"_Outlier_determining_genes.txt", index=False, header=False, sep="\t")

	#12. Creating a barplot with threshold line
	# Barplot for YoungAgeGroup	
	Genes = list(YoungAgeGroup.iloc[2:, 0]) #Selecting the ' Description' column with the gene names
	ExpressionDifference = list(YoungAgeGroup.iloc[2:, 3])
	barplot_expressiondifference(Genes, ExpressionDifference, "Young", YoungThreshold)
	# Barplot for Old
	Genes = list(OldAgeGroup.iloc[2:, 0]) #Selecting the ' Description' column with the gene names
	ExpressionDifference = list(OldAgeGroup.iloc[2:, 3])
	barplot_expressiondifference(Genes, ExpressionDifference, "Old", OldThreshold)



	##############################################################################
	######################## MACHINE LEARNING ####################################
	##############################################################################
	
	#_____________ Make folder for machinelearning results _______________________
	if "Outlier_ML_Tables" not in os.listdir(P.experiment_name+"/Compare_outliers"): #Making a folder for the machinelearning results
		os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Compare_outliers/Outlier_ML_Tables"))
	
	#______________Adding the description column__________________________________
	copyYoungAgeGroup['Description'] = ["-", "-"]+df_RNA_seq['Description'].tolist()
	copyOldAgeGroup['Description'] = ["-", "-"]+df_RNA_seq['Description'].tolist()
	for ThisAgeGroup in ["Young", "Old"]:
		if ThisAgeGroup == "Young": dfAgeGroup = copyYoungAgeGroup
		elif ThisAgeGroup == "Old": dfAgeGroup = copyOldAgeGroup
		#print(dfAgeGroup.name)

		print(dfAgeGroup)
		################## Transformationa and setting indexes #########################
		dfAgeGroup = dfAgeGroup.set_index('Description') #setting the "Desription" column (genenames) as index
		dfAgeGroup = dfAgeGroup.transpose() #transposing the dataframe so that the test-subjects can be used as instances
		dfAgeGroup = dfAgeGroup.iloc[0:, 1:]# .drop('AgeGroup')
		dfAgeGroup.rename(columns = {'-' : 'Inowngroup'}, inplace = True)

		dfAgeGroup = dfAgeGroup.fillna(0) #missing values become zero

		#dfAgeGroup.iloc[0:,1:] = log_adjustment(dfAgeGroup.iloc[0:,1:])

		#_____________Making file to write to_________________________________________
		if P.MODEL not in os.listdir(P.experiment_name+"/Compare_outliers/Outlier_ML_Tables"): #Making a folder for the machinelearning results
			os.mkdir(os.path.join(os.getcwd(), P.experiment_name+"/Compare_outliers/Outlier_ML_Tables/"+P.MODEL))

		################# Making the model #############################################

		mean_f1 = 0
		######## pREPARING THe train and test data ####################################
		print(st.CYAN)
		print("Method = ", P.METHOD, "\nModel = ", P.MODEL, "\nDataset = ", P.GENE_SELECTION, "\nGroup = ", ThisAgeGroup, st.RST)

		y = dfAgeGroup['Inowngroup'].values.copy() #Making y (prediction column)
		ycopy = [None]*len(y)
		if P.METHOD == "Classification":
			for i in range(0, len(y)): #converting strings to groups
				if y[i]: ycopy[i] = "InCorrectGroup"
				elif not y[i]: ycopy[i] = "InWrongGroup"
		"""elif P.METHOD == "Regression":
			for i in range(len(y)): #converting strings to integers (groups)
				if   y[i] == '20-29' or y[i] == '30-39' or y[i] == '40-49': y[i] = 1
				elif y[i] == '50-59'																			: y[i] = 2
				elif y[i] == '60-69' or y[i] == '70-79'										: y[i] = 3
		"""
		print(st.YELLOW, st.BOLD,"Number of examples", len(y), "(nr of people)", st.RST)
		y = ycopy

		X = dfAgeGroup.drop(['Inowngroup'],axis=1).values #Making X, on which the prediction has to be made
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.7, shuffle=False) #Splitting into train and test data
		
		print(len(y_test))
		print(len(X_test))
		print(len(X_train))
		print(len(y_train))

		#-------------------------------------------------------------------------------
		print("<-> Starting modeling process for {}...".format(P.MODEL))

		if P.METHOD == "Regression":
			if P.MODEL == "Support Vector Machine":	clf = svm.SVR(kernel="linear")
			elif P.MODEL == "RandomForest":					clf = RandomForestRegressor(n_estimators=15, max_depth=50, random_state=None) #, criterion='MSE', splitter='best')
																						#clf = BaggingRegressor(rf, n_estimators=45, max_samples=0.1, random_state=25)
			elif P.MODEL == "DecisionTree":					clf = DecisionTreeRegressor(random_state=None, max_depth=15) #, criterion='squared_error')
		elif P.METHOD == "Classification":
			if P.MODEL == "Support Vector Machine":	clf = svm.SVC(kernel='linear')
			elif P.MODEL == "RandomForest":					clf = RandomForestClassifier(n_estimators=15, max_depth=50, random_state=None) #, criterion='MSE', splitter='best')
			elif P.MODEL == "DecisionTree":					clf = DecisionTreeClassifier(criterion='entropy', max_depth=15, random_state=None)

		clf.fit(X_train, y_train)
		y_pred = clf.predict(X_test) #make the prediction on the test data

		print("Modeling done.")
		#-------------------------------------------------------------------------------
		################### EVALUATION OF THE MACHINELEARNING ##########################
		#-------------------------------------------------------------------------------
		#Evalueation of the preduction quality
		if P.METHOD =="Regression":
			for i in range(len(y_pred)):
				y_pred[i] = round(y_pred[i], 0)

		#---------------------------------------------------------------------------
		#Writing LaTeX table
		if P.WriteLaTeXTableforML:
			if P.random_baseline: help_name= help_name.replace(help_name.split("_")[0], "RandomBaseline")

			f = open(P.experiment_name+"/Compare_outliers/Outlier_ML_Tables/"+P.MODEL+"/"+ThisAgeGroup+"_"+P.GENE_SELECTION+"_"+P.MODEL, "w")
			f.write(help_name+":") #write the first line
			f.write("\n\n"+"\subsection*{}\n".format("{"+P.METHOD+"}"))
			f.write("\subsubsection*{}\n".format("{"+P.MODEL+"}"))

			write_latex_line(["Dataset","AgeGroups", "Accuracy", "Precision", "Recall", "F1", "Occ.Pred", "Occ.real", "Correct"], True, f)
			
			#ages = ["Young", "Middle", "Old"]
			if P.METHOD == "Classification": InOwnGroup = ["InCorrectGroup", "InWrongGroup"]
			if P.METHOD == "Regression": ages = [1,2,3]

			y_pred = list(y_pred)

			#_____________________MachineLearning Analysis____________________________
			Accuracy = sklearn.metrics.accuracy_score(y_test, y_pred, normalize=True, sample_weight=None)
			print("\033[33m\033[1m<-> Accuracy:", Accuracy, ' \033[0m')

			mean_f1 = 0
			for group in ["InCorrectGroup","InWrongGroup"]:
				Precision = sklearn.metrics.precision_score(y_test, y_pred, labels=None, pos_label=group, average='binary', sample_weight=None, zero_division='warn')
				Recall = sklearn.metrics.recall_score(y_test, y_pred, labels=None, pos_label=group, average='binary', sample_weight=None, zero_division='warn')
				F1 = sklearn.metrics.f1_score(y_test, y_pred, pos_label=group)
				Occurence = y_test.count(group)
				Predicted = y_pred.count(group)
				CorrectPredicted = sum([1 for i,j in zip(y_test,y_pred) if (i==j and j==group)])
				mean_f1 += F1

				print("\033[1m\033[34m<-> {}:".format(group))
				print("\033[33m\033[1m  -> Precision:\t",Precision,'\033[0m')
				print("\033[33m\033[1m  -> Recall:\t",Recall,'\033[0m')
				print("\033[33m\033[1m  -> F1:\t",F1,'\033[0m\n')
				print("\033[33m\033[1m  -> Database occurence:\t",Occurence,'\033[0m')
				print("\033[33m\033[1m  -> Occurence predicted:\t",Predicted,'\033[0m')
				print("\033[33m\033[1m  -> Correct predicted:\t",CorrectPredicted,'\033[0m\n')

				GenelistnameIsLong = False
				if len(help_name.split("_")[0].split("-")) > 2:
					GenelistnameIsLong = True
				word = " "
				if group in ["InCorrectGroup", 1]: #If "InCorrectGroup"
					word = help_name.split("_")[0].capitalize()
					if GenelistnameIsLong: word = "-".join(word.split("-")[0:2])
					write_latex_line([help_name.split("_")[0].capitalize() , group, Accuracy, Precision, Recall, F1, Predicted, Occurence, CorrectPredicted], False, f)
				else: #Group is InWrongGroup
					word = help_name.split("_")[1] 
					write_latex_line([ThisAgeGroup, group, " ", Precision, Recall, F1, Predicted, Occurence, CorrectPredicted], False, f)
			#####################################
			"""
			#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
				if group_predicted == 0: Precision = "not pred."
				else: Precision = correct/group_predicted #Calculating Precision
				if group_real == 0: Recall = "not occur."
				else: Recall = correct/group_real #Calculating Recall
				if (Precision != 0 or Recall != 0) and isinstance(Precision,(int,float)) and isinstance(Recall,(int,float)):
					F1 = (2*Precision*Recall)/(Precision+Recall) #F1 = sklearn.metrics.f1_score(y_test, y_pred, pos_label=group)
					group_mean_f1 += F1
				else:
					F1 = "Not poss."
					group_mean_f1 += 0
			"""

			f.write("\t\t\hline\n")
			f.write("\t\end{}".format("{tabular}\n"))
			f.write("\t\caption{}".format("{Evaluation of "+P.METHOD+" by "+P.MODEL+" using the "+help_name.replace("_","-")+" dataset.}\n"))
			f.write("\t\label{}".format("{tab:"+P.METHOD+P.MODEL+help_name+"}\n"))
			f.write("\end{}".format("{table}"))


			print(st.CYAN)
			print("Average F1 (over both groups (outlier or not)):	", mean_f1/2, st.RST )
			f.close()
			#-------------------------------------------------------------------------------



		#exit("JAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAPEEEEEEEEEEEEE")
		################################################################################
		########### EXTRACTING LIST OF MOST IMPORTANT GENES ############################
		
		if P.MODEL != "Support Vector Machine": importances = clf.feature_importances_ #creating a list with importance values for all the genes
		else: importances = clf.coef_[0]
		
		####### PLOTTING THE GENE IMPORTANCE FOR RANDOM FORREST #####################---
		if P.MODEL == "RandomForest":
			std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0).tolist() #calulating the standard deviation
			std.sort(reverse=True) #sorting the standard deviation
			std2 = std[:P.NrFoundRelevantGenes] #The X most imporatant genes

			forest_importances = pd.Series(importances, index=dfAgeGroup.iloc[0:,1:].columns) #linking importances to genes
			forest_importances = forest_importances.sort_values(ascending=False) #df.sort_values(by=['col1'])
			forest_importances2 = forest_importances[:P.NrFoundRelevantGenes] #The 50 most important genes

			#Making a plot that shows the gene importance of the many decision trees together with random forest with standarddeviation
			fig, ax = plt.subplots(figsize=(14,8))
			forest_importances2.plot.bar(yerr=std2, ax=ax)
			ax.set_title("Feature importances using MDI")
			ax.set_ylabel("Mean decrease in impurity")
			fig.tight_layout()
			plt.rcParams["figure.figsize"] = (10,6)
			
			plt.savefig(P.experiment_name+'/Compare_outliers/Outlier_Important_genes/'+P.GENE_SELECTION+'_'+ThisAgeGroup+'_RF_FeatureImportances.pdf')
			plt.close()
			#plt.show()


	###########################################################################-----

	#_____________________Selecting only the X best features______________________
	feature_dict = {}
	for feat, importance in zip(dfAgeGroup.iloc[0:,1:].columns, importances): #linking the genes and their importances together
		feature_dict[feat] = importance #filling a dictionary with genes as keys and their importance as values

	feature_dict = dict(sorted(feature_dict.items(), key=lambda item: item[1], reverse=True)) #Sorting the genes from important to not important
	best_features = list(islice(feature_dict.items(), P.NrFoundRelevantGenes)) #take the best Nr genes

	#____________________Writing the genelist_____________________________________
	if "Outlier_Important_genes" not in os.listdir('./'+P.experiment_name+"/Compare_outliers"):
		os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Compare_outliers/Outlier_Important_genes"))
	if P.GENE_SELECTION not in os.listdir('./'+P.experiment_name+"/Compare_outliers/Outlier_Important_genes"):
		os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Compare_outliers/Outlier_Important_genes/"+P.GENE_SELECTION))

	f = open("./"+P.experiment_name+"/Compare_outliers/Outlier_Important_genes/"+P.GENE_SELECTION+"/ImportantGenes_"+ThisAgeGroup"["+P.MODEL+"]["+P.GENE_SELECTION+"].txt", "w+")
	for i in best_features:
		f.write(i[0]+" "+str(i[1])+"\n")
	f.close

	#_____________________Plotting the Feature importance_________________________
	f_importances([x[1] for x in best_features], [x[0] for x in best_features]) #function for making a graph if gene importance, first argument is the importances, second arg is the featurenames
	
	#best_feat_list = [x[1] for x in best_features]
	#best_gene_list = [x[0] for x in best_features] 
	#print(best_gene_list)
	
	#print(best_feat_list)
	#exit()
	#best_feat_list = list(best_features.values())

	"""
	plt.bar([x for x in best_gene_list], best_feat_list, width=1.0)
	#plt.bar([x for x in range(len(feature_dict))], feature_dict)
	plt.show()
	exit("NOUNOU")
	"""

	"""
	######################## PLOTTING THE DECISION TREE ############################
	if P.MODEL == "DecisionTree":
		fn = list(df_RNA_seq)
		del fn[0]
		cn = ["Young", "Middle", "Old"]

		fig, axes = plt.subplots(nrows = 1,ncols = 1, dpi=500, figsize=(30,10))#,figsize = (30,28))
		tree.plot_tree(clf,
									feature_names = fn, 
									class_names=cn,
									#filled = True,
									fontsize=3
									)#;

		fig.savefig(P.GENE_SELECTION+' DecisionTree.png')
	#plt.figure(dpi=800) #,figsize=(10,10))  # set plot size (denoted in inches)
	#tree.plot_tree(clf, fontsize=2, filled=True, feature_names = fn)
	#plt.show()

	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	"""

	print('\a')
	################ LINKING MOST IMPORANT GENES TO PATHWAYS #######################
	#isabelle = {}
	#for gek in best_features:
	#	for i in temppathways: #loop over all pathwyays
	#		for tomke in temppathways[i]: #loop over all genes in a pathway
	#			if tomke == gek[0]:
	#				if i in isabelle: #hoog count op
	#					isabelle[i] += 1
	#				else:
	#					isabelle[i] = 1
	#
	#
	#for i in isabelle:
	#	isabelle[i] /= len(temppathways[i])
	#	isabelle[i] = round(isabelle[i], 3)
	#isabelle = dict(sorted(isabelle.items(), key=lambda item: item[1], reverse=True))
	#isabelle2 = list(islice(isabelle, 30))
	#print("&&&&&&&&&&&&&&&&&&&&&&&")
	#for k in isabelle2:
	#	print(isabelle[k],'\t', k)
	#print("&&&&&&&&&&&&&&&&&&&&&&&")

	#exit(0)






