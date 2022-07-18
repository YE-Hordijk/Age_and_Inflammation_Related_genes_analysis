#Machinelearning.py


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

#*******************************************************************************
def set_vars():
	global help_name
	global df_subjects
	global df_RNA_seq
	global help_name
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
	plt.savefig(P.experiment_name+'/Important_genes/'+P.GENE_SELECTION+"["+P.MODEL+"]"+'_FeatureImportances.pdf')
	plt.close()
	#plt.show()

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
	#print(w)
	#input("verder?")
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
#*******************************************************************************
#******************************* MAIN CODE *************************************
#*******************************************************************************
#*******************************************************************************

def machinelearning():
	set_vars()
	global df_subjects
	global df_RNA_seq
	global help_name
	
	#Reading SAMPLES dictionary and only selecting samples from a specific tissue
	df_samples = strip_data_file("Source/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 0, 0, '\t')
	lst_samples = df_samples['SAMPID'].tolist() #sample IDs from dataframe to list
	dict_samples = {lst_samples[i] for i in range(0, len(lst_samples), 1)} #list to dictionary

	#Readig the subject information 
	df_subjects = GF.readfile("Source/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", '\t', 0, "dataframe")
	df_subjects.columns = ['SUBJID','SEX','AGE','DTHHRDY'] #MISSCHIEN WEGHALEN, LIJKT OVERBODIG
	df_subjects = pd.DataFrame(df_subjects) #MISSCHIEN WEGHALEN, LIJKT OVERBODIG
	df_subjects = df_subjects.set_index('SUBJID') #Setting subjectID as index-column
	df_subjects.sort_values(by=['AGE','SEX'], inplace=True) #sort first by Age and then by Sex


	#_______Reading the normalized files_________________________
	df_RNA_seq = integrate_normalized_data(P.experiment_name+"/Files_from_R/"+help_name+"_NORMALIZED.txt")


	#_________________________Creating a random baseline____________________________
	if P.random_baseline: #fill the matrix with random numbers for machine learning baseline
		x = df_RNA_seq.iloc[0:, 2:].shape[0] # x = number of rows 
		y = df_RNA_seq.iloc[0:, 2:].shape[1] # y = number of columns
		df_RNA_seq.iloc[0:,2:]  = pd.DataFrame(np.random.randint(1,28,size=(x, y))) #creating random dataframe with same size
		print("Random matrix for baseline: ", df_RNA_seq)


	#________Adding new rows with subject information_______________________________
	print(df_RNA_seq)
	#Adding new rows with subject information
	new_row = make_subject_row('sex')
	df_RNA_seq = add_new_row(new_row)
	new_row = make_subject_row('DTHHRDY')
	df_RNA_seq = add_new_row(new_row)
	new_row = make_subject_row('age')
	df_RNA_seq = add_new_row(new_row)

	df_RNA_seq = df_RNA_seq.set_index('Name')
	print(df_RNA_seq)








	#############This code "cleans up" the data (optional ##########################
	if P.removing_outliers:
		young_middlepoint = []
		old_middlepoint = []
		young_num = 0
		old_num = 0


		#TODO 2) Calculating euclidian middle for "young" and "old"
		for i in df_RNA_seq.iloc[0:, 1:]: #looping over dataframe without the column "Descriptions"
			expr_List = df_RNA_seq[i].values[3:].tolist() #Making a list with of expressionlevels of all the genes for one sample
			if df_RNA_seq[i]['Subject age'] < "50-59": #this subject is young
				young_num += 1 #counting the number of young samples
				if (len(young_middlepoint) == 0): young_middlepoint = expr_List #first time
				else: young_middlepoint = np.add(young_middlepoint, expr_List) #add expression levels to the "young-list"
			elif df_RNA_seq[i]['Subject age'] > "50-59": #this subject is old
				old_num += 1 #counting the number of old samples
				if (len(old_middlepoint) == 0): old_middlepoint = expr_List #first time
				else: old_middlepoint = np.add(old_middlepoint, expr_List) #add expression levels to the "old-list"

		young_middlepoint = [counts / young_num for counts in young_middlepoint] #deviding all expressionlevels by number of counts creating the mean for every gene
		old_middlepoint = [counts / old_num for counts in old_middlepoint] ##deviding all expressionlevels by number of counts creating the mean for every gene

		#print(young_middlepoint[:7], "...")
		#print(old_middlepoint[:7], "...")
		#for i in range(len(young_middlepoint)):
		#	if ((young_middlepoint[i] - old_middlepoint[i]) > 1):
		#		print(young_middlepoint[i], " ", old_middlepoint[i])

		# 3) looping over instances and compare to which middle they are closer, if closer to the wrong middle remove them
		dist_to_young = 0
		dist_to_old = 0
		temp_new_data = {}
		temp_new_data['Description'] = df_RNA_seq['Description']
		outliers = {}
		outliers["Description"] = df_RNA_seq["Description"]


		for i in df_RNA_seq.iloc[0:, 1:]: #looping over dataframe without the column "Descriptions"
			expr_List = df_RNA_seq[i].values[3:].tolist() #Making a list with of expressionlevels of all the genes for one sample
			dist_to_young = np.linalg.norm(np.array(expr_List)-np.array(young_middlepoint)) # calculate dist_to_young
			dist_to_old = np.linalg.norm(np.array(expr_List)-np.array(old_middlepoint)) #calculate dist_to_old
			if (df_RNA_seq[i]['Subject age'] < "50-59") and not (dist_to_young > dist_to_old): #this subject is young and does not lie closer to old-middle
				temp_new_data[i] = df_RNA_seq[i] 		#i is not een outlier, so copy to temp_new_data
			elif df_RNA_seq[i]['Subject age'] > "50-59" and not (dist_to_young < dist_to_old): #this subject is old and does not lie closer to young-middle
				temp_new_data[i] = df_RNA_seq[i] #i is not een outlier, so copy to temp_new_data
			else: #dit is wel een outlier
				outliers[i] = df_RNA_seq[i] #wordt opgeslagen in een outlier bestand
				
		df_RNA_seq = pd.DataFrame.from_dict(temp_new_data) #, orient='index')
		df_outliers = pd.DataFrame.from_dict(outliers) #, orient='index')

		#Result
		print(P.GENE_SELECTION)
		print("clean data\n", df_RNA_seq)
		print("outliers\n", df_outliers)


	################## Transformationa and setting indexes #########################
	df_RNA_seq = df_RNA_seq.set_index('Description') #setting the "Desription" column (genenames) as index
	df_RNA_seq = df_RNA_seq.transpose() #transposing the dataframe so that the test-subjects can be used as instances


	for i in ['DTHHRDY','sex']:
		df_RNA_seq = df_RNA_seq.drop([i], axis=1) #remove the row "DTHHRDY" and "sex"
	df_RNA_seq = df_RNA_seq.fillna(0) #missing values become zero

	#print(df_RNA_seq.iloc[0:,1:])
	df_RNA_seq.iloc[0:,1:] = log_adjustment(df_RNA_seq.iloc[0:,1:])
	#print(df_RNA_seq.iloc[0:,1:])



	############## Make folder for machinelearning results #######################
	if "Machine_Learning_Results" not in os.listdir(P.experiment_name): #Making a folder for the machinelearning results
		os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Machine_Learning_Results"))

	if P.MODEL not in os.listdir(P.experiment_name+"/Machine_Learning_Results"): #Making a folder for the machinelearning results
		os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Machine_Learning_Results/"+P.MODEL))

	################# Making the model #############################################

	mean_f1 = 0
	######## pREPARING THe train and test data ####################################
	print(st.CYAN)
	print("Method = ", P.METHOD, "\nModel = ", P.MODEL, "\nDataset = ", P.GENE_SELECTION, st.RST)

	y = df_RNA_seq['age'].values.copy() #Making y (prediction column)

	if P.METHOD == "Classification":
		for i in range(len(y)): #converting strings to groups
			if   y[i] == '20-29' or y[i] == '30-39' or y[i] == '40-49': y[i] = "Young"
			elif y[i] == '50-59'																			: y[i] = "Middle"
			elif y[i] == '60-69' or y[i] == '70-79'										: y[i] = "Old"
	elif P.METHOD == "Regression":
		for i in range(len(y)): #converting strings to integers (groups)
			if   y[i] == '20-29' or y[i] == '30-39' or y[i] == '40-49': y[i] = 1
			elif y[i] == '50-59'																			: y[i] = 2
			elif y[i] == '60-69' or y[i] == '70-79'										: y[i] = 3

	print(st.YELLOW, st.BOLD,"Number of examples", len(y), "(nr of people)", st.RST)

	X = df_RNA_seq.drop(['age'],axis=1).values #Making X, on which the prediction has to be made
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, shuffle=False) #Splitting into train and test data

	#print(P.GENE_SELECTION)
	#print(P.MODEL) 	
	#print(len(X), " x ", len(X[0]))
	#print(len(y_test))
	#input("volgende dataset?")

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
	
	#print(y_train)

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

	#RMSE = math.sqrt(np.square(np.subtract(y_test,y_pred)).mean())
	#print("\n\033[33m\033[1m<-> P.Model = ", P.MODEL, "\n<-> RMSE:", round(RMSE, 4), '\a \033[0m')

	#Count the number of right and wrong guesses
	right = 0
	wrong = 0
	for i in range(len(y_test)):
		#print(y_test[i], "--", y_pred[i])
		if y_test[i] == y_pred[i]: right += 1
		else: wrong += 1
	Accuracy = round(((right/(right+wrong))*100), 1)
	print("\033[33m\033[1m<-> Accuracy:", Accuracy, ' \033[0m')

	#---------------------------------------------------------------------------

	#Writing LaTeX table
	if P.WriteLaTeXTableforML:
		if P.random_baseline: help_name= help_name.replace(help_name.split("_")[0], "RandomBaseline")
		
		
		#input(help_name)
		f = open(P.experiment_name+"/Machine_Learning_Results/"+P.MODEL+"/"+help_name+"_"+P.MODEL+"_ML_Results.txt", "w")
		f.write(help_name+":") #write the first line
		f.write("\n\n"+"\subsection*{}\n".format("{"+P.METHOD+"}"))
		f.write("\subsubsection*{}\n".format("{"+P.MODEL+"}"))

		write_latex_line(["Dataset","AgeGroups", "Accuracy", "Precision", "Recall", "F1", "Occ.Pred", "Occ.real", "Correct"], True, f)
		
		#ages = ["Young", "Middle", "Old"]
		if P.METHOD == "Classification": ages = ["Young", "Middle", "Old"]
		if P.METHOD == "Regression": ages = [1,2,3]

		group_mean_f1 = 0
		for group in ages:
			correct = 0 #Calculating Number of people correctly classified as <group>
			group_predicted = 0 #Total number of people in agegroup <group> that are predicted
			group_real = 0 #Total number of people in agegroup <group> that are really in the data
			for i in range(len(y_test)):
				correct += (y_test[i] == group and y_pred[i] == group) #counting the correct prediciotns of this specfic agegroup
				group_predicted += (y_pred[i] == group) #counting nr of times this agegroup was predicted (correct or false)
				group_real += (y_test[i] == group) #coutning nr of times this agegroup should have been predicted


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
			
			print("\033[1m\033[34m<-> {}:".format(group))
			print("\033[33m\033[1m  -> Precision:\t",Precision,'\033[0m')
			print("\033[33m\033[1m  -> Recall:\t",Recall,'\033[0m')
			print("\033[33m\033[1m  -> F1:\t",F1,'\033[0m\n')
			print("\033[33m\033[1m  -> Database occurance:\t",group_real,'\033[0m')
			print("\033[33m\033[1m  -> Occurance predicted:\t",group_predicted,'\033[0m')
			print("\033[33m\033[1m  -> Correct predicted:\t",correct,'\033[0m\n')
			
			GenelistnameIsLong = False
			if len(help_name.split("_")[0].split("-")) > 2:
				GenelistnameIsLong = True
			word = " "
			if group in ["Young", 1]: #If young group
				word = help_name.split("_")[0].capitalize()
				if GenelistnameIsLong: word = "-".join(word.split("-")[0:2])
				if P.use_middle_age: write_latex_line([word , group, " ", Precision, Recall, F1, group_predicted, group_real, correct], False, f)
				else: write_latex_line([help_name.split("_")[0], group, Accuracy, Precision, Recall, F1, group_predicted, group_real, correct], False, f)
			elif group in ["Middle", 2]: #If middle group
				if (GenelistnameIsLong): word = "-".join(help_name.split("_")[0].split("-")[2:])
				if P.use_middle_age: write_latex_line([word, group, Accuracy, Precision, Recall, F1, group_predicted, group_real, correct], False, f)
				#else: write_latex_line([word, " ", Accuracy, " ", " ", " ", " ", " ", " "], False)
			else: #Group is old
				word = help_name.split("_")[1] 
				write_latex_line([word, group, " ", Precision, Recall, F1, group_predicted, group_real, correct], False, f)
				

		f.write("\t\t\hline\n")
		f.write("\t\end{}".format("{tabular}\n"))
		f.write("\t\caption{}".format("{Evaluation of "+P.METHOD+" by "+P.MODEL+" using the "+help_name.replace("_","-")+" dataset}\n"))
		f.write("\t\label{}".format("{tab:"+P.METHOD+P.MODEL+help_name+"}\n"))
		f.write("\end{}".format("{table}"))

		if P.use_middle_age: numgroups = 3
		else: numgroups = 2
		group_mean_f1 /= numgroups
		mean_f1 += group_mean_f1
		#print("group_mean_f1:	",group_mean_f1)
		#print("mean_f1:	", mean_f1)

		print(st.CYAN)
		print("Average F1 (over all ",numgroups," age groups):	", mean_f1, st.RST )
		f.close()
		#-------------------------------------------------------------------------------

		


	################################################################################
	########### EXTRACTING LIST OF MOST IMPORTANT GENES ############################

	if P.MODEL != "Support Vector Machine": importances = clf.feature_importances_ #creating a list with importance values for all the genes
	else: importances = clf.coef_[0]
	####### PLOTTING THE GENE IMPORTANCE FOR RANDOM FORREST #####################---
	if P.MODEL == "RandomForest":
		std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0).tolist()
		std.sort(reverse=True)
		std2 = std[:P.NrFoundRelevantGenes] #The N most imporatant genes
		forest_importances = pd.Series(importances, index=df_RNA_seq.iloc[0:,1:].columns)
		forest_importances = forest_importances.sort_values(ascending=False) #df.sort_values(by=['col1'])
		forest_importances2 = forest_importances[:P.NrFoundRelevantGenes] #The N most important genes

		#Making a plot that shows the gene importance of the many decision trees together with random forest
		fig, ax = plt.subplots(figsize=(14,8))
		forest_importances2.plot.bar(yerr=std2, ax=ax)
		ax.set_title("Feature importances using MDI")
		ax.set_ylabel("Importance with standarddeviation")
		#XLABEL = str(P.NrFoundRelevantGenes)+"most important genes"
		ax.set_xlabel(str(P.NrFoundRelevantGenes)+" most important genes")
		fig.tight_layout()
		plt.axhline(y=0)
		plt.rcParams["figure.figsize"] = (10,6)
		plt.savefig(P.experiment_name+'/'+P.GENE_SELECTION+'.pdf')
		plt.close()
		#plt.show()


	###########################################################################-----


	#_____________________Selecting only the X best features______________________
	feature_dict = {}
	for feat, importance in zip(df_RNA_seq.iloc[0:,1:].columns, importances): #linking the genes and their importances together
		feature_dict[feat] = importance #filling a dictionary with genes as keys and their importance as values
	
	feature_dict = dict(sorted(feature_dict.items(), key=lambda item: item[1], reverse=True)) #Sorting the genes from important to not important
	best_features = list(islice(feature_dict.items(), P.NrFoundRelevantGenes)) #take the best 80 genes
	
	#____________________Writing the genelist_____________________________________
	if "Important_genes" not in os.listdir('./'+P.experiment_name):
		os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Important_genes"))
	if P.GENE_SELECTION not in os.listdir('./'+P.experiment_name+"/Important_genes"):
		os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Important_genes/"+P.GENE_SELECTION))
		
	f = open("./"+P.experiment_name+"/Important_genes/"+P.GENE_SELECTION+"/ImportantGenes["+P.MODEL+"]["+P.GENE_SELECTION+"].txt", "w+")
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





