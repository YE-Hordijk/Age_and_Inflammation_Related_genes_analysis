"""import dataframe_image as dfi
df = pd.DataFrame(np.random.randn(6, 6), columns=list('ABCDEF'))
df_styled = df.style.background_gradient()
print(df_styled)
dfi.export(df_styled, "mytable.png")"""

#Making_Gradient_Table.py

from Parameters import P
import os
import matplotlib.pyplot as plt
import numpy as np
import re
import csv
import pandas as pd
from matplotlib_venn import venn2, venn3, venn3_circles
import math

#For saving dataframe as png
import matplotlib.pyplot as plt
#from pandas.table.plotting import table # EDIT: see deprecation warnings below
from pandas.plotting import table 


import dataframe_image as dfi

import matplotlib
import seaborn as sns


#*******************************************************************************
def save_csv(d, filename):
	#df = pd.DataFrame.from_dict(d)
	df = pd.DataFrame(dict([ (k,pd.Series(v, dtype=pd.StringDtype())) for k,v in d.items() ]))
	df.to_csv(filename, index= False)
#*******************************************************************************

if "GenelistGradientTables" not in os.listdir(P.experiment_name): #Making a folder for the machinelearning results
	os.mkdir(os.path.join(os.getcwd(), P.experiment_name+"/GenelistGradientTables"))


#____________________Making the path____________________________________________

METHOD = "Classification"
#GENEDATASET = "senescence"
use_only_overlap = False

exp_type2 = ""
exp_type = input("TotalGeneDataset[1] or Compare_outliers[2]?")
if exp_type == "1": 
	exp_type1 = "Important_genes"
	foldername = "TotalGeneset_MLresults"
elif exp_type == "2": 
	exp_type2 = input("Outlier_Important_genes[1] or Outlier_determining_genes[2]?")
	if exp_type2 == "1": 
		exp_type1 = "Compare_outliers/Outlier_Important_genes"
		foldername = "CompareOutliers_MLresults"
		#xx = 
		if input("Use only overlapping genes for young and old? Yes[1] / No[2]") == '1': use_only_overlap = True
	elif exp_type2 == "2": 
		exp_type1 = "Compare_outliers/Outlier_determining_genes"
		foldername = "CompareOutliers_Outlierdeteremining"
		#exit("THIS OPTION IS NOT POSSIBLE YET AND MAYBE SHOULD BE BECAUSE WHAT ARE YOUGOING TO COMPARE IT TO?")
	else: exit("Invalid choice")
else: exit("Invalid choice")

#if exp_type2 != "2": #is machinelearning result genelist
lijst = ["senescence", "searchwords", "genes-from-papers","cell-age-signatures", "all"]





#if foldername not in os.listdir(P.experiment_name+"/VennDiagrams"): #Making a folder for the machinelearning results
#	os.mkdir(os.path.join(os.getcwd(), P.experiment_name+"/VennDiagrams/"+foldername))

DT_Allsets = {}
RF_Allsets = {}
SVM_Allsets = {}

#____________________Reading the lists_________________________________________
for GENEDATASET in lijst: 
	if exp_type2=="2": path = P.experiment_name+"/"+exp_type1+"/"+GENEDATASET
	else: path = P.experiment_name+"/"+exp_type1+"/"+METHOD+"/"+GENEDATASET

	print(path)
	#continue
	listoffiles = os.listdir(path)

	#__________________Filling the lists with genes_______________________________
	DT_list = []
	RF_list =[]
	SVM_list = []
	Young_list = []
	Old_list =[]


	#*****************************************************************************
	if exp_type2 == "2": #Outlier determining genes
		for genelist in listoffiles:
			res = genelist.split("_")[0]
			f = open(path+"/"+genelist, "r")
			for i in f: 
				if res == "Young": Young_list.append(i.split("	")[0])
				elif res == "Old": Old_list.append(i.split("	")[0])

		#for i in range(0,10):
		#	print(Young_list[i], "	", Old_list[i])

		A = set(Young_list)
		B = set(Old_list)


	#*****************************************************************************
	else: #Machine learning results (the lists are filled twice for the outlier importance because there is young and old for every classifier)
		for genelist in listoffiles:
			res = re.findall(r'\[.*?\]', genelist)[0] #finding the type of machinelearning in this file
			f = open(path+"/"+genelist, "r")
			for i in f: 
				if res == "[DecisionTree]": DT_list.append(i.split(" ")[0]+", "+str(round(float(i.split(" ")[1]), 3)))
				elif res == "[RandomForest]": RF_list.append(i.split(" ")[0]+", "+str(round(float(i.split(" ")[1]), 3)))
				elif res == "[Support Vector Machine]": SVM_list.append(i.split(" ")[0]+", "+str(round(float(i.split(" ")[1]), 3)))


		if exp_type2 == '1': # for ever model there are 2 lists (yound and old)
			if use_only_overlap:
				print(set(DT_list[:int(len(DT_list)/2)]))
				exit()
				DT_list = list(set(DT_list[:int(len(DT_list)/2)] ) & set(DT_list[int(len(DT_list)/2):]))
				RF_list = list(set(RF_list[:int(len(RF_list)/2)] ) & set(RF_list[int(len(RF_list)/2):]))
				SVM_list = list(set(SVM_list[:int(len(SVM_list)/2)] ) & set(SVM_list[int(len(SVM_list)/2):]))
		
	#*****************************************************************************
	print("\n***********************************************************************")
	print("DT_list:	", list(DT_list)[:5], "...", len(DT_list))
	print("RF_list:	", list(RF_list)[:5], "...", len(RF_list))
	print("SVM_list:	", list(SVM_list)[:5], "...", len(SVM_list))
	print("Young_list:	", list(Young_list)[:5], "...", len(Young_list))
	print("Old_list:	", list(Old_list)[:5], "...", len(Old_list))
	print("***********************************************************************\n")
	
	#*****************************************************************************

	DT_Allsets[GENEDATASET] = DT_list
	RF_Allsets[GENEDATASET] = RF_list
	SVM_Allsets[GENEDATASET] = SVM_list
	
	




	#******************** Making lists for table ***********************************
	

	#save_csv(intersect, P.experiment_name+"/VennDiagrams/"+foldername+'/'+GENEDATASET+".csv")



#ding = pd.from_dict(data, orient='columns', dtype=None, columns=None)
DT_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in DT_Allsets.items() ]))
RF_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in RF_Allsets.items() ]))
SVM_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in SVM_Allsets.items() ]))

#ding = pd.DataFrame.from_dict(DT_Allsets)

print("DecisionTree:\n", DT_df)
print("RandomForest:\n", RF_df)
print("SupportVectorMachine\n", SVM_df)
input("Start making the plots shall we?")


#global biggest
biggest = 0

#*******************************************************************************
def nu(x, p):
	#print(total)
	global biggest
	for i, s in enumerate(x):
		try: 
			if not math.isnan(x[i]): pass # x[i] = 0
		except: 
			if p == 0: 
				x[i] = s.split(", ")[p]
				if len(x[i]) > biggest: biggest = len(x[i])
				kar = 10 - len(x[i])
				x[i] += kar*" "
			else: x[i] = float(s.split(", ")[p])
	return x
#*******************************************************************************

naammodel = ""
tel = 0
for model in [DT_df, RF_df, SVM_df]:
	tel += 1
	if tel == 1: naammodel = "DT"
	elif tel == 2: naammodel = "RF"
	elif tel == 3: naammodel = "SVM"

	C__df = model.copy()
	copy1__df = C__df.apply(lambda x: nu(x, 0), axis=1)
	copy2__df = model.apply(lambda x: nu(x, 1), axis=1)

	print("biggest: ", biggest)

	print(copy1__df)
	print(copy2__df)

	def get_size(num_inst):
		param = 1
		return num_inst*0.25

	#results = np.random.rand(DT_df.shape[0], DT_df.shape[1])
	results = copy2__df.to_numpy()
	strings = copy1__df.to_numpy()
	labels = 	(np.asarray(["{0} {1:.3f}".format(string, value) for string, value in zip(strings.flatten(), results.flatten())]) ).reshape(model.shape[0], model.shape[1])
	

	from pylab import rcParams
	rcParams['figure.figsize'] = 13, get_size(model.shape[0])#4 #18
	fig, ax = plt.subplots()

	x_axis_labels = list(copy1__df.columns) # labels for x-axis
	s = sns.heatmap(results, annot=labels, fmt="", cmap="Blues", ax=ax, linewidths=.5, yticklabels=False, xticklabels=x_axis_labels, cbar_kws={'label': 'Importance'}) #vmin=0, vmax=1)
	sns.set(font_scale = 1)
	s.set_xlabel('Different gene sets', fontsize=10)
	#ax.yaxis.set_label_position("right")
	plt.xticks(rotation=0, fontsize=14)
	
	ax.figure.axes[-1].yaxis.label.set_size(25)


	#manager = plt.get_current_fig_manager()
	#manager.full_screen_toggle()
	plt.tight_layout()
	plt.savefig(P.experiment_name+'/GenelistGradientTables/'+naammodel+'.png', dpi=500)
	#plt.show()
	plt.close()

exit()





print("Run was succesfull!")
#**************************** The End ******************************************






