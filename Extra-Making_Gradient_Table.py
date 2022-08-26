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
DETERMINING_Allsets = {}


#____________________Reading the lists_________________________________________
for GENEDATASET in lijst: 
	if exp_type2=="2": path = P.experiment_name+"/"+exp_type1+"/"+GENEDATASET
	else: path = P.experiment_name+"/"+exp_type1+"/"+METHOD+"/"+GENEDATASET

	print(path)
	#continue
	listoffiles = os.listdir(path)

	#__________________Filling the lists with genes_______________________________
	DT_list = [[],[]]
	RF_list =[[], []]
	SVM_list = [[], []]
	Young_list = []
	Old_list =[]
	


	#*****************************************************************************
	if exp_type2 == "2": #Outlier determining genes
		for genelist in listoffiles:
			res = genelist.split("_")[0]
			f = open(path+"/"+genelist, "r")
			for i in f: 
				if res == "Young": Young_list.append(i.split("	")[0]+", "+str(round(float(i.split("	")[1]), 2)))
				elif res == "Old": Old_list.append(i.split("	")[0]+", "+str(round(float(i.split("	")[1]), 2)))

		#for i in range(0,10):
		#	print(Young_list[i], "	", Old_list[i])

		#A = set(Young_list)
		#B = set(Old_list)


	#*****************************************************************************
	else: #Machine learning results (the lists are filled twice for the outlier importance because there is young and old for every classifier)
		for genelist in listoffiles:
			print("-------------> ", genelist[14:18])
			res = re.findall(r'\[.*?\]', genelist)[0] #finding the type of machinelearning in this file
			f = open(path+"/"+genelist, "r")
			for i in f: 
				if res == "[DecisionTree]": 
					if genelist[14:18] == "_Old": DT_list[1].append(i.split(" ")[0]+", "+str(round(float(i.split(" ")[1]), 3))) # Old
					else: 												DT_list[0].append(i.split(" ")[0]+", "+str(round(float(i.split(" ")[1]), 3))) # Young

				elif res == "[RandomForest]": 
					if genelist[14:18] == "_Old": RF_list[1].append(i.split(" ")[0]+", "+str(round(float(i.split(" ")[1]), 3))) # Old
					else: 												RF_list[0].append(i.split(" ")[0]+", "+str(round(float(i.split(" ")[1]), 3))) # Young
					
				elif res == "[Support Vector Machine]": 
					if genelist[14:18] == "_Old": SVM_list[1].append(i.split(" ")[0]+", "+str(round(float(i.split(" ")[1]), 3))) # Old
					else: 												SVM_list[0].append(i.split(" ")[0]+", "+str(round(float(i.split(" ")[1]), 3))) # Young
		
		print(DT_list[0], "\n", DT_list[1])

		if exp_type2 == '1': # for ever model there are 2 lists (yound and old)
			if use_only_overlap:
				#print(set(DT_list[:int(len(DT_list)/2)]))
				DT_list[0] = list(set(DT_list[0]) & set(DT_list[1]))
				RF_list[0] = list(set(RF_list[0]) & set(RF_list[1]))
				SVM_list[0] = list(set(SVM_list[0] ) & set(SVM_list[1]))
		
	#*****************************************************************************
	print("\n***********************************************************************")
	print("DT_list:	", list(DT_list[0])[:5], "...", len(DT_list[0]))
	print("RF_list:	", list(RF_list[0])[:5], "...", len(RF_list[0]))
	print("SVM_list:	", list(SVM_list[0])[:5], "...", len(SVM_list[0]))
	print("Young_list:	", list(Young_list)[:5], "...", len(Young_list))
	print("Old_list:	", list(Old_list)[:5], "...", len(Old_list))
	print("***********************************************************************\n")
	
	#*****************************************************************************
	if GENEDATASET == "genes-from-papers": GENEDATASET = "papergenes"
	if GENEDATASET == "cell-age-signatures": GENEDATASET = "age-signa."
	if foldername == "TotalGeneset_MLresults": 
		DT_Allsets[GENEDATASET] = DT_list[0]
		RF_Allsets[GENEDATASET] = RF_list[0]
		SVM_Allsets[GENEDATASET] = SVM_list[0]

	elif foldername == "CompareOutliers_MLresults": 
		DT_Allsets[GENEDATASET] = {"Young": DT_list[0], "Old": DT_list[1]}
		RF_Allsets[GENEDATASET] = {"Young": RF_list[0], "Old": RF_list[1]}
		SVM_Allsets[GENEDATASET] = {"Young": SVM_list[0], "Old": SVM_list[1]}

	elif foldername == "CompareOutliers_Outlierdeteremining": 
		DETERMINING_Allsets[GENEDATASET] = {"Young": Young_list, "Old": Old_list}
		#exit("Nog niet geimplementeerd")


	else: exit("WTF")
	
#print(DT_Allsets)

	
	




	#******************** Making lists for table ***********************************
	

	#save_csv(intersect, P.experiment_name+"/VennDiagrams/"+foldername+'/'+GENEDATASET+".csv")


if foldername == "CompareOutliers_Outlierdeteremining":
	DET_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in DETERMINING_Allsets.items() ]))
	DET_dft = {} 
else:
	DT_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in DT_Allsets.items() ]))
	RF_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in RF_Allsets.items() ]))
	SVM_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in SVM_Allsets.items() ]))
	print("DecisionTree:\n", DT_df)
	print("RandomForest:\n", RF_df)
	print("SupportVectorMachine\n", SVM_df)
	DT_dft = {}
	RF_dft = {}
	SVM_dft = {}




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
				kar = 8 - len(x[i])
				
				if kar>0: x[i] += kar*" "
			else: x[i] = float(s.split(", ")[p])
	return x
#*******************************************************************************

naammodel = ""
tel = 0


if foldername == "TotalGeneset_MLresults": things = [DT_df, RF_df, SVM_df]
elif foldername == "CompareOutliers_MLresults": 
	for column in DT_df:
		DT_dft[column+"_Young"] = DT_df[column].values[0]
		DT_dft[column+"_Old"] = DT_df[column].values[1]
	for column in RF_df:
		RF_dft[column+"_Young"] = RF_df[column].values[0]
		RF_dft[column+"_Old"] = RF_df[column].values[1]
	for column in SVM_df:
		SVM_dft[column+"_Young"] = SVM_df[column].values[0]
		SVM_dft[column+"_Old"] = SVM_df[column].values[1]

	DT_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in DT_dft.items() ]))
	RF_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in RF_dft.items() ]))
	SVM_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in SVM_dft.items() ]))

	things = [DT_df, RF_df, SVM_df]

elif foldername == "CompareOutliers_Outlierdeteremining":
	for column in DET_df:
		DET_dft[column+"_Young"] = DET_df[column].values[0]
		DET_dft[column+"_Old"] = DET_df[column].values[1]
	DET_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in DET_dft.items() ]))
	things = [DET_df]

else: exit("WTF")


names = {}
names["TotalGeneset_MLresults"] = ["DT", "RF", "SVM"]
#names["CompareOutliers_MLresults"] = ["DTyoung", "DTold", "RFyoung", "RFold", "SVMyoung", "SVMold"]

for model in things: 
	
	if foldername == "CompareOutliers_Outlierdeteremining":
		naammodel = "MeanDifference"
	else: naammodel = names["TotalGeneset_MLresults"][tel]
	tel += 1
	
	model.index += 1 
	print(model)
	model.to_csv(P.experiment_name+'/GenelistGradientTables/'+foldername+'_'+naammodel+'.csv', index=False)


	C__df = model.copy()
	
	maxsize = 50
	if C__df.shape[0] > maxsize:
		C__df = C__df.iloc[0:maxsize, :]
	
	print(C__df)
	

	#print("*********************************************\n", C__df)
	copy1__df = C__df.apply(lambda x: nu(x, 0), axis=1)
	copy2__df = model.apply(lambda x: nu(x, 1), axis=1)
	if copy2__df.shape[0] > maxsize:
		copy2__df = copy2__df.iloc[0:maxsize, :]
	#exit("YEYEYEYEYYE")
	print("biggest: ", biggest)

	print(copy1__df)
	print(copy2__df)
	

	def get_size(num_inst):
		if num_inst<10:
			return 3
		return num_inst*0.4 #0.25

	#results = np.random.rand(DT_df.shape[0], DT_df.shape[1])
	results = copy2__df.to_numpy()
	strings = copy1__df.to_numpy()
	labels = 	(np.asarray(["{0} {1:.2f}".format(string, value) for string, value in zip(strings.flatten(), results.flatten())]) ).reshape(C__df.shape[0], C__df.shape[1])
	

	from pylab import rcParams
	if C__df.shape[1] > 5: rcParams['figure.figsize'] = 18, get_size(C__df.shape[0])#4 #18
	else: rcParams['figure.figsize'] = 13, get_size(C__df.shape[0])#4 #18
	fig, ax = plt.subplots()
	
	
	#plt.ylabel(r"My long label with unescaped {\LaTeX} $\Sigma_{C}$ math"
				#"\n"  # Newline: the backslash is interpreted as usual
				#r"continues here with $\pi$")
	
	#print(list(copy1__df.columns))
	newlist = [x.replace("_", "\n") for x in copy1__df.columns]
	#exit()
	x_axis_labels = newlist #list(copy1__df.columns) # labels for x-axis
	palette = "Blues"
	if foldername == "CompareOutliers_Outlierdeteremining": palette = "coolwarm"
	s = sns.heatmap(results, annot=labels, fmt="", cmap=palette, ax=ax, linewidths=.5, yticklabels=True, xticklabels=x_axis_labels, cbar_kws={'label': 'Importance'}) #vmin=0, vmax=1)
	sns.set(font_scale = 1)
	s.set_xlabel('Different gene sets', fontsize=20)
	#ax.yaxis.set_label_position("right")
	plt.xticks(rotation=0, fontsize=14)
	plt.yticks(np.arange(0.5, 0.5+C__df.shape[0]), np.arange(1, C__df.shape[0]+1), rotation=0, fontsize=14)
	
	ax.figure.axes[-1].yaxis.label.set_size(25)


	#manager = plt.get_current_fig_manager()
	#manager.full_screen_toggle()
	plt.tight_layout()
	plt.savefig(P.experiment_name+'/GenelistGradientTables/'+foldername+'_'+naammodel+'.png', dpi=120, bbox_inches = 'tight')
	#plt.show()
	plt.close()


exit()





print("Run was succesfull!")
#**************************** The End ******************************************





