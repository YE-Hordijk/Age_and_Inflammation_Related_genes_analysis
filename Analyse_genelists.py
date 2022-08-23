#Analyse_genelists.pyplot

from Parameters import P
import os
import matplotlib.pyplot as plt
import numpy as np
import re
import csv
import pandas as pd
from matplotlib_venn import venn2, venn3, venn3_circles

#*******************************************************************************
def save_csv(d, filename):
	#df = pd.DataFrame.from_dict(d)
	df = pd.DataFrame(dict([ (k,pd.Series(v, dtype=pd.StringDtype())) for k,v in d.items() ]))
	df.to_csv(filename, index= False)
#*******************************************************************************

if "VennDiagrams" not in os.listdir(P.experiment_name): #Making a folder for the machinelearning results
	os.mkdir(os.path.join(os.getcwd(), P.experiment_name+"/VennDiagrams"))


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
lijst = ["all", "senescence", "cell-age-signatures", "genes-from-papers", "searchwords"]



if foldername not in os.listdir(P.experiment_name+"/VennDiagrams"): #Making a folder for the machinelearning results
	os.mkdir(os.path.join(os.getcwd(), P.experiment_name+"/VennDiagrams/"+foldername))


#____________________Reading the lists_________________________________________
for GENEDATASET in lijst: 
	if exp_type2=="2": path = P.experiment_name+"/"+exp_type1+"/"+GENEDATASET
	else: path = P.experiment_name+"/"+exp_type1+"/"+METHOD+"/"+GENEDATASET


	listoffiles = os.listdir(path)

	#__________________Filling the lists with genes_______________________________
	if exp_type2 == "2":
		Young_list = []
		Old_list =[]
		for genelist in listoffiles:
			#res = re.findall(r'\[.*?\]', genelist)[0]
			res = genelist.split("_")[0]
			f = open(path+"/"+genelist, "r")
			for i in f: 
				if res == "Young": Young_list.append(i.split("	")[0])
				elif res == "Old": Old_list.append(i.split("	")[0])

		#for i in range(0,30):
		#	print(Young_list[i], "	", Old_list[i])
		
		#exit()
		A = set(Young_list)
		B = set(Old_list)

		v = venn2([A,B], ('Young', 'Old'))
		intersect = {}
		intersect["Young & Old"] = list(A&B)



	else:
		DT_list = []
		RF_list =[]
		SVM_list = []
		for genelist in listoffiles:
			res = re.findall(r'\[.*?\]', genelist)[0] #finding the type of machinelearning in this file
			f = open(path+"/"+genelist, "r")
			for i in f: 
				if res == "[DecisionTree]": DT_list.append(i.split(" ")[0])
				elif res == "[RandomForest]": RF_list.append(i.split(" ")[0])
				elif res == "[Support Vector Machine]": SVM_list.append(i.split(" ")[0])
	
		#print(SVM_list) #XXX
		#print(len(SVM_list)) #XXX
		
		
		if exp_type2 == '1': # for ever model there are 2 lists (yound and old)
			if use_only_overlap:
				DT_list = set(DT_list[:int(len(DT_list)/2)] ) & set(DT_list[int(len(DT_list)/2):])
				RF_list = set(RF_list[:int(len(RF_list)/2)] ) & set(RF_list[int(len(RF_list)/2):])
				SVM_list = set(SVM_list[:int(len(SVM_list)/2)] ) & set(SVM_list[int(len(SVM_list)/2):])
		
		
		A = set(DT_list)
		B = set(RF_list)
		C = set(SVM_list)
		
		#print(C) #XXX
		#print(len(C)) #XXX
		#exit()
		
		v = venn3([A,B,C], ('DecisionTree', 'Randomforest', 'Support Vector Machine'))
		
		intersect = {}
		intersect["DT & RF & SVM"] = list(A&B&C)
		intersect["DT & RF"] = list(A&B)
		intersect["DT & SVM"] = list(A&C)
		intersect["RF & SVM"] = list(B&C)




	plt.savefig(P.experiment_name+'/VennDiagrams/'+foldername+'/'+GENEDATASET+'.png')
	#plt.show()
	plt.close()

	#******************** Making lists for table ***********************************
	

	save_csv(intersect, P.experiment_name+"/VennDiagrams/"+foldername+'/'+GENEDATASET+".csv")

print("Run was succesfull!")
#**************************** The End ******************************************


"""
print("--------------------------------------\n")
AA = set([1,2,3,4])
BB = set([5,6,7,8])
CC = set([1,3,5,7])
D = AA-CC
print(D)

print("--------------------------------------\n")

#Setting fontsize for genes and labels
for text in v.set_labels: text.set_fontsize(15) #fontsize for labels
for text in v.subset_labels: text.set_fontsize(8) #fontsize for genes

def circle_content(content):
	tel = 0
	c2 = ""
	for i in content:
		if tel%2 == 0:
			c2 += i+" "
		else: c2 += i+"\n"
		tel += 1
	print(c2)
	return c2

ppp = v.get_label_by_id('100').set_text(circle_content(A-B-C))

#ppp = v.get_label_by_id('100').set_text(','.join(A-B-C))
v.get_label_by_id('110').set_text('\n'.join(A&B-C))
v.get_label_by_id('011').set_text('\n'.join(B&C-A))
v.get_label_by_id('001').set_text('\n'.join(C-A-B))
v.get_label_by_id('010').set_text('')

#plt.annotate(',\n'.join(B-A-C), xy=v.get_label_by_id('010').get_position() + np.array([0, 0.2]), xytext=(-20,40), ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1), arrowprops=dict(arrowstyle='->', connectionstyle='arc',color='gray'))
"""






