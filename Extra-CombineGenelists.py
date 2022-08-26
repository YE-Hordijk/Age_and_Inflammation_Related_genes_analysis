#Extra-CombiningGenelists.py

from Parameters import P
import os
import matplotlib.pyplot as plt
import numpy as np
import re
import csv
import pandas as pd
from matplotlib_venn import venn2, venn3, venn3_circles
import math

def save_csv(d, filename):
	#df = pd.DataFrame.from_dict(d)
	df = pd.DataFrame(dict([ (k,pd.Series(v, dtype=pd.StringDtype())) for k,v in d.items() ]))
	df.to_csv(filename, index= False)

#*******************************************************************************
GeneOccurence = {}
for FILE in os.listdir(P.experiment_name+'/Genelists'):
	print(FILE)
	f = open(P.experiment_name+'/Genelists/'+FILE, "r")
	for x in f:
		line = x[:-1]
		if line not in GeneOccurence: GeneOccurence[line] = 0
		GeneOccurence[line] += 1
		
		#print(line[:-1])
		#input()
#print(GeneOccurence)

#*******************************************************************************
path = P.experiment_name+'/GenelistGradientTables/'
threshold = 100 #P.NrFoundRelevantGenes

finallist = []
for i in os.listdir(path):
	if i[-4:] == ".csv": finallist.append(i)

print(finallist)
Total = []
MLOutliers = []
CompareOutliers = []

for i in finallist:
	if i[:18] == "TotalGeneset_MLres": Total.append(i)
	elif i[:18] == "CompareOutliers_ML": MLOutliers.append(i)
	elif i[:18] == "CompareOutliers_Ou": CompareOutliers.append(i)




if "CombinedGenelists!" not in os.listdir(P.experiment_name): #Making a folder for the machinelearning results
	os.mkdir(os.path.join(os.getcwd(), P.experiment_name+"/"+"CombinedGenelists!"))

#Make a combined list of important genes
combine_these_lists = MLOutliers + CompareOutliers
#combine_these_lists.remove("CompareOutliers_MLresults_DT.csv")
#combine_these_lists.remove("CompareOutliers_MLresults_RF.csv")
#combine_these_lists.remove("CompareOutliers_MLresults_SVM.csv")

for i in combine_these_lists:
	best = {}
	df = pd.read_csv(path+"/"+i)
	print(df)
	tel = 0
	for index, row in df.iterrows():
		for col in df.columns:
			tel += 1
			try: 
				if math.isnan(row[col]): pass
			except:
				if row[col].split(", ")[0] in best: best[row[col].split(", ")[0]] += abs(float(row[col].split(", ")[1])) #if item is already in best, add the value to this item to make it worth more
				elif not len(best)>=threshold or not all(value > abs(float(row[col].split(", ")[1])) for value in best.values()): #if best is not full or there are smaller items
					best[row[col].split(", ")[0]] = abs(float(row[col].split(", ")[1])) # Add to the best dictionary
					if len(best) > threshold: del best[min(best, key=best.get)] # remove smallest item if dict is full

	for q in best:
		best[q] = best[q]/(GeneOccurence[q])
	
	print(i)
	
	tellie = 0
	for k in best:
		tellie +=  1
		print("\item", tellie, k, "  ", round(best[k], 3))
	print("The above list is from ", i, "You can copy the list from the terminial, but it will also be saved in the folder CombinedGenelists!")
	input()
	P.experiment_name+"/"+"CombinedGenelists!"
	f = open(P.experiment_name+"/"+"CombinedGenelists!/"+i.split("_")[-1], "w+")
	for z in best:
		f.write(z+", "+str(best[z])+"\n")
	f.close


#Make Venn diagram if three lists were created above
if len(combine_these_lists)  <= 3:
	tell = -1
	A = []
	B = []
	C = []
	vennen = [A,B,C]
	vennamen = []
	for i in combine_these_lists:
		tell += 1
		f = open(P.experiment_name+"/"+"CombinedGenelists!/"+i.split("_")[-1], "r")
		for k in f: vennen[tell].append(k.split(", ")[0])
		f.close()
		vennamen.append(i.split("_")[-1].split(".")[0])
	
	A = set(A)
	B = set(B)
	C = set(C)
	if len(combine_these_lists) == 3: v = venn3([A,B,C], (vennamen))
	elif len(combine_these_lists) == 2: v = venn2([A,B], (vennamen))
	
	intersect = {}
	if len(combine_these_lists) == 3:
		intersect[vennamen[0]+" & "+vennamen[1]+" & "+vennamen[2]] = list(A&B&C)
	intersect[vennamen[0]+" & "+vennamen[1]] = list(A&B.symmetric_difference(list(A&B&C)))
	
	if len(combine_these_lists) == 3: 
		intersect[vennamen[0]+" & "+vennamen[2]] = list(A&C.symmetric_difference(list(A&B&C)))
		intersect["RF & "+vennamen[2]] = list(B&C.symmetric_difference(list(A&B&C)))


	if len(combine_these_lists) == 3: between = vennamen[0]+"-"+vennamen[1]+"-"+vennamen[2]
	elif len(combine_these_lists) == 2: between = vennamen[0]+"-"+vennamen[1]
	
	
	plt.savefig(P.experiment_name+"/"+"CombinedGenelists!/"+between+'_Venn.png', bbox_inches = 'tight')
	plt.show()
	plt.close()
	
	save_csv(intersect, P.experiment_name+"/"+"CombinedGenelists!/"+between+'_Venn.csv')




