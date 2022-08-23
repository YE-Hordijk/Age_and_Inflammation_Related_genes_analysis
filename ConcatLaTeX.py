#Add LaTeX tables

from Parameters import P
import os
thisfile = ""

#*******************************************************************************
def prepare_info_TGS(listoffiles):
	global thisfile
	thisfile = allfile = random = ""
	for k in listoffiles:
		if k == folder+"concatenated.tex": thisfile = k
		if k[0:3] == "all": allfile = k
		if k[0:7] == "RandomB": random = k
	if thisfile in listoffiles: listoffiles.remove(thisfile)
	listoffiles.remove(allfile)
	if random not in listoffiles: 
		cont = input("You forgot your baseline. Continue without? Yes[1]/No[2]")
		if cont == "2": exit()
	else: listoffiles.remove(random)
	listoffiles.sort(reverse = True)
	listoffiles.append(allfile)

	if random != "": listoffiles.append(random)
	nrfiles = len(listoffiles) #NUMBER OF FILES IN THIS FOLDER
	#print(folder, "\n", listoffiles)
	return listoffiles, nrfiles
#*******************************************************************************
def prepare_info_CO(listoffiles):
	
	thisfile, allfile, random = {}, {}, {} #setting everything to empty string
	ALLLISTS = [thisfile, allfile, random]
	for k in listoffiles:
		if k == folder+"concatenated.tex": thisfile[k] = None
		if k[0:3] == "all": allfile[k] = None
		if k[0:7] == "RandomB": random[k] = None

	for i in ALLLISTS:
		for j in i:
			listoffiles.remove(j)


	listoffiles.sort(reverse = True)
	for i in ALLLISTS:
		if i != thisfile:
			for j in i:
				listoffiles.append(j)

	#for i in listoffiles:
	#	print("--", i)
	
	nrfiles = len(listoffiles) #NUMBER OF FILES IN THIS FOLDER
	#print(folder, "\n", listoffiles)
	return listoffiles, nrfiles
#*******************************************************************************



regels = {"NoMiddleAge": [3,11,13,15,16,7],
					"WithMiddleAge": [3,11,14,16,17,7]
					}
color = 2
tel = 0
nrfiles = 0
regel = "NoMiddleAge"
m = input("Classification[1] or Regression[2]?")
z = input("TotalGeneDataset[1] or Compare_outliers[2]?")
if z == 1:
	l = input("No Middleage[1] or With Middleage[2]")
	if l == "2": regel = "WithMiddleAge" 
else: color = input("Use colors? Yes[1], No[2]")

if (m == "1"): method = "Classification"
elif (m == "2"): method = "Regression"

if z=="1": path = P.experiment_name+"/Machine_Learning_Results"
elif z=="2": path = P.experiment_name+"/Compare_outliers/Outlier_ML_Tables"



if not method in os.listdir(path): 
	exit(method+" folder does not exist")

else:
	for folder in os.listdir(path+"/"+method):
		if folder[0] == "#": continue
		tel = 0
		print(folder)
		c = open(path+"/"+method+"/"+folder+"/"+folder+"concatenated.tex", "w")
		
		#moving "all" to the last position
		listoffiles = os.listdir(path+"/"+method+"/"+folder)
		
		if z=="1": listoffiles, nrfiles = prepare_info_TGS(listoffiles)
		elif z=="2": listoffiles, nrfiles = prepare_info_CO(listoffiles)
		#else: exit("verkeerde waarde ingevoerd bij 2e vraag")




		input("volgende model?")
		colorflip = 0
		
		for File in listoffiles:
			f = open(path+"/"+method+"/"+folder+"/"+File, "r")
			subtel = 0
			
			if (tel%2)==0 and tel!=0 and z=="2": 
				c.write("\t\t\hline\n")
				if colorflip == 0: colorflip = 1
				elif colorflip == 1: colorflip = 0
				#print(colorflip)
			#print(colorflip,"\n")
		
			tel += 1
			for i in f: #looping over the readfile
				subtel += 1
				
				if color=="1" and subtel >= regels[regel][1] and subtel < regels[regel][2]:
					if colorflip==0: c.write("\t\t\\rowcolor{green!80!yellow!50}\n")
					else: c.write("\t\t\\rowcolor{green!40!yellow!40}\n")
					
					
				if tel == 1: #if the first file
					if subtel >= regels[regel][0] and subtel <= regels[regel][2]:
						c.write(i)

				elif tel == nrfiles: #the last file
					if subtel == regels[regel][3]: #caption
						if z=="1": c.write("\t\caption{Evaluation of "+method+" by "+folder+" using different datasets.}\n")
						elif z=="2": c.write("\t\caption{Evaluation of "+method+" by "+folder+" using different datasets, trying to predict if a sample is an outlier from its own agegroup.}\n")
					elif subtel == regels[regel][4]: #label
						if z == "1": c.write("\t\label{tab:"+method+folder+"No-MiddleAge}\n")
						elif z=="2": c.write("\t\label{tab:CompareOutliers"+method+folder+"}\n")
					elif subtel >= regels[regel][1]: 
						c.write(i)

				else: #a middle file
					if subtel >= regels[regel][1] and subtel <= regels[regel][2]:
						c.write(i)
				
				
				
	
				

			f.close()
		c.close()
exit()
