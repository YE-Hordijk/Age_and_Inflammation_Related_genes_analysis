#Add LaTeX tables for regression

from Parameters import P
import os
tel = 0
nrfiles = 0
if "Machine_Learning_Results" in os.listdir(P.experiment_name): 
	for folder in os.listdir(P.experiment_name+"/Machine_Learning_Results"):
		if folder[0] == "#": continue
		tel = 0
		print(folder)
		c = open(P.experiment_name+"/Machine_Learning_Results/"+folder+"/"+folder+"concatenated.txt", "w")
		
		##moving "all" to the last position
		listoffiles = os.listdir(P.experiment_name+"/Machine_Learning_Results/"+folder)
		thisfile = allfile = random = ""
		for k in listoffiles:
			if k == folder+"concatenated.txt": thisfile = k
			if k[0:3] == "all": allfile = k
			if k[0:6] == "Random": random = k
		listoffiles.remove(thisfile)
		listoffiles.remove(allfile)
		listoffiles.remove(random)
		listoffiles.sort(reverse = True)
		listoffiles.append(allfile)
		listoffiles.append(random)
		nrfiles = len(listoffiles) #NUMBER OF FILES IN THIS FOLDER
		#print(folder, "\n", listoffiles)
		input("volgende model?")
		
		for File in listoffiles:
			f = open(P.experiment_name+"/Machine_Learning_Results/"+folder+"/"+File, "r")
			subtel = 0
			tel += 1
			for i in f: 
				subtel += 1
				if tel == 1: #if the first file
					if subtel >= 3 and subtel <= 13:
						c.write(i)

				elif tel == nrfiles: #the last filename
					if subtel == 15: #caption
						c.write("\t\caption{Evaluation of Regression by "+folder+" using different datasets.}\n")
					elif subtel == 16: #label
						c.write("\t\label{tab:Regression"+folder+"No-MiddleAge}\n")
					elif subtel >= 11: 
						c.write(i)

				else: #a middle file
					if subtel >= 11 and subtel <= 13:
						c.write(i)

			f.close()
		c.close()
	exit()
