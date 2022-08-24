#Main.py
#Packages
import subprocess
import os
import re
from GeneralFunctions import st

#import Parameters as P
from Parameters import P


#*******************************************************************************
#******************************** FUNCTIONS ************************************
#*******************************************************************************
#*******************************************************************************
def arguments(*Args):
	temp = ""
	for arg in Args:
		temp += " "+str(arg)
	return temp
#*******************************************************************************
def setup_experiment(experiment_name):
	if experiment_name not in [f for f in os.listdir('.') if os.path.isdir(f)]:
		os.mkdir(os.path.join(os.getcwd(), experiment_name))
#*******************************************************************************
def calculate_filenames(Use_middle_age, Select_on_genes, Geneselection):
	#Calculating the names of the metadatafile and the countsfile
	if Select_on_genes: metafilename = Geneselection #variable for nameing files and folders
	else: metafilename = "NoGeneSelection"
	if Use_middle_age:	metafilename += "_With-MiddleAge_METADATA.txt"
	else: 								metafilename += "_No-MiddleAge_METADATA.txt"
	countsfilename = metafilename[:-12] + "unnormalized.txt"
	return countsfilename, metafilename
#*******************************************************************************
#************************************CODE***************************************
#*******************************************************************************

countsfilename, metafilename = calculate_filenames(P.use_middle_age, P.select_on_genes, P.GENE_SELECTION)
setup_experiment(P.experiment_name)


if P.random_baseline: datasets = ["senescence"]
else: datasets = ["senescence", "searchwords","cell-age-signatures", "genes-from-papers", "all"]
for g in datasets: #XXX
	P.GENE_SELECTION = g #Setting the parameter
	print("Dataset: ", g)


	#***CREATING A GENELIST****
	if (P.select_on_genes):
		print(st.GREEN, "\n*********** CREATE GENELIST **********", st.RST)
		import Create_genelist as CG
		genedict = CG.create_genelist()


	#****PREPROCESSING THE DATA for R******
	print(st.GREEN, "\n*********** PREPROCESSING DATA FOR R **********", st.RST)
	import Preprocessing as Pr
	Pr.preprocessing()


	#****USING R******
	countsfilename, metafilename = calculate_filenames(P.use_middle_age, True, g) #XXX
	print(st.GREEN, "\n*********** NORMALIZING AND VISUALIZING WITH R **********", st.RST)
	subprocess.call ("/usr/bin/Rscript --vanilla Normalize_and_visualize.R "+arguments(countsfilename, metafilename, P.experiment_name), shell=True) #Arguments: 1)countfile, 2)metadate, 3)project name
	
	
	for i in ["DecisionTree","RandomForest","Support Vector Machine"]:
		P.MODEL = i
		
		#****Machinelearning****
		print(st.GREEN, "\n*********** MACHINE LEARNING **********", st.RST)
		import Machinelearning as Ms
		Ms.machinelearning()
		
		
		#****Extracting important genes from PCfiles and finding outliers
		print(st.GREEN, "\n*********** Use PCs for MachineLearning2 **********", st.RST)
		import Use_PCs_for_ML2 as ML2
		ML2.use_pcs_for_ml2()

		



while True:
	print('\a\b\b', end='')

#****DO R again
#TODO
#subprocess.call ("/usr/bin/python3 Machinelearning.py "+arguments(), shell=True)




