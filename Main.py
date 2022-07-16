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
def update_parameter(par):
	#Giving a string quotationmarks
	if not (par.split(" = ")[1] in ["False", "True"]) and not isinstance(par.split(" = ")[1], int) and not par.split(" = ")[1].startswith(("{", "[")):
		par = par.split(" = ")[0] + " = " + "\""+par.split(" = ")[1]+"\"" #putting quotation marks around the value
	found = False
	f = open("Parameters.py", "r+")
	d = f.readlines()
	
	f.seek(0) #to the beginning of the file
	for i in d:
		if (par.split(" = ")[0] not in i): f.write(i) #if this line is not to be over writen by the new parameter, leave it by copying
		else:
			f.write(par+"\n")
			found = True
	if found == False:
		f.write(par+"\n")
	f.truncate()
	f.close
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
#************************************CODE***************************************
#*******************************************************************************







#*******************************************************************************

# 1) Setting parameters with which the experiment will be performed
"""
class S: #Settings
	experiment_name = "expTESTING_BEREN"
	update_files = False
	#⬇️For GENESELECTION: don't use underscores. It should be either: "all", "senescence", "searchwords","cell-age-signatures", "genes-from-papers"
	GENE_SELECTION = "cell-age-signatures" #"genes-from-papers" 
	inflam_synonyms = {'senescence','inflam', 'infection', ' cytokines', 'cytokine rush', 'immune response' , 'antibody'}
	tissue = "Whole Blood"
	select_on_genes = True #DONT PUT IT ON FALSE (r CANT HANDLE IT))
	use_middle_age = False
	YOUNG = "20-49"
	MIDDLE = "50-59"
	OLD = "60-79"
	random_baseline = False
	removing_outliers = False #True
	PredictionModel = "RandomForest"#"DecisionTree"#"RandomForest" #"DecisionTree" #"Support Vector Machine"
	PredictionMethod = "Classification" #"Regression"
"""
#*******************************************************************************
def calculate_filenames(Use_middle_age, Select_on_genes, Geneselection):
	#Calculating the names of the metadatafile and the countsfile
	if Select_on_genes: metafilename = Geneselection #variable for nameing files and folders
	else: metafilename = "NoGeneSelection"
	if Use_middle_age:	metafilename += "_With-MiddleAge_METADATA.txt"
	else: 								metafilename += "_No-MiddleAge_METADATA.txt"
	countsfilename = metafilename[:-12] + "unnormalized.txt"
	return countsfilename, metafilename



countsfilename, metafilename = calculate_filenames(P.use_middle_age, P.select_on_genes, P.GENE_SELECTION)




# 2)-------------------------------------------------------------------------------
setup_experiment(P.experiment_name)

# 3)-------------------------------------------------------------------------------
#Sending the parameters to extern script so that all scripts can use these parameters
#for attr in vars(S):
#	if not attr.startswith("__"):
#		str(attr+" = "+str(vars(S)[attr]))
#		update_parameter(str(attr+" = "+str(vars(S)[attr])))

#import test as tt
#tt.TEST("appel")
#tt.TEST("peer")

"""
for g in ["senescence", "searchwords","cell-age-signatures", "genes-from-papers", "all"]: #XXX
	P.GENE_SELECTION = g #Setting the parameter
	print(g, "-- ", P.GENE_SELECTION)
	#input("apen?")
	
	
	#***CREATING A GENELIST****
	if (P.select_on_genes):
		print(st.GREEN, "\n*********** CREATE GENELIST **********", st.RST)
		import Create_genelist as CG
		genedict = CG.create_genelist()
		#subprocess.call ("/usr/bin/python3 Create_genelist.py "+arguments(P.update_files, P.GENE_SELECTION, shell=True)


	#****PREPROCESSING THE DATA for R******
	print(st.GREEN, "\n*********** PREPROCESSING DATA FOR R **********", st.RST)
	#subprocess.call ("/usr/bin/python3 Preprocessing.py "+arguments(P.update_files, P.GENE_SELECTION), shell=True)
	import Preprocessing as Pr
	Pr.preprocessing()

	#****USING R******
	#countsfilename, metafilename = calculate_filenames(P.use_middle_age, True, g) #XXX
	print(st.GREEN, "\n*********** NORMALIZING AND VISUALIZING WITH R **********", st.RST)
	subprocess.call ("/usr/bin/Rscript --vanilla Normalize_and_visualize.R "+arguments(countsfilename, metafilename, P.experiment_name), shell=True) #Arguments: 1)countfile, 2)metadate, 3)project name



	for i in ["DecisionTree", "RandomForest", "Support Vector Machine"]: #XXX
		#update_parameter("PredictionModel = "+i) #XXX
		P.MODEL = i

		#****Machinelearning****
		print(st.GREEN, "\n*********** MACHINE LEARNING **********", st.RST)
		print(st.GREEN, "\b**", g, "\n**",i,st.RST)
		#subprocess.call ("/usr/bin/python3 Machinelearning.py "+arguments(), shell=True)
		import Machinelearning as Ms
		Ms.machinelearning()
		print("terug van Machine Learning 😃️")
"""


#****Extracting important genes from PCfiles and finding outliers
print(st.GREEN, "\n*********** Use PCs for MachineLearning2 **********", st.RST)
import Use_PCs_for_ML2 as ML2
ML2.use_pcs_for_ml2()
#subprocess.call ("/usr/bin/python3 ExtractingPCs.py "+arguments(), shell=True)


#****DO R again
#TODO

#******Machinelearning again
#TODO

#**** PRODUCE Gene list_x
#TODO



