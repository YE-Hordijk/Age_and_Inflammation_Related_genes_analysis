#Create_genelist.py

#Packages
import GeneralFunctions as GF
from GeneralFunctions import st
#import Parameters as P
from Parameters import P

import sys
import pandas as pd
import os
print("*Starting script Create_genelist.py")

#*******************************************************************************
#***********************Importing and setting variables*************************
#*******************************************************************************
#print(sys.argv[1])
#print(sys.argv[2])
#print sys.argv[1] 

class CG:
	cellAge = pd.DataFrame()
	cellSignatures = pd.DataFrame()




#*******************************************************************************
#MOST IMPORTANT GENES-LISTS COPIED FORM PAPERS
Asso_Inflam_Age = ['FCGR1A','GRAMD1C','PNP','CD274','DPEP2','TSTA3','KREMEN1',
										'SLC26A8','ALAS2','NDEL1','SAMD9L','RPS6KA5','SIAH1','GBP7',
										'SIAH2','KCNH7','DAPK2','SLC22A4','FAM20B','UBE2J1',
										'SLC25A37','HDGF','METTL9','GBP1','CD2AP'] #Transcriptome‐wide association study of inflammatory biologic age (Table 2)

sign_pathway = ['GABARAP','GBP4','GBP5','MAPK14','CTSB','NLRP1','GABARAPL1',
								'TAB3','GBP1','BIRC3','IFI16','IL1B','GBP7','MYD88','OAS1',
								'STAT1','CASP5','AIM2','UBE2F','UBOX5','NEDD4L','MDM2','UBE2J1',
								'UBE2O','SMURF2','SIAH1','ITCH','CUL4A','UBE2L6','CDC34','CLK1',
								'APAF1','ARF1','TLR5','CASP7','ARAP2','AP2M1','AP2S1','ZFYVE27',
								'AP2A1','IGF1R','IGF2R','LDLR','PRKCZ','PIP5K1B','SNX3',
								'DNAJC6','AKT3','AKT1','FLT3','HK1','SLC1A5','SLC2A1','SLC7A5',
								'KLC3','KLC1','PKN2','WASF1','RPS6KA5','PSME3','DDX58',
								'EIF2AK1','PPP2R2D','JAK3','IFIH1','CMPK2','DCK','PNP','RRM2B',
								'POLR1D','POLR2B','UPP1','ENTPD5'] #Transcriptome‐wide association study of inflammatory biologic age (Table 3)

GWAS_Loci = ['APOE','CHRNA3','CHRNA5','LPA','CDKN2A','CDKN2B','USP42','TMTC2',
							'IL6','ANKRD20A9P','LINC02227','FOXO3A146','RAD50','MC2R',
							'USP2-AS1','HLA-DQA1','HLA-DRB1','ATXN2','FURIN','EPHX2','PROX2',
							'CELSR2','PSRC1'] #Facing up to the global challenges of ageing (Table 1)

Age_asso_genes = ['CD248','LRRN3','NELL2','LEF1','CCR7','ABLIM1','GZMH',
											'MYC','CD27','FAM102A','SERPINE2','SLC16A10','FCGBP',
											'GPR56','BACH2','SYT11','PDE9A','NG','FLNB','NT5E',
											'FGFBP2','TGFBR3','ITM2C','ATF7IP2','CR2','FAIM3','PHGDH',
											'LDHB','SIRPG','FCRL6','PDE7A','NSIP','PAICS','BZW2',
											'OXNAD1','CX3CR1','SCML1','RPL22','LDLRAP1','RHOC','LTB',
											'FAM134B','LBH','PRSS23','SUSD3','PIK3IP1','MFGE8',
											'AGMAT','NKG7','PPP2R2B'] #The transcriptional landscape of age in human peripheral blood (Table 1)

genes_from_papers = sign_pathway + Age_asso_genes + GWAS_Loci + Asso_Inflam_Age
#*******************************************************************************
#*****************************Functions*****************************************
#*******************************************************************************

#Create a list of genes that are related to inflammatory pathways by using one of 3 methods, or all 3
def create_inflamrelated_gene_list(temppathways = {}):
	pathwaygenes = {}

	
	if P.GENE_SELECTION == "senescence": Tinflam_synonyms = {"senescence"}
	else: Tinflam_synonyms = P.inflam_synonyms
	
	#SEARCHWORDS	creating a list of related genes based on searchwords found in pathways desription
	if P.GENE_SELECTION in['searchwords', 'senescence', 'all']:
		countpathways = 0
		for pathway_description in temppathways: #looping over all the pathways (with lists of related genes)
			countpathways += 1
			for sw in Tinflam_synonyms: #looping over all the searchwords
				if sw.lower() in pathway_description.lower(): #if one of the searchwords occurs in the pathway description
					for gene in temppathways[pathway_description]: #looping over all genes in a pathway
						pathwaygenes[gene] = None #adding gene to pathwaygenes dictionary

	#CELL-AGE-SIGNATURES	creating list of realted genes by importing from cellAgen and Signatures
	if P.GENE_SELECTION=='cell-age-signatures' or P.GENE_SELECTION=='all':
		new_genes = pd.concat((CG.cellAge["gene_name"], CG.signatures["gene_symbol"]), axis=0, ignore_index=True, sort=False)
		for gene in new_genes:
			pathwaygenes[gene] = None #adding gene to pathwaygenes dictionary

	#GENES_FROM_PAPERS	creating list of related genes by reading predefined lists (based on genes from papers)
	if P.GENE_SELECTION=='genes-from-papers' or P.GENE_SELECTION=='all':
		for gene in genes_from_papers:
			pathwaygenes[gene] = None #adding gene to pathwaygenes dictionary
	
	try: del pathwaygenes[""]
	except: pass
	#Printing how manny genes where found
	if P.GENE_SELECTION == "searchwords":	woord = P.inflam_synonyms 
	else: 	woord = P.GENE_SELECTION
	print(st.YELLOW,st.BOLD,len(pathwaygenes.keys()),"genes are related to", woord, st.RST)
	
	return pathwaygenes

################################################################################
#*******************************************************************************
#********************************Script*****************************************
#*******************************************************************************
################################################################################

def create_genelist():
	print(P.GENE_SELECTION)

	
	#***********Importing nessesary files*****************************************
	GF.ensure_file('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Human', 'WikiPathways_2019_Human', P.update_files)
	GF.ensure_file('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human', 'KEGG_2019_Human', P.update_files)
	GF.ensure_file('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2017b', 'GO_Biological_Process_2017b', P.update_files)
	GF.ensure_file('https://genomics.senescence.info/cells/cellAge.zip', 'cellAge1.zip', P.update_files)
	GF.ensure_file('https://genomics.senescence.info/cells/cellSignatures.zip', 'signatures1.zip', P.update_files)


	#*********Reading the gene-pathway linking files (if nessesary)***************
	temppathways = {}
	if P.GENE_SELECTION in ["searchwords", "senescence", 'all']: #if P.GENE_SELECTION is one of these 3 options
		wikipath = GF.readfile("Source/WikiPathways_2019_Human", '\t', 0, "dictionary")
		kegg = GF.readfile("Source/KEGG_2019_Human", '\t', 0, "dictionary")
		go = GF.readfile("Source/GO_Biological_Process_2017b", '\t', 0, "dictionary")
		temppathways = {**wikipath, **kegg, **go}

	if P.GENE_SELECTION == "cell-age-signatures" or P.GENE_SELECTION=='all':
		CG.cellAge = GF.readfile('Source/cellAge1.csv', ';', 0, "dataframe")
		CG.signatures = GF.readfile('Source/signatures1.csv', ';', 0, "dataframe")
	
	#**********Create dictionary with selected genes******************************
	#if select_on_genes:
	pathwaygenes = create_inflamrelated_gene_list(temppathways) #retrieving dictionary of genes related to inflammatory pathwyas (gens are keys without values)
	

	#*********Writing the genelist to a file**************************************
	if "Genelists" not in os.listdir('./'+P.experiment_name):
		os.mkdir(os.path.join(os.getcwd(), "./"+P.experiment_name+"/Genelists"))

	f = open("./"+P.experiment_name+"/Genelists/"+P.GENE_SELECTION+"_Genelist.txt", "w+")
	for i in pathwaygenes.keys():
		f.write(i+"\n")
	f.close
	#*****************************************************************************
	
	return(pathwaygenes)



#print"*Closing script Create_genelist.py")








