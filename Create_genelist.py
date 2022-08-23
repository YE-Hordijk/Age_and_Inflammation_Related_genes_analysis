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
#InflamBiologicAge = ["FCGR1A", "GRAMD1C", "PNP", "CD274", "DPEP2", "TSTA3", "KREMEN1", "SLC26A8", "ALAS2", "NDEL1", "SAMD9L", "RPS6KA5", "SIAH1", "GBP7", "SIAH2", "KCNH7", "DAPK2", "SLC22A4", "FAM20B", "UBE2J1", "SLC25A37", "HDGF", "METTL9", "GBP1", "CD2AP", "MCOLN3", "CTNNAL1", "GNAS", "STOM", "UPP1", "TIFA", "IFI27", "ACSL6", "VNN3", "ATL2", "RAP2A", "GBP5", "ABCG2", "AHSP", "ALDH1A1", "HPS1", "TNIP1", "DAAM2", "FLT3", "ZAK", "IGF2BP2", "DOCK4", "RAP1GAP", "SELENBP1", "PBX1", "TMEM49", "ANK1", "ANKRD22", "C20orf11", "SPTB", "ATP1B2", "THRA", "CARM1", "GPR146", "GRINA", "GMPR", "EPB42", "ARSG", "NFKBIZ", "IGF1R", "ITGA1", "NRG1", "PTGER2", "C18orf10", "RAD23A", "DNAJC6", "RPGRIP1", "TARS", "TMEM71", "PKD1L3", "CCNI", "TBC1D22B", "PARP9", "NEDD4L", "SLC4A1", "SLC38A5", "DAPP1", "HEATR5A", "CD1C", "SQLE", "MAF1", "FBXL20", "CACNA1E", "MPZL1", "ZER1", "ARL6IP5", "SLC4A10", "TRIM22", "APAF1", "PDE3B", "ZDHHC2", "CYB5R3", "C9orf156", "ABHD5", "C17orf39", "CASS4", "STAT1", "SGIP1", "ANKH", "CDK19", "PSME3", "CHI3L1", "UBE2L6", "IKBIP", "PARP14", "ANGPT1", "TP53INP2", "CISD2", "GSTK1", "AP2M1", "CNPY3", "SLC14A1", "PLCL2", "NCOA7", "TRIM58", "FLVCR2", "SLC2A1", "DCK", "HBM", "JARID2", "ITCH", "KLC3", "TNFAIP6", "SC4MOL", "CRLF3", "EIF2AK1", "GOLPH3L", "GLRX5", "SPAG9", "STXBP1", "CPPED1", "IFIH1", "IL1RN", "PLSCR4", "MCOLN1", "BANP", "C9orf84", "CUL4A", "CD300LF", "VTI1B", "GPRIN3", "LAP3", "DPCD", "RNF10", "MYL4", "MAPK14", "CDC34", "PLSCR1", "SPATA6", "JAZF1", "SLC1A5", "AP2S1", "EGLN2", "ZNF701", "GYPC", "PA2G4", "HK1", "DDX60", "SLC25A37", "ABHD2", "WSB2", "BSG", "C10orf54", "CEACAM4", "ADIPOR1", "P2RY14", "RNF123", "AIM2", "GCLC", "ZDHHC3", "OSBP2", "DEF8", "SLC43A2", "FIS1", "TM6SF1", "WNK1", "TMCC3", "SERPING1", "TMOD1", "MBD5", "ARRDC3", "XK", "SMURF2", "SNX3", "CCNG2", "ACP1", "FAM53B", "DSC2", "AGPAT3", "TPST1", "NTAN1", "NFE2", "SLC6A9", "POLR1D", "EPB49", "RNF175", "HBD", "EPSTI1", "ANXA3", "HAL", "PDCD1LG2", "TGM2", "NFIA", "C5orf32", "PINK1", "CIZ1", "PKN2", "SIN3B", "PRR11", "POLR2B", "IL1B", "VSIG10", "MXI1", "RSAD2", "BLVRB", "DUSP3", "SLC7A5", "RRM2B", "TAB3", "TMEM144", "CD93", "DPM1", "PI3", "KLRB1", "DPM2", "ARF1", "PLEK2", "WDR26", "CCDC125", "TOPORS", "GPR84", "CMTM1", "CA1", "CDYL", "NFIX", "CCNA2", "H1F0", "EIF1B", "EMR3", "GABARAP", "AP2A1", "DCAF12", "GBP4", "KCNJ15", "HSDL2", "POP7", "HEPACAM2", "UBE2F", "RGS2", "DOCK5", "UBE2O", "MARCH3", "GUK1", "PPP2R5B", "IGF1R", "CSDA", "TFDP1", "MDM2", "LUC7L", "PPP2R2D", "SP100", "MAP4K5", "GNA12", "MGEA5", "DHX29", "PCGF5", "DNAJB2", "C15orf54", "ZNF586", "RBM38", "CLK1", "AFF1", "ASCC2", "KLHL22", "UBXN6", "C2orf24", "IFRD1", "HP", "SLC45A4", "PSME4", "C1orf183", "UROD", "NPRL3", "IFITM1", "ZNF143", "C15orf39", "ERGIC1", "ZFAND5", "TRAK2", "PPP4R1", "ARNTL2", "ZFYVE27", "CMPK2", "CDH1", "LDLR", "PGM2L1", "CA2", "UBOX5", "YPEL3", "RANBP10", "ARAP2", "SLC6A8", "DYRK3", "AMPD2", "ARHGEF12", "AKT3", "ZNF776", "NARF", "TREM1", "LGALS3", "GPRASP1", "EAPP", "PIP4K2A", "CTSB", "CARHSP1", "SLC25A39", "WASF1", "NCAPH2", "GBP6", "SIRT1", "JAK3", "LSP1", "PROK2", "KIAA1257", "C10orf46", "TSPAN17", "BAT3", "TNS1", "GABARAPL1", "UBN1", "USP46", "KCNE3", "PPFIA1", "RALGDS", "TSPAN5", "HERPUD1", "CASP5", "YTHDC1", "RBMS1", "HSD17B11", "KIAA0513", "IFI16", "ANKRA2", "PIP5K1B", "OAS1", "DCAF10", "WDR89", "TFPI", "PHF17", "PRPF18", "ERN1", "LEMD3", "GCLM", "DEDD", "RBM12B", "LNX2", "KEL", "PARP12", "NLRP1", "PHF20", "C13orf18", "TBC1D24", "EVI5", "ZNF721", "SLC39A8", "MICALCL", "C20orf108", "SEC14L4", "SRRD", "TLR5", "C20orf197", "CASP7", "FZD5", "IGF2R", "HAGH", "EEPD1", "ENTPD5", "MYD88", "SP140", "ASGR2", "OPTN", "FBXO6", "GCH1", "RFXANK", "CMAH", "DCAF11", "TNFSF14", "C22orf25", "TMEM156", "FLJ38717", "KCNJ2", "FAM198B", "RELL1", "ZNF44", "BIRC3", "RAB37", "DYNLL2", "SFRP2", "DCTN4", "C17orf49", "SULF2", "FAM174A", "AMD1", "ASCC3", "PRDX2", "SLPI", "CKAP2L", "PAPOLG", "RASSF2", "TUT1", "SH3BGRL3", "GSK3A", "ZNF160", "KLC1", "RAPGEF2", "PGAP1", "STAU1", "C3AR1", "EIF3B", "FRMD3", "AKT1", "ICAM3", "DDX58", "AMFR", "LMBRD1", "ANXA4", "C18orf8", "CCNDBP1", "MTHFS", "EPS8", "ARSB", "ZNF658", "C2orf88", "SEC14L1", "EZH2", "ENST00000331346", "PPDPF", "CYB5R4", "DUSP11", "XPO7", "XPO7", "PTTG1", "ST6GALNAC4", "PRKCZ"]


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
							'IL6','ANKRD20A9P','LINC02227','FOXO3A','RAD50','MC2R',
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


#print(len(Asso_Inflam_Age))
#print(len(set(Asso_Inflam_Age)&set(sign_pathway)))
#print(set(Asso_Inflam_Age)&set(sign_pathway))
#exit()
genes_from_papers = sign_pathway + Age_asso_genes + GWAS_Loci + Asso_Inflam_Age
#*******************************************************************************
#*****************************Functions*****************************************
#*******************************************************************************

#Create a list of genes that are related to inflammatory pathways by using one of 3 methods, or all 3
def create_inflamrelated_gene_list(temppathways = {}):
	pathwaygenes = {}

	
	if P.GENE_SELECTION == "senescence": Tinflam_synonyms = {"senescence"}
	elif P.GENE_SELECTION == "searchwords": Tinflam_synonyms = P.inflam_synonyms
	elif P.GENE_SELECTION == "all": Tinflam_synonyms = set(list(P.inflam_synonyms)+["senescence"])
	
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









