class P:
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
	MODEL = "RandomForest"#"DecisionTree"#"RandomForest" #"DecisionTree" #"Support Vector Machine"
	METHOD = "Classification" #"Regression"
	WriteLaTeXTableforML = True
