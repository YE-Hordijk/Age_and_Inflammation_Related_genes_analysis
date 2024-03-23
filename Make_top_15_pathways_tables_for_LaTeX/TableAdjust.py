#TableAdjustment.py
import requests
import gzip
from zipfile import ZipFile
import shutil
import os
import re
from datetime import datetime
import psutil #for checking if enough diskspace
import operator #for sorting



import numpy as np
from numpy import log as ln
import pandas as pd
import math as mt
import statistics
from itertools import islice #for selecting dictionary items
#import sys
import matplotlib.pyplot as plt



import sys
import pandas as pd
import os

name = input("Give the name of the file: ")
name = "KEGG_2021_Human_table.txt"
name = "BioPlanet_2019_table.txt"
name = "SVM_BioPlanet_2019_table.txt"
name = "MD-BioPlanet_2019_table.txt"
File = pd.read_csv(name, sep='\t')
#File = pd.read_csv('KEGG_2021_Human_table.txt/KEGG_2021_Human_table.csv')

print(File)
ding = {}
File = File.loc[:, ["Term", "Overlap", "Adjusted P-value", "Genes"]]

File= File.sort_values(by=['Adjusted P-value'], ascending=True)
print(File)


for i in File.loc[:,:]:
	#input(i)
	if i == "Genes":
		for k in File.loc[:, i]:
			if "Genes" not in ding: 
				ding["Genes"] = []
				ding["Genes"].append(str(k.replace(";", "; ")))
			else: ding["Genes"].append((k.replace(";", "; ")))
	
	elif i == "Adjusted P-value": 
		for k in File.loc[:, i]:
			#input(k)
			if "Adjusted P-value" not in ding: 
				ding["Adjusted P-value"] = []
			ding["Adjusted P-value"].append(f'{k:.2e}')

			
	else: ding[i] = [k for k in File.loc[:, i]]

ding = pd.DataFrame.from_dict(ding, orient='columns')# columns=['A', 'B', 'C', 'D'])



ding.index += 1 
print(ding)
ding.to_csv(name[:-4]+"NEW.csv", index=True, sep=",")  




