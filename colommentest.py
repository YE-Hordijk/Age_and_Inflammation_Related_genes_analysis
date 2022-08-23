#import GeneralFunctions as GF
#GeneralFunctions.py

#Automatic updating
import requests
import gzip
from zipfile import ZipFile
import shutil
import os
import re
from datetime import datetime
import psutil #for checking if enough diskspace

#General packages
import numpy as np
import pandas as pd
import math as mt
import statistics
from numpy import log as ln
from itertools import islice #for selecting dictionary items
import csv
import matplotlib
from matplotlib import pyplot as plt
from collections import OrderedDict


#*******************************************************************************
def readfile(filename, separator, heading, datatype):
	print("\nReading", filename ,"...")

	if (datatype == "dictionary"):
		data = {}
		with open(filename, mode='r') as f:
			reader = csv.reader(f, delimiter = separator)
			#data = {row[0] for row in reader}
			for row in reader:
				if len(row[1:]) > 1: data["tada"] = row[0:] #the key is the first item in the row row[0] and the value is fileld with the other items row[1:]
				else: data[row[0]] = row[1]

	elif (datatype == "dataframe"):
		data = pd.read_csv(filename, encoding = "latin-1", sep=separator, header=heading)
	else:
		print("Error: No valid datatype!")
	print('Reading complete')
	return data

banaan = readfile("RNAseq_colommen.txt", "\t", 0, "dictionary")
#*******************************************************************************

print(banaan)
print(len(banaan['tada']))
exit()
banaan = list(banaan)
print(banaan)
