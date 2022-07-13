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



class st():
	YELLOW = '\033[33m'
	BLUE = '\033[34m'
	MAGENTA = '\033[35m'
	CYAN = '\033[36m'
	BOLD = '\033[1m'
	ITA = '\x1B[3m' #append it with ' \x1B[0m
	RST = '\033[0m'
	RED = '\033[38;2;{};{};{}m'. format(255,69,0)
	ORANGE = '\033[38;2;{};{};{}m'. format(255,165,0)
	GREEN = '\033[38;2;{};{};{}m'. format(144,238,144)



#*******************************************************************************
#**************************Downloading and updating*****************************
#*******************************************************************************
def update_line(line):
	found = False
	date = str(datetime.today().date())

	f = open("Source/Updating_date.txt", "r+")
	d = f.readlines()
	f.seek(0)
	for i in d:
		if (line not in i): f.write(i)
		else: 
			f.write(date + i[10:])
			found = True
	if found == False:
		f.write(date+": "+line+"\n")
	f.truncate()
	f.close
#*******************************************************************************
#Function that downloads the required files, stores it in the Source folder and removes an older version
def download_file(url_adress, filename, check_existing):
	if (check_existing): #checking if file already exists, ifso: not downloading. Force downloading by setting check_existing to False
		#Removing extentions from filename
		TFN = filename #TN = TempFileName
		for ext in [".html", ".zip", ".gz", ".gct", "txt", ".csv"]:
			TFN = TFN.lower().replace(ext, "")

		#Making a list of all the files in 'Source' without extesions
		FILELIST = os.listdir('Source')
		for filenum in range(len(FILELIST)):
			for ext in [".html", ".zip", ".gz", ".gct", "txt", ".csv"]:
				FILELIST[filenum] = FILELIST[filenum].lower().replace(ext, "")

		#Checking if filename occurs in the folder 'Source'
		if TFN in FILELIST:
			print(st.ITA,"\'"+TFN+"\' exists in Source folder, not downloading this file.",st.RST)
			return()
		else: print(st.ITA,"\'"+TFN+"\' Does not yet exist.", st.RST)

	update_line(filename) #Updating file with downloading dates
	print("Downloading <{}>...".format(filename))
	r = requests.get(url_adress, allow_redirects=True)
	open('Source/'+filename, 'wb').write(r.content)
	print("Download finished.")
	#Unzipping if nessesary
	for item in os.listdir("Source"):
		if (item.endswith(".html") or item.endswith(".zip") or item.endswith(".gz")): 
			extens = item.split('.')[-1]
			if (extens != 'html'):
				print("Unpacking .{}...".format(extens))
				if (extens == "gz"): #if a .gz file needs to be extracted
					with gzip.open('Source/'+filename, 'rb') as f_in:
						with open('Source/'+filename[:-3], 'wb') as f_out:
							shutil.copyfileobj(f_in, f_out)
				if (extens == "zip"): #If a zip file needs to be extracted
					with ZipFile('Source/'+filename, 'r') as zipObj:
						zipObj.extractall('Source')
				print("Unpacking done.")
			os.remove("Source/"+item)
			print("Original .{} removed!".format(extens))
#*******************************************************************************

def check_diskspace(minimum, printing):
	total, used, free = shutil.disk_usage("/")
	if (printing):
		print("Total: %d GB" % (total // (2**30)))
		print("Used: %d GB" % (used // (2**30)))
		print("Free: %d GB" % (free // (2**30)))
	return (free // (2**30)) > minimum
#*******************************************************************************
#Updating of downloading source files
def ensure_file(url_adress, filename, update):
	#If Source folder does not exist
	if "Source" not in [f for f in os.listdir('.') if os.path.isdir(f)]:
		os.mkdir(os.path.join(os.getcwd(), "Source"))
		print("New Source-files directory")
		update = True

	#If update_date file does not exist
	if "Updating_date.txt" not in os.listdir('./Source'):
		f = open("Source/Updating_date.txt", "w")
		f.write("Updated			Filename\n")
		f.close()

	#Updating the source files if the variable "update" is set to true
	CFE = True # CFE = Check if File Exists, if the file already exist it will not be downloaded
	if (update): check_fileexist = False
	download_file(url_adress, filename, CFE)
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#******************************************************************************* ->preprocessing
#Reading datafiles and converting to specified datatype
def readfile(filename, separator, heading, datatype):
	print("\nReading", filename ,"...")

	if (datatype == "dictionary"):
		data = {}
		with open(filename, mode='r') as f:
			reader = csv.reader(f, delimiter = separator)
			#data = {row[0] for row in reader}
			for row in reader:
				if len(row[1:]) > 1: data[row[0]] = row[1:] #the key is the first item in the row row[0] and the value is fileld with the other items row[1:]
				else: data[row[0]] = row[1]

	elif (datatype == "dataframe"):
		data = pd.read_csv(filename, encoding = "latin-1", sep=separator, header=heading)
	else:
		print("Error: No valid datatype!")
	print('Reading complete')
	return data

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************




#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#OUde functies
"""
#*******************************************************************************
def check_available_disk_space():
	mb_free = data._get_mb_free_diskspace(video_dir)
	if mb_free < 10:
			self.message_handler.set_message('INFO', 'insufficient space on disk')
			return False
	else:
			return True 

#*******************************************************************************
if(GF.check_diskspace(3, False)): #Checking if there is enough diskspace
else: exit("Not enough diskspace for downloading. Need 3 GB free storage space.")

#*******************************************************************************





"""

