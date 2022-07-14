#ExtractingPCs.py


#General packages
import GeneralFunctions as GF
from GeneralFunctions import st
import Parameters as P

import numpy as np
from numpy import log as ln
import os
import pandas as pd
import math as mt
import statistics
from itertools import islice #for selecting dictionary items
import matplotlib
from matplotlib import pyplot as plt
from collections import OrderedDict
#import sys


if select_on_genes: help_name = GENE_SELECTION #variable for nameing files and folders
else: help_name = "NoGeneSelection"
if use_middle_age:	help_name += "_With-MiddleAge"
else: 								help_name += "_No-MiddleAge"

print(help_name)
