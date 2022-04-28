#!/usr/bin/env python
"""Script to run holdout experiments
MIST scirpts should be placed at the path '../MIST' to get this work.
"""
import sys
sys.path.append("../MIST/")
import Data
import os
from time import time
import utils
from imputers import *
import pandas as pd


__author__ = "Linhua Wang"
__maintainer__ = "Linhua Wang"
__email__ = "linhuaw@bcm.edu"

if __name__ == '__main__':
	folder = sys.argv[1] # Path to the folder containing hold-out data
	ep = 0.7 # initial epsilon value

	#count, meta = utils.read_ST_data(os.path.join(folder, "norm.csv"))
	for i in range(5):
		# Read and organize hold-out masks
		mask = pd.read_csv(os.path.join(folder, "ho_mask_%d.csv" %i), index_col=0)
		#mask = mask.loc[count.index, count.columns]
		# Removing values to be held out
		fn = os.path.join(folder, "ho_data_%d.csv" %i)
		data = Data.Data(countpath=fn)
		# Run each imputaiton method
		for imputer_name in ["MAGIC", "knnSmooth", "mcImpute", "spKNN", "MIST"]:
			norm_outF = os.path.join(folder, "%s_%d.csv" %(imputer_name, i))
			st_time = time()
			imputer = Imputer(imputer_name, data)
			imputed_data = imputer.fit_transform()
			# Get the selected epsilon value from MIST
			if imputer_name == "MIST":
				ep = imputer.data.epsilon
				norm_outF = os.path.join(folder, "%s_%d_%d.csv" %(imputer_name, ep*100, i))
			imputed_data.to_csv(norm_outF) # save imputed gene expression to path
			print("[%d, %s] elapsed %.1f seconds.." %(i, imputer_name, time() - st_time))

