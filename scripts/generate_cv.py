#!/usr/bin/env python
"""Generate data for hold-out experiments
MIST scirpts should be placed at the path '../MIST' to get this work.
"""

import pandas as pd
import numpy as np
import sys
sys.path.append("../MIST/")
import utils
from time import time
from os.path import exists, join
from os import mkdir

__author__ = "Linhua Wang"
__maintainer__ = "Linhua Wang"
__email__ = "linhuaw@bcm.edu"

def main(data_folder, norm="cpm", kFold=5):
	"""Method to generate 5-fold CV for hold-out experiments
	
	Parameters:
	---------
	data_folder: path to the folder containing the raw data,
				generated hold-out data will be saved here.
	norm: normalization method, cpm is used for the manuscript.
	kFold: #CV, 5-fold CV is used in the manuscript.

	"""
	norm_fn = join(data_folder, "norm.csv")
	raw_fn = join(data_folder, "raw.csv")

	raw, meta = utils.read_ST_data(raw_fn)
	normed, libsize =  utils.data_norm(raw, norm)
	if norm in ["median", "logMed", "logCPM"]:
		normed = normed.round(2)

	### Write ground truth
	normed.to_csv(join(data_folder,"norm.csv"))
	### Write the library size for each spot
	libsize.to_csv(join(data_folder, "libsize.csv"))

	### Use high-coverage genes for holdout experiments
	normFilt = utils.filterGene_sparsity(normed,0.5)
	genes = normFilt.columns.tolist()
	print("Holdouut %s genes." %(normFilt.shape[1]))
	# generate 5 fold cross validation datasets
	ho_dsets, ho_masks = utils.generate_cv_masks(normFilt, genes, kFold)

	for fd in range(kFold):
		ho_data, ho_mask = ho_dsets[fd], ho_masks[fd]
		ho_mask.to_csv(data_folder + "/ho_mask_%d.csv" %fd)
		ho_data.to_csv(data_folder + "/ho_data_%d.csv" %fd)

if __name__ == "__main__":
	folder = sys.argv[1]
	norm = 'cpm'
	assert norm in ["cpm", "logCPM", "median", "logMed"]
	if "slideseq" in folder:
		main(folder, norm=norm)
	else:
		main(folder, norm=norm)
