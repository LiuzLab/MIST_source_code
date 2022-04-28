#!/usr/bin/env python
"""Script to evaluate all imputation methods for holdout experiments.
"""
import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import r2_score, roc_auc_score, recall_score, precision_score
import sys
import utils
from  tqdm import trange
from time import time
from os.path import join
import Data
from neighbors import spatialCCs
from multiprocessing import Pool
from glob import glob

__author__ = "Linhua Wang"
__maintainer__ = "Linhua Wang"
__email__ = "linhuaw@bcm.edu"

def evalSpot(ori, mask, meta, model_data, model_name):
	"""Method to evaluate spot level performance for holdout experiments

	Parameters:
	----------
	ori: data frame, spot by gene,
		original gene expression matrix, with ground truth of holdout values
	mask: data frame, spot by gene,
		1 indicates a holdout event, 0 otherwise
	meta: data frame, spot by (x, y) coordinates
	model_data: data frame, spot by gene,
				imputed gene expression
	model_name: str, imputation method name
	"""
	spots = mask.index[(mask == 1).any(axis=1)] # Find spots that have at least one holdout
	meta = meta.loc[spots,:]
	rmses, pccs_all, snrs = [], [], []
	for i in trange(len(spots)):
		spot = spots[i] # evaluate every spot
		genes = mask.columns[mask.loc[spot,:] == 1].tolist() # get holdout genes for this spot

		try:
			tru = ori.loc[spot, genes].to_numpy() # get ground truth of holdout values
			imp = model_data.loc[spot, genes].to_numpy() # get estimated holdout values
			rmses.append(np.sqrt(np.mean(np.square(imp - tru)))) # calculate rmse
			pccs_all.append(pearsonr(tru,imp)[0]) # calculate PCC
			snrs.append(np.log2((np.sum(imp) +1) /
				(1+np.sum(np.absolute(tru-imp))))) # calculate SNR
		except:
			# scipy PCC fails if there is only 1 unique value or less than 3 values,
			# here we assign nan for such cases 
			rmses.append(np.nan)
			pccs_all.append(np.nan)
			snrs.append(np.nan)
	# organize results to output
	spot_perf = pd.DataFrame({"spot":spots, "x":meta.iloc[:,0],
				"y":meta.iloc[:,1],"rmse":rmses,
				"pcc": pccs_all,"snr": snrs,
				"model":model_name})
	return spot_perf

## Evaluate gene level performance for holdout test
def evalGene(ori, mask, ho, meta, model_data, model_name):
	"""Method to evaluate gene level performance for holdout experiments

	Parameters:
	----------
	ori: data frame, spot by gene,
		original gene expression matrix, with ground truth of holdout values
	mask: data frame, spot by gene,
		1 indicates a holdout event, 0 otherwise
	ho: data frame, spot by gene,
		gene expression matrix after removing holdout values
	meta: data frame, spot by (x, y) coordinates
	model_data: data frame, spot by gene,
				imputed gene expression
	model_name: str, imputation method name
	"""
	genes = mask.columns[(mask == 1).any(axis=0)] # Find spots that have at least one holdout

	rmses, pccs_all, snrs, mapes = [], [], [], []
	mrs = []
	genes_ho = []

	for i in trange(len(genes)):
		gene = genes[i] # evaluate every spot
		spots = mask.index[mask.loc[:, gene] == 1].tolist()# get holdout spots for this gene

		if len(spots) == 0: # skip if no spots were held out
			continue

		genes_ho.append(gene)
		mr = (ho.loc[:, gene] == 0).sum()/float(ho.shape[0]) # get non-zero proportion of hold out values
		mrs.append(mr)
		tru = ori.loc[spots, gene].to_numpy() # get ground truth of holdout values
		imp = model_data.loc[spots, gene].to_numpy() # get estimated holdout values
		rmses.append(np.sqrt(np.mean(np.square(imp-tru)))) # calculate rmse
		pccs_all.append(pearsonr(tru,imp)[0]) # calculate PCC

		snrs.append(np.log2((np.sum(imp) +1) /
			(1+np.sum(np.absolute(tru-imp))))) # calculate SNR
		mapes.append(np.mean(np.divide(np.absolute(imp - tru), tru))) # calculate MAPE
	# integrate results for output
	gene_perf = pd.DataFrame({"gene": genes_ho,
				"rmse":rmses,
				"pcc": pccs_all,
				"snr": snrs,
				"mape":mapes,
				"model":model_name,
				"mr": mrs})

	return gene_perf

def evalSlide(ori, mask, ho, model_data, model_name, spots=None):
	"""Method to evaluate whole-slide level performance for holdout experiments

	Parameters:
	----------
	ori: data frame, spot by gene,
		original gene expression matrix, with ground truth of holdout values
	mask: data frame, spot by gene,
		1 indicates a holdout event, 0 otherwise
	ho: data frame, spot by gene,
		gene expression matrix after removing holdout values
	meta: data frame, spot by (x, y) coordinates
	model_data: data frame, spot by gene,
				imputed gene expression
	model_name: str, imputation method name
	spots: evaluate performance on certain spots. Useful when evaluating on LCN spots
			only.
	"""
	if spots != None:
		# LCN spots data
		observed = ori.loc[spots,:]
		ho_mask = mask.loc[spots,:]
		ho_data = ho.loc[spots,:]
		imputed_data = model_data.loc[spots,:]
	else:
		observed = ori.copy()
		ho_data = ho.copy()
		imputed_data = model_data.copy()
		ho_mask = mask.copy()
	# Get holdout indices
	M = np.ravel(ho_mask) 
	inds = np.where(M)
	# Get ground truth and estimated holdout values
	tru = np.ravel(observed.values)[inds]
	imp = np.ravel(imputed_data.values)[inds]
	# Calculate SNR, RMSE, MAPE and PCC
	snr = np.log2(np.sum(imp) / np.sum(np.absolute(tru-imp)))
	rmse = np.sqrt(np.mean(np.square(imp - tru)))
	mape = np.mean(np.divide(np.absolute(imp - tru), tru))
	try:
		pcc = pearsonr(tru, imp)[0]
	except:
		pcc = None
	# Calculate non-zero proportion before and after imputation
	MR1 = float((ho_data == 0).sum().sum()) / np.prod(ho_data.shape)
	MR2 = float((imputed_data == 0).sum().sum()) / np.prod(imputed_data.shape)
	# Integrate results for output
	perf_df = pd.DataFrame(data=[[rmse, mape, snr, pcc, model_name, MR1, MR2, MR1-MR2]],
			 columns= ['RMSE', 'MAPE', 'SNR', 'PCC', 'ModelName', 'hoMR', 'impMR', 'redMR'])
	print(perf_df)
	return perf_df

def evalAll(data_folder, model_names, slideonly=True, cvFold=5):
	"""Method to evaluate all levels of holdout experiments 

	Parameters:
	----------
	data_folder: path to the folder storing holdout related files.
	model_names: names of imputation methods to be compared.
	slideonly: True -- evaluate slide level only. False -- evaluate all levels.
	cvFold: int, the CV fold of holdout experiments.
	"""

	model_perf_dfs = []
	spot_perf_dfs = []
	gene_perf_dfs = []
	for seed in range(cvFold):
		# Performance for each fold
		st = time()
		mask = pd.read_csv("%s/ho_mask_%d.csv" %(data_folder, seed), index_col=0)
		genes = mask.columns.tolist()
		ori, meta = utils.read_ST_data("%s/norm.csv" %data_folder) # read ground truth
		t1 = time()
		print("[Fold %d] Ground truth data loading elapsed %.1f seconds." %(seed, t1 - st))
		ori = ori.loc[mask.index, genes]
		ho = ori.copy()
		ho[mask==1] = 0
		ori = np.log2(ori + 1) #CPM to logCPM (log scale is variance stablized)
		for model_name in model_names:
			# Read imputed gene expression by every method
			t2 = time()
			fn = os.path.join(data_folder, "imputed/%s_%d.csv" %(model_name, seed))
			# if model_name == "MIST":
			# 	fn = glob(os.path.join(data_folder, "%s_*_%d.csv" %(model_name, seed)))[0]
			model_data = pd.read_csv(fn, index_col=0)

			## SAVER imputed by R language, - became . in the header
			if model_name == "SAVER":
				model_data.columns = ho.columns.tolist()
				model_data.index = ho.index.tolist()

			model_data = model_data.loc[ori.index, genes]
			model_data = np.log2(model_data + 1) # log transform
			t3 = time()
			print("[Fold %d, %s] Model data loading elapsed %.1f seconds." %(seed, model_name, t3-t2))
			model_perf_df = evalSlide(ori, mask, ho, model_data, model_name) # evaluate slide-level performance
			t4 = time()
			print("[Fold %d, %s] Slide-level performance evaluation elapsed %.1f seconds." %(seed, model_name, t4-t3))
			model_perf_df['cvFold'] = seed
			model_perf_dfs.append(model_perf_df)

#			spot_perf_df = evalSpot(ori, mask, meta, model_data, model_name)# evaluate spot-level performance
			t5 = time()
			print("[Fold %d, %s] Spot-level  performance evaluation elapsed %.1f seconds." %(seed, model_name, t5-t4))
			gene_perf_df = evalGene(ori, mask, ho,  meta, model_data, model_name)
			t6 = time()
			print("[Fold %d, %s] Gene-level  performance evaluation elapsed %.1f seconds." %(seed, model_name, t6-t5))
# #			spot_perf_df['cvFold'] = seed
#			spot_perf_dfs.append(spot_perf_df)
			gene_perf_df['cvFold'] = seed
			gene_perf_dfs.append(gene_perf_df)
	# Integrate and save all results
	model_perf_dfs = pd.concat(model_perf_dfs)
	print(model_perf_dfs)
	model_perf_dfs.to_csv(os.path.join(data_folder, "performance/mcImpute_slide_level_results.csv"))
# #	spot_perf_dfs = pd.concat(spot_perf_dfs)
	gene_perf_dfs = pd.concat(gene_perf_dfs)
# #	spot_perf_dfs.to_csv(os.path.join(data_folder, "performance/spot_level_results.csv"))
	gene_perf_dfs.to_csv(os.path.join(data_folder, "performance/mcImpute_gene_level_results.csv"))


def eval_LCN_runner(param):
	"""Method to get the holdout performance for only LCN spots
	"""
	projDir, dn, fd = param
	models = ["MSIT", "mcImpute", "MAGIC", "knnSmooth","spKNN", "SAVER"]
	LCN_spots = LCN_captured_spots(join(projDir, dn), fd) # Get LCN spots

	with open(join(join(projDir, dn), "LCN_spots_%d.csv" %fd), "w") as f:
		f.write(",".join(LCN_spots))

	# Start to evaluate
	ho = pd.read_csv(join(join(projDir, dn), "ho_data_%d.csv" %fd), index_col=0)
	mask = pd.read_csv(join(join(projDir, dn), "ho_mask_%d.csv" %fd), index_col=0)
	observed = pd.read_csv(join(join(projDir, dn), "norm.csv"), index_col=0)
	observed = np.log2(observed + 1)
	results = []
	for model in models:
		model_df = pd.read_csv(join(join(projDir, dn), "%s_%d.csv" %(model, fd)), index_col=0)
		model_df = np.log2(model_df + 1)
		model_perf = evalSlide(observed, mask, ho, model_df, model, spots=LCN_spots)
		results.append(model_perf)
	results = pd.concat(results)
	results["cvFold"] = fd
	results["data"] = dn
	return results

def eval_LCNs(projDir):
	"""Method to get the holdout performance for only LCN spots
		This function was used to get the LCN level improvement of
		MIST as compared to other methods.
	"""

	data_names = ["MouseWT", "MouseAD", "Melanoma", "Prostate"]
	folds = range(5)

	# Make parameters for parallel computing
	params = []
	for dn in data_names:
		for fd in folds:
			params.append([projDir, dn, fd])
	# Start parallel computing
	p = Pool(10)
	model_perfs = p.map(eval_LCN_runner, params)
	# Integrate results and save results
	model_perfs = pd.concat(model_perfs)
	model_perfs.to_csv(join(projDir, "LCNspots_slide_level_results.csv"))


def LCN_captured_spots(folder, fd):
	"""Method to get the LCN captured spots for this fold of holdout"""

	cp = join(folder, "ho_data_%d.csv" %fd)
	ho_data_obj = Data.Data(countpath=cp)
	mist_fn = glob(join(folder, "MIST_*_%d.csv" %fd))[0]
	# Get the epsilon value selected for this fold
	ep = float(mist_fn.split("/")[-1].split('_')[1]) / 100
	ho_data_obj.update_ep(ep)
	# Get CCs
	ccs = spatialCCs(ho_data_obj.nodes,ho_data_obj.cormat,
		ho_data_obj.epsilon, merge=0)
	# Get LCN captured CCs
	spots = []
	for cc in ccs:
		if len(cc) > 5:
			for node in cc:
				spots.append(node.name)
	return spots

def main(data_folder, cvFold):
	"""Main method to evaluate all-level performance"""
	perf_folder = os.path.join(data_folder, "performance")
	if not os.path.exists(perf_folder):
		os.mkdir(perf_folder)
	# model_names = ["MIST", "mcImpute","MAGIC", "spKNN", "knnSmooth"]
	model_names = ['minipatch_mcImpute']
	evalAll(data_folder, model_names, cvFold=cvFold)


if __name__ == "__main__":
	dataDir = sys.argv[1]
	data_names = os.listdir(dataDir)
	data_names = [d for d in data_names if ('Human' in d) or ('Mouse' in d) or ('Melanoma' == d)]
	cvFold=5
	for i, data_name in enumerate(data_names):
		wkd = os.path.abspath(f"{dataDir}/{data_name}")
		main(wkd, cvFold)
		print(f'[{data_name}] {i} / {len(data_names)}')