#!/usr/bin/env python
"""The main script to run MIST
"""

from neighbors import *
import pandas as pd
import numpy as np
from multiprocessing import Pool
from time import time
import Data

__author__ = "Linhua Wang"
__license__ = "GNU General Public License v3.0"
__maintainer__ = "Linhua Wang"
__email__ = "linhuaw@bcm.edu"

def MIST(rdata, nExperts=10, ncores=1, verbose=1):
	"""
	The function that takes ST data object and return the imputed gene expression and \
	region membership for all ST spots.

	Parameters
	----------
	"""
	start_time = time()

	## Get the region assignment
	member_df = rdata.adata.obs
	isolated_spots = member_df.loc[member_df.region_ind == 'isolated'].index.tolist()
	core_regions = set(member_df.region_ind)
	core_regions.remove("isolated")
	core_regions = list(core_regions)

	## Construct a data frame for denoising
	data = pd.DataFrame(data=rdata.adata.layers['CPM'].toarray().astype(int),
		index = rdata.adata.obs.new_idx, columns = rdata.adata.var_names)

	spots = data.index.tolist()
	imputed_whole = rankMinImpute(data)
	
	# Indices of non-zero gene expression values in the matrix
	known_idx = np.where(data.values)
	# Make template to save results
	imputed = data.copy()

	imputed.loc[isolated_spots,:] = imputed_whole.loc[isolated_spots, :]
	nspots = 0
	# Impute each CC
	for i1 in range(len(core_regions)):
		region = core_regions[i1]
		region_spots = member_df.loc[member_df.region_ind == region, :].index.tolist()
		nspots += len(region_spots)
		print(f"[Region] {i1} / {len(core_regions)} | [Spots] {nspots} / {member_df.shape[0] - len(isolated_spots)}.")
		# QC
		other_spots = [s for s in spots if s not in region_spots]
		m = len(region_spots)
		s = np.min([int(len(other_spots) / 10), m]) # sampling is based on current cluster size
		# Matrix to save multiple runs of rank minimization
		values = np.zeros(shape=(m, data.shape[1], nExperts))
			
		### Get parameters for rank minimization
		# Parallel computing
		ensemble_inputs = []
		for i2 in range(nExperts):
			np.random.seed(i2)
			sampled_spot = list(np.random.choice(other_spots, s, replace=False))
			ri_spots = region_spots + sampled_spot
			ensemble_inputs.append(data.loc[ri_spots,:])
		p = Pool(ncores)
		ensemble_outputs = p.map(rankMinImpute, ensemble_inputs)
		p.close()

		# Reorganize and average multiple runs of rank minimization results for CC
		for i2 in range(nExperts):
			values[:,:,i2] = ensemble_outputs[i2].loc[region_spots,:].values
		imputed.loc[region_spots,:] = np.mean(values, axis=2) 

	## Make sure original non-zero values are not perturbed
	imputed.values[known_idx] = data.values[known_idx] 
	return imputed

def rankMinImpute(data):
	"""
	The function that impute each LCN. Modified based on the\
	 MATLAB code from Mongia, A., et. al., (2019).

	Input: expression matrix with rows as samples and columns as genes
	Output: imputed expression matrix of the same size.
	"""
	t1 = time()
	D = np.ravel(data.values) # flattened count matrix
	idx = np.where(D) # nonzero indices of D
	y = D[idx] # nonzero observed values
	n = np.prod(data.shape)
	err= 1E-12
	x_initial = np.zeros(np.prod(data.shape))
	tol= 1E-4
	decfac = 0.7
	alpha = 1.1
	x = x_initial
	lamIni = decfac * np.max(np.absolute(y))
	lam = lamIni

	f1 = np.linalg.norm(y - x[idx], 2) + np.linalg.norm(x, 1) * lam

	while lam > lamIni * tol:
		for i in range(20):
			f0 = f1
			z = np.zeros(n)
			z[idx] = y - x[idx]
			b = x + (1/alpha) * z
			B = np.reshape(b, data.shape)
			U, S, V = np.linalg.svd(B,full_matrices=False)
			s = softThreshold(S, lam/(2*alpha))
			X = U @ np.diag(s) @ V
			X[X<0] = 0
			x = X.ravel()
			f1 = np.linalg.norm(y - x[idx], 2) + np.linalg.norm(x, 1) * lam
			e1 = np.linalg.norm(f1-f0)/np.linalg.norm(f1+f0)
			if e1 < tol:
				break
		e2 = np.linalg.norm(y-x[idx])
		if e2 < err:
			break
		lam = decfac*lam
	imputed = pd.DataFrame(data=X, index=data.index, columns=data.columns)
	assert not imputed.isna().any().any()
	t2 = time()
	return imputed

def softThreshold(s, l):
	"""
	Helper function to compute the L1 gradient for function rankMinImpute.
	"""
	return np.multiply(np.sign(s), np.absolute(s - l))