import pandas as pd
import numpy as np
import ReST
import scipy as sp
import scanpy as sc
from anndata import AnnData as ad
from time import time
import os
import joblib
from multiprocessing import set_start_method

def createReSTobject(ho_fn):
	ho_df = pd.read_csv(ho_fn, index_col=0)
	acs, ars = [int(i.split("x")[0]) for i in ho_df.index],  [int(i.split("x")[1]) for i in ho_df.index]
	obs_df = pd.DataFrame({'array_col': acs, 'array_row': ars}, index=ho_df.index)
	var_df = pd.DataFrame({'gene': ho_df.columns}, index=ho_df.columns)
	ho_rd = ReST.ReST(counts=ho_df, coordinates=obs_df, gene_df=var_df)
	ho_rd.adata.obs['new_idx'] = ho_rd.adata.obs[['array_col', 'array_row']].apply(lambda x: 'x'.join(x.astype(str)), axis=1)
	return ho_rd

def preprocess(adata0, hvg_prop=0.8, n_pcs=30):
	### IMPORTANT TO NOTE:
	### Holdout data is already CPM normalized, therefore, do not normalize it here again.
	adata = adata0.copy()
	adata.var_names_make_unique()
	adata.obs['new_idx'] = adata.obs[['array_col', 'array_row']].apply(lambda x: 'x'.join(x.astype(str)), axis=1)
	## Data normalization 
	sc.pp.filter_genes(adata, min_cells=10)
	sc.pp.normalize_total(adata, inplace=True, target_sum=10**6)
	# LOG2CPM normalization
	sc.pp.log1p(adata, base=2)
	sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=int(hvg_prop * adata.shape[1]))
	## Construct a knn graph based on PC similarity (top 50, nneighbors=100)
	sc.pp.pca(adata, n_comps=n_pcs)
	adata.obsp['raw_weights'] = pd.DataFrame(data=adata.obsm['X_pca'], 
		index=adata.obs_names).T.corr().loc[adata.obs_names, adata.obs_names].values
	return adata

def main(dataDir, methods=["MIST", "MAGIC", "knnSmooth", "spKNN", "mcImpute"], k=5):
	for i in range(k):
		ho_fn = f'{dataDir}/ho_data_{i}.csv'
		ho_rd = createReSTobject(ho_fn)
		ho_rd.adata.layers['CPM'] = ho_rd.adata.X

		# Region extraction
		region_fn = f'{dataDir}/region_inds.csv'
		node_fn = f'{dataDir}/nodes.job'

		raw_fn = f'{dataDir}/raw.csv'
		raw_rd = createReSTobject(raw_fn)
		raw_rd.adata = preprocess(raw_rd.adata)

		if raw_rd.shape[0] < 1000:
			min_region = 10
		else:
			min_region = 40

		raw_rd.extract_regions(min_sim=0.1, min_region = min_region)
		ho_rd.nodes = raw_rd.nodes

		ho_rd.adata.obs['region_ind'] = raw_rd.adata.obs.loc[ho_rd.adata.obs.index,'region_ind']

		outDir = f'{dataDir}/imputed'
		if not os.path.exists(outDir):
			os.mkdir(outDir)
		# ho_rd.extract_regions(min_sim=0.1, gap=0.05, min_region = min_region)
		print(f"{len(set(ho_rd.adata.obs.region_ind)) - 1} regions detected.")

		for method in methods:
			outFn =  f'{outDir}/{method}_{i}.csv'
			if not os.path.exists(outFn):
				ts = time()
				imputed_data = ho_rd.impute(method=method, n_neighbors=4, ncores=10, nExperts=10)
				imputed_data.to_csv(outFn)
				te = time()
				print(f"[{method}] {i}/{k} imputation done in {(te-ts):.2f} seconds.")
			else:
				print(f"[{method}] {i}/{k} imputation exists.")

if __name__ == '__main__':
	import sys
	dataDir = sys.argv[1]
	set_start_method("spawn")
	main(dataDir)
