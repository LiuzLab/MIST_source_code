import sys
import pandas as pd
import numpy as np
import ReST
from time import time
import os
from multiprocessing import set_start_method, get_context
from  MIST2 import rankMinImpute as mcImpute

def run_mcImpute_minipatch(folder, data_name):

    for i in range(5):
        out_fn = f'{folder}/{data_name}/imputed/minipatch_mcImpute_{i}.csv'

        if not os.path.exists(f"{folder}/{data_name}/imputed"):
            os.mkdir(f"{folder}/{data_name}/imputed")
        if os.path.exists(out_fn):
            continue

        ho_fn = f'{folder}/{data_name}/ho_data_{i}.csv'
        ho_data = pd.read_csv(ho_fn, index_col=0)
        nspots = ho_data.shape[0]
        spots = ho_data.index.to_numpy()

        new_params = []
        n_patches = 200
        n_samples = int(np.sqrt(nspots))

        for k in range(n_patches):
            np.random.seed(k)
            sampled_spots = np.random.choice(spots, n_samples, replace=False)
            df = ho_data.loc[sampled_spots,:]
            new_params.append([df, f'Patch {k} | {n_patches}'])

        print(f"[{data_name} | fold {i}] Start running mcImpute with {len(new_params)} mini patches.")

        t0 = time()
        with get_context("spawn").Pool(20) as pool:
            patch_results = pool.map(mcImpute, new_params)
        t1 = time()
        elapsed = t1 - t0

        print(f"[{data_name} | fold {i}] Finished in {elapsed:.1f} seconds.")
        patch_results = pd.concat(patch_results)
        patch_results = patch_results.groupby(patch_results.index).mean()
        return_df = ho_data.copy()
        return_df.loc[patch_results.index, patch_results.columns] = patch_results.values      
        return_df.to_csv(out_fn)

if __name__ == '__main__':
    folder = sys.argv[1] # Path to the folder containing hold-out data
    set_start_method("spawn")
    if len(sys.argv) < 3:
        data_names = os.listdir(folder)
        data_names  = [d for d in data_names if ('Human' in d) or ('Mouse' in d) or (d == 'Melanoma')]
        for data_name in data_names:
            run_mcImpute_minipatch(folder, data_name)
    else:
        data_name = sys.argv[2]
        print(data_name)
        run_mcImpute_minipatch(folder, data_name)