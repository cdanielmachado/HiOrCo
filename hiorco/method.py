import os
import numpy as np
import pandas as pd
from multiprocessing import Pool
from numpy.random import choice, seed
from datetime import datetime
from scipy.stats import binom
from statsmodels.stats.multitest import fdrcorrection
from functools import partial


def evaluate(rows, data):
    x = data.loc[rows,:]
    total = x.min(axis=0).sum()
    n_samples = data.shape[1]
    probability = (x.sum(axis=1) / n_samples).prod()
    expected = n_samples * probability
    p_value = 1 - binom.cdf(total, n_samples, probability)
    return total, expected, p_value


def compute(data, k=2, n=10, p=100, part_size=100, abundance_cutoff=0.001, min_samples=10, 
    fold_cutoff=2, fdr_cutoff=0.05, output_folder=".", random_selection=True, parallel=True):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    print(f"Inital dataset: {data.shape[0]} species, {data.shape[1]} samples.", flush=True)

    data = data / data.sum(axis=0) 
    data = (data > abundance_cutoff).astype(int)
    data = data.loc[data.sum(axis=1) >= min_samples, data.sum(axis=0) > 1]

    print(f"Processed data: {data.shape[0]} species, {data.shape[1]} samples.", flush=True)

    combinations = []
    species = list(data.index.values)
    pop_sample = [(x,) for x in species]

    eval_func = partial(evaluate, data=data)

    if parallel:
        pool = Pool()
        map_func = pool.map
    else:
        map_func = map

    for i in range(2, k+1):

        print(f"\nComputing co-occurring species sets of size {i}...", flush=True)

        # combine all tuples with one more species (without repetition)
        combinations = (tuple(sorted(set(x)|{y})) for x in pop_sample
                        for y in species if y not in set(x))
        combinations = list(set(combinations))

        # compute scores for all combinations
        scores = map_func(eval_func, combinations)
        scores = pd.DataFrame(scores, columns=["total", "expected", "p_value"])
        scores["community"] = combinations
        
        # compute FDR
        scores["significant"], _ = fdrcorrection(scores["p_value"], alpha=fdr_cutoff)
        scores["final"] = np.sqrt(scores["total"])

        #filter FDR and min occurrence threshold
        selected = scores.query("significant and total >= @min_samples and total >= @fold_cutoff * expected")

        # if no results, stop at this point
        if len(selected) == 0:
            print("No significant hits found.", flush=True)
            break
        else:
            print(f"Selected: {len(selected)} sets.", flush=True)

        top_size = min(n, len(selected))
        pop_size = min(p, len(selected))

        # get the top results for this iteration
        selected = selected.sort_values("final", ascending=False)
        top = selected["community"].head(top_size)
        
        # select population for next iteration
        if random_selection:
            weights = selected["final"] / selected["final"].sum()
            pop_sample = choice(selected["community"], size=pop_size, replace=False, p=weights)
        else:
            pop_sample = selected["community"].head(pop_size)
        
        # save communities to splited files

        save_output(top, output_folder, i, part_size)


def save_output(results, output_folder, i, part_size):

    comms = []
    for j, org_ids in enumerate(results):
        comm_id = 'comm_{}_{}'.format(i, j)
        entries = [(comm_id, org_id) for org_id in org_ids]
        comms.extend(entries)

        k, l = (j+1) % part_size, j // part_size
        if k == 0 or j == len(results) - 1:
            df_k = pd.DataFrame(comms, columns=['community', 'species'])
            filename = f"{output_folder}/size_{i}_part_{l+1}.tsv"
            df_k.to_csv(filename, sep='\t', index=False, header=False)
            comms = []
