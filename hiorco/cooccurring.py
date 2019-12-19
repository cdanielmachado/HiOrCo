#!/usr/bin/env python

#SBATCH -t 14-0:00:00
#SBATCH --job-name longcooc
#SBATCH --mem=256000
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --array=1-12
#SBATCH -o ../log/out_%A_%a.txt
#SBATCH -e ../log/err_%A_%a.txt


from __future__ import division

import os
import numpy as np
import pandas as pd
from multiprocessing import Pool
from numpy.random import choice, seed
from datetime import datetime
from scipy.stats import binom
from statsmodels.stats.multitest import fdrcorrection


cases = [(x, y) for x in [0.01, 0.001, 0.0001, 0.00001] for y in [1, 2, 3]]
x, y = cases[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]

min_rel_abundance = x
min_samples = 10
min_fold_change = 2
fdr_cutoff = 0.05
random_selection = True
count_studies = False
min_studies = 0
MAX_SIZE = 50
TOP_SIZE = 1000
POP_SIZE = 10000
PART_SIZE = 50
env_filter = None

#seed(0)

output_folder = f"../communities/repl_10000_{x}_{y}/"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#metadata = pd.read_csv('../data/emp_qiime_mapping_qc_filtered.tsv', sep='\t')
#metadata.rename(columns={'#SampleID': 'sample'}, inplace=True)
#host_samples = metadata.query("empo_1 == 'Host-associated'")["sample"]
#free_samples = metadata.query("empo_1 == 'Free-living'")["sample"]

df = pd.read_csv('../data/emp_150bp_filtered.tsv', sep='\t')

#if env_filter == "host":
#    df = df[df["sample"].isin(host_samples)]

#if env_filter == "free":
#    df = df[df["sample"].isin(free_samples)]

df = df.query("value > @min_rel_abundance")

data = pd.pivot_table(df, index='org_id', columns='sample', values='value', aggfunc=sum, fill_value=0)
bin_data = np.sign(data)
bin_data = bin_data[(bin_data.sum(axis=1) >= min_samples)]

#df2 = pd.read_csv('../data/emp_groupby_study.tsv', sep='\t')
#data2 = pd.pivot_table(df2, index='org_id', columns='study_id', values='value', aggfunc=sum, fill_value=0)
#bin_data2 = np.sign(data2)

def evaluate(rows):
    
    x = bin_data.loc[rows,:]
    total = x.min(axis=0).sum()
    n_samples = bin_data.shape[1]
    probability = (x.sum(axis=1) / n_samples).prod()
    expected = n_samples * probability
    p_value = 1 - binom.cdf(total, n_samples, probability)

#    studies = bin_data2.loc[rows,:].min(axis=0).sum()

    return total, expected, p_value#, studies

combinations = []
species = list(bin_data.index.values)
pop_sample = [(x,) for x in species]
pool = Pool()

for i in range(2, MAX_SIZE+1):

    print("Size:", i, "at", datetime.now(), flush=True)

    # combine all tuples with one more species (without repetition)
    combinations = (tuple(sorted(set(x)|{y})) for x in pop_sample
                    for y in species if y not in set(x))
    combinations = list(set(combinations))

    # compute scores for all combinations
    scores = pool.map(evaluate, combinations)
    scores = pd.DataFrame(scores, columns=["total", "expected", "p_value"])#, "studies"])
    scores["community"] = combinations
    
    # compute FDR
    scores["significant"], _ = fdrcorrection(scores["p_value"], alpha=fdr_cutoff)

    if count_studies:
        scores["final"] = np.sqrt(scores["total"] * scores["studies"])
    else:
        scores["final"] = np.sqrt(scores["total"])

    #filter FDR and min occurrence threshold
    selected = scores.query("significant and total >= @min_samples and total >= @min_fold_change * expected")# and studies >= @min_studies")
    
    # if no results, stop at this point
    if len(selected) == 0:
        print("no significant found.", flush=True)
        break
    else:
        print("selected:", len(selected), flush=True)

    top_size = min(TOP_SIZE, len(selected))
    pop_size = min(POP_SIZE, len(selected))

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

    comms = []
    for j, org_ids in enumerate(top):
        comm_id = 'comm_{}_{}'.format(i, j)
        entries = [(comm_id, org_id) for org_id in org_ids]
        comms.extend(entries)

        k, l = (j+1) % PART_SIZE, j // PART_SIZE
        if k == 0 or j == len(top) - 1:
            df_k = pd.DataFrame(comms, columns=['community', 'species'])
            df_k.to_csv(output_folder + 'cooccurring_{}_{}.tsv'.format(i, l),
                        sep='\t', index=False, header=False)
            comms = []
