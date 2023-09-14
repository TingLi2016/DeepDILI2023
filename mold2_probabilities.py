import os
from os import path
import pandas as pd
from functools import reduce

from mold2_base_classifier import classifier_names

def combine_validation_probabilities(base_path, mcc, probability_path, run_name):
  def read_validation_class_df(cl, base_path, run_name):
    suffix = "_K_7" if cl == "knn" else ""
    return pd.read_csv(base_path+'/'+cl+'/validation_class/validation_'+cl+'_paras_'+run_name+suffix +'.csv')

  base_cl_dfs = {cl : read_validation_class_df(cl, base_path, run_name) for cl in classifier_names}

  data = combine_probabilities(base_cl_dfs, mcc, include_ytrue=True)

  os.makedirs(probability_path, exist_ok=True)
  data.to_csv(path.join(probability_path +'/validation_probabilities_'+run_name+'.csv'))

def combine_probabilities(base_cl_dfs, mcc, include_ytrue):
  seeds = {
    cl : mcc[(mcc[cl+"_MCC"] >= 0.128) & (mcc[cl+"_MCC"] <= 0.351)].seed.unique()
    for cl in classifier_names
  }

  def make_derived_cl_df(cl, incl_ytrue):
    orig_df = base_cl_dfs[cl]
    df = orig_df[['id','y_true'] if incl_ytrue else ['id']]
    for _, seed in enumerate(seeds[cl]):
      cols = [col for col in orig_df.columns if 'prob_'+cl+'_seed_'+seed in col]
      df[cl+'_seed_'+seed] = orig_df[[*cols]]
    return df

  derived_cl_dfs = [
    make_derived_cl_df(cl, incl_ytrue=(include_ytrue and cl == 'knn'))
    for cl in classifier_names
  ]

  return reduce(lambda acc_df, df: pd.merge(acc_df, df, on='id', how='left'), derived_cl_dfs)