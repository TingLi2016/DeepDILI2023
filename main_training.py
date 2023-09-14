import sys
import warnings
warnings.filterwarnings('ignore'); warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import pickle

import verification
import mold2_training

if len(sys.argv) != 4 and len(sys.argv) != 5:
  raise ValueError('Expected arguments: static-data-dir results-dir output-training-objects-file [#-sig-features]')

static_data_dir = sys.argv[1]
res_dir = sys.argv[2]
tr_objs_fpath = sys.argv[3]
n_sig_feats = int(sys.argv[4], 10) if len(sys.argv) == 5 else 135

np.random.seed(1)

sig_feats = sorted(pd.read_csv(static_data_dir+'/important_features_order.csv').feature[n_sig_feats:])

tr_objs = mold2_training.train(static_data_dir, res_dir, sig_feats)

print("Training results vs expected:")
verification.compare_test_results(tr_objs, static_data_dir)

print("Writing training objects...")

with open(tr_objs_fpath, 'wb') as f:
  pickle.dump(tr_objs, f)

print("Training objects have been written to \"" + tr_objs_fpath + "\".")
