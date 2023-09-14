import pandas as pd
import keras

import mold2_base_classifier
from mold2_training_objs import TrainingObjects
from mold2_probabilities import combine_validation_probabilities
from mold2_base_classifier import classifier_names
import mold2_DeepDILI

def train(static_data_dir, res_dir, sig_feats):
  qsar = pd.read_csv(static_data_dir + '/QSAR_year_338_pearson_0.9.csv', low_memory=False)
  data = qsar.filter(items=qsar.columns[:5].tolist() + sig_feats, axis=1)
  data_split = pd.read_csv(static_data_dir + '/data_split.csv')
  best_model = keras.models.load_model(static_data_dir + '/mold2_best_model.h5')
  mcc = pd.read_csv(static_data_dir + '/combined_score.csv')
  res_name = 'train'

  cl_objs = {
    cl :  mold2_base_classifier.train(cl, data, data_split, res_name, res_dir + '/' + cl)
    for cl in classifier_names
  }

  probs_dir = res_dir + '/probabilities_output'

  combine_validation_probabilities(res_dir, mcc, probs_dir, res_name)

  scaler = mold2_DeepDILI.validation_predict(probs_dir, res_name, best_model, res_dir + '/performance')

  return TrainingObjects(sig_feats, best_model, mcc, scaler, cl_objs)