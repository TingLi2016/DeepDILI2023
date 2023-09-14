import os
import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.metrics import confusion_matrix, f1_score, roc_auc_score, balanced_accuracy_score
from sklearn.preprocessing import StandardScaler

import mold2_base_classifier
from mold2_base_classifier import classifier_names
from mold2_probabilities import combine_probabilities
from mold2_training_objs import TrainingObjects

def model_predict(X, y, model, col_name):
  y_pred = model.predict(X)
  y_pred_class = np.where(y_pred > 0.5, 1, 0)
  pred_result = pd.DataFrame()
  pred_result['id'] = y.index
  pred_result['y_true'] = y.values
  pred_result['prob_'+col_name] = y_pred
  pred_result['class_'+col_name] = y_pred_class

  result=measurements(y, y_pred_class, y_pred)
  return pred_result, result

def measurements(y_test, y_pred, y_pred_prob):
  acc = metrics.accuracy_score(y_test, y_pred)
  precision = metrics.precision_score(y_test, y_pred)
  f1 = metrics.f1_score(y_test, y_pred)
  mcc = metrics.matthews_corrcoef(y_test, y_pred)
  auc = roc_auc_score(y_test, y_pred_prob)
  sensitivity = metrics.recall_score(y_test, y_pred)
  balanced_accuracy = balanced_accuracy_score(y_test, y_pred)

  TN, FP, FN, TP = confusion_matrix(y_test, y_pred).ravel()
  specificity = TN/(TN+FP)
  npv = TN/(TN+FN)
  return [TN, FP, FN, TP, acc, auc, sensitivity, specificity, precision, npv, f1, mcc, balanced_accuracy]

def dim_reduce(df, col_name, best_model):
  X = df.iloc[:, 3:]
  y = df.loc[:, 'y_true']

  sc = StandardScaler()
  sc.fit(X)
  X = sc.transform(X)

  train_class, train_result = model_predict(X, y, best_model, col_name)

  return train_class, train_result, sc

def sep_performance(df):
  cols = ['TN', 'FP', 'FN', 'TP', 'Accuracy', 'AUC', 'Sensitivity', 'Specificity', 'PPV', 'NPV', 'F1', 'MCC', 'Balanced_accuracy']
  for i, col in enumerate(cols):
    if i == 0:
      df[col] = df.value.str.split(',').str[i].str.split('[').str[1].values
    elif i == len(cols)-1:
      df[col] = df.value.str.split(',').str[i].str.split(']').str[0].values
    else:
      df[col] = df.value.str.split(',').str[i].values

  for i, col in enumerate(cols):
    if i < 4:
      df[col] = df[col].astype(int)
    else:
      df[col] = df[col].astype(float)
      df[col] = round(df[col], 3)
  del df['value']

  return df

def reform_result(results):
  df = pd.DataFrame(data=results.items())
  ###reform the result data format into single colum
  df = df.rename(columns={0:'name', 1:'value'})
  df['name'] = df['name'].astype('str')
  df['value'] = df['value'].astype('str')
  df = sep_performance(df)
  return df


def validation_predict(probability_path, var, best_model, result_path):
  data = pd.read_csv(probability_path+'/validation_probabilities_' + var + '.csv')

  train_results={}
  col_name = 'feature_' + var

  train_class, train_result, scaler = dim_reduce(data, col_name, best_model)

  train_results[col_name] = train_result
  vcPath = result_path + '/validation_class'
  os.makedirs(vcPath, exist_ok=True)
  train_class.to_csv(vcPath+'/validation_'+col_name+'.csv')

  vpPath = result_path + '/validation_performance'
  os.makedirs(vpPath, exist_ok=True)
  reform_result(train_results).to_csv(vpPath+'/validation_'+col_name+'.csv')

  return scaler

def predict(data, tr_objs: TrainingObjects):

  base_cl_dfs = {
    cl: mold2_base_classifier.predict(tr_objs.base_cl_objs[cl], data, tr_objs.features)
    for cl in classifier_names
  }

  base_cls_df = combine_probabilities(base_cl_dfs, tr_objs.mcc, include_ytrue=False)

  X_test = tr_objs.scaler.transform(base_cls_df.iloc[:, 1:])
  y_pred = tr_objs.best_model.predict(X_test)

  res = pd.DataFrame()
  res['mol_id'] = base_cls_df['id']
  res['prob_pred'] = y_pred
  res['class_pred'] = np.where(y_pred > 0.5, 1, 0)
  return res