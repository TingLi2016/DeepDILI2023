import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn import metrics
from sklearn.metrics import confusion_matrix, f1_score, roc_curve, roc_auc_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from xgboost import XGBClassifier

from mold2_training_objs import BaseClassifierTrainingObjects

classifier_names = ['knn', 'lr', 'svm', 'rf', 'xgboost']

def create_classifier(cl_name):
  if cl_name == 'knn':
    return KNeighborsClassifier(n_neighbors=7)
  elif cl_name == 'lr':
    return LogisticRegression(C=0.1, max_iter=300, class_weight=None)
  elif cl_name == 'rf':
    return RandomForestClassifier(random_state=1, n_estimators=700, max_depth=11,  min_samples_leaf=5, class_weight='balanced', bootstrap = True, max_features='log2')
  elif cl_name == 'svm':
    return SVC(kernel='rbf', C=1, gamma='scale', probability=True,  class_weight = 'balanced', random_state=1)
  elif cl_name == 'xgboost':
    return XGBClassifier(learning_rate=0.01, n_estimators=700, max_depth=11, subsample=0.7, scale_pos_weight=0.66, eval_metric="logloss")
  else:
    raise ValueError('classifier name "' + cl_name + '" not recognized.')

def measurements(y_test, y_pred, y_pred_prob):
    acc = metrics.accuracy_score(y_test, y_pred)
    sensitivity = metrics.recall_score(y_test, y_pred)
    TN, FP, FN, TP = confusion_matrix(y_test, y_pred).ravel()
    specificity = TN/(TN+FP)
    precision = metrics.precision_score(y_test, y_pred)
    f1 = metrics.f1_score(y_test, y_pred)
    mcc = metrics.matthews_corrcoef(y_test, y_pred)
    auc = roc_auc_score(y_test, y_pred_prob)
    npv = TN/(TN+FN)
    return [TN, FP, FN, TP, acc, auc, sensitivity, specificity, precision, npv, f1, mcc]

def model_predict(X, y, model, col_name):
  y_pred_prob = model.predict_proba(X)
  # keep probabilities for the positive outcome only
  y_pred_prob = y_pred_prob[:, 1]
  y_pred_class = np.where(y_pred_prob > 0.5, 1, 0)

  ###create dataframe
  pred_result = pd.DataFrame()
  pred_result['prob_' + col_name] = y_pred_prob
  pred_result['class_' + col_name] = y_pred_class

  if y is not None:
    pred_result['id'] = y.index
    pred_result['y_true'] = y.values

  performance = measurements(y, y_pred_class, y_pred_prob) if y is not None else None

  return pred_result, performance

def make_col_name(classifier, seed, runName):
  colNameSuffix = '_K_7' if classifier == 'knn' else ''
  return classifier + '_seed_' + seed + '_paras_'+ runName + colNameSuffix

def train(cl_name, orig_data, split, run_name, output_base_dir):
    data = orig_data[['DILI_label', 'final_year', *orig_data.columns[5:]]] # feature columns start at index 5

    # split the data into training, validation
    X_org,  y_org = data[data.final_year<1997].iloc[:,2:], data[data.final_year<1997]['DILI_label']
    X, X_val, y, y_val = train_test_split(X_org,  y_org, test_size=0.2, stratify=y_org, random_state=7)

    tpPath = output_base_dir + '/training_performance'
    vpPath = output_base_dir + '/validation_performance'
    tcPath = output_base_dir + '/training_class'
    vcPath = output_base_dir + '/validation_class'

    os.makedirs(tpPath, exist_ok=True)
    os.makedirs(vpPath, exist_ok=True)
    os.makedirs(tcPath, exist_ok=True)
    os.makedirs(vcPath, exist_ok=True)

    train_results = {}
    validation_results = {}
    pred_val_df = pd.DataFrame()

    seeds = [['' for _ in range(5)] for _ in range(20)]
    col_names = [['' for _ in range(5)] for _ in range(20)]
    scalers = [[None for _ in range(5)] for _ in range(20)]
    classifiers = [[None for _ in range(5)] for _ in range(20)]

    for i in range(20):
      for j in range(5):
        seed = str(i)+'_skf_'+str(j)
        train_index = split[split[seed + '_status'] == 'train'][seed].unique()
        validation_index = split[split[seed + '_status'] == 'validation'][seed].unique()

        ###get train, validation dataset
        X_train, X_validation = X.iloc[train_index,:], X.iloc[validation_index,:]
        y_train, y_validation = y.iloc[train_index], y.iloc[validation_index]

        ### scale the input
        sc = MinMaxScaler()
        sc.fit(X_train)
        X_train = sc.transform(X_train)
        X_validation = sc.transform(X_validation)
        X_val_s = sc.transform(X_val)

        col_name = make_col_name(cl_name, seed, run_name)

        cl = create_classifier(cl_name)
        cl.fit(X_train, y_train)

        ### predict training results
        train_class, train_results[col_name] = model_predict(X_validation, y_validation, cl, col_name)

        ### predict validation results
        validation_class, validation_results[col_name] = model_predict(X_val_s, y_val, cl, col_name)

        pred_val_df = pd.concat([pred_val_df, validation_class], axis=1, sort=False)

        train_class.to_csv(tcPath+'/train_'+col_name+'.csv')

        seeds[i][j] = seed
        col_names[i][j] = col_name
        scalers[i][j] = sc
        classifiers[i][j] = cl

    ###save the result of validation results
    col_name_suffix = '_K_7' if cl_name == 'knn' else ''
    file_name_suffix = cl_name + '_paras_' + run_name + col_name_suffix
    pd.DataFrame(data=train_results.items()).to_csv(tpPath+'/train_'+file_name_suffix+'.csv')
    pred_val_df.to_csv(vcPath+'/validation_'+file_name_suffix+'.csv')
    pd.DataFrame(data=validation_results.items()).to_csv(vpPath+'/validation_'+file_name_suffix+'.csv')

    return BaseClassifierTrainingObjects(cl_name, classifiers, seeds, scalers)

def predict(cl_tr_objs: BaseClassifierTrainingObjects, data, sig_feats):
  feature_cols = data.filter(items=sig_feats, axis=1).columns.tolist()
  X_test = data[feature_cols]

  res_df = pd.DataFrame()
  res_df['id'] = data['id'] if 'id' in data.columns else data.index

  for i in range(20):
    for j in range(5):
      X_test_s = cl_tr_objs.scalers[i][j].transform(X_test)
      col_name = make_col_name(cl_tr_objs.classifier_name, cl_tr_objs.seeds[i][j], "predict")
      test_class, _ = model_predict(X_test_s, None, cl_tr_objs.classifiers[i][j], col_name)
      res_df = pd.concat([res_df, test_class], axis=1, sort=False)

  return res_df