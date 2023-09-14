import pandas as pd
import mold2_DeepDILI

def compare_test_results(tr_objs, static_data_dir):
  req_data = pd.read_csv(static_data_dir+'/QSAR_year_338_pearson_0.9.csv', low_memory=False)
  req_data = req_data[req_data.final_year >= 1997]
  pred_result = mold2_DeepDILI.predict(req_data, tr_objs)
  print("Performing comparison vs. expected values for QSAR_... data subset:")
  pred_probs = pred_result['prob_pred'].tolist()
  expected_probs = pd.read_csv(static_data_dir+'/expected_ge1997_test_results.csv')['prob_feature_test'].tolist()
  assert len(pred_probs) == len(expected_probs)
  for i in range(len(pred_probs)):
    print(str(pred_probs[i]) + "   " + str(expected_probs[i]))
