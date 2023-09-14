from io import StringIO
import subprocess
import sys
import warnings
warnings.filterwarnings('ignore'); warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
from flask import Flask, request, Response, send_from_directory, jsonify, make_response
from rdkit import Chem

import mold2_training
import mold2_DeepDILI

if len(sys.argv) != 4 and len(sys.argv) != 5:
  raise ValueError(
    'Expected arguments: static-data-dir results-dir mold2exe [#-sig-features]'
  )

static_data_dir = sys.argv[1]
res_dir = sys.argv[2]
mold2exe = sys.argv[3]
n_sig_feats = int(sys.argv[4], 10) if len(sys.argv) == 5 else 135

np.random.seed(1)

sig_feats = sorted(pd.read_csv(static_data_dir+'/important_features_order.csv').feature[n_sig_feats:])

print("Training...")

tr_objs = mold2_training.train(static_data_dir, res_dir, sig_feats)

print("Training completed.")
print("Starting web application.")

app = Flask('deep-dili', static_folder='webapp', static_url_path='')

@app.route('/', methods=['GET'])
def index():
  return send_from_directory('webapp', 'index.html')

@app.route('/deep-dili/predict-sdf', methods=['POST'])
def predict_sdf_body():
  file = request.files['sdf']
  if not file:
    raise ValueError('SD file must be provided.')
  sdfContent = file.read().decode()
  return predict_sdf(sdfContent)

@app.route('/deep-dili/predict-smiles', methods=['GET'])
def predict_smiles():
  smiles = request.args["s"]
  mol = Chem.MolFromSmiles(smiles)
  sio = StringIO()
  with Chem.SDWriter(sio) as w:
    w.write(mol)
  sdfContent = sio.getvalue()
  return predict_sdf(sdfContent)

def predict_sdf(sdfContent):
  mold2Proc = subprocess.run(
    [mold2exe],
    input = sdfContent, capture_output=True, encoding="utf-8"
  )
  if mold2Proc.returncode != 0:
    return make_response(jsonify(mold2Proc.stderr), 400)

  mold2_df = pd.read_table(
    StringIO(mold2Proc.stdout),
    header=0, sep='\t', dtype={'id': str}, keep_default_na=False
  )

  # Fill in any missing molecule ids with defaults based on row index.
  mold2_df['id'] = mold2_df.apply(
    lambda row: row['id'] or "anon-mol-{}".format(int(row.name)+1),
    axis=1
  );

  # Balk if final molecule ids are not unique.
  if len(mold2_df.id.unique()) != len(mold2_df.id):
    return make_response(jsonify("Duplicate molecule identifiers in SD file."), 400)

  res = mold2_DeepDILI.predict(mold2_df, tr_objs)

  return Response(res.to_json(orient='records'), mimetype="application/json")

app.run(host="0.0.0.0")