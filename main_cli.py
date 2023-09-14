from io import StringIO
import subprocess
import pickle
import sys
import warnings
warnings.filterwarnings('ignore'); warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

import mold2_DeepDILI

def predict_sdf(istream, tr_objs):
  mold2Proc = subprocess.run(
    ["/deep-dili/mold2", "-i", "/dev/stdin", "-o", "/dev/stdout"],
    stdin = istream, encoding="utf-8", capture_output=True, check=True
  )

  mold2_df = pd.read_table(
    StringIO(mold2Proc.stdout),
    header=0, sep='\t', dtype={'id': str}, keep_default_na=False
  )

  # Fill in any missing molecule ids with defaults based on row index.
  mold2_df['id'] = mold2_df.apply(
    lambda row: row['id'] or "anon-mol-{}".format(int(row.name)+1),
    axis=1
  )

  # Abort if final molecule ids are not unique.
  if len(mold2_df.id.unique()) != len(mold2_df.id):
    print("Duplicate molecule ids were found in SDF input (on first lines of entries).", file=sys.stderr)
    sys.exit(1)

  res = mold2_DeepDILI.predict(mold2_df, tr_objs)

  res.to_csv(sys.stdout, sep="\t", index=False)

if __name__ == "__main__":
  print(sys.argv)
  print(len(sys.argv))
  if len(sys.argv) != 3:
    raise ValueError('Expected arguments: <training-objects-file> <input-file>')

  with open(sys.argv[1], 'rb') as tr_file, open(sys.argv[2], 'rb') as sdf_file:
    tr_objs = pickle.load(tr_file)
    predict_sdf(sdf_file, tr_objs)
