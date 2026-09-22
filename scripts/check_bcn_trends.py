#!/usr/bin/env python3
"""Check the optional BCN R trend demonstration against study estimates."""
from pathlib import Path
import argparse
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);a=p.parse_args()
x=pd.read_csv(a.output/'trends.csv').sort_values('Phase').reset_index(drop=True)
y=pd.read_csv(ROOT/'examples/expected/BCN_UB_trends.csv').sort_values('Phase').reset_index(drop=True)
pd.testing.assert_frame_equal(x,y,check_dtype=False,rtol=1e-10,atol=1e-10)
print('PASS: BCN_UB trends and confidence limits agree with the study outputs, including unavailable estimates.')
