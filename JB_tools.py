__author__ = 'jordanburke'

import pandas as pd
import sys
#sys.path.append('./SP_ANALYSIS/SPTools/')
#import SPTools


## Open JuncBase output
def read_junc_out(juncbase_output):
    with open(juncbase_output, 'r') as fin:
        junc_df = pd.read_csv(fin, sep='\t')
    return junc_df

##