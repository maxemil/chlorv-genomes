import pandas as pd
from itertools import pairwise
import glob 

for f in glob.glob('pyani/blastn_output/ChlorV-*vs_ChlorV-*.blast_tab.dataframe'):
    df = pd.read_csv(f, sep='\t')
    print(f)
    lgap = 0
    gaploc = (0,0)
    for i, j in pairwise(df['0'].tolist()):
        i = int(i.replace('frag', '')) * 500
        j = int(j.replace('frag', '')) * 500
        gap = j-i
        if gap > lgap:
            lgap = gap
            gaploc = (i, j)
    print(lgap, gaploc)