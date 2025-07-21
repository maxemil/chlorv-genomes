#! /usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df_perc = pd.read_csv('pyani/ANIb_percentage_identity.tab', index_col=0, sep='\t')
df_perc = df_perc *100

df_cov = pd.read_csv('pyani/ANIb_alignment_coverage.tab', index_col=0, sep='\t')
df_cov = df_cov *100

g = sns.clustermap(df_perc, figsize=(5,5), cmap="crest", annot=df_perc, fmt='.2f', dendrogram_ratio=(0, .1), vmin=0, vmax=100)
g.ax_col_dendrogram.set_title('ANIb percentage identity')
plt.savefig('ANIb_percentage_identity.pdf', dpi=300, bbox_inches='tight')
plt.cla()

g = sns.clustermap(df_cov, figsize=(5,5), cmap="crest", annot=df_cov, fmt='.2f', dendrogram_ratio=(0, .1), vmin=0, vmax=100)
g.ax_col_dendrogram.set_title('ANIb alignment coverage')
plt.savefig('ANIb_alignment_coverage.pdf', dpi=300, bbox_inches='tight')
plt.cla()
