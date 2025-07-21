#! /usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
from Bio import SeqIO
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def coverage_window(covs, windowsize):
    cv = pd.Series()
    i = windowsize
    while i <= len(covs):
        cv[i-windowsize/2] = covs[i-windowsize:i].mean()
        i += windowsize
    return cv

bed = 'ChlorV-1.6mA.context.bed'
threshold = 60
min_threshold = 10
infile = 'ChlorV-1.fna'
contig = 'ChlorV-1'
gff_file = 'ChlorV-1.gff'

df = pd.read_csv(bed, sep='\t', header=None)
bedMethylHeader = ['chrom', 'start', 'end', 'mod_base', 'score', 'strand',
                'start_pos', 'end_pos', 'color', 'N_valid_cov', 'fraction_modified',
                'N_mod', 'N_canonical', 'N_other_mod', 'N_delete', 'N_fail', 'N_diff', 'N_nocall', 'context']
df.columns = bedMethylHeader

modification_labels = {'a':'6mA', 'm':'5mC', '21839':'4mC'}
df = df[df.N_valid_cov > df.N_diff]
df.mod_base = df.mod_base.apply(lambda x: modification_labels[x])

# df = df[df.context.apply(lambda x: True if len(x) == 4 else False)]
df['context2'] = df.context.apply(lambda x: x[1:3])

gff = pd.read_csv(gff_file, sep='\t', comment='#', header=None)
# genes = pd.Series(np.zeros(gff[4].max() + 1))
# for f in gff[gff[6] == '-'].iterrows():
#     genes.loc[f[1][3]:f[1][4]] = 1



df5perc = df[df.fraction_modified >= threshold]
df5perc = df5perc[df5perc.mod_base == '5mC'] 

recdict = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
seqlen = len(recdict[contig].seq)

figsize = (24,4)
fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=300)
sns.scatterplot(data=df5perc, x='start', y="fraction_modified", hue='mod_base', style='strand', ax=ax)

dfp = df5perc[df5perc['strand'] == '-']
adfm = pd.Series(np.zeros(seqlen + 1))
adfm.loc[dfp['start']] = dfp['fraction_modified'].to_list()

windowsize = 1000
ax2 = ax.twinx()
sns.lineplot(coverage_window(adfm.astype(bool).astype(int)*100, windowsize), ax=ax2, color=(0.8392156862745098, 0.15294117647058825, 0.1568627450980392))
# ax2.set_ylim(0, 2)
ax2.tick_params(axis='y', colors=(0.8392156862745098, 0.15294117647058825, 0.1568627450980392))

dfp = df5perc[df5perc['strand'] == '+']
adfp = pd.Series(np.zeros(seqlen + 1))
adfp.loc[dfp['start']] = dfp['fraction_modified'].to_list()

ax3 = ax.twinx()
sns.lineplot(coverage_window(adfp.astype(bool).astype(int)*100, windowsize), ax=ax3, color=(1.0, 0.4980392156862745, 0.054901960784313725))
ax3.set_ylim(ax2.get_ylim())
ax3.tick_params(axis='y', colors=(1.0, 0.4980392156862745, 0.054901960784313725))

# ax3 = ax.twinx()
# sns.lineplot(genes, ax=ax3, color='goldenrod')
# ax3.tick_params(axis='y', colors='goldenrod')

# sns.move_legend(ax, "center left", bbox_to_anchor=(1, .5))
lns = ax.get_lines() + ax2.get_lines() + ax3.get_lines()
labs = [l.get_label() for l in ax.get_lines()] + ['- rolling', '+ rolling']
ax.legend(lns, labs)
ax.set_xlabel('Genomic position')

fig.tight_layout()
outfile = f'{contig}_6mA_fraction_modified_scatter.pdf'
fig.savefig(outfile)
plt.close()

df5perc = df[df.fraction_modified >= min_threshold]
df5perc.sort_values('mod_base')

def annotate(data, **kws):
    n = len(data)
    ax = plt.gca()
    ax.text(.8, .8, f"N = {(data.fraction_modified >= 60).sum()}", transform=ax.transAxes)
    ax.text(.1, .8, f"N = {(data.fraction_modified < 60).sum()}", transform=ax.transAxes)

g = sns.FacetGrid(df5perc, col="mod_base", row='chrom', sharey=False, hue='mod_base', height=4, aspect=1)
g.map(sns.histplot, "fraction_modified", binwidth=2)
g.map_dataframe(annotate)
g.refline(x=60, color='gray')
plt.tight_layout()
outfile = f'{contig}_fraction_modified_hist.pdf'
plt.savefig(outfile)
plt.close()



context = {'mAT':'[ATCG]AT[ATCG]', 
            'mAA':'[ATCG]AA[ATCG]',
            'mAG':'[ATCG]AG[ATCG]',
            'mAC':'[ATCG]AC[ATCG]',
            'WmATW':'[AT]AT[AT]',
            'SmATS':'[CG]AT[CG]'}
count = {k:0 for k in context.keys()}
for k,v in context.items():
    count[k] = df5perc.context.apply(lambda x: True if re.match(v, x) else False).sum()

figsize = (6,4)
fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=300)
sns.barplot(count, ax=ax, color='darkgrey')
ax.set_xlabel('context')
ax.set_ylabel(f'count 6mA modified > {threshold}%')
fig.tight_layout()
outfile = f'{contig}_6mA_fraction_modified_context.pdf'
fig.savefig(outfile)
plt.close()


import logomaker
from Bio import SeqIO
import re

recdict = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))

# pattern = 'CA[ATGC][ATGC][ATGC][ATGC][ATGC][ATGC]TG'
pattern = '[AG]CG[CT]'
occurences = []
matches = []
position = 1
for m in re.finditer(pattern, str(recdict[contig].seq)):
    occurences.append(m.group(0))
    if adfp[m.start() + position] > threshold:
        matches.append(m.group(0))
adfm = np.flip(adfm).reset_index(drop=True)
for m in re.finditer(pattern, str(recdict[contig].seq.reverse_complement())):
    occurences.append(m.group(0))
    if adfm[m.start() + position] > threshold:
        matches.append(m.group(0))

lm = logomaker.alignment_to_matrix(matches, to_type='probability')
# df5perc = df[df.fraction_modified >= threshold]
# df5perc = df5perc[df5perc.mod_base == '6mA'] 
# lm = logomaker.alignment_to_matrix(df5perc.context, to_type='probability')

color_scheme = {
    'A' : (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    'T' : (1.0, 0.4980392156862745, 0.054901960784313725),
    'G' : (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    'C' : (0.8392156862745098, 0.15294117647058825, 0.1568627450980392)}

figsize = (10,2)
fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=300)
logomaker.Logo(lm, ax=ax, color_scheme=color_scheme)
# fig.savefig(f'{contig}_6mA_all_modified_logo.pdf')
fig.savefig('RCGY_logo.pdf')
plt.close()