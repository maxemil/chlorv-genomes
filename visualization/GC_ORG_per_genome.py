
from Bio import SeqIO, SeqUtils
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob

gen2gc = {}
for f in glob.glob('ffn/*'):
    base = f.split('/')[1].split('.')[0]
    genegc = []
    for rec in SeqIO.parse(f, 'fasta'):
        genegc.append(SeqUtils.GC(rec.seq))
    gen2gc[base] = pd.DataFrame(genegc)
    gen2gc[base]['genome'] = base

df = pd.concat(list(gen2gc.values()))
df = df.sort_values('genome')
sns.boxplot(data=df, y=0, x='genome')
plt.savefig('gene_GC.pdf')

# plot ORF lengths
genome2orf = []
for f in glob.glob('Aliimimivirinae_wt_isolates/*'):
    base = f.split('/')[1].split('.')[0]
    genome2orf += [(base,len(rec.seq)) for rec in SeqIO.parse(f, 'fasta')]

df = pd.DataFrame(genome2orf)
order = df.groupby(0).median().sort_values(1).index

fig, ax = plt.subplots(1, 1, figsize=(10,10))
sns.boxplot(data=df, y=1, x=0, ax=ax, order=order)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,
        horizontalalignment='right')
fig.tight_layout()
fig.savefig('Aliimimivirinae_ORF_len.pdf')

# plot ORF numbers
df = df.groupby(0).count()
df.sort_values(1, inplace=True)

fig, ax = plt.subplots(1, 1, figsize=(10,10))
sns.scatterplot(data=df, y=1, x=df.index, ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,
        horizontalalignment='right')
fig.tight_layout()
fig.savefig('Aliimimivirinae_ORF_count.pdf')

# plot genome GC
genome2GC = []
for f in glob.glob('Aliimimivirinae_wt_isolates_genomes/*'):
    base = f.split('/')[1].split('.')[0]
    genome2GC += [(base,SeqUtils.GC(rec.seq)) for rec in SeqIO.parse(f, 'fasta')]

df = pd.DataFrame(genome2GC)
df = df.groupby(0).mean()
df.sort_values(1, inplace=True)

fig, ax = plt.subplots(1, 1, figsize=(10,10))
sns.scatterplot(data=df, y=1, x=df.index, ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,
        horizontalalignment='right')
fig.tight_layout()
fig.savefig('Aliimimivirinae_GC_genome.pdf')

# plot ORF GC
genome2GC = []
for f in glob.glob('Aliimimivirinae_wt_isolates_ffn/*'):
    base = f.split('/')[1].split('.')[0]
    genome2GC += [(base,SeqUtils.GC(rec.seq)) for rec in SeqIO.parse(f, 'fasta')]

df = pd.DataFrame(genome2GC)
order = df.groupby(0).median().sort_values(1).index

fig, ax = plt.subplots(1, 1, figsize=(10,10))
sns.boxplot(data=df, y=1, x=0, ax=ax, order=order)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,
        horizontalalignment='right')
fig.tight_layout()
fig.savefig('Aliimimivirinae_GC_ORF.pdf')