from Bio import SeqIO
import glob
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt

genome_marker = defaultdict(dict)
markers = []
for f in glob.glob('alignments_manual/*.faa'):
    base = f.split('/')[1].split('.')[0]
    markers.append(base)
    for rec in SeqIO.parse(f, 'fasta'):
        genome_marker[rec.id.split('..')[0]][base] = 1

df = pd.DataFrame.from_dict(genome_marker)
df.fillna(0, inplace=True)


aliimimivirinae = [g.strip() for g in open('GVDB_references.txt')]
df.columns = pd.Series(df.columns).apply(lambda x:x.replace('Megaviricetes_Imitervirales_Mimiviridae_', ''))
aliimimivirinae += ['ChlorV-1', 'ChlorV-2', 'ChlorV-3', 'ChlorV-4', 'ChlorellaV', 'crov']
df = df[aliimimivirinae]

plt.rcParams['svg.fonttype'] = 'none'
fig, axs = plt.subplots(4,int(df.shape[1]/4), figsize=(40,5))
for genome, ax in zip(df.T.iterrows(), axs.flat):
    wedges, text = ax.pie(df['crov']/7, colors=[(0, 0, 0, 0)]*7, normalize=False, wedgeprops = {'linewidth': 3, 'edgecolor':'black'}, labeldistance=1.3)
    ax.set_title(genome[0])
    for i, (index, marker) in enumerate(genome[1].items()):
        if marker:
            wedges[i].set_facecolor('green')
        if genome[0] == 'crov':
            text[i].set_text(index)
fig.savefig('Coulson_plots.svg')
plt.cla()