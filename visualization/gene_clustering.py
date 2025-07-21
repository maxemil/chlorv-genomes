import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

orthofile = 'ChlorVs_orthos/Results_ChlorVs_orthos/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv'
df = pd.read_csv(orthofile, sep='\t', index_col=0)

fig, ax = plt.subplots(1, 1, figsize=(30,22))
sns.heatmap(df, cmap=sns.color_palette("Blues", as_cmap=True), annot=True, fmt=".1f", vmin=100, vmax=350)
fig.tight_layout()
fig.savefig('Orthogroup_sharing.pdf')

orthofile = 'ChlorVs_orthos/Results_ChlorVs_orthos/Orthogroups/Orthogroups.tsv'
df = pd.read_csv(orthofile, sep='\t', index_col=0)
df.fillna('', inplace=True)
# df = df.applymap(lambda x: len(x.split(',')) if x else 0)
df = df.applymap(lambda x: 1 if x else 0)
#remove nodes of degree 1
df = df[-(df.sum(axis=1) == 1)]
#drop duplcated ogs with same taxonomic distribution
# df = df.drop_duplicates()
#select orthogroups present in at least 10 genomes
df = df[df.sum(axis=1) > 10]

df.columns = pd.Series(df.columns).apply(lambda x:'_'.join(x.split('_')[1:]) if x.startswith('GV') else x)
aliimimivirinae = [g.strip() for g in open('GVDB_references.txt')]
genome2group = {}
for f in aliimimivirinae[0:25]:
    genome2group[f] = 'group I'
genome2group[aliimimivirinae[25]] = 'undetermined'
for f in aliimimivirinae[26:]:
    genome2group[f] = 'group II'
for f in ['ChlorV-1', 'ChlorV-2', 'ChlorV-3', 'ChlorV-4', 'ChlorellaV', 'CroV']:
    genome2group[f] = 'group I'
lut = dict(zip(set(genome2group.values()), "rbg"))
groups = pd.Series(df.columns).apply(lambda x: genome2group[x])
groups.name = 'groups'
col_colors = groups.map(lut)
col_colors.index = df.columns

clusmap = sns.clustermap(df, figsize=(30,25), cbar_pos=None, col_colors=col_colors, 
                        method='average', metric='jaccard', 
                        cmap=sns.light_palette("darkkhaki", as_cmap=True), yticklabels=False,
                        dendrogram_ratio=(0, .2))
ordering = [t.get_text() for t in heatmap.get_xticklabels()]
clusmap.ax_row_dendrogram.set_visible(False)                        
plt.tight_layout()
plt.savefig('Orthogroup_sharing_cluster.pdf')
plt.cla()

dfm = df.melt(ignore_index=False)
dfm = dfm[dfm['value'] > 0]
dfm['variable'].to_csv('ortho2genome.tsv', sep='\t', index=True, header=False)

with open('ortho2genome_types.tsv', 'w') as out:
    for i in set(df.index):
        print(i, 'geneFamily', sep='\t', file=out)
with open('ortho2genome_types.tsv', 'a') as out:
    for i in set(df.columns):
        print(i, 'genome', sep='\t', file=out)


from Bio import SeqIO, SeqUtils
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob

gen2len = {}
for f in glob.glob('GVDB_references_fna/*') + glob.glob('genomes/*'):
    base = f.split('/')[1].split('.')[0]
    length = 0
    for rec in SeqIO.parse(f, 'fasta'):
        length += len(rec.seq)
    gen2len[base] = length
df = pd.DataFrame.from_dict(gen2len, orient='index')
df.index = pd.Series(df.index).apply(lambda x:'_'.join(x.split('_')[1:]) if x.startswith('GV') else x)
df.columns = ['Genome Length']
df['Group'] = pd.Series(df.index, index=df.index).apply(lambda x: genome2group[x])
df = df.loc[ordering]

fig, ax = plt.subplots(1, 1, figsize=(10,5))
sns.barplot(data=df, x=df.index, y='Genome Length', hue='Group', dodge=False, ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,
        horizontalalignment='right')
fig.tight_layout()
fig.savefig('genome_sizes.pdf')