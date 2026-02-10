import pandas as pd
import glob
from Bio import SeqIO

# load all .prots.o6 and .add.o6 files into a single dataframe
df = pd.concat([pd.read_csv(f, sep="\t", header=None) for f in glob.glob('ChlorV-*.o6')])
header = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'e-value', 'bitscore']
df.columns = header

df['e-value'] = df['e-value'].astype(float)
df['annotation'] = ''    

annot = pd.read_excel('gene_clusters_annotations.xlsx')

for spec in ['ChlorV-1', 'ChlorV-2', 'ChlorV-3', 'ChlorV-4']:
    recdict = SeqIO.to_dict(SeqIO.parse(f'{spec}.faa', 'fasta'))

    for row in annot.iterrows():
        a = f'putative {row[1]["annot_consensus"]}' if row[1]['annot_consensus'] != 'hypothetical protein' else row[1]['annot_consensus']
        if pd.isna(row[1][spec]): 
            continue
        for rec in row[1][spec].split(','):
            if rec in recdict:
                recdict[rec].description = a
            else:
                print(f'Warning: {rec} not found in recdict')
        df.loc[df['target'] == rec, 'annotation'] = a


df.to_excel('ChlorVs_vs_Bigelowiella_GEVEs.xlsx', index=False)
