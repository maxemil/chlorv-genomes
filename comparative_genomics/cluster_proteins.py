import pandas as pd
import networkx as nx
from collections import defaultdict
from Bio import SeqIO

seqdict = SeqIO.to_dict(SeqIO.parse('Aliimimivirinae_only_isolates.faa', 'fasta'))

df = pd.read_csv('Aliimimivirinae_only_isolates.o6', sep='\t', header=None)
df = df[df[0].apply(lambda x: 'ChlorV' in x)]
df = df[df[1].apply(lambda x: 'ChlorV' in x)]

prots = set(list(df[0]) + list(df[1]))

df = df[df.apply(lambda x: x[0] != x[1], axis=1)]
df = df[df[2] > 0.7]

df['qlen'] = df[0].apply(lambda x: len(seqdict[x].seq))
df['tlen'] = df[1].apply(lambda x: len(seqdict[x].seq))
df = df[df.apply(lambda x: True if x[3] >= x['qlen']*0.7 or x[3] >= x['tlen']*0.7 else False, axis=1)]

sets = [set(x) for x in df[[0,1]].itertuples(index=False, name=None)]

edges = []
for i, s1 in enumerate(sets):
    for j, s2 in enumerate(sets):
        if i != j:
            if set(s1).intersection(set(s2)):
                edges.append((i,j))

G = nx.Graph()
G.add_nodes_from(range(len(sets)))
G.add_edges_from(edges)

merged_sets2 = []
for c in nx.connected_components(G):
    combined_lists = [sets[i] for i in c]
    flat_list = set([item for sublist in combined_lists for item in sublist])
    merged_sets2.append(flat_list)


covered_prots = [i for s in merged_sets2 for i in s]
for p in prots:
    if p not in covered_prots:
        merged_sets2.append(set([p]))

all_prots = []
for f in merged_sets:
    prots = {s:'' for s in specs}
    for s in specs:
        prots[s] = ','.join([p for p in f if p.startswith(s)])
    prots['annot'] = ','.join(set([annots[a] for a in f]))
    all_prots.append(prots)
pd.DataFrame(all_prots).to_excel('all_prots.xlsx', index=False)