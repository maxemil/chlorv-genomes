import pandas as pd
from itertools import combinations, product

df = pd.read_csv('Aliimimivirinae_only_isolates.o6', sep='\t', header=None)
df.columns = ['query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits']

gff = pd.read_csv('Aliimimivirinae_only_isolates.gff', sep='\t', header=None, comment='#')
gff.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
gff.index = gff['attributes'].apply(lambda x: x.split(';')[0].split('=')[1])

df['qstart'] = df['query'].apply(lambda x: gff.loc[x, 'start'])
df['qend'] = df['query'].apply(lambda x: gff.loc[x, 'end'])
df['tstart'] = df['target'].apply(lambda x: gff.loc[x, 'start'])
df['tend'] = df['target'].apply(lambda x: gff.loc[x, 'end'])
df['tstrand'] = df['target'].apply(lambda x: gff.loc[x, 'strand'])
df['qstrand'] = df['query'].apply(lambda x: gff.loc[x, 'strand'])
df['qgenome'] = df['query'].apply(lambda x: x.split('..')[0])
df['tgenome'] = df['target'].apply(lambda x: x.split('..')[0])
df = df[df.fident > 0.5]
df = df.drop(['fident', 'alnlen', 'mismatch', 'gapopen', 'bits', 'evalue'],axis=1)

def subset_df(df, genome1, genome2):
    sdf = df[(df['qgenome'] == genome1) & (df['tgenome'] == genome2)]
    sdf = sdf.sort_values('qstart')
    return sdf

def get_pplus(pair, step1, step2):
    pp1 = f"{pair[0].split('..')[0]}..{(int(pair[0].split('..')[1])+step1):003d}"
    pp2 = f"{pair[1].split('..')[0]}..{(int(pair[1].split('..')[1])+step2):003d}"
    return (pp1, pp2)

def find_next_match(p, pairwise):
    for step1, step2 in product([0,1], [-1,0,1]):
        pp = get_pplus(p, step1, step2)
        if pp in pairwise:
            return pp
    return None

def get_ranges(df):
    pairwise = df.apply(lambda x: (x['query'], x['target']), axis=1).to_list()

    p = pairwise.pop(0)
    synrange = [p]
    all_ranges = []
    while pairwise:
        p = find_next_match(p, pairwise)
        if p:
            pairwise.remove(p)
            synrange.append(p)
        else:
            if len(synrange) > 1:
                all_ranges.append(synrange)
            p = pairwise.pop(0)
            synrange = [p]
    return all_ranges


def get_strand(f, gff):
    strand = 0
    for g in f:
        if gff.loc[g[0], 'strand'] == gff.loc[g[1], 'strand']:
            strand += 1
        else:
            strand -= 1
    if strand > 0:
        return '+'
    else:
        return '-'

def write_ranges(ranges, gff):
    with open('Aliimimivirinae_only_isolates.filtered.links', 'w') as out:
        for f in ranges:
            min1 = min([gff.loc[g[0], 'start'] for g in f])
            max1 = max([gff.loc[g[0], 'end'] for g in f])
            min2 = min([gff.loc[g[1], 'start'] for g in f])
            max2 = max([gff.loc[g[1], 'end'] for g in f])
            genome1 = gff.loc[f[0][0], 'seqname']
            genome2 = gff.loc[f[0][1], 'seqname']
            strand = get_strand(f, gff)
            print(genome1, genome2, min1, max1, min2, max2, strand, sep='\t', file=out)

ranges = []
for g, c in combinations(df['qgenome'].unique(), 2):
    sdf = subset_df(df, g, c)
    ranges += get_ranges(sdf)
write_ranges(ranges, gff)
