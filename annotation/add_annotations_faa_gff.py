import pandas as pd 
from Bio import SeqIO

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
    SeqIO.write(recdict.values(), f'{spec}.annot.faa', 'fasta')

    gff_file = '{}.gff'.format(spec)
    commentlines = []
    for line in open(gff_file):
        if line.startswith('#'):
            commentlines.append(line)
    gff = pd.read_csv(gff_file, sep='\t', comment='#', header=None)
    gff['id'] = gff[8].apply(lambda s: dict(item.split('=') for item in s.split(';'))['ID'])
    # add locus_tag to gff
    gff['locus_tag'] = gff['id'].apply(lambda s: re.sub(r'[^a-zA-Z0-9]', '', s.split('..')[0]).upper() + '_' + s.split('..')[1])
    gff[8] = gff.apply(lambda x: x[8] + ';' + 'locus_tag=' + x['locus_tag'], axis=1)
    # add product to gff
    gff['product'] = gff['id'].apply(lambda s: recdict[s].description)
    gff[8] = gff.apply(lambda x: x[8] + ';' + 'product=' + x['product'], axis=1)
    gff.drop(columns=['id', 'product', 'locus_tag'], inplace=True)
    with open('{}.annot.gff'.format(spec), 'w') as out:
        print(''.join(commentlines), file=out, end='')
        gff.to_csv(out, index=False, header=False, sep='\t')