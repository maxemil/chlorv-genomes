#! /usr/bin/env python3
import ete3
import pandas as pd
ncbi = ete3.NCBITaxa()
import rich_click as click

@click.command()
@click.option('-i', '--diamond_tax', type=str)
@click.option('-o', '--outfile', type=str)

def main(diamond_tax, outfile):
    df = parse_diamond_tax(diamond_tax)
    df.to_csv(outfile, sep='\t', index=False)

def check_lineage(t, tax_level):
    if ncbi.get_rank([t])[t] == tax_level:
        return ncbi.get_taxid_translator([t])[t]

def get_tax_level(t, tax_level):
    tax = None
    if t == 0:
        return tax
    try:
        for l in ncbi.get_lineage_translator([t])[t]:
            tax = check_lineage(l, tax_level)
            if tax:
                return tax
    except:
        print('failed lineage')
    return tax

def parse_diamond_tax(infile):
    df = pd.read_csv(infile,  sep='\t', header=None)
    df.columns = ['prot_id', 'taxid', 'e-value']
    for rank in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        df[rank] = df['taxid'].apply(lambda x: get_tax_level(x, rank))
    return df


if __name__ == '__main__':
    main()

# df.loc[df['taxid'] > 0, 'tax'] = 'root'
# df.loc[(df['tax'] == 'root') & (df['superkingdom']), 'tax'] = df.loc[(df['tax'] == 'root') & (df['superkingdom']), 'superkingdom']
# df.loc[(df['superkingdom'] == 'Viruses') & (df['class']), 'tax'] = df.loc[(df['superkingdom'] == 'Viruses') & (df['class']), 'class']
# df.loc[(df['class'] == 'Megaviricetes') & (df['order']), 'tax'] = df.loc[(df['class'] == 'Megaviricetes') & (df['order']), 'order']
# df['tax'].value_counts()
