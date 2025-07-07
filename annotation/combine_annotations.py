#! /usr/bin/env python3

from Bio import SeqIO
from Bio.SearchIO import HmmerIO
import pandas as pd
import click

@click.command()
@click.option('-f', '--faa', type=str, required=True)
@click.option('-ips', '--interproscan', type=str)
@click.option('-g', '--gvogs', type=str)
@click.option('-ga', '--gvog_annotations', type=str)
@click.option('-e', '--emapper', type=str)
@click.option('-t', '--tax_nr', type=str)
@click.option('-o', '--outfile', type=str, required=True)
@click.option('--xlsx', is_flag=True, default=False)

def main(faa, gvogs, gvog_annotations, emapper, tax_nr, interproscan, outfile, xlsx):
    recdict = SeqIO.to_dict(SeqIO.parse(faa,'fasta'))
    
    columns = ['seqid', 'database', 'accession', 'description', 'tag', 'start', 'stop', 'evalue']
    annots = [pd.DataFrame(columns=columns)]

    if emapper:
        annots.append(parse_emapper(emapper))
    if interproscan:
        annots.append(parse_interproscan(interproscan))
    if tax_nr:
        annots.append(parse_tax_nr(tax_nr))
    if gvogs:
        annots.append(parse_gvog(gvogs, gvog_annotations))

    annotations = pd.concat(annots, axis=0)
    annotations = annotations[columns]

    annotations.loc[annotations['stop'].isna(), 'stop'] = annotations[annotations['start'].isna()].apply(lambda x: len(recdict[x['seqid']].seq), axis=1)
    annotations.loc[annotations['start'].isna(), 'start'] = 0
    
    annotations = annotations.astype({'start':int, 'stop':int, 'seqid':str})
    annotations.sort_values(by='seqid', axis=0, inplace=True)

    if xlsx:
        annotations.to_excel(outfile, index=False, header=True)
    else:
        annotations.to_csv(outfile, index=False, header=True, sep='\t')

    for f in recdict.keys():
        if not f in annotations['seqid'].values:
            print(F'found no annotation for {f}')

def parse_interproscan(interproscan):
    interproscan = pd.read_csv(interproscan, sep='\t', header=None)
    cols = ['seqid', 'hash', 'len', 'database', 'accession', 'description', 'start', 
            'stop', 'evalue', 'status', 'date', 'ips accession', 'ips description']
    interproscan.columns = cols
    return interproscan


def parse_emapper(emapper):
    emapper = pd.read_csv(emapper, skiprows=2, skipfooter=3, header=2, sep='\t', engine='python')
    emapper.rename(columns={'#query':'seqid', 'Description':'description', 'Preferred_name':'tag',
                            'eggNOG_OGs':'accession'}, inplace=True)
    emapper['database'] = 'eggNOG'
    return emapper


def parse_tax_nr(tax_nr):
    tax_nr = pd.read_csv(tax_nr, sep='\t')
    tax_nr.fillna('', inplace=True)
    tax_nr['description'] = tax_nr.apply(lambda x: ','.join([l for l in x[['superkingdom', 
                'phylum', 'class','order', 'family', 'genus', 'species']] if l]), axis=1)
    tax_nr.rename(columns={'prot_id':'seqid', 'taxid':'accession', 'e-value':'evalue'},
                  inplace=True)
    tax_nr['database'] = 'NCBI NR/Taxonomy'
    return tax_nr


def parse_gvog(gvog, gvog_annotations):
    gvog_df = pd.DataFrame(columns=['seqid', 'evalue', 'accession', 'description'])
    gvog_annotations = pd.read_csv(gvog_annotations, sep='\t')
    gvog_annotations.index = gvog_annotations['GVOG']
    gvog_annotations.fillna('', inplace=True)
    with open(gvog) as infile:
        for rec in HmmerIO.Hmmer3TabParser(infile):
            for hit in rec.hits:
                description = ','.join(set(gvog_annotations.loc[hit.query_id, 'NCVOG_descs'].split(' | ')))
                ncvog_accessions = ','.join(set(gvog_annotations.loc[hit.query_id, 'consensus_NCVOG'].split(' | ')))
                if ncvog_accessions:
                    description += f';{ncvog_accessions}'
                gvog_df = pd.concat([gvog_df, pd.DataFrame([[hit.id, hit.evalue, hit.query_id, description]], 
                                                           columns=list(gvog_df.columns))], ignore_index=True)
    gvog_df['database'] = 'GVOG/NCVOG'
    return gvog_df


if __name__ == '__main__':
    main()