#! /usr/bin/env python3

import pandas as pd
import numpy as np
import click
import re
from Bio import SeqIO, SeqRecord

@click.group()
def main():
    pass

def parse_gff(gffin):
    header = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df = pd.read_csv(gffin, comment='#', sep='\t', header=None)
    df.columns = header
    return df

@main.command()
@click.option('-l', '--upstream_length', type=int, default=150)
@click.option('-f', '--fastain', required=True, type=str)
@click.option('-g', '--gffin', required=True, type=str)
@click.option('-o', '--outfile', required=True, type=str)
def upstream(fastain, gffin, outfile, upstream_length):
    seqs = {rec.id:rec for rec in SeqIO.parse(fastain, 'fasta')}
    df = parse_gff(gffin)
    with open(outfile, 'w') as out:
        for row in df.iterrows():
            prot_id = re.search('ID=([A-Za-z0-9\._-]+);?', row[1]['attributes']).group(1)
            if row[1]['strand'] == '+':
                start = row[1]['start'] - upstream_length if row[1]['start'] >= upstream_length else 0
                end = row[1]['start']
                up = seqs[row[1]['seqname']].seq[start:end]
                id = "{}:{}:{}".format(row[1]['seqname'], start, end)
                up = SeqRecord.SeqRecord(id=id, seq=up, description="strand='+';")
            elif row[1]['strand'] == '-':
                start = row[1]['end']
                end = row[1]['end'] + upstream_length
                up = seqs[row[1]['seqname']].seq[start:end].reverse_complement()
                id = "{}:{}:{}".format(row[1]['seqname'], start, end)
                up = SeqRecord.SeqRecord(id=id, seq=up, description="strand='-';")
            up.description += f'downstream_prot_id={prot_id};'
            if len(up.seq) > 5:
                SeqIO.write(up, out, 'fasta')

@main.command()
@click.option('-f', '--fastain', required=True, type=str)
@click.option('-g', '--gffin', required=True, type=str)
@click.option('-o', '--outfile', required=True, type=str)
@click.option('--coding/--noncoding', default=False)
def noncoding(fastain, gffin, outfile, coding):
    seqs = [rec for rec in SeqIO.parse(fastain, 'fasta')]
    df = parse_gff(gffin)
    desc = 'noncoding region'
    with open(outfile, 'w') as out:
        for seq in seqs:
            sub = df[df['seqname'] == seq.id]
            mask = np.ones(len(seq.seq), bool)
            for row in sub.iterrows():
                mask[int(row[1]['start']):int(row[1]['end'])] = False
            if coding:
                mask = np.invert(mask)
                desc = 'coding region'
            z = np.concatenate(([False], mask, [False]))
            start = np.flatnonzero(~z[:-1] & z[1:])
            end = np.flatnonzero(z[:-1] & ~z[1:])
            intergenic = np.column_stack((start, end-1))
            for a,b in intergenic:
                id = "{}:{}:{}".format(seq.id, a, b)
                inter = SeqRecord.SeqRecord(id=id, seq=seq.seq[a:b], description=desc)
                SeqIO.write(inter, out, 'fasta')

if __name__ == '__main__':
    main()