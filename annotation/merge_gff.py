#! /usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import click
from BCBio import GFF

@click.command()
@click.option('-p', '--prodigal', type=str)
@click.option('-g', '--genemark', type=str)
@click.option('-f', '--fasta', type=str)
@click.option('-s', '--sep', type=str, default='_')
@click.option('-l', '--length', type=int, default=100)


def main(prodigal, genemark, fasta, sep, length):
    recdict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    assert len(recdict) == 1
    feats_selected = compare_select_feats(recdict, prodigal, genemark)
    feats_filtered = filter_length(feats_selected, length)
    rename_write_gff(recdict, feats_filtered, sep)

def compare_select_feats(recdict, prodigal, genemark):
    feats_prodigal = parse_gff_feats(prodigal, recdict)
    feats_genemark = parse_gff_feats(genemark, recdict)

    rm_prodigal = set()
    rm_genemark = set()
    feats_selected = set()

    for f in feats_prodigal:
        for g in feats_genemark:
            if g.location.strand != f.location.strand:
                continue
            if f.location == g.location:
                f.qualifiers['source'][0] += f';{g.qualifiers["source"][0]}'
                feats_selected.add(f)
                rm_prodigal.add(f)
                rm_genemark.add(g)
            elif f.location.start == g.location.start:
                if f.location.end > g.location.end:
                    feats_selected.add(f)
                    rm_prodigal.add(f)
                    rm_genemark.add(g)
                if f.location.end < g.location.end:
                    rm_prodigal.add(f)
                    rm_genemark.add(g)
                    feats_selected.add(g)
            elif f.location.end == g.location.end:
                if f.location.start < g.location.start:
                    rm_prodigal.add(f)
                    rm_genemark.add(g)
                    feats_selected.add(f)
                if f.location.start > g.location.start:
                    feats_selected.add(g)
                    rm_prodigal.add(f)
                    rm_genemark.add(g)
            elif f.location.start < g.location.start and f.location.end > g.location.end:
                    feats_selected.add(f)
                    rm_prodigal.add(f)
                    rm_genemark.add(g)
            elif f.location.start > g.location.start and f.location.end < g.location.end:
                    feats_selected.add(g)
                    rm_prodigal.add(f)
                    rm_genemark.add(g)
    [feats_prodigal.remove(r) for r in rm_prodigal]
    [feats_genemark.remove(r) for r in rm_genemark]
    feats_selected = sorted(list(feats_selected) + feats_prodigal + feats_genemark, key=lambda x: x.location.start)
    return feats_selected

def filter_length(feats_selected, length):
    feats_filtered = []
    for f in feats_selected:
        if f.location.end - f.location.start >= length * 3:
            feats_filtered.append(f)
    return feats_filtered

def rename_write_gff(recdict, feats, sep):
    prefix, rec = recdict.popitem()
    for i, f in enumerate(feats):
        f.id = "{0}{2}{1:03d}".format(prefix, i+1, sep)
        f.qualifiers['ID'] = f.id
        if 'Parent' in f.qualifiers:
            del f.qualifiers['Parent']
    rec.features = feats
    with open(f'{prefix}.merge.gff', 'w') as out:
        GFF.write([rec], out)

def parse_gff_feats(gff, recdict):
    gff_iter = GFF.parse(gff, recdict, limit_info=dict(gff_type=["CDS"]))
    feats = []
    for rec in gff_iter:
        for f in rec.features:
            feats.append(f)
    return feats

if __name__ == '__main__':
    main()