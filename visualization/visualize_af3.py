import pandas as pd
import glob
from pymol import cmd,CmdException
from Bio import PDB
import os
from Bio import SeqIO

with open('header.md', 'w') as out:
    initialize_md(out)

seqdict = SeqIO.to_dict(SeqIO.parse('ChlorV-1.annotated.faa', 'fasta'))
# seqdict = SeqIO.to_dict(SeqIO.parse('ChlorV-234.faa', 'fasta'))

diff_index = {'chlorv-1..071':2, 'chlorv-1..119':1, 'chlorv-1..133':1, 'chlorv-1..163':1, 'chlorv-1..202':2, 
              'chlorv-1..378':1, 'chlorv-1..379':2}
manual = {'chlorv-1..190':0, 'chlorv-1..292':0, 'chlorv-1..359':0}

for foldout in sorted(glob.glob('*/*.foldout')):
    if glob.glob(f"{os.path.basename(foldout).replace('.foldout', '')}*.png"):
        print(f"Skipping {foldout}, already processed.")
        continue
    index = 0
    # Manually set the index for the foldseek file in some cases where the best hit doesn't have a good structure
    if foldout.split('/')[0] in diff_index.keys():
        index = diff_index[foldout.split('/')[0]]
    elif foldout.split('/')[0] in manual.keys():
        continue
    rec = seqdict[foldout.split('/')[0].replace('chlorv', 'ChlorV')]
    seq_annot = rec.description.replace(f'{rec.id}', '')
    df = pd.read_csv(foldout, sep='\t')
    df.sort_values(by=['evalue'], axis=0, ascending=True, inplace=True, ignore_index=True)
    print(foldout)
    if not df.size == 0 and df.loc[index]['evalue'] < 0.0001:
        pdb_file = df.loc[index]['target'].split('_')[0].split('-')[0]
        chain = df.loc[index]['target'].split('_')[1].split('-')[0]
        chain_description, png_file = visualize_alignment(foldout.replace('.foldout', '_model.cif'), pdb_file, chain, df.loc[index])
        with open(png_file.replace('png', 'md'), 'w') as out:
            write_md(os.path.basename(foldout).replace('.foldout', ''), df.loc[0:2], pdb_file, chain, chain_description, seq_annot, png_file, out)
    else:
        png_file = visualize_structure(foldout.replace('.foldout', '_model.cif'))
        with open(png_file.replace('png', 'md'), 'w') as out:
            write_md_nohit(os.path.basename(foldout).replace('.foldout', ''), seq_annot, png_file, out)

def initialize_md(md_file):
    print("---", file=md_file)
    print("geometry: margin=20mm,landscape", file=md_file)
    print("...", file=md_file)

def write_md_nohit(model, seq_annot, png_file,  md_file):
    print('\\newpage', file=md_file)
    print(f"# {model}", file=md_file)
    print(f"* Sequence-based annotation for {model} is {seq_annot} \n", file=md_file)
    print(f"* No significant structural hit found \n", file=md_file)
    print(f'![predicted structure of {model}]({png_file} "{png_file}"){{ width=50% }}\n', file=md_file)

def write_md(model, df, pdb_file, chain, chain_description, seq_annot, png_file, md_file):
    print('\\newpage', file=md_file)
    print(f"# {model}", file=md_file)
    print(f"* Sequence-based annotation for {model} is {seq_annot}\n", file=md_file)
    print(f"* Best hit was {pdb_file} chain {chain}: {chain_description}\n", file=md_file)
    sub_df = df[['target', 'prob', 'fident', 'alnlen', 'evalue', 'theader']]
    sub_df.loc[:,'theader'] = sub_df.theader.apply(lambda x: x.partition(' ')[2])
    print(f"{sub_df.to_markdown(buf=None, index=False)}", file=md_file)
    print(f'![left: reference structure of {pdb_file} chain {chain}. right: predicted structure of {model}, unaligned sequences are shown as transparent]({png_file} "{png_file}"){{ width=80% }}\n', file=md_file)

def color_plddt(selection="all"):
	# Define the pLDDT bind boundaries
	bin_lower = [0,	50,	70,	90]

	bin_upper = [50, 70, 90, 100]

	n_colors = 4
	colors = [[0.8784,0.5059,0.3333],
	[0.9765,0.8627,0.3020],
	[0.4980,0.7843,0.9294],
	[0.2157,0.3373,0.6157]]
     
	# Loop through color intervals
	for i in range(n_colors):
		lower = bin_lower[i]
		upper = bin_upper[i]
		color = colors[i]

		# Define a unique name for the atoms which fall into this group
		group = selection + "_plDDT_" + str(lower) + "_to_" + str(upper)

		# Compose a selection command which will select all atoms which are
		#	a) in the original selection, AND
		#	b) have B factor in range lower <= b < upper
		sel_string = selection + " & ! b < " + str(lower)

		if(i < n_colors - 1):
			sel_string += " & b < " + str(upper)
		else:
			sel_string += " & ! b > " + str(upper)

		# Select the atoms in bin
		cmd.select(group, sel_string)

		# Create a name for the color
		color_name = "color_" + str(i+1)
		cmd.set_color(color_name, color)

		# Color the atoms in the selected group
		cmd.color(color_name, group)

def visualize_structure(model):
    cmd.reinitialize()

    modelname = os.path.basename(model).replace('_model.cif', '')

    modelpath = os.path.abspath(model)
    cmd.load(modelpath, modelname)

    color_plddt(selection=f'{modelname}')

    cmd.set('cartoon_discrete_colors', 1)
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_trace_color', 'black')
    cmd.set('ray_shadow', 0)
    cmd.set('surface_quality', 1)
    cmd.set('ambient', 0.5)
    cmd.set('ray_opaque_background', 'off')
    cmd.set('opaque_background', 'off')
    cmd.set('cartoon_transparency', 0)

    cmd.set('grid_mode', 1)
    
    cmd.hide('all')
    cmd.show('cartoon', f'model {modelname}')
    cmd.png(f'{modelname}.png', width='10cm', dpi=300)
    return f'{modelname}.png'

def visualize_alignment(model, ref, chain, aln):
    cmd.reinitialize()

    modelname = os.path.basename(model).replace('_model.cif', '')

    cmd.fetch(ref)
    refname = ref
    refchain = f"{refname}_{chain}"

    cmd.split_chains(refname)

    for f in cmd.get_object_list():
        if not f.lower() == refchain.lower() and not f.lower() == modelname.lower():
            cmd.delete(name=f)

    if refchain.lower() == cmd.get_object_list()[0].lower():
        refchain = cmd.get_object_list()[0]

    cif = PDB.MMCIF2Dict.MMCIF2Dict(f'{ref}.cif')
    for i, s in enumerate(cif['_entity_poly.pdbx_strand_id']):
        if chain in s.split(','):
            chainindex = i
    chain_description = cif['_entity.pdbx_description'][chainindex]

    cmd.remove('solvent or organic or inorganic')

    modelpath = os.path.abspath(model)
    cmd.load(modelpath, modelname)

    cmd.select('modsel', selection=f'i. {aln.qstart}-{aln.qend} and model {modelname}')
    selected_atoms = cmd.select('refsel', selection=f'resi {aln.tstart}-{aln.tend} and model {refchain} and name CA')

    if selected_atoms < aln['tend']-aln['tstart']+1:
        print('residue numbering is likely off')
        myspace = {'firstresi': []}
        cmd.iterate(f"first {refchain} and name CA", "firstresi.append(int(resi))", space=myspace)
        cmd.alter(refchain, f'resi=int(resi)+1-{myspace['firstresi'][0]}')
        # cmd.iterate(f"{refchain} and name CA", "int(resi), resi, resn")
        selected_atoms = cmd.select('refsel', selection=f'i. {aln.tstart}-{aln.tend} and model {refchain}')
    cmd.rebuild()

    try:
        cmd.alignto(modelname, method='cealign', object='aln')
    except pymol.CmdException:
        cmd.alignto(modelname, method='cealign', object='aln', window=3)

    color_plddt(selection=f'{modelname}')

    color_scheme = ['0xE5F3FF', '0xf8d7b2', '0xebecf0']
    cmd.color(color_scheme[0], f'ss h and model {refchain}')
    cmd.color(color_scheme[1], f'ss s and model {refchain}')
    cmd.color(color_scheme[2], f"ss l+'' and model {refchain}")


    cmd.set('cartoon_discrete_colors', 1)
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_trace_color', 'black')
    cmd.set('ray_shadow', 0)
    cmd.set('surface_quality', 1)
    cmd.set('ambient', 0.5)
    cmd.set('ray_opaque_background', 'off')
    cmd.set('opaque_background', 'off')
    cmd.set('cartoon_transparency', 0.75)
    cmd.set('cartoon_transparency', 0, 'refsel')
    # cmd.set('cartoon_transparency', 0, f'model {refchain}')
    cmd.set('cartoon_transparency', 0, 'modsel')
    
    cmd.set('grid_mode', 1)
    cmd.reset()    

    cmd.hide('all')
    cmd.show('cartoon', f'model {modelname}')
    cmd.show('ribbon', f'model {refchain}')
    cmd.show('cartoon', f'model {refchain}')
    cmd.png(f'{modelname}_{refchain}.png', width='10cm', height='5cm', dpi=300)
    return chain_description, f'{modelname}_{refchain}.png'