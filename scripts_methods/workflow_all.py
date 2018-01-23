#Same analysis using gwf:

# This workflow is made usinf gwf and has as objective the execution of Rhizobium analysis. 
# The directory where this workflow is placed is ./NChain/faststorage/rhizobium/new_assemblies/genomes. 
# The outputs of all pipelines are transfered to the folder called 'data', so we won't have problem to address to these files.
# You also should initate this analyses with a folder for each strain containing the .fna files 

from gwf import Workflow
import glob
import os
import shutil

gwf = Workflow()


def prokka_run(directory, ID):
	outdir = f'{directory}/{ID}'
	infile = f'{directory}/{ID}.contigs.fna'
	inputs = [f'{directory}/{ID}.contigs.fna']
	outputs = [f'{outdir}/{ID}{x}' for x in ['.gff', '.gbk', '.fna', '.faa']]
	#print('\n'.join([outdir, infile, "".join(outfiles)]))
	options = {
		'memory': '4g',
		'cores': '4',
		'walltime': '10:00:00',
		'account': 'NChain'
	}

	spec = '''/project/NChain/faststorage/tools/prokka-1.12/bin/prokka --outdir {outdir} --prefix {ID} {infile}'''.format(ID=ID, directory=directory, infile=infile, outdir=outdir)
	
	print(spec)
	return inputs, outputs, options, spec

def moving_file(directory, ID, out_file):        
	outdir = f'{directory}/{ID}'
	inputs = [f'{outdir}/{ID}.faa', f'{outdir}/{ID}.gff']
	outputs = [f'{out_file}/{ID}{x}' for x in ['.faa', '.gff']]
        #print('\n'.join([outdir, infile, "".join(outfiles)]))
	
	options = {'memory': '2g',
		   'cores': '2',
		   'walltime': '00:00:10',
		   'account': 'NChain'}
	spec = '''cp {directory}/{ID}/{ID}.faa {directory}/{ID}/{ID}.gff {out_file}'''.format(ID=ID, directory=directory, out_file=out_file)
	
	print(spec)
	return inputs, outputs, options, spec

def proteinortho_run_step1(directory, inputs, project_name):
	outdir = f'{directory}'
	print(inputs)
	
	outputs_list = list()
	for ID in inputs:
		if ID[-4::] == '.faa':
			outputs = [f'{ID}{x}' for x in ['.phr', '.pin', '.psq']]
			outputs_list.append(outputs)
	options = {
		'memory': '8g',
		'cores': '16',
		'walltime': '00:30:00',
		'account': 'NChain'
	}

	spec = f'/project/NChain/faststorage/tools/proteinortho_v5.16_final/proteinortho5.pl -project={project_name} -cpus=16 -synteny -step=1 {directory}/*.faa'
	
	outputs = sum(outputs_list, [])
	print(outputs)
	return inputs, outputs, options, spec


def proteinortho_run_step2(directory, inputs, project_name, N, M = 200):

	# Creating the input files for 
	inputs = inputs
	outputs = [f'{project_name}.ffadj-graph_{N}_{M}']
	print(outputs)
	options = {
		'memory': '8g',
		'cores': '8',
		'walltime': '24:00:00',
		'account': 'NChain'
	}
	spec = f'/project/NChain/faststorage/tools/proteinortho_v5.16_final/proteinortho5.pl -project={project_name} -cpus=8 -synteny -step=2 -jobs={N}/{M} {directory}/*.faa'

	print(spec)
	return inputs, outputs, options, spec


def proteinortho_run_step3(directory, project_name, M):
	inputs1 = [f'{project_name}.ffadj-graph_{x}_{M}' for x in range(1,M,1)]
	inputs2 = [f'{project_name}.blast-graph_{x}_{M}' for x in range(1,M,1)]
	inputs = inputs1 + inputs2
	print(inputs)
	outputs = [f'{project_name}.proteinortho', f'{project_name}.proteinortho-graph']
	options = {
		'memory': '8g',
		'cores': '8',
		'walltime': '5:00:00',
		'account': 'NChain'
	}
	spec = f'/project/NChain/faststorage/tools/proteinortho_v5.16_final/proteinortho5.pl -project={project_name} -cpus=8 -synteny -step=3 {directory}/*.faa'	
	print(inputs)
	print(spec)
	return inputs, outputs, options, spec


def camous_pipeline():
	pass

def gene_alignments(directory, out_directory, gene, gene_name):
	inputs = [(f) for f in glob.glob(f'{directory}/*.fna')]
	outputs = [(f) for f in glob.glob(f'{out_directory}/*.fna')]
	
	options = {
                'memory': '1g',
                'cores': '4',
                'walltime': '00:10:00',
                'account': 'NChain'
	}
	spec = f'/home/mica16/anaconda2/envs/py36/bin/python codon_aware_clustal.py {gene} /project/NChain/faststorage/tools/clustalo > {out_directory}/{gene_name}'
	print(spec)
	return inputs, outputs, options, spec

# For each fasta file run prokka, many output files are generated: .gff, .gbk, .fna, .faa and so on.
# Find all the subdirectories existent in that directory
d = '.'
folders = [os.path.join(d, o) for o in os.listdir(d) 
					if os.path.isdir(os.path.join(d,o))]

inputs_next = list()
outputs_next = list()
for strain in folders: 
	print(strain)
	if strain[2] == '3':
		workflow_name = 'prokka{}'.format(strain[2:6])
		gwf.target_from_template(workflow_name, prokka_run(directory = strain, ID = strain[2::])) 
		inputs_next.append([f'{strain}/{strain[2::]}/{strain[2::]}{x}' for x in ['.gff', '.faa']])
		outputs_next.append([f'./data/{strain[2::]}{x}' for x in ['.gff', '.faa']]) # outputs are adressed to one unique file
		
		gwf.target_from_template(workflow_name+'movingt',moving_file(directory = strain, ID = strain[2::], out_file = './data/'))
flatted_input = sum(inputs_next, [])
flatted_output = sum(outputs_next, [])

## Running proteinortho in steps

# Running step 1
gwf.target_from_template('Index_creation', proteinortho_run_step1(directory = './data/', inputs=flatted_output, project_name = '201_strains_proteinortho'))
inputs = proteinortho_run_step1(directory = './data/', inputs=flatted_output, project_name = '201_strains_proteinortho')[1]
print(inputs)

# Running step 2
for part in range(1,10,1):
	gwf.target_from_template('proteinortho_step2_{}'.format(part), proteinortho_run_step2(directory = './data', inputs=inputs, project_name = '201_strains_proteinortho', N = part, M = 10))

# Running step 3
gwf.target_from_template('Graph_compilation_step3', proteinortho_run_step3(directory = './data', project_name = '201_strains_proteinortho', M = 10))


# Running Camous pipeline
# This pipeline must generate .fna files for each gene group
# Make the aligns directory as well (group_alns) 

genes = [(f) for f in glob.glob(f'./group_fnas/*.fna')]
## Running gene-alignments
for gene in genes:
	workflow_name = gene.split('/')[2]
	print(workflow_name)
	gwf.target_from_template('alignments_{}'.format(workflow_name), gene_alignments(directory = './group_fnas', gene = gene, out_directory = './group_alns', gene_name = workflow_name))













