#Same analysis using gwf:

# This workflow is made usinf gwf and has as objective the execution of Rhizobium analysis. 
# The directory where this workflow is placed is ./NChain/faststorage/rhizobium/new_assemblies/genomes. 
# The outputs of all pipelines are transfered to the folder called 'data', so we won't have problem to address to these files.


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

	spec = f"/project/NChain/faststorage/tools/prokka-1.12/bin/prokka --outdir {outdir} --prefix {ID}  {infile}"
	print(spec)

	return inputs, outputs, options, spec

def proteinortho_run_step1(directory, strain_list, project_name):
	outdir = f'{directory}/'

	for ID in strain_list:
		infile = [f'{outdir}/{ID[2::]}{x}' for x in ['.faa', '.gff']]
		outputs = [f'{outdir}/{ID[2::]}{x}' for x in ['.phr', '.pin', '.psq']]
	options = {
		'memory': '8g',
		'cores': '16',
		'walltime': '00:30:00',
		'account': 'NChain'
	}

	spec = f'/project/NChain/faststorage/tools/proteinortho_v5.16_final/proteinortho5.pl -project={project_name} -cpus=16 -synteny -step=1 {directory}/*.faa'

	return inputs, outputs, options, spec


def proteinortho_run_step2(directory, strain_list, project_name, N, M = 200):

	for ID in strain_list:
		inputs = [f'{outdir}/{ID[2::]}{x}' for x in ['.phr', '.pin', '.psq']]
		outputs = []
	options = {
		'memory': '8g',
		'cores': '8',
		'walltime': '24:00:00',
		'account': 'NChain'
	}
	spec = "/project/NChain/faststorage/tools/proteinortho_v5.16_final/proteinortho5.pl -project={project_name} -cpus=8 -synteny -step=2 -jobs={N}/200 {directory}/*.faa"


	return inputs, outputs, options, spec


def proteinortho_run_step3(directory, strain_list, project_name, N, M):
	myproject.blast-graph_186_200
	inputs = [f'{directory}/{project_name}.blast-graph_{N}_{M}' for x in range(1,M,1)]
	outputs = [f'{project_name}.proteinortho', f'{project_name}.proteinortho-graph']
	options = {
		'memory': '8g',
		'cores': '8',
		'walltime': '4:00:00',
		'account': 'NChain'
	}
	spec = "/project/NChain/faststorage/tools/proteinortho_v5.16_final/proteinortho5.pl -project={project_name} -cpus=8 -synteny -step=3 {directory}/*.faa"	

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
		outputs_next.append([f'./data2/{strain[2::]}{x}' for x in ['.gff', '.faa']]) # outputs are adressed to one unique file

flatted_input = sum(inputs_next, [])
flatted_output = sum(outputs_next, [])

print(flatted_input)
print(flatted_output)
gwf.target('Moving_files_test', inputs= flatted_input, outputs = []) << """ ./manipulating_folders.sh """


# Running proteinortho in steps
gwf.target_from_template('Index_creation', proteinortho_step_1(directory = './data/' strain_list = d, project_name = '201_strains_proteinortho'))


for part in range(1,200,1):
	gwf.target_from_template('proteinortho_step2_{}'.format(part), proteinortho_step_2(directory = './data/' strain_list = d, project_name = '201_strains_proteinortho', N = part, M = 200))


gwf.target_from_template('Graph_compilation_step3', proteinortho_step_3(directory = './data/' strain_list = d, project_name = '201_strains_proteinortho', N = part, M = 200))






