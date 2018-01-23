#Same analysis using gwf:

# This workflow is made usinf gwf and has as objective the execution of Rhizobium analysis. 
# The directory where this workflow is placed is ./NChain/faststorage/rhizobium/new_assemblies/genomes. 
# The outputs of all pipelines are transfered to the folder called 'data', so we won't have problem to address to these files.


from gwf import Workflow
import glob
import os

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


def proteinortho_run():

# For each fasta file run prokka, many output files are generated: .gff, .gbk, .fna, .faa and so on.
# Find all the subdirectories existent in that directory
d = '.'
folders = [os.path.join(d, o) for o in os.listdir(d) 
					if os.path.isdir(os.path.join(d,o))]

for strain in folders: 
	print(strain)
	if strain[2] == '3':
		workflow_name = 'prokka{}'.format(strain[2:6])
		gwf.target_from_template(workflow_name, prokka_run(directory = strain, ID = strain[2::])) 

		inputs_next = [f'{strain}/{strain[2::]}{x}' for x in ['.gff', '.faa']
		print inputs_next


#gwf.target_from_template('Moving_files', move_faa_files(directory = strain, ID = strain[2::], new_dir = './data2'))


