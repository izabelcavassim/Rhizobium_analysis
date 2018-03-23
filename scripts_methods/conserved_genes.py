
#blastn -db concatenated_core_genes.fasta -query conserved_genes.fasta -evalue 0.1 -outfmt 6 -max_target_seqs 1 -max_hsps 1 > conserved_genes_conversion.txt

# Concatenate conserved genes, estimate ANI across species
import glob
import pandas as pd
import collections
from sys import argv
import pandas as pd
import matplotlib.pyplot as plt

def parse_fasta(filename):
	file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
	file_separe = file.split('>') #spliting each entry by the > 
	file_separe.remove('')
	parse_dict = collections.OrderedDict()
	header = []
	for entry in file_separe:
		seq = entry.splitlines()
		header = seq[0] #these are the first elements of the list 
		#header = header.split('|')
		header = header.split('|')[1]
		seq = ''.join(seq[1:]) #joining the sequences 
		parse_dict[header] = seq
	return parse_dict


blast_results_conserved_genes = pd.read_table(argv[2], header=None)

#print blast_results_conserved_genes

default_outfmt6_cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.strip().split(' ')
blast_results_conserved_genes.columns = default_outfmt6_cols


list_conserved_genes = blast_results_conserved_genes['sseqid'].tolist()
#print len(list_conserved_genes)
conserved_core_genes  = list()
not_core = list()
for i in list_conserved_genes:
	gene_name = i.split('|')[0]
	#print gene_name
	gene_temp = parse_fasta(argv[1]+"/{}.fna".format(gene_name))
	if len(gene_temp.keys()) == 201:
		conserved_core_genes.append(gene_temp)
	else:
		not_core.append(gene_temp.keys())

print len(conserved_core_genes)

dd = collections.defaultdict(list)
for d in conserved_core_genes: # you can list as many input dicts as you want here
      for key, value in d.iteritems():
         dd[key].append(value)

with open('{}'.format(argv[3]), 'w') as data:
	for i in sorted(dd.keys()):
	  	data.write('>{}\n'.format(i))
	  	sequences = ''.join(dd[i])
		print len(sequences)
	  	chunks = [sequences[i:i+60] for i in range(0, len(sequences), 60)]
	  	for i in chunks:
	  		data.write('{}\n'.format(i))

histogram_strains_missing = dict()
for i in not_core:
	temp_missing = list(set(conserved_core_genes[0].keys()) - set(i))
	for strain in temp_missing:
		if strain not in histogram_strains_missing:
			histogram_strains_missing[strain] = 1
		else:
			histogram_strains_missing[strain] += 1

print histogram_strains_missing
for key, values in histogram_strains_missing.items():
	print "{} {}\n".format(key, values)
