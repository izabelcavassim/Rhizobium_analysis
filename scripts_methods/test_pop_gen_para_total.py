import dendropy
from dendropy.calculate import popgenstat
import glob as glob
import math
import pandas as pd
import numpy as np
from GFFReader import gff_iter 
from collections import Counter

def parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
	from itertools import izip   
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, origin, country, origin2 in izip(t['Seq ID'], t['Genospecies'], t['Country'], t['Origin2']):
		pop_map[str(strain_id)]={'genospecies':origin, 'country':country, 'origin':origin2}
	return pop_map

def invert_dict(d): 
	inverse = dict() 
	for key in d: 
		# Go through the list that is saved in the dict:
		for item in d[key]:
			# Check if in the inverted dict the key exists
			if item not in inverse: 
				# If not create a new list
				inverse[item] = [key] 
			else: 
				inverse[item].append(key) 
	return inverse

def gene_groups_maker():
	proteinortho_directory = '/Users/PM/Desktop/New_data/final_project.poff_disambiguated.groups'
	gene_groups = {}
	with open(proteinortho_directory) as f:
		for line in f:
			temp = line.strip('\n')
			temp = temp.split(':')
			group_name = temp[0]
			group_name = group_name.split('|')[0]
			members = temp[1].split(' ')
			members = members[1::]
			new_members = list()
			for i in members:
				try:
					new_members.append(i.split('|')[1]) # changing here to accept the new strains
				except:
					continue
			if group_name not in gene_groups:
				gene_groups[group_name] = new_members # creating a list

	gene_groups_inv = invert_dict(gene_groups)
	np.save("gene_groups_dict_inverted.npy", gene_groups_inv)

#gene_groups_maker()
# seqs = dendropy.DnaCharacterMatrix.get_from_string(seqstr, 'nexus')
# taxon_namespace = seqs.taxon_namespace
# print taxon_namespace

# print dict_map['3404']
# tax_pop_label_map = {}
# for t in taxon_namespace:
#     tax_pop_label_map[t] = dict_map[t.label]['genospecies']

# print tax_pop_label_map
# tax_parts = taxon_namespace.partition(membership_dict=tax_pop_label_map)

# for s in tax_parts.subsets():
#     print(s.description())

def tajimas_d_hand(num_sequences, avg_num_pairwise_differences, num_segregating_sites):

	### VERIFICATION ###
	###
	### Given: num_sequences = 10, num_pairwise_differences = 3.888889, num_segregating_sites = 16
	###  i.e.: tajimas_d(10, 3.888889, 16)  == -1.44617198561
	###  Then:    a1 == 2.82896825397
	###           a2 == 1.53976773117
	###           b1 == 0.407407407407
	###           b2 == 0.279012345679
	###           c1 == 0.0539216450284
	###           c2 == 0.0472267720013
	###           e1 == 0.0190605338016
	###           e2 == 0.0049489277699
	###           D ==  -1.44617198561

	a1 = sum([1.0/i for i in range(1, num_sequences)])
	a2 = sum([1.0/(i**2) for i in range(1, num_sequences)])
	b1 = float(num_sequences+1)/(3*(num_sequences-1))
	b2 = float(2 * ( (num_sequences**2) + num_sequences + 3 )) / (9*num_sequences*(num_sequences-1))
	c1 = b1 - 1.0/a1
	c2 = b2 - float(num_sequences+2)/(a1 * num_sequences) + float(a2)/(a1 ** 2)
	e1 = float(c1) / a1
	e2 = float(c2) / ( (a1**2) + a2 )
	try:
		D = (
			float(avg_num_pairwise_differences - (float(num_segregating_sites)/a1))
			/ math.sqrt(
				(e1 * num_segregating_sites )
			+ ((e2 * num_segregating_sites) * (num_segregating_sites - 1) ))
			)
		return D
	except ZeroDivisionError:
		return "nan"
		#return 0

def wattersons_theta_hand(char_matrix, num_segregating_sites):
	"""
	Returns Watterson's Theta (per sequence)
	"""
	sequences = char_matrix.sequences()
	a1 = sum([1.0/i for i in range(1, len(sequences))])

	if num_segregating_sites != 0:
		return float(num_segregating_sites) / a1
	else:
		return 'nan'


def do_basic_popgen(seqs):
	sequences = seqs.sequences()
	num_sequences = len(sequences)
	num_seg_sites = popgenstat.num_segregating_sites(seqs)
	avg_pair = popgenstat.average_number_of_pairwise_differences(seqs)
	nuc_div = popgenstat.nucleotide_diversity(seqs)
	nuc_watterson = wattersons_theta_hand(char_matrix = seqs, num_segregating_sites = num_seg_sites)
	tajimas_D = tajimas_d_hand(num_sequences = num_sequences, avg_num_pairwise_differences = avg_pair, num_segregating_sites = num_seg_sites)
	
	return([num_sequences, num_seg_sites, avg_pair, nuc_div, nuc_watterson, tajimas_D])

def pop_stats():

	gene_alignments = [f for f in glob.glob("/Users/PM/Desktop/New_data/group_alns/*.fna")]
	
	#strain_158 = pd.read_table('gene_locations_sorted_sm158.txt', sep = '\t', header = 0)
	#strain_158_genes = strain_158['qseqid'].tolist()
	
	results = []
	count = 0 
	for gene in gene_alignments:
		gene_name = gene.split('/')[6]
		gene_name = gene_name.split('.fna')[0]
		
		count += 1
		print 'Gene, %d' % count
		seqs = dendropy.DnaCharacterMatrix.get(path=gene,schema="fasta")
		stats = do_basic_popgen(seqs)
		print stats
		stats.append(gene_name)
		#print stats
		results.append(stats)
	return results

def plasmid_origin(gff, save_directory, gene_groups, dict_gene, save = False):
	'''Look at the most often plasmid origin, if the plasmids are equally probable, take the first'''
	

	strain_name = gff.split('/')[6]
	#print strain_name

	# Making a pandas dataframe
	plasmid = []
	start = []
	end = []
	group = []
	strand = []
	gene_product = []

	for entry in gff_iter(gff, type_filter=["CDS", "gene", "product"]):
		try:
			group.append(gene_groups[entry.attributes["ID"]][0])

			#  Splitting plasmid group:
			#plasmid =  entry.sequence.split('_')[3]

			if gene_groups[entry.attributes["ID"]][0] not in dict_gene:

				dict_gene[gene_groups[entry.attributes["ID"]][0]] = [entry.attributes["product"]]
			else:
				dict_gene[gene_groups[entry.attributes["ID"]][0]].append(entry.attributes["product"])

			plasmid.append(entry.sequence)
			start.append(entry.start)
			end.append(entry.end)
			strand.append(entry.strand)
			gene_product.append(entry.attributes["product"])

		
		except:
			pass


	if save == True:
		pandas_gff_ortholous_group = pd.DataFrame(
			{'plasmid': plasmid,
			'start': start,
			'end': end,
			'gene_group': group,
			'strand': strand,
			'product': gene_product
			})

		column_order = ['gene_group', 'plasmid', 'start', 'end', 'strand', 'product']
		pandas_gff_ortholous_group[column_order].to_csv(save_directory+strain_name)
	#sequence id, start and end positions, and ID of CDS' and genes in file.gff

	return(dict_gene)
	# Parsing the gff files
gff_files = [(f) for f in glob.glob('/Users/PM/Desktop/New_data/gff_files/*.gff')]

# Parsing the genegroup
dict_gene = {}
gene_groups = np.load("gene_groups_dict_inverted.npy").item()
for f in gff_files:
	results = plasmid_origin(gff = f, save_directory = '/Users/PM/Desktop/New_data/gwas_gff/', gene_groups = gene_groups, dict_gene = dict_gene, save = False)

for gene_group in results:
	print '{}, {}'.format(gene_group, Counter(results[gene_group]).most_common(1)[0][0]) # If there is more than one mode this returns an arbitrary one

#test = pop_stats()
#header = ['num_sequences', 'num_seg_sites', 'avg_pair', 'nuc_div', 'theta', 'tajimas_D', 'gene_group' ] #, 'plasmid', 'function']
#df = pd.DataFrame(test, columns=header)
#print df

#df.to_csv('total_statistics_all_new_fasta_files.csv', sep='\t')
