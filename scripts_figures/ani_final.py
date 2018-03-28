import numpy as np
import pandas as pd

def parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
	from itertools import izip
	
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, sara_id, origin, country in izip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Origin2']):
		pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country, 'origin2':origin}
	return pop_map

def parse_fasta(filename):
	file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
	file_separe = file.split('>') #spliting each entry by the > 
	file_separe.remove('')
	parse_dict = {}
	header = []
	for entry in file_separe:
		seq = entry.splitlines()
		header = seq[0] #these are the first elements of the list 
		seq = ''.join(seq[1:]) #joining the sequences 
		parse_dict[header] = seq
	return parse_dict


def compare_equal(array1, array2):
	
	N = len(array1)
	count = 0
	for i in xrange(N):
		if array1[i]==array2[i] and array1[i] != 5 and array2[i] != 5:
			count += 1.0
	ani = count/float(N)
	return ani


fasta = parse_fasta('/Users/PM/Desktop/New_data/concatenated_genes.fasta')

pop = parse_pop_map()
pop_map = pop.keys()
# Names
strains_names = sorted(pop_map, key=lambda x: pop[x]['genospecies'])
print 'These are the strains evaluated', strains_names
strains_names.remove('3260')
strains_names.remove('3381')
strains_names.remove('3339')
strains_names.remove('3211')


# Making the matrix to update the ani 
ani_matrix = np.zeros((196,196))
ani_matrix = pd.DataFrame(ani_matrix, index = strains_names, columns = strains_names)

for ind1 in fasta.keys(): # sequences are in the first entry
	for ind2 in fasta.keys():
		ani_matrix[ind1[0:4]][ind2[0:4]] += compare_equal(fasta[ind1], fasta[ind2])


print ani_matrix

ani_matrix.to_csv('ani_sorted_by_genospecies_test.csv', header = True)


