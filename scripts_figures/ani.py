import scipy as sp
import h5py 
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('PDF')
import pylab as pl
import time
import matplotlib.pyplot as plt
from scipy import stats
#import seaborn as sns

def parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
	from itertools import izip   
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, origin, country in izip(t['Seq ID'], t['Genospecies'], t['Country']):
		pop_map[str(strain_id)]={'genospecies':origin, 'country':country}
	return pop_map


def compare_equal(array1, array2):
	N = len(array1)
	count = 0
	for i in xrange(N):
            if array1[i]==array2[i] and array1[i] != 5:
                count += 1
	ani = count/float(N)

	return ani

def average_nucleotide_identity(fig_dir = '/Users/PM/Desktop/New_Data/final_figures/',
								gt_hdf5_file= '/Users/PM/Desktop/New_data/final_snps.hdf5',
								conserved_genes = '/Users/PM/Desktop/New_data/conserved_core_genes.txt'):

	pop = parse_pop_map()
	pop_map = pop.keys()
	ct_array = pop.values()

	text_file = open(conserved_genes, "r").read() 
	file_separe = text_file.split('\n')
	concatenated_genes = file_separe


	#print concatenated_genes 

	h5f = h5py.File(gt_hdf5_file)
	print h5f.keys()
	ag = h5f['alignments']
	gene_big_groups = sorted(ag.keys())
	gene_groups = list()

	gene_cleaned = list()

	# Names
	strains_names = sorted(pop_map, key=lambda x: pop[x]['genospecies'])
	print 'These are the strains evaluated', strains_names
	strains_names.remove('3260')
	strains_names.remove('3381')
	strains_names.remove('3339')
	strains_names.remove('3211')

	# Making the matrix to update the ani 
	ani_matrix = np.zeros((len(strains_names),len(strains_names)))
	ani_matrix = pd.DataFrame(ani_matrix, index = strains_names, columns = strains_names)
		

	counting = 0
	all_results = []
	print len(gene_cleaned)
	for gg1 in concatenated_genes:
			all_comparisons = []
			start_time = time.time()
			if counting <= len(concatenated_genes[0:15]):
				print 'Working on gene group: %s'%gg1
				g1 = ag[gg1[5::]]
				print 'Number of strains is %d' % len(g1['strains'][...])

				filenames = g1['strains'][...]

				for idx, fname1 in enumerate(filenames[:-1]):
					for fname2 in filenames[idx+1:]:
						result = compare_equal(g1['nsequences'][...][idx], g1['nsequences'][...][idx+1])
						ani_matrix[fname1[:4]][fname2[:4]] += result
						ani_matrix[fname2[:4]][fname1[:4]] += result
			else:
				break
			print "The time used was %6f" % (time.time()-start_time)
			counting += 1
			print counting	

	ani_matrix = ani_matrix/counting

	#print ani_matrix.diagonal(offset=0, axis1=0, axis2=1)
	print ani_matrix
	
	# Saving the matrix 
	#header = ['gene', 'min', 'max', 'mean', 'std', 'ste']
	#df = pd.DataFrame(all_results, columns=header)
	#print df

	#df.to_csv('ani_all.csv', sep='\t')
	return ani_matrix
				  
ani = average_nucleotide_identity()
ani.to_csv('ani_sorted_by_genospecies_282.csv', header = True)

# A heatmap of the matrix
pl.pcolor(ani)
pl.colorbar()
pl.xlim([0,ani.shape[1]])
pl.ylim([0,ani.shape[0]])
pl.title('Average nucleotide identity - 198 strains')
pl.show()
#pl.savefig('/faststorage/project/NChain/rhizobium/ANI/figures/heat_map_all_coregenes.png')
