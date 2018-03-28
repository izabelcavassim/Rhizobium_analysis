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

def average_nucleotide_identity(fig_dir = '/Users/PM/Desktop/',
								gt_hdf5_file= '/Users/PM/Desktop/New_data/newsnps_100.hdf5'):

	pop = parse_pop_map()
	print pop
	pop_map = pop.keys()
	ct_array = pop.values()

	h5f = h5py.File(gt_hdf5_file)
	ag = h5f
	gene_big_groups = sorted(ag.keys())

	print gene_big_groups
	g1 = h5f[gene_big_groups[0]]

	print g1.keys()

	gene_cleaned = list()

	# Names
	strains_names = sorted(pop_map, key=lambda x: pop[x]['genospecies'])
	print 'These are the strains evaluated', strains_names
	strains_names.remove('3260')
	strains_names.remove('3381')
	strains_names.remove('3339')
	strains_names.remove('3211')

	strain_index = pd.Index(strains_names)

	# Making the matrix to update the ani 
	ani_matrix = np.zeros((len(strains_names),len(strains_names)))
	ani_matrix = pd.DataFrame(ani_matrix, index = strains_names, columns = strains_names)
		
	# Taking just the core genes
	for gene in gene_big_groups:
			strains_total = ag[gene]['strains'][...]
			if len(strains_total) == 196:
					gene_cleaned.append(gene)
	snps_number = 0
	counting = 0
	ani = None

	ani = sp.zeros((len(strains_names), len(strains_names)))
	all_results = []
	print 'Number of core genes: %f' % len(gene_cleaned)
	for gg1 in gene_cleaned:
			all_comparisons = []
			start_time = time.time()
			if counting <= len(gene_cleaned):
				print 'Working on gene group: %s'%gg1
				g1 = ag[gg1]
				
				strains = g1['strains'][...]
				snps = g1['snps'][...]

				# all major allele:
				ani += sp.dot(np.float64(snps.T), np.float64(snps))

				# all minor allele:
				snps[snps == 1.0] = np.ma.masked

				snps[snps == 0] = 1.0
				
				ani += sp.dot(np.float64(snps.T), np.float64(snps))

				ani += np.divide(ani, np.diag(ani))

	ani = np.divide(ani, np.diag(ani))

	ani_dataframe = pd.DataFrame(data=ani,    # values
            index=strains,    # 1st column as index
            columns=strains)  # 1st row as the column names

	return(ani_dataframe)

ani = average_nucleotide_identity()
print  np.diag(ani)


ani.to_csv('ani_sorted_by_genospecies_snps_new_data.csv', header = True)
pl.pcolor(ani)
pl.colorbar()
pl.xlim([0,ani.shape[1]])
pl.ylim([0,ani.shape[0]])
pl.title('test')
#pl.show()
pl.savefig('/Users/PM/Desktop/scripts_Asger/Rhizobium_analysis/scripts_figures/heat_map_all_coregenes.png')