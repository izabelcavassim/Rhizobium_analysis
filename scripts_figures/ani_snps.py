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
								gt_hdf5_file= '/Users/PM/Desktop/newsnps_100.hdf5'):

	pop = parse_pop_map()
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

	# Making the matrix to update the ani 
	ani_matrix = np.zeros((len(strains_names),len(strains_names)))
	ani_matrix = pd.DataFrame(ani_matrix, index = strains_names, columns = strains_names)
		
	# Taking just the core genes
	deleted = ['3859', '3260']
	for gene in gene_big_groups:
			strains_total = ag[gene]['strains'][...]
			if len(strains_total) == len(strains_names):
					gene_cleaned.append(gene)


	snps_number = 0
	counting = 0
	ani = None

	ani = sp.zeros((len(strains_names), len(strains_names)))
	all_results = []
	print len(gene_cleaned)
	for gg1 in gene_cleaned[0:100]:
			all_comparisons = []
			start_time = time.time()
			if counting <= len(gene_cleaned):
				print 'Working on gene group: %s'%gg1
				g1 = ag[gg1]
				print g1.keys()
				snps = g1['snps'][...]
				print snps.shape[0]
				print np.sum(snps, axis = 0)

				ani += sp.dot(snps.T, snps)
				print ani
				ani = np.divide(ani, np.diag(ani))
	print ani