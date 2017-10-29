import h5py
import numpy as np
import scipy as sp
import pandas as pd
import os
import glob 
import math
from scipy.stats.stats import pearsonr
import pylab as pl
from numpy import linalg
import collections
from collections import Counter
from collections import OrderedDict
import itertools
from matplotlib import pyplot as plt


def parse_pop_map(file_name):
    from itertools import izip
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, sara_id, origin, country in izip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Origin2']):
        pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country}
    return pop_map

def one_by_one_orthologous_groups():
	t = pd.read_table('/Users/PM/Desktop/PHD_incomplete/Methods/Topology/proteinortho_to_orthomcl_groupmap.txt', sep = ':', names=['proteinortho', 'orthomcl'])
	print t
	t = t.dropna(axis=0, how='any', thresh=None, subset=None, inplace=False)

	print 'The number of genes that are present in both Proteinortho and Orthomcl is: %d' % len(t['proteinortho'].tolist())


def ANI_figure():
	import pandas as pd
	import numpy as np
	import pylab as pl
	import time
	from matplotlib import cm as cm
	from matplotlib import pyplot as plt
	import seaborn as sns
	import matplotlib

	working_folder = '/Users/PM/Dropbox/NCHAIN_share/Pop_gene_paper/Figures/figure1_B_ANI/'


	ani = pd.DataFrame.from_csv(working_folder + 'ani_sorted_by_genospecies_194.csv')

	gs2colour = {
	'gsA':'#386CB0',
	'gsB':'#FB8072',
	'gsC':'#1B9E77',
	'gsD':'#984EA3',
	'gsE':'#F0027F'
	}
	country2colour = {
	'UK':'green',
	'F':'blue',
	'DK':'red',
	}  
	origin22colour = {
	'UK':'green',
	'F':'blue',
	'DK':'red',
	'DKO':'orange'
	} 

	maps = parse_pop_map(file_name = working_folder + 'Rhizobium_soiltypes_new.txt')

	gs_list = [maps[str(i)]['genospecies'] for i in ani.index]
	gs_colours = [gs2colour[gs] for gs in gs_list]

	country_list = [maps[str(i)]['country'] for i in ani.index]
	country_colours = [country2colour[country] for country in country_list]

	origin2_list = [maps[str(i)]['origin2'] for i in ani.index]
	origin2_colours = [origin22colour[origin2] for origin2 in origin2_list]

	cmap = cm.get_cmap('jet', 30)

	ax = plt.figure(figsize=(10, 10))

	ax = sns.clustermap( ani, 
	row_colors = origin2_colours,
	col_colors = gs_colours,
	cmap = cmap, 
	xticklabels=False,
	yticklabels=False
	)

	leg1 = ax.ax_heatmap
	labels = gs2colour.keys()
	patches1 = [
	matplotlib.patches.Patch(color=color1, label=label1)
	for label1, color1 in zip(gs2colour.keys(),gs2colour.values())]

	leg2 = ax.ax_heatmap
	labels = origin22colour.keys()
	patches2 = [
	matplotlib.patches.Patch(color=color2, label=label2)
	for label2, color2 in zip(origin22colour.keys(),origin22colour.values())]
	leg2.legend(patches2, labels, loc=(-.3,0.9), frameon=False)

	plt.legend(handles = patches1, loc = (3,0.1))
	plt.savefig(working_folder + 'ani_sorted_by_genospecies_194_PY9.pdf')

#ANI_figure()


def kinship_all_genes(snps_file= '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/new_snps.HDF5',
				 meta_data = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt',
                 min_maf=0.10,
                 max_strain_num=198,
                 save = False):
    """
    Calculates the kinship. 
    """
    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    all_strains = set()
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains) == max_strain_num:
            all_strains = set(strains).union(all_strains)
    num_strains = len(all_strains)
    print 'Found %d "distinct" strains' % num_strains

    ordered_strains = sorted(list(all_strains))

    strain_index = pd.Index(ordered_strains)
    K_snps = sp.zeros((num_strains, num_strains))
    counts_mat_snps = sp.zeros((num_strains, num_strains))

    for i, gg in enumerate(gene_groups):
        if i % 100 == 0:
            print 'Working on gene nr. %d' % i
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains) < max_strain_num:
            strain_mask = strain_index.get_indexer(strains)

            # Already normalized snps
            snps = data_g['norm_snps'][...]
            freqs = data_g['freqs'][...]
            mafs = sp.minimum(freqs, 1 - freqs)
            maf_mask = mafs > min_maf
            snps = snps[maf_mask]

            if len(snps) == 0:
                continue
            K_snps_slice = K_snps[strain_mask]
            K_snps_slice[:, strain_mask] += sp.dot(snps.T, snps)
            K_snps[strain_mask] = K_snps_slice
            counts_mat_snps_slice = counts_mat_snps[strain_mask]
            counts_mat_snps_slice[:, strain_mask] += len(snps)
            counts_mat_snps[strain_mask] = counts_mat_snps_slice

    K_snps = K_snps / counts_mat_snps  # element-wise division
    print 'The mean of the GRM diagonal is %f' % np.mean(np.diag(K_snps))

    tuple_index_kinship = (K_snps, ordered_strains)

    headers = list()
    maps = parse_pop_map(meta_data)
    for i in ordered_strains:
        headers.append( 'SM' +maps[i]['sara_id'])

    np.savetxt('strains_order.csv', headers,  delimiter=",", fmt='%s')
    np.savetxt('kinship_maf_01.csv', K_snps, fmt='%.18e', delimiter=',', header = ','.join(headers), comments="")
    
    # Saving as a compressed numpy file:
    if save == True:
    	np.savez_compressed("{}/{}".format(kinship_matrix), matrix=snps, strains=strains, maf = maf) 

    return(tuple_index_kinship)

def plot_dirty_PCA(figure_fn = 'pca.png', k_figure_fn = 'kinship_heatmap.pdf', title=None,
                   figure_dir = '/project/NChain/faststorage/rhizobium/ld/figures',strains=None,
                   meta_data = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
    from scipy import linalg
    
    kinship_mat, strains = (kinship_all_genes())

    evals, evecs = linalg.eig(kinship_mat)  #PCA via eigen decomp
    evals[evals<0]=0
    sort_indices = sp.argsort(evals,)
    ordered_evals = evals[sort_indices]
    print ordered_evals[-10:]/sp.sum(ordered_evals)
    pc1,pc2 = evecs[:,sort_indices[-1]],evecs[:,sort_indices[-2]]
    pl.clf()
    

    if strains is not None:    
        ct_marker_map = {'DK':'*', 'DKO': '^','UK':'v', 'F':'o'}
        gs_color_map = {'gsA':'#386CB0','gsB':'#FB8072', 'gsC':'#1B9E77', 'gsE': '#F0027F', 'gsD': '#984EA3'}
        pop_map = parse_pop_map(meta_data)

        for i, strain in enumerate(strains):
            d = pop_map.get(strain,'NA')
            if d=='NA':
                gs = 'NA'
                country = 'NA'
            else:
                gs = d['genospecies']
                country = d['country']
            pl.scatter(pc1[i],pc2[i], marker=ct_marker_map[country], c=gs_color_map[gs], alpha=0.2, s=100, edgecolor='none')
        for gs in gs_color_map:
            pl.scatter([], [], color=gs_color_map[gs], marker = 's', label=gs, s=100, edgecolor='none')
        for country in ct_marker_map:
            if country !='NA':
                pl.scatter([], [], color='k', marker = ct_marker_map[country], label=country, s=100, facecolors='none')

        
        pl.legend(scatterpoints=1)
    
    # Ploting Principal components    
    else:
        pl.plot(pc1,pc2,'k.')
    if title is not None:
        pl.title(title)

    # Ploting the variance explained by each PC 
    tot = sum(evals)
    var_exp = [(i / tot)*100 for i in sorted(evals, reverse=True)]
    cum_var_exp = np.cumsum(var_exp)

    pl.xlabel('PC1 (35.03%)')
    pl.ylabel('PC2 (22.34%)')
    pl.tight_layout()
    pl.savefig('Final_PCA_SNP_level.pdf',format='pdf')

    # Ploting the cumulative variance explained
    with pl.style.context('seaborn-whitegrid'):
        pl.figure(figsize=(6, 6))
        pl.bar(range(198), var_exp, alpha= 1, align='center', label='individual explained variance')
        pl.step(range(198), cum_var_exp, where='mid', label='cumulative explained variance')
        pl.ylabel('Explained variance ratio')
        pl.xlabel('Principal components')
        pl.legend(loc='best')
        #plt.tight_layout()
        pl.savefig('variance_explained_PC1234')
        #pl.show()

    pl.savefig(figure_dir+'/'+k_figure_fn)
#plot_dirty_PCA()


def simple_intergenic_ld():
	pass

def findInstances(list1, list2):
    """For each item in list1, return a list of offsets to its occurences in list2"""
    for i in list1:
        yield [pos for pos,j in enumerate(list2) if i==j]


def presence_absence_plot(figure = False):
	# Opening the list of genes and its respective strains
	#genes_strains = open('C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/PCA_gene_level_all_strains/presence_absence_headers.txt', 'r').read().split()

	# Using ordered dictionary to maintain the gene structure
	gene_groups = collections.OrderedDict()

	with open('/Users/PM/Desktop/PHD_incomplete/Methods/PCA_gene_level_all_strains/proteinortho1to1.groups.txt') as f:
		for line in f:
			temp = line.split(':')
			group_name = temp[0]
			members = temp[1].split(' ')
			members = members[1::]
			new_members = list()
			for i in members:
				new_members.append(int(i[0:4]))
			if group_name not in gene_groups:
				gene_groups[group_name] = new_members # creating a list

	#print len(gene_groups.keys())
	print sorted(gene_groups['group10'])

	# Filling up first the core genes
	keys_sorted = sorted(gene_groups, key=lambda k: len(gene_groups[k]))

	# Open the data frame with the origin information
	sorted_strain_names = pd.read_csv('/Users/PM/Desktop/PHD_incomplete/nchain/strains_id.txt', sep='\t')

	strain_names = list(sorted_strain_names['Seq ID'])
	Matrix_counts = np.zeros((200,len(gene_groups)))
	lands = list(sorted_strain_names['Genospecies'])
	countries = {'gsA': 1, 'gsB':2, 'gsC':3, 'gsD':4, 'gsE': 5}

	paralogous = 0
	histogram = list()
	count = 0
	gene_group_wo_paralogous = []
	for gene in keys_sorted:
		strains = gene_groups[gene]
		
		# Avoind paralogous:
		if len(strains) == len(set(strains)):
			gene_group_wo_paralogous.append(gene)
  			index_match = list(findInstances(strains, strain_names))
  			histogram.append(len(index_match))
  			merged = list(itertools.chain(*index_match))
  			for m in merged:      
   				Matrix_counts[m,count] = countries[lands[m]]
   			count += 1
   		else:
   			#print 'It is paralogous'
   			for m in merged:      
   				Matrix_counts[m,count] = countries[lands[m]]
   			paralogous += 1
   			count += 1

   	print 'Number of genes that are paralagous', paralogous
	print Counter(histogram)

	if figure == True:
		arr = plt.hist(histogram, 50, facecolor='darkblue', alpha=0.75)
		for i in range(50):
			plt.text(arr[1][i],arr[0][i],str(arr[0][i]), fontsize=7, rotation = 'vertical', va = 'top')

		plt.xlabel('Number of strains contain in a gene', fontsize= 8)
		plt.ylabel('Frequency of genes', fontsize= 8)
		#plt.title('Gene distribution')
		plt.savefig('Histogram_genes_shared.pdf')

	print Matrix_counts
	print 'The shape of the matrix of counts is:', Matrix_counts.shape

	# Saving the matrix
	np.savetxt("presence_absence_matrix_by_genospecies.csv", Matrix_counts, delimiter=",")  

	return(gene_groups)  


#print presence_absence_plot(figure = False)

def define_combinations(genospecies):
	n = len(genospecies)
	d = {gs:0 for gs in genospecies}

	for k in xrange(2,n+1):
		combinations = itertools.combinations(genospecies, k)
		for comb in combinations:
			d["".join(comb)] = 0
	return d

def venn_diagram():
	pop = parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt')

	gene_groups = presence_absence_plot()

	genospecies = ['gsA', 'gsB', 'gsC', 'gsD', 'gsE']

	venn_dictionary = define_combinations(genospecies)
	total_geno = {'gsA': 0, 'gsB':0, 'gsC': 0,'gsD': 0, 'gsE': 0}
	hist_per_geno = {'gsA': {}, 'gsB':{}, 'gsC': {},'gsD': {}, 'gsE': {}}

	for gene in gene_groups:
  		strains = gene_groups[gene]
  		print strains
  		gs_list = sp.array([pop[str(strain)]['genospecies'] for strain in strains])
  		print gs_list
  		counts_dict = Counter(gs_list)

  		# How many genes in total are present in each genospecies?
  		gs_list = sp.unique(gs_list)
  		for i in gs_list:
  			total_geno[i] += 1

  		# How many genes contain a given amount of strains from a certain genospecies?:
  		for geno, counts in counts_dict.items():
  			tup = (geno, counts)
  			if tup[1] not in hist_per_geno[tup[0]]:
  				hist_per_geno[tup[0]][tup[1]] = 1
    		else:
      			hist_per_geno[tup[0]][tup[1]] += 1

  		# How many genes are shared by each possible combination of genospecies?:    
  		gs_list_joined = "".join(sorted(sp.unique(gs_list)))
  		venn_dictionary[gs_list_joined] += 1

	print venn_dictionary

	# The barplot figure is written in r, script: barplot_genes_genospecies.R 
#print venn_diagram()

def len_forced(obj):
	ob_len = len(obj)
	if ob_len > 5:
		ob_len = 5
	return(ob_len)

def colorful_heatmap(Matrix_counts):
	import seaborn as sns
	import sys
	sys.setrecursionlimit(100000)

	pop = parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt')

	# Open the data frame with the origin information
	sorted_strain_names = pd.read_csv('/Users/PM/Desktop/PHD_incomplete/nchain/strains_id.txt', sep='\t')
	strain_names = list(sorted_strain_names['Seq ID'])
	#strain_names.remove(3260)

	counts = sp.array([pop[str(strain)]['genospecies'] for strain in strain_names])
	print sum(counts == 'gsA')
	print sum(counts == 'gsB')
	print sum(counts == 'gsC')
	print sum(counts == 'gsD')
	print sum(counts == 'gsE')
	
	# Removing one empty row
	Matrix_counts = Matrix_counts[~(Matrix_counts==0).all(1)]
	b = np.sum(Matrix_counts, axis = 0)
	idx = b.argsort()

	idx = sorted(range(0,Matrix_counts.shape[1]), key=lambda x: list(np.unique(Matrix_counts[:,x])))
	idx = sorted(idx, key=lambda x: len_forced(list(np.unique(Matrix_counts[:,x]))))
	Matrix_counts = Matrix_counts[:,idx]

	print 'The shape of the matrix of counts is:', Matrix_counts.shape

	# Trying clustermap 
	gsA = ['#386CB0']
	gsB = ['#FB8072']
	gsC = ['#1B9E77']
	gsD = ['#984EA3']
	gsE = ['#F0027F']

	colors = 33*gsA + 34*gsB + 116*gsC + 5*gsD + 11*gsE

	#method='average', metric= 'hamming'
	pandas_data_frame = pd.DataFrame(Matrix_counts[0:199, 0:Matrix_counts.shape[1]])

	g = sns.clustermap(pandas_data_frame, row_cluster=True, row_colors = colors, xticklabels=True)
	plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
	plt.savefig('presence_absence_genospecies_total_rowcluster.pdf')

Matrix_counts = np.genfromtxt("presence_absence_matrix_by_genospecies.csv", delimiter=",")
colorful_heatmap(Matrix_counts = Matrix_counts) 

def pca_gene_level(Matrix_counts, meta_data = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):

	# Removing one empty row
	Matrix_counts = Matrix_counts[~(Matrix_counts==0).all(1)]

	strains = pd.read_csv('/Users/PM/Desktop/PHD_incomplete/nchain/strains_id.txt', sep='\t')
	strain_names = list(strains['Seq ID'])
	strain_names.remove(3260)

	cov = np.dot(Matrix_counts, Matrix_counts.T)
	evals, evecs = linalg.eig(cov)  #PCA via eigen decomp
	evals[evals<0]=0
	sort_indices = sp.argsort(evals,)
	ordered_evals = evals[sort_indices]
	print ordered_evals[-10:]/sp.sum(ordered_evals)
	pc1,pc2 = evecs[:,sort_indices[-1]],evecs[:,sort_indices[-2]]
	pl.clf()


	if strain_names is not None:    
		ct_marker_map = {'DK':'*', 'DKO': '^','UK':'v', 'F':'o'}
		gs_color_map = {'gsA':'#386CB0','gsB':'#FB8072', 'gsC':'#1B9E77', 'gsE': '#F0027F', 'gsD': '#984EA3'}
		pop_map = parse_pop_map(meta_data)

		for i, strain in enumerate(strain_names):
			d = pop_map.get(str(strain),'NA')
			print d
			if d=='NA':
				gs = 'NA'
				country = 'NA'
			else:
				gs = d['genospecies']
				country = d['country']
			pl.scatter(pc1[i],pc2[i], marker=ct_marker_map[country], c=gs_color_map[gs], alpha=0.2, s=100, edgecolor='none')
		for gs in gs_color_map:
			pl.scatter([], [], color=gs_color_map[gs], marker = 's', label=gs, s=100, edgecolor='none')
		for country in ct_marker_map:
			if country !='NA':
					pl.scatter([], [], color='k', marker = ct_marker_map[country], label=country, s=100, facecolors='none')


		pl.legend(scatterpoints=1)

	# Ploting the variance explained by each PC 
	tot = sum(evals)
	var_exp = [(i / tot)*100 for i in sorted(evals, reverse=True)]

	pl.xlabel('PC1 (86.25%)')
	pl.ylabel('PC2 (3.34%)')
	pl.tight_layout()
	pl.savefig('Final_PCA_gene_level.pdf',format='pdf')

#Matrix_counts = np.genfromtxt("presence_absence_matrix_zeros_ones.csv", delimiter=",")
#pca_gene_level(Matrix_counts)
    



