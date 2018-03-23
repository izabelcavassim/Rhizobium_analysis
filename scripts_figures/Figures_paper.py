import h5py
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
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
import matplotlib
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
from matplotlib import pyplot as plt
from sys import argv
figure_dir = '/Users/PM/Desktop/New_data/final_figures/'
tables_directory = '/Users/PM/Desktop/New_data/final_results/'

def parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
	from itertools import izip
	
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, sara_id, origin, country in izip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Origin2']):
		pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country, 'origin2':origin}
	return pop_map


def ANI_figure(figure_dir = figure_dir, level = 'snp_level'):
	import time
	from matplotlib import cm as cm


	#working_folder = '/Users/PM/Dropbox/NCHAIN_share/Pop_gene_paper/Figures/figure1_B_ANI/'
	working_folder = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/'

	if level == 'snp_level':
		ani = pd.DataFrame.from_csv('/Users/PM/Desktop/scripts_Asger/Rhizobium_analysis/scripts_figures/ani_sorted_by_genospecies_snps_new_data.csv', sep = ';')
		ani.values[[np.arange(ani.shape[0])]*2] = 1
		strain_names = ani.index
		#ani = pd.DataFrame.from_csv('/Users/PM/Desktop/scripts_Asger/Rhizobium_analysis/scripts_figures/ani_sorted_by_genospecies_snps_new_data_maf_0.1.csv')

	if level == 'gene_level':
		presence_absence = pd.read_csv('/Users/PM/Desktop/scripts_Asger/Rhizobium_analysis/scripts_figures/presence_absence_matrix_by_genospecies_new_dataset_ones_final.csv', sep = ',', header = 0, index_col = 0)
		strain_names = presence_absence.index

		ani = sp.dot(presence_absence, presence_absence.T)

	if level == 'ANI_level':
		ani = pd.DataFrame.from_csv('/Users/PM/Desktop/scripts_Asger/Rhizobium_analysis/scripts_figures/ani_sorted_by_genospecies_282.csv', sep = ',')
		
		strain_names = ani.index	
		print ani	
		ani=ani.values
		print ani
		print np.fill_diagonal(ani, 1)

		print ani

	#ani = pd.DataFrame.from_csv('/Users/PM/Dropbox/BACKUPWINDOWS/MASTER/PHD/Bjarnicode/scripts/ani_sorted_by_genospecies_500_genes.csv', sep =';')

	gs2colour = {
	'gsA':'#386CB0',
	'gsB':'#FB8072',
	'gsC':'#1B9E77',
	'gsD':'#984EA3',
	'gsE':'#F0027F'
	}
	country2colour = {
	'UK': '#8dd3c7',
	'F': '#ffffb3',
	'DK': '#bebada',	
	'Ref':'#80b1d3'
	}  
	origin22colour = {
	'UK': '#8dd3c7',
	'F': '#ffffb3',
	'DK': '#bebada',	
	'DKO': '#fb8072',
	'Ref':'#80b1d3'
	} 

	maps = parse_pop_map(file_name = working_folder + 'Rhizobium_soiltypes_new.txt')

	gs_list = [maps[str(i)]['genospecies'] for i in strain_names]
	gs_colours = [gs2colour[gs] for gs in gs_list]
	origin2_list = [maps[str(i)]['country'] for i in strain_names]


	origin2_colours = [origin22colour[origin2] for origin2 in origin2_list]

	cmap = cm.get_cmap('jet', 30)

	ax = plt.figure(figsize=(10, 10))
	ax = sns.clustermap(ani, 
	row_colors = origin2_colours,
	col_colors = gs_colours,
	cmap = cmap, 
	#method="median",
	#metric = 'cosine',
	xticklabels=False,
	yticklabels=False,
	row_cluster=True, col_cluster=True,
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
	plt.savefig(figure_dir + 'ani_sorted_by_genospecies_{}.pdf'.format(level))

ANI_figure(level = 'ANI_level')


def kinship_all_genes(snps_file= '/Users/PM/Desktop/New_data/newsnps_100.hdf5',
				 meta_data = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt',
				 min_maf=0.05,
				 max_strain_num=201,
				 save = False):
	"""
	Calculates the kinship. 
	"""
	h5f = h5py.File(snps_file)
	gene_groups = h5f.keys()
	all_strains = set()

	data_g = h5f['1258']
	all_strains = data_g['strains'][...]

	print all_strains
	num_strains = len(all_strains)
	print 'Found %d "distinct" strains'%num_strains

	ordered_strains = sorted(list(all_strains))

	strain_index = pd.Index(ordered_strains)
	K_snps = sp.zeros((num_strains, num_strains))
	counts_mat_snps = sp.zeros((num_strains, num_strains))

	for i, gg in enumerate(gene_groups):
		if i % 100 == 0:
			print 'Working on gene nr. %d' % i
		data_g = h5f[gg]
		strains = data_g['strains'][...]

		print len(strains)
		if len(strains) <= max_strain_num:
			strain_mask = strain_index.get_indexer(strains)

			#print len(strain_mask)

			# Already normalized snps
			snps = data_g['snps'][...]
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

	print np.diag(K_snps)
	print 'The mean of the GRM diagonal is %f' % np.mean(np.diag(K_snps))

	K_snps = K_snps / np.diag(K_snps)

	tuple_index_kinship = (K_snps, ordered_strains)

	K_snps = pd.DataFrame(data=K_snps,    # values
			index=ordered_strains,    # 1st column as index
			columns=ordered_strains)  # 1st row as the column names


	K_snps.to_csv('ani_sorted_by_genospecies_snps_new_data_maf_0.05.csv', header = True)

	#headers = list()
	#maps = parse_pop_map(meta_data)
	#for i in ordered_strains:
	#	headers.append( 'SM' +maps[i]['sara_id'])

	#np.savetxt('strains_order.csv', headers,  delimiter=",", fmt='%s')
	#np.savetxt('kinship_maf_01.csv', K_snps, fmt='%.18e', delimiter=',', header = ','.join(headers), comments="")
	
	# Saving as a compressed numpy file:
	#if save == True:
	#	np.savez_compressed("{}/{}".format(kinship_matrix), matrix=snps, strains=strains, maf = maf) 

	return(tuple_index_kinship)

#print kinship_all_genes()

def plot_dirty_PCA(figure_fn = 'pca.png', k_figure_fn = 'kinship_heatmap.pdf', title=None,
				   figure_dir = '/Users/PM/Desktop/New_data/final_figures',strains=None,
				   meta_data = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt', pc_1 = 5, pc_2 = 6):
	from scipy import linalg

	kinship_mat, strains = (kinship_all_genes())

	evals, evecs = linalg.eig(kinship_mat)  #PCA via eigen decomp

	print evecs
	evals[evals<0]=0
	sort_indices = sp.argsort(evals,)
	ordered_evals = evals[sort_indices]
	print ordered_evals[-10:]/sp.sum(ordered_evals)
	pc1,pc2 = evecs[:,sort_indices[-pc_1]],evecs[:,sort_indices[-pc_2]]
	pl.clf()
	
	with pl.style.context('default'):    

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

	print var_exp

	#pl.xlim(-0.165,0.10)
	#pl.ylim(-0.10,0.15)
	pl.xlabel('PC5')
	pl.ylabel('PC6')
	pl.tight_layout()
	pl.savefig(figure_dir+'/''Final_PCA_SNP_level_PC5_PC6.pdf',format='pdf')

	# Ploting the cumulative variance explained
	with pl.style.context('default'):
		pl.figure(figsize=(6, 6))
		pl.bar(range(len(var_exp)), var_exp, alpha= 1, align='center', label='individual explained variance')
		pl.step(range(len(var_exp)), cum_var_exp, where='mid', label='cumulative explained variance')
		pl.ylabel('Explained variance ratio')
		pl.xlabel('Principal components')
		pl.grid(True)
		pl.legend(loc='best')
		#plt.tight_layout()
		pl.savefig(figure_dir+'/variance_explained.pdf')
		#pl.show()

#plot_dirty_PCA()


def simple_intergenic_ld():
	pass

def findInstances(list1, list2):
	"""For each item in list1, return a list of offsets to its occurences in list2"""
	for i in list1:
		yield [pos for pos,j in enumerate(list2) if i==j]

def shared_orthologous_genes(figure_dir = '/Users/PM/Desktop/New_data/final_figures/'):
	import seaborn as sns
	

	presence_absence = pd.read_csv('/Users/PM/Desktop/scripts_Asger/Rhizobium_analysis/scripts_figures/presence_absence_matrix_by_genospecies_new_dataset_ones_final.csv', sep = ',', header = 0, index_col = 0)
	strain_names = presence_absence.index

	dot_product = sp.dot(presence_absence, presence_absence.T)
	print dot_product

	
	ax = sns.clustermap(dot_product)
	plt.savefig(figure_dir+'shared_genes.pdf')

	np.fill_diagonal(dot_product, 0)
	mean_value = np.mean(dot_product, axis = 0)
	print mean_value


#print shared_orthologous_genes()

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


def presence_absence_plot(figure = True, proteinortho_directory = '/Users/PM/Desktop/New_data/final_project.poff_disambiguated.groups', figure_dir = '/Users/PM/Desktop/New_data/final_figures/', 
	tables_directory = '/Users/PM/Desktop/New_data/final_results/'):
	# Opening the list of genes and its respective strains
	#genes_strains = open('C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/PCA_gene_level_all_strains/presence_absence_headers.txt', 'r').read().split()
	#import seaborn as sns
	import sys
	from matplotlib import cm as cm
	import matplotlib
	plt.style.use('seaborn-dark-palette')
	sys.setrecursionlimit(100000)
	# Using ordered dictionary to maintain the gene structure
	gene_groups = collections.OrderedDict()
	pop = parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt')


	#'/Users/PM/Desktop/PHD_incomplete/Methods/PCA_gene_level_all_strains/proteinortho1to1.groups.txt'
	with open(proteinortho_directory) as f:
		for line in f:
			temp = line.split(':')
			group_name = temp[0]
			group_name = group_name.split('|')[0]
			members = temp[1].split(' ')
			members = members[1::]
			new_members = list()
			for i in members:
				try:
					new_members.append(i[0:4]) # changing here to accept the new strains
				except:
					continue
			if group_name not in gene_groups:
				gene_groups[group_name] = new_members # creating a list

	#print len(gene_groups.keys())
	a = invert_dict(gene_groups)
	
	for i, j in a.items():
		print('{} {}'.format(i,len(j)))

	strain_names = sorted(gene_groups['group10'])
	print strain_names

	# Filling up first the core genes
	keys_sorted = sorted(gene_groups, key=lambda k: len(gene_groups[k]))

	# Open the data frame with the origin information
	sorted_strain_names = pd.read_csv('/Users/PM/Desktop/PHD_incomplete/nchain/strains_id.txt', sep='\t')

	Matrix_counts = np.zeros((len(strain_names),len(gene_groups)))
	countries = {'gsA': 1, 'gsB':2, 'gsC':3, 'gsD':4, 'gsE': 5}

	paralogous = 0
	histogram = list()
	count = 0
	headers = []
	for gene in keys_sorted:
		headers.append(gene)
		strains = sorted(gene_groups[gene])

		index_match = list(findInstances(strains, strain_names))

		histogram.append(len(index_match))
		merged = list(itertools.chain(*index_match))

		for m in merged:      
			Matrix_counts[m,count] = countries[pop[strain_names[m]]['genospecies']]
		count += 1

	stats = Counter(histogram)

	for values, keys in stats.items():
		print keys, values

	if figure == True:
		with plt.style.context('default'):    
			arr = plt.hist(histogram, 50, facecolor='black', alpha=0.75)
			#for i in range(50):
			#	plt.text(arr[1][i],arr[0][i],str(arr[0][i]), fontsize=7, rotation = 'vertical', va = 'top')


			plt.xlabel('Number of genomes', fontsize= 8)
			plt.ylabel('Number of orthologous genes', fontsize= 8)
		
			plt.savefig(figure_dir + 'Histogram_genes_shared_genomes_final.pdf')

	#print Matrix_counts
	#print 'The shape of the matrix of counts is:', Matrix_counts.shape
	
		#print headers
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
		'DKO': 'yellow',
		'Ref':'grey'
		}

		gs_list = [pop[str(strain)]['genospecies'] for strain in strain_names]
		gs_colours = [gs2colour[gs] for gs in gs_list]

		print Counter(gs_colours)

		country_list = [pop[str(strain)]['country'] for strain in strain_names]
		country_colours = [country2colour[country] for country in country_list]


		pandas_data_frame = pd.DataFrame(Matrix_counts)


		# Initializing figure

		ax = plt.figure(figsize=(10, 10))
		ax = sns.clustermap(pandas_data_frame, 
		row_colors = gs_colours,
		xticklabels=False,
		yticklabels=False
		)

		columns = ax.linkage.dendrogram_col()
		sns.heatmap(arr2, ax=axes[1,0], yticklabels = ['SNPs', 'Members'], xticklabels = ' ')

		leg1 = ax.ax_heatmap
		labels = gs2colour.keys()

		## Genospecies colors
		patches1 = [
		matplotlib.patches.Patch(color=color1, label=label1)
		for label1, color1 in zip(gs2colour.keys(),gs2colour.values())]
	
	## Country colors 
	#patches2 = [
	#matplotlib.patches.Patch(color=color2, label=label2)
	#for label2, color2 in zip(gs2colour.keys(),gs2colour.values())]

		plt.legend(handles = patches1, loc = (3,0.1))
		plt.savefig('Final_figures/Presence_absence_plot_final.pdf')



	df = pd.DataFrame(Matrix_counts, index=strain_names, columns= headers)
	df.to_csv(tables_directory + "presence_absence_matrix_by_genospecies_new_dataset_final.csv", sep=',', index= strain_names)

	return(Matrix_counts, headers, strain_names, gene_groups)
	# Saving the matrix
	#np.savetxt("presence_absence_matrix_by_genospecies_new_dataset_all.csv", Matrix_counts, delimiter=",", fmt='%.4e')  

#print presence_absence_plot(figure = True)

def define_combinations(genospecies = ['gsA', 'gsB', 'gsC', 'gsD', 'gsE']):
	n = len(genospecies)
	d = {gs:0 for gs in genospecies}

	for k in xrange(2,n+1):
		combinations = itertools.combinations(genospecies, k)
		for comb in combinations:
			d["".join(comb)] = 0
	return d
import csv

def make_histograms(dicts, xlab = None, ylab = None, save_name = None, fig_dir = figure_dir):
	import pylab as pl
	X = np.arange(len(dicts))


	with open('%s/figure_%s.csv'%(fig_dir,save_name), 'wb') as f:  # Just use 'w' mode in 3.x
		w = csv.DictWriter(f, dicts.keys())
		w.writeheader()
		w.writerow(dicts)
   	pl.bar(X, dicts.keys(), align='center', width=0.5)
	pl.xticks(X, dicts.keys())
  	ymax = max(dicts.values())*1.05
  	xmax = len(dicts)
  	pl.xlabel(xlab)
  	pl.ylabel(ylab)
  	pl.ylim(0, ymax)
  	pl.xlim(-1,xmax)
  	pl.grid(True)
  	pl.savefig('%s/figure_%s.pdf'%(fig_dir,save_name))
   	pl.show()

def venn_diagram():
	pop = parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt')

	gene_groups = presence_absence_plot(figure = False)[3]

	#genospecies = ['gsA', 'gsB', 'gsC', 'gsD', 'gsE']

	venn_dictionary = define_combinations()
	total_geno = {'gsA': 0, 'gsB':0, 'gsC': 0,'gsD': 0, 'gsE': 0}
	hist_per_geno = {'gsA': {}, 'gsB':{}, 'gsC': {},'gsD': {}, 'gsE': {}}

	for gene in gene_groups:
		strains = gene_groups[gene]

		gs_list = sp.array([pop[str(strain)]['genospecies'] for strain in strains])

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

	#print venn_dictionary
	print hist_per_geno
	print total_geno

	print sum(venn_dictionary.values())

	strains_geno = {'gsC': 116, 'gsB': 32, 'gsA': 32, 'gsE': 11, 'gsD': 5}


	print '\nAll the possible combinations is:', venn_dictionary
	print '\nTotal amount of genes present in each genospecies', total_geno
	print '\nDistribution of genes per genospecies', hist_per_geno

	'''Total number of genes present in each genospecies'''

	make_histograms(total_geno, xlab = 'Genospecies', ylab = 'Genes', save_name = 'genes_genospecies.pdf')

	make_histograms(strains_geno, xlab = 'Genospecies', ylab = 'Strains', save_name = 'strains_genospecies.pdf')

	for i in strains_geno.keys():
		make_histograms(hist_per_geno[i], xlab= 'Strains', ylab = 'Genes', save_name = i)

	# The barplot figure is written in r, script: barplot_genes_genospecies.R 
#print venn_diagram()

def len_forced(obj):
	ob_len = len(obj)
	if ob_len > 5:
		ob_len = 5
	return(ob_len)

def colorful_heatmap():
	import seaborn as sns
	import sys
	from matplotlib import cm as cm
	import matplotlib
	sys.setrecursionlimit(100000)

	pop = parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt')

	Matrix_counts = pd.DataFrame.from_csv('/Users/PM/Desktop/scripts_Asger/Rhizobium_analysis/scripts_figures/presence_absence_matrix_by_genospecies_new_dataset_final.csv')
	
	Info_genes = pd.read_table('/Users/PM/Desktop/New_data/Final_Assembly_scores_statistics/final_project.poff_disambiguated.groups.scores', sep = '\t', 	index_col=False )

	print Info_genes
	print Info_genes.columns.values
	# Tranforming the pandas dataframe in a numpy array
	#results =  presence_absence_plot()
	#Matrix_counts = results[0]
	#headers = results[1]
	strain_names = Matrix_counts.index	
	gene_names = Matrix_counts.columns.values 

	print 'The shape of the matrix of counts is:', Matrix_counts.shape


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
	'DKO': 'yellow',
	'Ref':'grey'
	}

	gs_list = [pop[str(strain)]['genospecies'] for strain in strain_names]
	gs_colours = [gs2colour[gs] for gs in gs_list]

	country_list = [pop[str(strain)]['country'] for strain in strain_names]
	country_colours = [country2colour[country] for country in country_list]


	pandas_data_frame = pd.DataFrame(Matrix_counts)

	#ax = plt.figure(figsize=(10, 10))
	g = sns.clustermap(pandas_data_frame.iloc[:,0:10], 
	row_colors = gs_colours,
	xticklabels=False,
	yticklabels=False
	)

	columns = g.dendrogram_col.reordered_ind

	sorted_genes = list(gene_names[columns])

	print sorted_genes
	sub = Info_genes[Info_genes['Group'].isin(sorted_genes)]
	sub = sub.fillna(0)
	sub = sub[['Connectivity', 'Synteny']]

	#sns.heatmap(sub, yticklabels = ['Connectivity', 'Synteny'], xticklabels = ' ')

	leg1 = g.ax_heatmap
	labels = gs2colour.keys()

	## Genospecies colors
	patches1 = [matplotlib.patches.Patch(color=color1, label=label1)
	for label1, color1 in zip(gs2colour.keys(),gs2colour.values())]

	plt.legend(handles = patches1, loc = (3,0.1))


	# Compute and plot first dendrogram.
	fig = pl.figure(figsize=(8,8))
	ax1 = fig.add_axes([0.09,0.1,0.2,0.6])

	# Tell pointplot to plot on ax1 with the ax argument
	#sns.clustermap(x="x", y="y", data=data, ax=ax1)
	ax1.set_xticks([])
	ax1.set_yticks([])

	im = sns.clustermap(pandas_data_frame.iloc[:,0:10], 
	row_colors = gs_colours,
	xticklabels=False,
	yticklabels=False)
	im.set_xticks([])
	im.set_yticks([])

	# Tell the factorplot to plot on ax2 with the ax argument
	# Also store the FacetGrid in 'g'
	#g=sns.factorplot(x="x", y="y", data=data, ax=ax2)

	# Compute and plot second dendrogram.
	ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
	sns.heatmap(sub, yticklabels = ['Connectivity', 'Synteny'], xticklabels = ' ', ax = ax2)
	ax2.set_xticks([])
	ax2.set_yticks([])
	# Close the FacetGrid figure which we don't need (g.fig)
	#plt.close(g.fig)

	plt.show()


	plt.savefig('Final_figures/Presence_absence_plot_final_test.png')



	#print g.dendrogram_col.reordered_ind
	#print g.dendrogram_row.reordered_ind

	#g = sns.clustermap(pandas_data_frame, row_cluster=True, row_colors = gs_colours[g.dendrogram_row.reordered_ind], xticklabels=False)
	#plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

#colorful_heatmap()


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
		ct_marker_map = {'DK':'*', 'DKO': '^','UK':'v', 'F':'o', 'Ref': 'o'}
		gs_color_map = {'gsA':'#386CB0','gsB':'#FB8072', 'gsC':'#1B9E77', 'gsE': '#F0027F', 'gsD': '#984EA3', 'Ref': '#1B9E77'}
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
			print country
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


def gc_content(sequence):
	#1: GC content per genes, split by core and accessory genes
	#2: GC content per genospecies, split by core and accessory genes

	sequence = np.reshape(sequence, sequence.shape[0]*sequence.shape[1])

	sequence = sequence.tolist()

	GC = 0
	for i in sequence:
		if i == 2 or i == 3:
			GC += 1
	return GC/float(len(sequence))


def GC_content_per_genospecies(	fig_dir = figure_dir,
								 geno = 'gsA',
								 #geno_species=[argv[1], argv[1]], bin_size=0.2,
								 gt_hdf5_file='/Users/PM/Desktop/New_data/final_snps.hdf5'):
	pop = parse_pop_map()
	pop_map = pop.keys()
	ct_array = pop.values()

	h5f = h5py.File(gt_hdf5_file)
	ag = h5f['alignments']
	gene_big_groups = sorted(ag.keys())
	gene_groups_accessory = list()

	geno_species = [geno, geno]

	# Names
	strains_names = sorted(pop_map, key=lambda x: pop[x]['genospecies'])
	print 'These are the strains evaluated', strains_names

	strains_names.remove('3260')
	strains_names.remove('3381')
	strains_names.remove('3339')
	strains_names.remove('3211')
	strains_list = strains_names

	gs_content_results = []
	number_of_strains = []

	# Taking just the accesory genes genes
	for gene in gene_big_groups:
		if len(ag[gene]['strains']) < 196:
			gene_groups_accessory.append(gene)

	# Taking just the core genes
	for gene in gene_groups_accessory:

		# Gene strain list
		strains_list = ag[gene]['strains'][...]

		# Looking at specific genospecies
		gs_list = []
		for strain in strains_list:
			gs_list.append(pop[strain[:4]]['genospecies'])
	
		# Transforming the strain list in array
		strains_list = np.asarray(strains_list)
		gs_filter1, gs_filter2 = [sp.in1d(gs_list, gs) for gs in geno_species]

		# Extracting the nucleotide sequences
		g1 = ag[gene]
		g1 = g1['nsequences'][...]

		# Filtering by genospecies
		g1_gno = g1[gs_filter1,:]

		if g1_gno.shape[0] != 0:
			gs_content_results.append(gc_content(g1_gno))
			number_of_strains.append(g1_gno.shape[0])
		#print ' These are the strains belonging to this genospecies', g1_gno
		name = 'GC_' + geno + '.csv'
		np.savetxt(tables_directory+ name, np.c_[gs_content_results, number_of_strains], delimiter=",")

#GC_content_per_genospecies(geno = 'gsA')	
#GC_content_per_genospecies(geno = 'gsB')
#GC_content_per_genospecies(geno = 'gsC')
#GC_content_per_genospecies(geno = 'gsD')
#GC_content_per_genospecies(geno = 'gsE')

def GC3s_content(sequence):

	codontable = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
	'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
	'---':'X'
	}

	nt_decode_map = {1:'A', 2:'C', 3:'G', 4:'T', 5:'-', 6:'N'}
	#1: GC content per genes, split by core and accessory genes
	#2: GC content per genospecies, split by core and accessory genes

	sequence = np.reshape(sequence, sequence.shape[0]*sequence.shape[1])

	sequence = sequence.tolist()

	GC3 = 0
	codons = [sequence[x:x+3] for x in xrange(0, len(sequence), 3)]
	for codon in codons:
		if codon != [4, 1, 3] and codon != [4, 3, 3]:
			if codon[2] == 2 or codon[2] == 3: 
				GC3 += 1
	print GC3
	return GC3/float(len(codons))


def GC_content_per_gene(fig_dir = figure_dir,
						geno_species='gsA', bin_size=0.2,
						meta_data = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt',
						gt_hdf5_file='/Users/PM/Desktop/New_data/final_snps.hdf5'):

	pop = parse_pop_map(meta_data)
	pop_map = pop.keys()
	
	ct_array = pop.values()

	h5f = h5py.File(gt_hdf5_file)
	ag = h5f['alignments']
	gene_big_groups = sorted(ag.keys())

	#print gene_big_groups
	core_gene_groups = list()
	accessory_gene_groups = list()

	# Names
	strains_names = sorted(pop_map, key=lambda x: pop[x]['genospecies'])
	print 'These are the strains evaluated', strains_names
	strains_names.remove('3260')
	strains_names.remove('3381')
	strains_names.remove('3339')
	strains_names.remove('3211')
	
	# Taking just the core genes and paralogous
	for gene in gene_big_groups:
			
		strains_gene = ag[gene]['strains'][...]

		if len(strains_gene) == 196:
			core_gene_groups.append(gene)
		else:
			accessory_gene_groups.append(gene)

	print 'Starting core gene groups', len(core_gene_groups)
	results_core = []
	gene_name_core = []
	gene_members = []
	results_core_gc3 = []

	for gg1 in core_gene_groups:
		print 'Working on gene group: %s'%gg1
		g1 = ag[gg1]
		strains_gene = g1['strains'][...]

		gene_members.append(len(strains_gene))
		g1 = g1['nsequences'][...]
		print g1.shape

		gene_name_core.append(gg1)
		results_core.append(gc_content(g1))
		results_core_gc3.append(GC3s_content(g1))

	data1 = np.array([gene_name_core, results_core, results_core_gc3, gene_members])
	data1 = data1.T
	np.savetxt("GC_core_genes_test.csv", data1, fmt= "%s", delimiter=",")

	print 'Starting accessory gene groups', len(accessory_gene_groups)
	results_accessory = []
	gene_name_accessory = []
	gene_members = []
	results_accessory_gc3 = []
	for gg1 in accessory_gene_groups:
		print 'Working on gene group: %s'%gg1

		if gg1 != 10705:
			g1 = ag[gg1]

			strains_gene = g1['strains'][...]


			gene_members.append(len(strains_gene))
			g1 = g1['nsequences'][...]

			gene_name_accessory.append(gg1)
			results_accessory.append(gc_content(g1))
			results_accessory_gc3.append(GC3s_content(g1))

	data2 = np.array([gene_name_accessory, results_accessory, results_accessory_gc3, gene_members])
	data2 = data2.T

	np.savetxt("GC_accessory_genes_test.csv", data2, fmt= "%s", delimiter = ",")

	
#print GC_content_per_genospecies()
#print GC_content_per_gene()



