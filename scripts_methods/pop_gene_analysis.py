# Site Frequency Spectrum
# 1. Define the group of genes to be analysed
# 2. Segment it in genospecies
# 3. Calculate the SFS for each genospecies
import h5py 
import pandas as pd
import numpy as np
import scipy as sp
from sys import argv
import matplotlib
from collections import Counter
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
matplotlib.rcParams.update({'font.size': 10})
import matplotlib.pyplot as plt
from decimal import *
import pylab as pl

def parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
	from itertools import izip   
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, origin, country, origin2 in izip(t['Seq ID'], t['Genospecies'], t['Country'], t['Origin2']):
		pop_map[str(strain_id)]={'genospecies':origin, 'country':country, 'origin':origin2}
	return pop_map

def num_segregating_sites(gene_matrix):
	"""
	Input snp matrix
	Returns the raw number of segregating sites (polymorphic sites).
	Sum over collumns, if sum != 0 or sum != nrow(matrix) : segregating site
	"""

	#for i in len(gene_matrix.shape[0] - 1):
	#gene_matrix = numpy.delete(gene_matrix, (0), axis=0)

	freqs = sp.mean(gene_matrix, 0)
	mafs = sp.minimum(freqs, 1 - freqs)

	maf_filter = mafs > 0.001
	mafs = mafs[maf_filter]

	sum_list = mafs * gene_matrix.shape[0]
	data = [float(Decimal("%.2f" % e)) for e in sum_list]
	
	SFS = Counter(data)

	total = sum(SFS.itervalues(), 0.0)
	SFS = {k: v / total for k, v in SFS.iteritems()}
	return SFS

def SFS_per_geno(geno_species=[argv[1], argv[1]], bin_size=0.2,
				 gt_hdf5_file= '/Users/PM/Desktop/New_data/newsnps_100.hdf5'):

	#[u'aacids', u'blosum62_scores', u'codon_snp_freqs', u'codon_snp_positions', u'codon_snps', u'codons', u'diversity',
	# u'dn_ds_ratio', u'freqs', u'is_synonimous_snp', u'norm_codon_snps', u'norm_snps', u'nts', u'num_non_syn_sites', 
	# u'num_syn_sites', u'num_vars', u'raw_snp_positions', u'raw_snps', u'snp_positions', u'snps', u'strains']

	pop = parse_pop_map()
	pop_map = pop.keys()
	ct_array = pop.values()

	h5f = h5py.File(gt_hdf5_file)

	ag = h5f

	gene_big_groups = sorted(ag.keys())
	gene_groups = list()

	# Taking just the core genes
	for gene in gene_big_groups:
		if len(ag[gene]['strains']) == 196:
			gene_groups.append(gene)

	print 'Number of genes analysed: %f' % len(gene_groups)
	# Names
	# sorted strains by genospecies

	strains_names = sorted(pop_map, key=lambda x: pop[x]['genospecies'])
	#print 'These are the strains evaluated', strains_names
	strains_names.remove('3260')
	strains_names.remove('3381')
	strains_names.remove('3339')
	strains_names.remove('3211')
	strains_list = strains_names

	tajimas_D = []
	tajimas_D_syn = []
	tajimas_D_non = []
	g1_list = []
	g1_syn = []
	g1_non = []
	

	for gene in gene_groups:

		#print 'Working on gene group: %s'%gene
		# Gene strain list
		strains_list = ag[gene]['strains'][...]

		# Looking at specific genospecies
		gs_list = []
		for strain in strains_list:
			gs_list.append((pop[strain]['genospecies'])) # +pop[strain]['origin'])

		# Transforming the strain list in array
		strains_list = np.asarray(strains_list)
		gs_filter1, gs_filter2 = [sp.in1d(gs_list,[gs]) for gs in geno_species]
			
		gs1 = strains_list[gs_filter1]
		gs2 = strains_list[gs_filter2]
		total_gs = np.append(gs1,gs2)

		# Transforming the strain list in array
		strains_list = np.asarray(strains_list)
		gs_filter1, gs_filter2 = [sp.in1d(gs_list, gs) for gs in geno_species]

		# Extracting the nucleotide sequences
		g1 = ag[gene]
		g1 = g1['codon_snps'][...].T

		try:
			g1_geno = g1[gs_filter1, :]
			g1_list.append(g1_geno)

			syn_index = ag[gene]['is_synonimous_snp'][...]
			g1_syn.append(g1_geno[:,syn_index])

			#non_index = np.negative(syn_index) # it doesnt work anymore
			g1_non.append(g1_geno[:,~syn_index])
		except:
			pass

	g1_conc = np.concatenate(g1_list, axis = 1) 
	g1_syn_conc = np.concatenate(g1_syn, axis = 1)
	g1_non_conc = np.concatenate(g1_non, axis = 1)
	
	print 'all %f' % g1_conc.shape[1]	
	print 'synonymous %f' % g1_syn_conc.shape[1]
	print 'non-synonymous %f' % g1_non_conc.shape[1]

	sfs_codon =  num_segregating_sites(g1_conc)
	sfs_codon = pd.DataFrame.from_dict(sfs_codon, orient='index')

	sfs_syn =  num_segregating_sites(g1_syn_conc)
	sfs_syn = pd.DataFrame.from_dict(sfs_syn, orient='index')

	sfs_non =  num_segregating_sites(g1_non_conc)
	sfs_non = pd.DataFrame.from_dict(sfs_non, orient='index')

	
	df = pd.concat([sfs_codon, sfs_syn, sfs_non], axis=1)
	df.columns = ['Codon', 'Syn', 'Non-syn']

	print df
	t = 'SFS for genospecies %s organic field' % argv[1]
	ax = df.plot(kind='bar',  title = t)
	ax.set_xlabel('Minor allele count')
	ax.set_ylabel('Fraction of SNPs')
	ax.set_title('SFS of %s' % argv[1])
	plt.tight_layout()
	plt.savefig('SFS_%s.pdf' % argv[1])
	#plt.show()

#SFS_per_geno()

def gene_genospecies_corr(snps_hdf5_file = '/Users/PM/Desktop/New_data/newsnps_100.hdf5',
						  min_maf = 0.10, min_num_snps = 20):
	from itertools import izip
	h5f = h5py.File(snps_hdf5_file)
	gene_groups = sorted(h5f.keys())

	pop_map = parse_pop_map()
	unique_gs = ['gsA', 'gsB', 'gsC', 'gsD', 'gsE']

	avg_gene_genosp_ld_dict = {'all': {}, 'nonsyn': {}, 'syn': {}}
		
	for i, gg in enumerate(gene_groups):
		if i%100==0:
			print '%d: Gene %s'%(i,gg)  
		g = h5f[gg]

		#Filtering SNPs with small MAFs
		freqs = g['codon_snp_freqs'][...]
		mafs = sp.minimum(freqs,1-freqs)
		maf_filter = mafs>min_maf
		if sp.sum(maf_filter)>min_num_snps:            
			is_synonimous_snp = g['is_synonimous_snp'][...]
			is_nonsynonimous_snp = ~is_synonimous_snp
			syn_snp_filter = is_synonimous_snp*maf_filter
			nonsyn_snp_filter = is_nonsynonimous_snp*maf_filter

			if sp.sum(syn_snp_filter)>0:
				all_norm_snps = g['norm_codon_snps'][...]
				norm_snps = all_norm_snps[maf_filter]
				M,N = norm_snps.shape
				strains = g['strains']
				gs_list = sp.array([pop_map[strain]['genospecies'] for strain in strains])
				
				d = {}
				for gs in unique_gs:          
					gs_snp = sp.array(sp.in1d(gs_list,[gs]),dtype='single')
					gs_snp = (gs_snp - sp.mean(gs_snp))/sp.std(gs_snp)
					r_list = sp.dot(gs_snp,norm_snps.T)/float(N)
					r2_list = r_list**2
					assert M==len(r2_list), 'A bug detected.'
					d[gs] = {'mean_r2': sp.mean(r2_list), 'num_snps': M, 'r2s': r2_list}
				avg_gene_genosp_ld_dict['all'][gg] = d
					
				syn_snp_filter = is_synonimous_snp[maf_filter]
				nonsyn_snp_filter = ~syn_snp_filter
				M = sp.sum(nonsyn_snp_filter)
				if M>5:
					d = {}
					for gs in unique_gs:          
						gs_snp = sp.array(sp.in1d(gs_list,[gs]),dtype='single')
						gs_snp = (gs_snp - sp.mean(gs_snp))/sp.std(gs_snp)
						r2_list = avg_gene_genosp_ld_dict['all'][gg][gs]['r2s'][nonsyn_snp_filter]
						assert M==len(r2_list), 'A bug detected.'
						d[gs] = {'mean_r2': sp.mean(r2_list), 'num_snps': M, 'r2s': r2_list}
					avg_gene_genosp_ld_dict['nonsyn'][gg] = d
					
				M = sp.sum(syn_snp_filter)
				if sp.sum(syn_snp_filter)>5:
					d = {}
					for gs in unique_gs:          
						gs_snp = sp.array(sp.in1d(gs_list,[gs]),dtype='single')
						gs_snp = (gs_snp - sp.mean(gs_snp))/sp.std(gs_snp)
						r2_list = avg_gene_genosp_ld_dict['all'][gg][gs]['r2s'][syn_snp_filter]
						assert M==len(r2_list), 'A bug detected.'
						d[gs] = {'mean_r2': sp.mean(r2_list), 'num_snps': M, 'r2s': r2_list}
					avg_gene_genosp_ld_dict['syn'][gg] = d
	return avg_gene_genosp_ld_dict

def summarize_genospecies_correlations(snps_hdf5_file = '/Users/PM/Desktop/New_data/newsnps_100.hdf5', 
								 seq_file = '/Users/PM/Desktop/New_data/final_snps.hdf5',
								 fig_dir = '/Users/PM/Desktop/New_data/final_figures',
								 geno_species='gsA', bin_size=0.2):
	h5f = h5py.File(snps_hdf5_file)
	sh5f = h5py.File(seq_file)
	gene_groups = sorted(h5f.keys())
	avg_gene_genosp_ld_dict = gene_genospecies_corr()
	#Now plot things per genospecies...
	dn_ds_ratios = {'all':[], 'nonsyn':[], 'syn':[]}
	num_seg_sites = {'all':[], 'nonsyn':[], 'syn':[]}
	pi_diversity = {'all':[], 'nonsyn':[], 'syn':[]}
	mean_r2s = {'all':[], 'nonsyn':[], 'syn':[]}
	for gg in gene_groups:
		g = h5f[gg]
		print len(g['nts'][...])
		#sg = sh5f['snps'][gg]
		for snp_type in ['all','nonsyn','syn']:
			ok_genes = set(avg_gene_genosp_ld_dict[snp_type].keys())
			if gg in ok_genes:
				mean_r2s[snp_type].append(float(avg_gene_genosp_ld_dict[snp_type][gg][geno_species]['mean_r2']))
				dn_ds_ratio = g['dn_ds_ratio'][...]
				dn_ds_ratios[snp_type].append(float(dn_ds_ratio))
				
				raw_snp_positions = g['raw_snp_positions'][...]
				num_seg_sites_per_base = len(raw_snp_positions)/len(g['nts'][...])
				num_seg_sites[snp_type].append(float(num_seg_sites_per_base))
				 
				diversity = g['diversity'][...]
				pi_diversity[snp_type].append(float(diversity))

	
	for snp_type in ['all','nonsyn','syn']:
		avg_r2s = sp.array(mean_r2s[snp_type])
		dn_ds_list = sp.array(dn_ds_ratios[snp_type])
		
		pl.clf()
#         bins = sp.arange(0,1+bin_size,bin_size)
		bins = [0.0,0.4,1.0]
		digitize = sp.digitize(avg_r2s, bins)    
		xs = []
		ys = []
		for bin_i in range(len(bins)-1):
			bin_filter = digitize==(bin_i+1)
			if len(dn_ds_list[bin_filter])>0:
#                 xs.append(sp.mean(avg_r2s[bin_filter]))
#                 ys.append(sp.mean(dn_ds_list[bin_filter]))
				xs.append(bins[bin_i]+0.5*(bins[bin_i+1]-bins[bin_i]))
				ys.append(dn_ds_list[bin_filter])

#         pl.plot(xs, ys, 'k.', alpha=0.3)
		fig = pl.figure(1, figsize=(9, 6))
		ax = fig.add_subplot(111)
		ax.boxplot(ys)
		ax.set_xticklabels(['%0.1f'%i for i in xs])
		pl.xlabel('Mean r2 between SNPs within a gene and %s'%geno_species)
		pl.ylabel(r'$\frac{K_a}{K_s}$')
		pl.savefig(fig_dir+'/Ka_Ks_vs_%s_corr_%s.png'%(geno_species,snp_type))
		
		nss_list = sp.array(num_seg_sites[snp_type])
		xs = []
		ys = []
		for bin_i in range(len(bins)):
			bin_filter = digitize==(bin_i+1)
			if len(nss_list[bin_filter])>0:
				xs.append(sp.mean(avg_r2s[bin_filter]))
				ys.append(sp.mean(nss_list[bin_filter]))

		pl.clf()
		pl.plot(xs, ys, 'k.', alpha=0.3)
		pl.xlabel('Mean r2 between SNPs within a gene and %s'%geno_species)
		pl.ylabel(r'Number of segregating sites per nucleotide ($S$)')
		pl.savefig(fig_dir+'/Num_seg_sites_vs_%s_corr_%s.png'%(geno_species,snp_type))
		
		pi_div_list = sp.array(pi_diversity[snp_type])
		xs = []
		ys = []
		for bin_i in range(len(bins)):
			bin_filter = digitize==(bin_i+1)
			if len(pi_div_list[bin_filter])>0:
				xs.append(sp.mean(avg_r2s[bin_filter]))
				ys.append(sp.mean(pi_div_list[bin_filter]))

		pl.clf()
		pl.plot(xs, ys, 'k.', alpha=0.3)
		pl.xlabel('Mean r2 between SNPs within a gene and %s'%geno_species)
		pl.ylabel(r'Nucleotide diversity ($\pi$)')
		pl.savefig(fig_dir+'/nt_diversity_vs_%s_corr_%s.png'%(geno_species,snp_type))

#summarize_genospecies_correlations()   


def gen_ld_plots(snps_hdf5_file = '/Users/PM/Desktop/New_data/newsnps_100.hdf5', 
				 max_dist=3000, min_maf=0.1, bin_size=30,
				 fig_dir = '/Users/PM/Desktop/New_data/final_figures', filter_pop=None,
				 fix_syn_nonsyn_ratio=True):
	
	pop_map = parse_pop_map()

	from itertools import izip
	h5f = h5py.File(snps_hdf5_file)
	gene_groups = sorted(h5f.keys())
	ld_dist_dict = {'all':{}, 'nonsyn':{}, 'syn':{}}
	distances = range(1,max_dist)
	
	for dist in distances:
		ld_dist_dict['all'][dist]={'r2_sum':0.0, 'snp_count':0.0}
		ld_dist_dict['nonsyn'][dist]={'r2_sum':0.0, 'snp_count':0.0}
		ld_dist_dict['syn'][dist]={'r2_sum':0.0, 'snp_count':0.0}
	
	for i, gg in enumerate(gene_groups):
		if i%100==0:
			print '%d: Gene %s'%(i,gg)  
		g = h5f[gg]

		if g['codon_snp_freqs'].size>1:
			if filter_pop is not None:
				strains = g['strains']
				indiv_filter = sp.zeros((len(strains)),dtype='bool8')
				for s_i, s in enumerate(strains):
					try:
						if pop_map[s]['genospecies']==filter_pop:
							indiv_filter[s_i]=True
					except:
						continue
				if sp.sum(indiv_filter)<20:
					continue
	
			
			if filter_pop is not None:
				codon_snps = g['codon_snps'][...]
				codon_snps = codon_snps[:,indiv_filter]
				norm_codon_snps = sp.transpose(codon_snps)
				freqs = sp.mean(norm_codon_snps,0)
				norm_codon_snps = (norm_codon_snps-freqs)/sp.sqrt(freqs*(1-freqs))
				norm_codon_snps = sp.transpose(norm_codon_snps)
				mafs = sp.minimum(freqs,1-freqs)
				maf_filter = mafs>min_maf
				if sp.sum(maf_filter)>1:
					is_synonimous_snp = g['is_synonimous_snp'][...]
					is_nonsynonimous_snp = ~is_synonimous_snp
					syn_snp_filter = is_synonimous_snp*maf_filter
					nonsyn_snp_filter = is_nonsynonimous_snp*maf_filter
	
					if sp.sum(syn_snp_filter)>sp.sum(nonsyn_snp_filter):
						all_norm_snps = norm_codon_snps
						all_positions = g['codon_snp_positions'][...]
						norm_snps = all_norm_snps[maf_filter]
						positions = all_positions[maf_filter]
						M,N = norm_snps.shape
						
						ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
						assert M==len(positions), 'A bug detected.'
						for i in range(M-1):
							for j in range(i+1,M):
								dist = positions[j] - positions[i]
								if dist<max_dist:
									ld_dist_dict['all'][dist]['r2_sum']+=ld_mat[i,j]**2
									ld_dist_dict['all'][dist]['snp_count']+=1.0
		
						
						if sp.sum(nonsyn_snp_filter)>1:
							norm_snps = all_norm_snps[nonsyn_snp_filter]
							positions = all_positions[nonsyn_snp_filter]
							M,N = norm_snps.shape
							
							ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
							assert M==len(positions), 'A bug detected.'
							for i in range(M-1):
								for j in range(i+1,M):
									dist = positions[j] - positions[i]
									if dist<max_dist:
										ld_dist_dict['nonsyn'][dist]['r2_sum']+=ld_mat[i,j]**2
										ld_dist_dict['nonsyn'][dist]['snp_count']+=1.0
									 
						snps_sample_size = sp.sum(nonsyn_snp_filter)
						   
						if sp.sum(syn_snp_filter)>1:
							norm_snps = all_norm_snps[syn_snp_filter]
							positions = all_positions[syn_snp_filter]
							
							sample_indices = sorted(sp.random.choice(sp.arange(len(positions)), snps_sample_size, replace=False))
							norm_snps = norm_snps[sample_indices]
							positions = positions[sample_indices]
							M,N = norm_snps.shape
							
							ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
							assert M==len(positions), 'A bug detected.'
							for i in range(M-1):
								for j in range(i+1,M):
									dist = positions[j] - positions[i]
									if dist<max_dist:
										ld_dist_dict['syn'][dist]['r2_sum']+=ld_mat[i,j]**2
										ld_dist_dict['syn'][dist]['snp_count']+=1.0            
			else:
				#Filtering SNPs with small MAFs
				freqs = g['codon_snp_freqs'][...]
				mafs = sp.minimum(freqs,1-freqs)
				maf_filter = mafs>min_maf
				if sp.sum(maf_filter)>1:
					is_synonimous_snp = g['is_synonimous_snp'][...]
					is_nonsynonimous_snp = ~is_synonimous_snp
					syn_snp_filter = is_synonimous_snp*maf_filter
					nonsyn_snp_filter = is_nonsynonimous_snp*maf_filter
	
					if sp.sum(syn_snp_filter)>sp.sum(nonsyn_snp_filter):
						all_norm_snps = g['norm_codon_snps'][...]
						all_positions = g['codon_snp_positions'][...]
						norm_snps = all_norm_snps[maf_filter]
						positions = all_positions[maf_filter]
						M,N = norm_snps.shape
						
						ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
						assert M==len(positions), 'A bug detected.'
						for i in range(M-1):
							for j in range(i+1,M):
								dist = positions[j] - positions[i]
								if dist<max_dist:
									ld_dist_dict['all'][dist]['r2_sum']+=ld_mat[i,j]**2
									ld_dist_dict['all'][dist]['snp_count']+=1.0
		
						
						if sp.sum(nonsyn_snp_filter)>1:
							norm_snps = all_norm_snps[nonsyn_snp_filter]
							positions = all_positions[nonsyn_snp_filter]
							M,N = norm_snps.shape
							
							ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
							assert M==len(positions), 'A bug detected.'
							for i in range(M-1):
								for j in range(i+1,M):
									dist = positions[j] - positions[i]
									if dist<max_dist:
										ld_dist_dict['nonsyn'][dist]['r2_sum']+=ld_mat[i,j]**2
										ld_dist_dict['nonsyn'][dist]['snp_count']+=1.0
									 
						snps_sample_size = sp.sum(nonsyn_snp_filter)
						   
						if sp.sum(syn_snp_filter)>1:
							norm_snps = all_norm_snps[syn_snp_filter]
							positions = all_positions[syn_snp_filter]
							
							sample_indices = sorted(sp.random.choice(sp.arange(len(positions)), snps_sample_size, replace=False))
							norm_snps = norm_snps[sample_indices]
							positions = positions[sample_indices]
							M,N = norm_snps.shape
							
							ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
							assert M==len(positions), 'A bug detected.'
							for i in range(M-1):
								for j in range(i+1,M):
									dist = positions[j] - positions[i]
									if dist<max_dist:
										ld_dist_dict['syn'][dist]['r2_sum']+=ld_mat[i,j]**2
										ld_dist_dict['syn'][dist]['snp_count']+=1.0
	
	
	
	for plot_type in ld_dist_dict.keys():
		avg_r2s = []
		dist_0_r2s = []
		dist_1_r2s = []
		dist_2_r2s = []
		plot_distances = []
		dist_0s = []
		dist_1s = []
		dist_2s = []
		for dist in distances:
			if ld_dist_dict[plot_type][dist]['snp_count']>0:
				avg_r2 = ld_dist_dict[plot_type][dist]['r2_sum']/float(ld_dist_dict[plot_type][dist]['snp_count'])
				avg_r2s.append(avg_r2)
				plot_distances.append(dist)
				if dist%3==0:
					dist_0_r2s.append(avg_r2)
					dist_0s.append(dist)
				elif dist%3==1:
					dist_1_r2s.append(avg_r2)
					dist_1s.append(dist)
				elif dist%3==2:
					dist_2_r2s.append(avg_r2)
					dist_2s.append(dist)
			
		plot_distances = sp.array(plot_distances)
		dist_0s = sp.array(dist_0s)
		dist_1s = sp.array(dist_1s)
		dist_2s = sp.array(dist_2s)
		avg_r2s = sp.array(avg_r2s)
		dist_0_r2s = sp.array(dist_0_r2s)
		dist_1_r2s = sp.array(dist_1_r2s)
		dist_2_r2s = sp.array(dist_2_r2s)
	
		
		pl.clf()
		print plot_distances 
		bins = sp.arange(0,max(plot_distances),bin_size)
		digitize = sp.digitize(plot_distances, bins)    
		xs = []
		ys = []
		for bin_i in range(len(bins)):
			bin_filter = digitize==(bin_i+1)
			if len(plot_distances[bin_filter])>0:
				xs.append(sp.mean(plot_distances[bin_filter]))
				ys.append(sp.mean(avg_r2s[bin_filter]))
		
		pl.plot(xs, ys, color='k', linestyle='None', marker='.', alpha=0.5)
		pl.xlabel(r'Pairwise distance ($d$)')
		pl.ylabel(r'Squared correlation ($r^2$)')
		if filter_pop is not None:
			pl.title('LD decay of %s' % argv[1])
			pl.savefig('%s/ld_%s_codons_%s.pdf'%(fig_dir,plot_type,filter_pop))
		else:
			pl.title('Overall LD decay')
			pl.savefig('%s/ld_%s_codons.pdf'%(fig_dir,plot_type))
	
		pl.clf()
		bins = sp.arange(0,max(dist_1s),bin_size)
		digitize = sp.digitize(dist_1s, bins)    
		xs = []
		ys = []
		for bin_i in range(len(bins)):
			bin_filter = digitize==(bin_i+1)
			if len(dist_1s[bin_filter])>0:
				xs.append(sp.mean(dist_1s[bin_filter]))
				ys.append(sp.mean(dist_1_r2s[bin_filter]))
		pl.plot(xs,ys, linestyle='None', marker='.', color='green', alpha=0.5, label=r'$d$ mod $3 = 1$')
	
		bins = sp.arange(0,max(dist_2s),bin_size)
		digitize = sp.digitize(dist_2s, bins)    
		xs = []
		ys = []
		for bin_i in range(len(bins)):
			bin_filter = digitize==(bin_i+1)
			if len(dist_2s[bin_filter])>0:
				xs.append(sp.mean(dist_2s[bin_filter]))
				ys.append(sp.mean(dist_2_r2s[bin_filter]))    
		pl.plot(dist_2s,dist_2_r2s, linestyle='None', marker='.', color='red', alpha=0.5, label=r'$d$ mod $3 = 2$')
	
		bins = sp.arange(0,max(dist_0s),bin_size)
		digitize = sp.digitize(dist_0s, bins)    
		xs = []
		ys = []
		for bin_i in range(len(bins)):
			bin_filter = digitize==(bin_i+1)
			if len(dist_0s[bin_filter])>0:
				xs.append(sp.mean(dist_0s[bin_filter]))
				ys.append(sp.mean(dist_0_r2s[bin_filter]))    
		pl.plot(dist_0s, dist_0_r2s, linestyle='None', marker='.', color='blue', alpha=0.5, label=r'$d$ mod $3 = 0$')
		pl.xlabel(r'Pairwise distance ($d$)')
		pl.ylabel(r'Squared correlation ($r^2$)')
		pl.legend()
		if filter_pop is not None:
			pl.title('LD decay of %s' % argv[1])
			pl.savefig('%s/part_ld_%s_codons_%s.pdf'%(fig_dir,plot_type,filter_pop))
		else:
			pl.title('Overall LD decay')
			pl.savefig(fig_dir+'/part_ld_%s_codons2.pdf'%(plot_type))
gen_ld_plots(filter_pop=argv[1])
