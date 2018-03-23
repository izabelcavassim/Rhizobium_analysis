import scipy as sp
import h5py 
import os
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import pylab
# pylab.rcParams['legend.numpoints'] = 1

import collections
import pandas as pd
#import cPickle
#import gzip

nt_map = {'A':1, 'C':2, 'G':3, 'T':4, '-':5, 'N':6}
nt_decode_map = {1:'A', 2:'C', 3:'G', 4:'T', 5:'-', 6:'N'}

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

def translate_dna(sequence):

	proteinsequence = ''
#     start = sequence.find('ATG')
#     sequencestart = sequence[int(start):]
#     stop = sequencestart.find('TAA')
#     cds = str(sequencestart[:int(stop)+3])
	cds=sequence
	
	for i in range(0,len(cds),3):
		codon = cds[i:i+3]

		proteinsequence += codontable.get(codon,'X')
	return proteinsequence


def parse_fasta(filename):
	file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
	file_separe = file.split('>') #spliting each entry by the > 
	file_separe.remove('')
	data_dict = {'iids':[], 'sequences':[], 'psequences':[], 'nsequences':[]}
	header = []

	for entry in file_separe:
		seq = entry.splitlines()
		header = seq[0] #these are the first elements of the list
		l = header.split('|')

		data_dict['iids'].append(l[2])

		sequence = ''.join(seq[1:]) #joining the sequences 
		data_dict['sequences'].append(sequence)
		psequence = translate_dna(sequence)
		data_dict['psequences'].append(psequence)
		nsequence = list(map(lambda x: nt_map[x], sequence))
		data_dict['nsequences'].append(nsequence)
	return data_dict 

def parse_fasta_file(filename):
	
	header = None
	sequence = ''
	data_dict = {'iids':[], 'sequences':[], 'psequences':[], 'nsequences':[]}
	with open(filename) as f:
		for line in f:
			print line
			if line[0] == ">":
				print line
				#if header is not None:
				assert len(sequence)%3==0, 'Sequence length is not divisible by 3.'
				l = header.split('|')
				print l
				data_dict['iids'].append(l[2])
				data_dict['sequences'].append(sequence)
				psequence = translate_dna(sequence)
				data_dict['psequences'].append(psequence)
				nsequence = list(map(lambda x: nt_map[x], sequence))
				data_dict['nsequences'].append(nsequence)
				header = line.strip()
				sequence = ''
			else:
				sequence += line.strip()
	return data_dict

def get_codon_syn_map():
	all_codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 
				  'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 
				  'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 
				  'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 
				  'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 
				  'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 
				  'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 
				  'TTT']
	
	all_bases = {'A','T','C','G'}

	ret_dict = {}
	
	for codon in all_codons:
		syn_list = []
		for i in range(3):
			base = codon[i]
			other_bases = all_bases - {base}
			syn = 0
			for new_base in other_bases:
				new_codon = codon[:i] + new_base + codon[i + 1:]
				syn += int(codontable[codon]==codontable[new_codon])
			syn_list.append(syn/3.0)
		ret_dict[codon]=sp.sum(syn_list)
	return ret_dict

def parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
	#from itertools import zip
	
	pop_map = {}
	t = pd.read_table(file_name)
	t = t.rename(columns=lambda x: x.strip())
	for strain_id, sara_id, origin, country in zip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Origin2']):
		pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country}
	return pop_map


def gen_genotype_hdf5file(out_hdf5_file ='/Users/PM/Desktop/New_data/final2_snps.hdf5', 
						  snps_directory='/Users/PM/Desktop/New_data/group_snps/',
						  fna_files_directory='/Users/PM/Desktop/New_data/group_alns/',
						  #fna_files_directory= '/Users/PM/Desktop/PHD_incomplete/Methods/group_alns/',
						  snps_wo_struct_directory='/Users/PM/Desktop/New_data/group_snps_corrected_structure/'):

	h5f = h5py.File(out_hdf5_file)
	snps_g = h5f.create_group('snps')
	snps_wostruct_g = h5f.create_group('snps_wo_struct')
	align_g = h5f.create_group('alignments')
	
	print("Parsing alignments")
	aln_files = os.listdir(fna_files_directory)
	print(aln_files)
	for i, a_f in enumerate(aln_files[1::]):
		if i%100==0:
			print i
		l = a_f.split('.')
		group_num = l[0][5:]

		aln_dict = parse_fasta(fna_files_directory+'/'+a_f)

		g = align_g.create_group('%s'%group_num)
		g.create_dataset('strains', data=aln_dict['iids'])
		g.create_dataset('sequences', data=aln_dict['sequences'], compression='lzf')
		g.create_dataset('psequences', data=aln_dict['psequences'], compression='lzf')
		g.create_dataset('nsequences', data=sp.array(aln_dict['nsequences'], dtype='int8'))


	print "Raw SNPs files"
	snps_files = os.listdir(snps_directory)
	for i, snps_f in enumerate(snps_files):
		if i%100==0:
			print i
		l = snps_f.split('.')
		group_num = l[0][5:]
		snps_data = sp.load(snps_directory+'/'+snps_f)

		#print snps_data['strains']
		g = snps_g.create_group('%s'%group_num)
		g.create_dataset('snps', data=snps_data['matrix'], compression='lzf')
		g.create_dataset('alignment_length',data=snps_data['alignment_length'])
		g.create_dataset('minor_frequencies',data=snps_data['minor_frequencies'])
		g.create_dataset('strains',data=[a.encode('utf8') for a in snps_data['strains']])
		g.create_dataset('positions',data=snps_data['positions'])
		h5f.flush()

	print "Now SNPs wo structure"
	snps_files = os.listdir(snps_wo_struct_directory)
	for i, snps_f in enumerate(snps_files):
		if i%100==0:
			print(i)
		l = snps_f.split('.')
		group_num = l[0][5:]
		snps_data = sp.load(snps_wo_struct_directory+'/'+snps_f)
		g = snps_wostruct_g.create_group('%s'%group_num)
		g.create_dataset('snps', data=snps_data['matrix'], compression='lzf')        
		g.create_dataset('alignment_length',data=snps_data['alignment_length'])
		g.create_dataset('minor_frequencies',data=snps_data['minor_frequencies'])
		g.create_dataset('strains',data=[a.encode('utf8') for a in snps_data['strains']])
		g.create_dataset('positions',data=snps_data['positions'])

	
	h5f.close()

def call_good_snps(sequence, ok_snps, snp_positions, codon_syn_map=None, ok_seq_filter=None, seq_num_vars=None):
	
	M,N = ok_snps.shape
	snps = []
	nts = []
	
	codon_snps = []
	codon_snp_positions = []
	codons = []
	aacids = []
	is_synonimous_snp =  []
	tot_num_syn_sites = 0
	tot_num_non_syn_sites = 0    
	for ok_snp, snp_pos in zip(ok_snps, snp_positions):                    
		mean_snp = sp.mean(ok_snp)
		snp = sp.zeros(N)
		snp[ok_snp>mean_snp]=1
		snps.append(snp)
		
		#Get nucleotides 
		nt0 = nt_decode_map[ok_snp.min()]
		nt1 = nt_decode_map[ok_snp.max()]
		nts.append([nt0,nt1])
		
		#5. Check codon position
		codon_pos = snp_pos%3
		
		#6. Identify non-ambiguous codon changes.
		#Check if there is a preceding/succeeding SNP within the codon.
		if codon_pos==0:
			if not (ok_seq_filter[snp_pos+1] and ok_seq_filter[snp_pos+2]):
				continue
			if not(seq_num_vars[snp_pos+1]==1 and seq_num_vars[snp_pos+2]==1):
				continue
			cdseq12 = sequence[snp_pos+1:snp_pos+3]
			codon0 = nt0+cdseq12
			codon1 = nt1+cdseq12
			
		elif codon_pos==1:
			if not (ok_seq_filter[snp_pos-1] and ok_seq_filter[snp_pos+1]):
				continue
			if not(seq_num_vars[snp_pos-1]==1 and seq_num_vars[snp_pos+1]==1):
				continue
			cdseq0 = sequence[snp_pos-1]
			cdseq2 = sequence[snp_pos+1]
			codon0 = cdseq0+nt0+cdseq2
			codon1 = cdseq0+nt1+cdseq2
		
		elif codon_pos==2:
			if not (ok_seq_filter[snp_pos-1] and ok_seq_filter[snp_pos-2]):
				continue
			if not(seq_num_vars[snp_pos-1]==1 and seq_num_vars[snp_pos-2]==1):
				continue
			cdseq01 = sequence[snp_pos-2:snp_pos]
			codon0 = cdseq01+nt0
			codon1 = cdseq01+nt1
		
		assert codon0!=codon1, 'Codons are identical?'
		
		#This appears to be a unique codon change with a dimorphic SNP.
		codons.append([codon0,codon1])
		freq = sp.mean(snp,0)
		
		#7. Check non-synonimous/synonimous
		num_syn_sites = freq*codon_syn_map[codon0]+(1-freq)*codon_syn_map[codon1]
		num_non_syn_sites = 3-num_syn_sites
		tot_num_syn_sites += num_syn_sites
		tot_num_non_syn_sites += num_non_syn_sites
		
		aa0 = codontable[codon0]
		aa1 = codontable[codon1]
		aacids.append([aa0,aa1])
		is_synon = aa0==aa1
		is_synonimous_snp.append(is_synon)
				
		codon_snps.append(snp)
		codon_snp_positions.append(snp_pos)
		
	return {'snps':snps, 'nts':nts, 'codon_snps':codon_snps, 'codon_snp_positions':codon_snp_positions, 
			'codons':codons, 'aacids': aacids, 'is_synonimous_snp':is_synonimous_snp,
			'num_syn_sites':tot_num_syn_sites, 'num_non_syn_sites':tot_num_non_syn_sites, }

def call_variants(gt_hdf5_file='/Users/PM/Desktop/New_data/final_snps.hdf5', 
				  out_file='/Users/PM/Desktop/New_data/newsnps_100.hdf5', min_num_strains=100):
				  #blosum62_file='/project/NChain/faststorage/rhizobium/ld/blosum62.txt'):
	"""
	Generate a new set of SNPs to look at.
	
	For all nts:
		if it is a SNP
			count # of variants. 
			check AA changes
			quantify AA change severity    
	
	"""
	pop_map = parse_pop_map()
	print pop_map
	
	from itertools import izip
#     blosum62_matrix, blosum62_dict = parse_blosum62(blosum62_file)
	codon_syn_map = get_codon_syn_map()
	h5f = h5py.File(gt_hdf5_file)
	ag = h5f['alignments']
	oh5f = h5py.File(out_file)
	gene_groups = sorted(ag.keys())
	num_parsed_genes = 0
	for gg in gene_groups:
		g = ag[gg]
		
		#0. Check if there is evidence for CNVs/paralogs?
		seq_ids = g['strains']
		strains_list = map(lambda x: x.split('-')[0], seq_ids)
		strains, strain_counts = sp.unique(strains_list, return_counts=True)
		if len(strains)<len(strains_list):
			print 'Evidence for paralogs/CNVs'
			print strain_counts
			print '%d strains have unique gene copies'.format(len(strains))
		elif len(seq_ids)>=min_num_strains:
			strains = map(lambda x: x[0:4], seq_ids)
						
			#1. Filter indel/bad rows
			nt_mat = g['nsequences'][...]
			num_vars = sp.apply_along_axis(lambda x: len(sp.unique(x)), 0, nt_mat)
			no_gaps_no_missing = sp.all(nt_mat<5,0)
			nt_mat = sp.transpose(nt_mat)
			bad_rows_filter = (num_vars<5)*no_gaps_no_missing
			if sp.sum(bad_rows_filter)>0:
				print 'passed bad filter control'
				raw_snps = nt_mat[bad_rows_filter]
				
				#Calculate nucleotide diversity and ani
				M,N = raw_snps.shape
				diversity = 0.0
				ani=0.0
				for i in range(N-1):
					for j in range(i+1,N):
						diversity+=sp.sum(raw_snps[:,i]!=raw_snps[:,j])
						ani+=sp.sum(raw_snps[:,i]==raw_snps[:,j])
				
				diversity = diversity/len(raw_snps)
				diversity = 2*diversity/(N*(N-1.0))
				ani = ani/len(raw_snps)
				ani = 2*ani/(N*(N-1.0))
				
				#2. Filter non-variable rows
				ok_num_vars = num_vars[bad_rows_filter]
				var_filter = ok_num_vars>1                
				num_raw_snps = sp.sum(var_filter)
				if num_raw_snps>0:
					print 'Working on gene group: %s'%gg
					
					M,N = nt_mat.shape
					non_gap_positions = sp.arange(M)[bad_rows_filter]
					all_snps = raw_snps[var_filter]
					all_snp_positions = non_gap_positions[var_filter]
					
					#3. Identify good SNPs (dimorphic SNPs)
					good_snp_filter = ok_num_vars==2
					ok_snps = raw_snps[good_snp_filter]
					snp_positions = non_gap_positions[good_snp_filter]
					assert len(ok_snps)==len(snp_positions), 'A bug detected!'
					
					#4. Call good SNPs
					good_snps_dict = call_good_snps(g['sequences'][0], ok_snps, snp_positions, codon_syn_map=codon_syn_map,
									ok_seq_filter = no_gaps_no_missing, seq_num_vars=num_vars)
					
					snps = good_snps_dict['snps']
					nts = good_snps_dict['nts']
					codon_snps = good_snps_dict['codon_snps']
					codon_snp_positions = good_snps_dict['codon_snp_positions']
					codons = good_snps_dict['codons']
					aacids = good_snps_dict['aacids']
					is_synonimous_snp = good_snps_dict['is_synonimous_snp']
					num_syn_sites = good_snps_dict['num_syn_sites']
					num_non_syn_sites = good_snps_dict['num_non_syn_sites']

					
					#Normalize SNPs
					norm_snps = sp.transpose(snps)
					freqs = sp.mean(norm_snps,0)
					norm_snps = (norm_snps-freqs)/sp.sqrt(freqs*(1-freqs))
					norm_snps = sp.transpose(norm_snps)
					
					
					norm_codon_snps = sp.transpose(codon_snps)
					codon_snp_freqs = sp.mean(norm_codon_snps,0)
					norm_codon_snps = (norm_codon_snps-codon_snp_freqs)/sp.sqrt(codon_snp_freqs*(1-codon_snp_freqs))
					norm_codon_snps = sp.transpose(norm_codon_snps)
					
					#Calculate dn/ds ratios
					num_syn_subt = sp.sum(is_synonimous_snp)
					num_non_syn_subt = len(is_synonimous_snp)-num_syn_subt
					if num_syn_subt>0:
						dn_ds_ratio = (num_non_syn_subt/num_non_syn_sites)/(num_syn_subt/num_syn_sites)
					else:
						dn_ds_ratio=-1


					#Calculate McDonald-Kreitman Statistics..
					
					#Store everything to a HDF5 file
					og = oh5f.create_group(gg)   
					og.create_dataset('num_vars', data=num_vars)
					og.create_dataset('raw_snps', data=sp.array(all_snps,dtype='int8'), compression='lzf')
					og.create_dataset('raw_snp_positions', data=all_snp_positions)
					og.create_dataset('snps', data=sp.array(snps,dtype='int8'), compression='lzf')
					og.create_dataset('norm_snps', data=sp.array(norm_snps,dtype='single'), compression='lzf')
					og.create_dataset('freqs', data=sp.array(freqs,dtype='single'))
					og.create_dataset('snp_positions', data=snp_positions)
					og.create_dataset('codon_snps', data=sp.array(codon_snps,dtype='single'), compression='lzf')
					og.create_dataset('norm_codon_snps', data=sp.array(norm_codon_snps,dtype='single'), compression='lzf')
					og.create_dataset('codon_snp_freqs', data=sp.array(codon_snp_freqs,dtype='single'))
					og.create_dataset('is_synonimous_snp', data=is_synonimous_snp)
					og.create_dataset('strains', data=strains)
					og.create_dataset('codon_snp_positions', data=codon_snp_positions)
#                     og.create_dataset('blosum62_scores', data=blosum62_scores)
					og.create_dataset('aacids', data=sp.array(aacids))
					og.create_dataset('nts', data=sp.array(nts))
					og.create_dataset('codons', data=sp.array(codons))
					og.create_dataset('num_syn_sites', data=num_syn_sites)
					og.create_dataset('num_non_syn_sites', data=num_non_syn_sites)
					og.create_dataset('dn_ds_ratio', data=dn_ds_ratio)
					og.create_dataset('diversity', data=diversity)
					og.create_dataset('ani', data=ani)
					oh5f.flush()
					num_parsed_genes +=1
		else:
			print 'Too few strains..'
	print 'Parsed %d'%num_parsed_genes


gen_genotype_hdf5file()
#call_variants()
