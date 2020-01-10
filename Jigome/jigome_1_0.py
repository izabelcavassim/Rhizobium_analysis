#!/usr/bin/env python
from __future__ import division
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import sys
import subprocess
import copy
import re
import os
import pickle
import datetime
import numpy as np

script_file = os.path.basename(__file__)
version = script_file[7:-3]

"""
GENERAL DESCRIPTION
This script takes sets of contigs assembled by SPAdes and carries out further assembly to 
create larger scaffolds.  Assembly is guided by the connections in the contig graph and by
homologous links found a set of reference genomes.  The chromosomal contigs are oriented
and ordered using a set of core genes that have the same order in all members of a set of
fully-assembled genomes. Plasmid replication genes are identified and used to label the
corresponding scaffolds.

PRINCIPAL OUTPUTS
scaffolds: the principal fasta files with the final scaffolds
scaffold_descriptions: indicates how input contigs are placed in the scaffolds
contig_counts: for each input contig, shows whether used, and if so how many times
logfile.txt: includes warnings of unexpected events

OTHER OUTPUTS 
clean_contigs: fasta files of the contigs after discarding contaminants
discarded_contigs: fasta file of contigs rejected as contaminant (low coverage)
contig_overlaps: linkages between contigs (overlapping ends, usually 127 bp)
marker_contigs: lists of contigs identified as chromosomal or plasmid
blast_files: results of blast searches for matches to reference genomes and genes

PROJECT-SPECIFIC FEATURES
This version was specifically written for genomes of Rhizobium leguminosarum sequenced by 
MicrobesNG on Illumina HiSeq using a 250bp paired-end protocol and assembled by SPAdes 3.7
with mostly 127 bp k-mers. Species-specific reference files are:
(1) a set of public genomes (not all completely assembled)
(2) a set of RepA protein sequences representing the known plasmid classes
(3) a set of chromosomal gene sequences (starting with dnaA) that are in a conserved order
(4) a representative sequence of the dnaA gene

The input contig names have a specific format, e.g. '>NODE_1_length_911879_cov_22.4383_ID_5807'
If a coverage threshold is selected (in order to remove low-coverage contaminant contigs),
the script uses the length and coverage figures from this header. If this information is
not present, select threshold = 0 to omit this stage, and remove the contaminants some other way.

DEPENDENCIES
Besides the libraries imported at the start, this script calls NCBI BLAST+. Developed in 
Python 3.6.5.

VERSION 
This version is functionally identical to version 0_1p, which was used to process the genomes
described by Cavassim et al. (2019), except that the block on using the strain's own PacBio
assembly as reference, introduced for fair testing, has been removed.

AUTHOR
Peter Young, December 2018

"""

#######
# FILES
#######

#Input folders and files
########################
#These are written into the script, since the script was developing but the data were constant!

#Folders containing the input fasta files and the outputs
source_folder = "/Users/peter/Documents/Sequence_data/NCHAIN_genomes/Illumina_assembly/Contigs_from_MicrobesNG/"
working_folder = "/Users/peter/Documents/Sequence_data/NCHAIN_genomes/Illumina_assembly/Jigome_output_" + version + "/"

#Reference sequences (these are specific to Rhizobium leguminosarum)
reference_db = "/Users/peter/Dropbox/NCHAIN_share/Rhizobium_data/reference_files/Rleg_47_genomes.fas"
rep_file = "/Users/Peter/Dropbox/NCHAIN_share/Rhizobium_data/reference_files/RepA_new20_standards.fas"
chr_file = "/Users/peter/Dropbox/NCHAIN_share/Rhizobium_data/reference_files/ordered_orthocore4.fas"
dnaA_file = "/Users/peter/Dropbox/NCHAIN_share/Rhizobium_data/reference_files/dnaA.fas"


#Output folders and files
#########################

if not os.path.isdir(working_folder): os.mkdir(working_folder)
os.mkdir(working_folder + "clean_contigs_" + version)
os.mkdir(working_folder + "discarded_contigs_" + version)
os.mkdir(working_folder + "contig_overlaps_" + version)
os.mkdir(working_folder + "scaffolds_" + version)
os.mkdir(working_folder + "contig_counts_" + version)
os.mkdir(working_folder + "marker_contigs_" + version)
os.mkdir(working_folder + "scaffold_descriptions_" + version)
os.mkdir(working_folder + "blast_files_" + version)


dnaA_blast_out = working_folder + "dnaA_blast.tab"
data_table_file = working_folder + "summary_" + version + ".tab"
logfile = working_folder + "logfile_" + version + ".txt"
blast_files_folder = working_folder + "blast_files_" + version + "/"


#the following command creates fasta_file_list.txt with a list of all .fasta files in the 
#source folder.  To use only certain files, comment this out and create the list manually

subprocess.call("ls " + source_folder + "*.fasta > " + working_folder +"fasta_file_list.txt", shell=True)

strain_list = []
with open(working_folder +"fasta_file_list.txt") as list_file:
    for strain_file in list_file:
        strain = strain_file.split("/")[-1].split(".")[0]
        strain_list.append(strain)

#strain_list = ["strain-1", "strain-2"]
strain_list = ["3206-3"]


############
# PARAMETERS
############

coverage_threshold = 0.3 #as a fraction of the median coverage of contigs >10kb 
#Low-coverage contigs will be removed at the outset unless coverage_threshold = 0
file_extension = ".fasta" #for the input contig files
contig_num_pos = 1 #Position of the contig number in the name (i.e. after the nth "_")
overlap = 127
tag_length = 500
max_mismatches = 36 #changed from 4 to allow 95-base overlaps
l_min_score = 850
r_min_score = 600
scaff_l_min_score = 400
scaff_r_min_score = 400

verbose = False
#Setting verbose = True creates additional outputs for test runs:
#files of unique, repeat, etc. contigs
#files of intermediate stage0, 1, 2 and final scaffolds
   

###########
# FUNCTIONS
###########

def paf(print_string):
    """
    Print and file: everything printed to stdout is also written in logfile.
    """
    print(print_string)
    with open(logfile, "a") as outfile:
        outfile.write(str(print_string) + "\n")

def get_median_coverage(contigs_dict):
    """
    Returns median coverage of the large contigs (>10kb)
    """
    import statistics
    coverage_list = []
    for contig in contigs_dict:
        if int(contig.split("_")[3]) >= 10000:
            coverage_list.append(float(contig.split("_")[5]))
    med_cov = statistics.median(coverage_list)
    return med_cov
        
    
def end_list(seqs, from_base, to_base):
    """
    Returns a short sequence from each end of each contig
    """
    list_of_ends = []
    for seq in seqs:
        fiveprime_end = seq[from_base:to_base]
        fiveprime_end.id = "fiveprime_"+seq.id
        list_of_ends.append(fiveprime_end)
        threeprime_end = seq.reverse_complement()[from_base:to_base]
        threeprime_end.id = "threeprime_"+seq.id
        list_of_ends.append(threeprime_end)
    return list_of_ends


def matching_ends(list_of_ends, max_mis):
    """
    Finds matches between contig ends and returns a dict of the matches
    """
    paf("starting matching_ends at " + str(datetime.datetime.now()))
    
    dict_of_matches = {}
    for n in range(len(list_of_ends)):
        seq1 = list_of_ends[n]
        neighbours = {}
        for m in range(n+1, len(list_of_ends)):
            seq2 = list_of_ends[m]
            alignment_worthwhile = ((seq2.seq[:30] in seq1.reverse_complement().seq) 
                                or (seq2.seq[30:60] in seq1.reverse_complement().seq)) 
            if alignment_worthwhile:
                alignments = pairwise2.align.localms(seq1.reverse_complement().seq, seq2.seq, 1, -100, -100, -100)
                gaps = alignments[0][0].count("-")
                if alignments[0][2] >= len(seq1.seq) - max_mis:
                    neighbours[seq2.id] = len(seq1.seq) - gaps
        dict_of_matches[seq1.id] = neighbours
                
    for end1 in dict_of_matches:
            for end2 in dict_of_matches[end1]:        
                dict_of_matches[end2][end1] = dict_of_matches[end1][end2]
        
    return dict_of_matches
    
    
def other_end(end1):
    """
    Changes scaffold end name from 'fiveprime' to 'threeprime' and vice versa.
    """
    name_parts = end1.split("prime_")
    if name_parts[0] == "three":
        end2 = "fiveprime_" + name_parts[1]
    elif name_parts[0] == "five":
        end2 = "threeprime_" + name_parts[1]
    else:
        end2 = "other_end_error_" + end1       
    return end2
    
    
def other_side(end1):
    """
    Changes scaffold end name from 'left' to 'right' and vice versa.
    """
    name_parts = end1.split("_")
    if name_parts[0] == "left":
        end2 = "right_" + name_parts[1]
    elif name_parts[0] == "right":
        end2 = "left_" + name_parts[1]
    else:
        end2 = "other_side_error_" + end1       
    return end2

 
def order_chromosomal_contigs(chr_blast_output):
    """
    Creates an ordered list of chromosomal contigs based on hits in Blast output 
    """
    ordered_chr_contigs = []
    current_contig = "null"
    current_contig_direction = 0
    current_contig_hits = 0

    with open(chr_blast_output) as blast_matches:
        for hit in blast_matches:
            hit_data = hit.rstrip("\n").split("\t")
            core_gene_dir = int(hit_data[0].split("|")[1])
            if float(hit_data[2]) >= 90.0:
                new_contig = hit_data[1]
                new_contig_direction = core_gene_dir*np.sign(int(hit_data[9])-int(hit_data[8]))
            
                if new_contig == current_contig and new_contig_direction == current_contig_direction:
                    current_contig_hits += 1
                else:            
                    contig_tuple = (current_contig, current_contig_direction, current_contig_hits)
                    ordered_chr_contigs.append(contig_tuple)
                    current_contig = new_contig
                    current_contig_direction = new_contig_direction
                    current_contig_hits = 1

        contig_tuple = (current_contig, current_contig_direction, current_contig_hits)
        ordered_chr_contigs.append(contig_tuple)
        ordered_chr_contigs.pop(0)

    #If hits to a contig are not contiguous, keep only the longest run    
    chr_contig_dict = {}    #stores the longest run for each contig
    remove_list = []    #stores the shorter runs for deletion
    n = -1
    for entry in ordered_chr_contigs:
        n += 1
        contig = entry[0]
        hits = entry[2]
        if contig not in chr_contig_dict:
            chr_contig_dict[contig] = (n, entry)
        elif hits > chr_contig_dict[contig][1][2]:
            remove_list.append(chr_contig_dict[contig])
            chr_contig_dict[contig] = (n, entry)
        else:
            remove_list.append((n, entry))

    #The first contig will usually also be the last - both should be kept        
    for item in remove_list:
        
        if int(item[0]) == 0 or int(item[0]) == len(ordered_chr_contigs)-1:
            remove_list.remove(item)
            
    remove_list.sort(reverse = True)
    for item in remove_list:
        position = item[0]
        ordered_chr_contigs.pop(position)
    
    return ordered_chr_contigs


def censor_contig(contig_end, u_contigs, o_dict):
    """
    removes the entries for both ends of a contig in a contig list and overlap dict
    """
    for c_e in [contig_end, other_end(contig_end)]:
        if c_e in u_contigs:
            u_contigs.remove(c_e)
        if c_e in o_dict:
            o_dic = o_dict[c_e]
            if o_dic != {}:
                overlapped_contig = list(o_dic.keys())[0]
                if overlapped_contig in o_dict: del o_dict[overlapped_contig][c_e]
            del o_dict[c_e]
    return 

    
def get_next_contig(contig_list):
    """
    Removes and returns the right-hand end of the first contig in a list of contig tuples 
    (e.g. from order_chromosomal_contigs()), adding a fiveprime_ or threeprime_ prefix 
    depending on the direction given in the second item of the tuple
    """
    next_contig_tuple = contig_list.pop(0)
    next_contig = next_contig_tuple[0]
    if next_contig_tuple[1] == 1: next_contig = "threeprime_" + next_contig
    else: next_contig = "fiveprime_" + next_contig
    return next_contig
        
    
def grow_scaffold(scaffold, o_dict, u_cont, r_cont):
    """
    Extends a scaffold by adding contigs that are connected by matching end sequences
    using the overlap dict, and choosing the best alternatives using trio_hits()
    Calls trio_hits(), other_end()
    """
    free_end = True
    current_end = scaffold.pop()
    
    
    while free_end:
    
        next_contig_dict = o_dict[current_end]
        right_contig = "link_null"
        
        if len(next_contig_dict) == 1:
            next_contig = list(next_contig_dict.keys())[0]
            next_end = other_end(next_contig)
            scaffold.extend([current_end, next_contig])
                                        
            if len(scaffold) > 1 and next_end == scaffold[1]:
                    right_contig = "link_to_start_(2)"
                    free_end = False
            elif next_end in u_cont:

                if next_end in scaffold:
                    right_contig = "link_loop"
                    free_end = False
                else:
                    current_end = next_end
                    right_contig = other_end(next_end)

            else:              
                if next_end in r_cont and next_end in o_dict:           
                    right_contig = trio_hits(current_end, next_end, hit_list, o_dict, contigs_dict)                    
                    if right_contig in unique_contigs_master and right_contig in scaffold:                    
                        right_contig = "link_duplicate"
                else:                
                    right_contig = "link_not_found"

                if right_contig[0:4] != "link":
                    candidate = other_end(right_contig)
                    if not candidate in scaffold:
                        scaffold.extend([next_end, right_contig])

                        if candidate in o_dict and len(o_dict[candidate]) > 0:
                            current_end = candidate
                        else: 
                            right_contig = "link_not_in_overlaps_(2)"
                    else:
                        right_contig = "link_in_scaffold"
        
            if right_contig[0:4] == "link":
                scaffold.append(right_contig)
                free_end = False
     
        else:
            if len(next_contig_dict) == 0:
                scaffold.append("link_not_in_overlaps")
            elif len(next_contig_dict) >1:
                scaffold.append("link_ambiguous_(2)")
            else: scaffold.append("link_other") #should never happen!    
            free_end = False  
                        
    return scaffold


def trio_hits(l_contig, mid_contig_end, blast_hits, olap_dict, cont_dict):
    """
    From a unique contig (l_contig), steps across a repeat (mid_contig) to find the next unique contig,
    based on best case of adjacent hits in a reference genome (product of the two blast scores).
    Called by grow_scaffold()
    """
    mid_contig = mid_contig_end.split("prime_")[1]
    length_mid = len(cont_dict[mid_contig].seq)    
    right_links = []
    right_links = list(olap_dict[mid_contig_end].keys())
    
    #If contigs are chromosomal, ensure they are adjacent
    if chr_links:
        if l_contig in chr_links:        
            new_right_links= []
            for r_l in right_links:
                if r_l not in chr_links or r_l in chr_links[l_contig]:
                    new_right_links.append(r_l)
            right_links = new_right_links        
        
    if len(right_links) == 1:
        outcome = right_links[0]
    elif len(right_links) == 0:
        outcome = "link_listempty"
    else:
        left_matches = []
        for hit in blast_hits:
            if (hit[0] == l_contig) and (int(hit[11]) >= l_min_score):
                left_matches.append(hit)
        link_count = {}
        for link in right_links:
      
            right_matches = []

            for hit in blast_hits:
                if (hit[0] == link) and (int(hit[11]) >= r_min_score):                   
                    right_matches.append(hit)
                
            for lhit in left_matches:
                for rhit in right_matches:
                    if lhit[1] == rhit[1]:
                        lh_start = int(lhit[8])
                        lh_end = int(lhit[9])
                        rh_start = int(rhit[8])
                        rh_end = int(rhit[9])
                        if abs(lh_start - rh_start) < length_mid + 100:
                            if (lh_end - lh_start)/(rh_end - rh_start) < 0:
                                if abs(lh_end - rh_end) > abs(lh_start - rh_start):
                                    link_score = int(lhit[11]) * int(rhit[11])
                                    if not link in link_count: 
                                        link_count[link] = link_score
                                    elif link_score > link_count[link]:
                                        link_count[link] = link_score
                                                 
        number_of_matches = len(link_count)
        if number_of_matches == 1:
            outcome = list(link_count.keys())[0]
        if number_of_matches == 0:
            outcome = "link_unmatched"
        if number_of_matches > 1:
            outcome = max(link_count, key = link_count.get)

    return outcome


def unique_contigs_are_unique(scaffold_list, unique_contigs_list):
    """
    If unique contigs occur more than once, resolve iteratively until unique
    Calls resolve_unique_contigs()
    """
    i= 0
    old_scaffold_list = copy.deepcopy(scaffold_list)
    old_scaffold_list = purge_redundancy(old_scaffold_list)
    new_scaffold_list = []
    while new_scaffold_list != old_scaffold_list and i < 20:
        
        i += 1
        if i != 1: 
            old_scaffold_list = copy.deepcopy(new_scaffold_list)
            #new list is now old list
        new_scaffold_list = new_resolve_unique_contigs(old_scaffold_list, unique_contigs_list)    
        new_scaffold_list = purge_redundancy(new_scaffold_list)

    return new_scaffold_list
    

def new_resolve_unique_contigs(scaffold_list, unique_contigs_list):
    """    
    #find unique contigs used more than once and merge the scaffolds containing them
    #or, if conflicting, split them at points of conflict
    Called by unique_contigs_are_unique()
    Calls split_siamese(), find_unique_contig(), make_scaff_overlap_dict(), combine_overlapping_contigs()
    
    The original version failed when there was more than one point of conflict in a 
    scaffold because the unresolved conflicts were returned to the scaffold list.
    This new version deals with one conflicted contig at a time and uses the resulting
    scaffold list as the starting point for the next.
    """
    
    contig_location = {}
    s_l = copy.deepcopy(scaffold_list)
    
    #first deal with any scaffolds that have more than one copy of a unique contig
    to_remove = []
    for scaf in s_l:   
        for contig in unique_contigs_list:
            if scaf.count(contig) > 1:
                scaffold_parts = split_siamese(contig, scaf)
                to_remove.append(scaf)
                s_l.extend(scaffold_parts)
                break 
    for scaf in to_remove:
        s_l.remove(scaf) 


    for contig in unique_contigs_list:
        #if contig[:4] == "five": 
            finds = find_unique_contig(contig, s_l)

            if len(finds) > 1:
                contig_location[contig] = finds

    sc_ov = {}
    sc_ov = make_scaff_overlap_dict(contig_location)

    #This is the new bit that takes just the first conflicted contig   
    first_k = list(sc_ov.items())[0:1]
    first_sc_ov = dict(first_k)
    new_scaffold_list = combine_overlapping_contigs(first_sc_ov, s_l)

    #Split off unique scaffolds attached by their 3' ends to multiple scaffolds
    
    for contig in contig_location:
        if contig[:5] == "three":
            for scaf in contig_location[contig]:
                conflict = False
                if scaf.index(contig) == 1:
                    conflict = True
                    new_left_scaf = scaf[:3]
                    new_right_scaf = scaf[3:]
                if scaf.index(contig) == len(scaf) - 2:
                    conflict = True
                    new_left_scaf = scaf[:-3]
                    new_right_scaf = scaf[-3:]
                if conflict:
                    new_left_scaf.append("link_conflict6")
                    new_right_scaf.insert(0,"link_conflict6")
                    if len(new_left_scaf) >= 4: 
                        new_scaffold_list.append(new_left_scaf)
                    if len(new_right_scaf) >= 4:
                        new_scaffold_list.append(new_right_scaf)
                    if scaf in new_scaffold_list:
                        new_scaffold_list.remove(scaf)

    return new_scaffold_list

def split_siamese(contig, scaf):
    """
    Splits a contig that has more than one copy of a unique contig     
    Called by new_resolve_unique_contigs()
    """
    pos1 = scaf.index(contig)
    pos2 = len(scaf) - list(reversed(scaf)).index(contig) -1
    linker = "link_siamese"
    offset1, offset2 = 1,1
    if other_end(contig) in scaf:
        opos1 = scaf.index(other_end(contig))
        opos2 = len(scaf) - list(reversed(scaf)).index(other_end(contig)) -1
        if opos1 == pos1 - 1:
            offset1 = 0
        elif opos1 == pos1 + 1:
            offset1 = 1
        if opos2 == pos2 - 1:
            offset2 = 0
        elif  opos2 == pos2 + 1:
            offset2 = 1
    if pos1 == 1 and scaf[2] != other_end(contig):
        offset1 = 2
    left_scaf = scaf[:pos1+offset1] + [linker]
    mid_scaf = [linker] + scaf[pos1+offset1:pos2+offset2] + [linker]
    right_scaf = [linker] + scaf[pos2+offset2:]
    out_scaffold_list = [left_scaf, mid_scaf, right_scaf]
    
    return out_scaffold_list   
    

def find_unique_contig(contig, s_l):
    """
    Returns a list of scaffolds that contain a certain contig
    Called by new_resolve_unique_contigs()
    """    
    in_scaffold = []
    for scaffold in s_l:
        if contig in scaffold:
            in_scaffold.append(scaffold)
            
    return in_scaffold


def make_scaff_overlap_dict(contig_location):
    """
    Creates a dict of supposedly unique contigs that occur in more than one scaffold (keys)
    with the scaffolds that contain them (values)
    Called by new_resolve_unique_contigs()
    """
    scaffold_overlaps = []
    sc_ov = {}
    for contig in contig_location:
        
        if contig[:4] == "five": 

            if not contig_location[contig] in scaffold_overlaps:
                scaffold_overlaps.append(contig_location[contig])
                sc_ov[contig] = copy.deepcopy(contig_location[contig])
 
    #orient each scaffold so that contig k is fiveprime-threeprime
    #unless it is the first link in the scaffold
    # *** this will fail if the 'unique' contig occurs >1 time in the scaffold!
    # - but split_siamese should have taken care of that
    for k, v in sc_ov.items():
        for scaf in v:
    
            if scaf[1] == k or (other_end(k) in scaf and scaf.index(k) - scaf.index(other_end(k)) == 1):
                if k[:4] == "five": scaf.reverse()  

    return sc_ov

    
def combine_overlapping_contigs(sc_ov, scaffold_list):
    """
    Combines scaffolds that share unique contigs
    Called by new_resolve_unique_contigs()
    Calls split_at_conflict(), purge_redundancy()
    """         
    for k in sc_ov:
       
        conflict = False
        nos = len(sc_ov[k])
        sca_lis = []
        l_length = {}
        r_length = {}
        for n in range(nos):
           
            sca_lis.append(sc_ov[k][n])
            p = sca_lis[n].index(k)
            l_length[n] = p+1
            r_length[n] = len(sca_lis[n]) - p-1
        
        l_longest = max(l_length, key=l_length.get)
        r_longest = max(r_length, key=r_length.get)        
        new_scaff =   sca_lis[l_longest][:l_length[l_longest]] + sca_lis[r_longest][-r_length[r_longest]:]
        
        alt_scaff = []
        for n in range(nos):
            if str(sca_lis[n][1:-1])[1:-1] not in str(new_scaff): 
                conflict = True                
                split_scaffs = split_at_conflict(new_scaff, sca_lis[n], k)
                for scaff in split_scaffs:
                    if scaff not in alt_scaff:
                        alt_scaff.append(scaff)

        if not conflict:
            scaffold_list.append(new_scaff)
        else:                
            alt_scaff2 = purge_redundancy(alt_scaff) 
            for new_scaff in alt_scaff2:
                if len(new_scaff) > 2: #exclude empty scaffolds
                    scaffold_list.append(new_scaff)
                                
        for scaff in sca_lis:
            if scaff in scaffold_list:
                scaffold_list.remove(scaff)
            else:
                scaff.reverse()
                if scaff in scaffold_list:
                    scaffold_list.remove(scaff)
    
    return scaffold_list   


def split_at_conflict(scaff1, scaff2, k):
    """
    Takes two conflicting scaffolds (first is longer) sharing a contig k
    Splits them at the points of conflict and returns the split parts
    Called by combine_overlapping_contigs()
    """
    linker = "link_conflict5"
    p1 = scaff1.index(k)
    p2 = scaff2.index(k)
    
    l_split1=0
    r_split1=0
    l_split2=0
    r_split2=0
    out1=[]
    out2=[]
    out3=[]
    out_scaffs = []

    for q in range(len(scaff2)-p2-1):
        if scaff1[p1+q] != scaff2[p2+q]:
            r_split1 = p1+q-1
            r_split2 = p2+q-1
            break
    for q in range(p2):
        if scaff1[p1-q] != scaff2[p2-q]:
            l_split1 = p1-q+2
            l_split2 = p2-q+2
            break

    if l_split1 and r_split1:    
        out1 = scaff1[:l_split1]
        out1.append(linker)
        out2 = scaff1[l_split1:r_split1]
        out2.insert(0,linker)
        out2.append(linker)
        out3 = scaff1[r_split1:]
        out3.insert(0,linker)
        out_scaffs.extend([out1,out2,out3])
       
    elif l_split1:
        out1 = scaff1[:l_split1]
        out1.append(linker)
        out3 = scaff1[l_split1:]
        out3.insert(0,linker)
        out_scaffs.extend([out1,out3])    
    
    elif r_split1:
        out1 = scaff1[:r_split1]
        out1.append(linker)
        out3 = scaff1[r_split1:]
        out3.insert(0,linker)
        out_scaffs.extend([out1,out3])
        
    else:
        paf("NB! Expected conflict not found")

    if l_split2 and r_split2:
        out1 = scaff2[:l_split2]
        out1.append(linker)
        out2 = scaff2[l_split2:r_split2]
        out2.insert(0,linker)
        out2.append(linker)
        out3 = scaff2[r_split2:]
        out3.insert(0,linker)
        out_scaffs.extend([out1,out2,out3])
        
    elif l_split2:
        out2 = scaff2[:l_split2]
        out2.append(linker)
        out3 = scaff2[l_split2:]
        out3.insert(0,linker)
        out_scaffs.extend([out2,out3])
    
    elif r_split2:
        out2 = scaff2[:r_split2]
        out2.append(linker)
        out3 = scaff2[r_split2:]
        out3.insert(0,linker)
        out_scaffs.extend([out2,out3])
        
    else:
        paf("NB! Expected conflict not found")

    return out_scaffs


def purge_redundancy(scaff_list):
    """
    Removes scaffolds that are empty or entirely contained within another
    Called by unique_contigs_are_unique(), combine_overlapping_contigs(), MAIN
    Calls list_in_list()
    """
    for scaff in list(scaff_list):
        if len(scaff) < 4:
            scaff_list.remove(scaff)

    to_delete = ["deleted"] #place-marker for deleted scaffolds
        
    for n in range(0,(len(scaff_list)-1)):

        if scaff_list[n] != to_delete:    
            n_core = scaff_list[n][1:-1]
            for m in range((n+1),len(scaff_list)):
                if scaff_list[m] != to_delete:
                    m_core = scaff_list[m][1:-1]
                    if list_in_list(m_core, scaff_list[n]):
                        scaff_list[m] = to_delete
                    elif list_in_list(n_core, scaff_list[m]):
                        scaff_list[n] = to_delete
                    
                    if "dummy" in m_core[0]:
                        if list_in_list([m_core[1]], scaff_list[n]) or list_in_list([m_core[2]], scaff_list[n]):
                            scaff_list[m] = to_delete
                    elif "dummy" in n_core[0]:
                        if list_in_list([n_core[1]], scaff_list[m]) or list_in_list([n_core[2]], scaff_list[m]):
                            scaff_list[n] = to_delete
                
    while to_delete in scaff_list:
        scaff_list.remove(to_delete)
            
    return scaff_list
 
 
def list_in_list(a,b):
    """
    Determines whether list a is part of list b, either forwards or backwards
    """
    if any(a == b[offset:offset+len(a)] for offset in range(len(b)-len(a)+1)):
        return True
    else: 
        a.reverse()
        if any(a == b[offset:offset+len(a)] for offset in range(len(b)-len(a)+1)):
            return True
        else: return False
       
    
def unique_scaffold_ends(scaffold_list, r_cont, cont_dict, overlap):
    """
    Takes a list of scaffolds and returns the names of the ends of the first and last 
    unique contigs and their distance from the scaffold end (the flange)
    Calls get_unique_end()
    """
    scaff_end_dict = {}
    lsl = len(scaffold_list)
    for n in range(lsl):
        end_name = "left_" + str(n)
        scaff_end_dict[end_name] = get_unique_end(scaffold_list[n], r_cont, cont_dict, overlap)
        end_name = "right_" + str(n)
        scaff_end_dict[end_name] = get_unique_end(reversed(scaffold_list[n]), r_cont, cont_dict, overlap)
    n = lsl
    return scaff_end_dict

    
def get_unique_end(scaffold, r_cont, cont_dict, overlap):
    """
    Returns name of first unique contig end and its distance from the start (flange)
    Called by unique_scaffold_ends()
    Calls other_end()
    """
    flange_len = overlap
    temp_scaffold = list(scaffold)    
    contig_end = temp_scaffold.pop(1)
    contig = contig_end.split("prime_")
    if contig_end in r_cont:
        #add the length of the terminal contig if it is a repeat       
        #(originally just if attached by its 5' end, but rRNA was a problem)
        flange_len += (len(cont_dict[contig[1]].seq) - overlap)
        while len(temp_scaffold) > 1:
            if temp_scaffold[1] in r_cont:
                #if the next contig is a repeat, add its length to the flange and move on
                contig_name = temp_scaffold[1].split("prime_")[1]
                flange_len += (len(cont_dict[contig_name].seq) - overlap)
                contig_end = temp_scaffold[1]
                temp_scaffold.pop(1)
                temp_scaffold.pop(1)
            else:   
                contig_end = temp_scaffold[1]
                temp_scaffold = ["done"]
    else: contig_end = other_end(contig_end)
    return contig_end, flange_len
    
        
def best_pairing(current_end, end_dict, inverse_dict, blast_hits, l_min_score, r_min_score):
    """
    Returns a dict of possible connections to a scaffold end, based on adjacency in reference genomes.
    """
    #this duplicates part of trio_hits - should try to rewrite that to use this function
    
    l_flange = int(end_dict[current_end][1])
    l_contig = end_dict[current_end][0]
    
    #first find blast hits for the target scaffold end
    left_matches = []
    for hit in blast_hits:
        if hit[0] == l_contig and int(hit[11]) >= l_min_score:
            left_matches.append(hit)
                        
    link_count = {}
    
    #then find other ends with correctly oriented hits adjacent to the target hits
    for slink in end_dict:
        link = end_dict[slink][0]
  
        right_matches = []

        for hit in blast_hits:
            if hit[0] == link and int(hit[11]) >= r_min_score:                   
                right_matches.append(hit)
            
        for lhit in left_matches:
            for rhit in right_matches:
                srhit = inverse_dict[rhit[0]]
                r_flange = end_dict[srhit][1]
                joint_flange = l_flange + r_flange
                
                if lhit[1] == rhit[1]:
                    lh_start = int(lhit[8])
                    lh_end = int(lhit[9])
                    rh_start = int(rhit[8])
                    rh_end = int(rhit[9])

                    if abs(lh_start - rh_start) < joint_flange + 3000:
                        if (lh_end - lh_start)/(rh_end - rh_start) < 0:
                            if abs(lh_end - rh_end) > abs(lh_start - rh_start):
                                link_score = int(lhit[11]) * int(rhit[11])
                                if not link in link_count: 
                                    link_count[link] = link_score
                                elif link_score > link_count[link]:
                                    link_count[link] = link_score
    return link_count


def make_pairs_dict(scaffold_list, scaff_end_dict, hit_list):
    """
    Creates a dict that lists scaffold ends that could be adjacent to each scaffold end
     - the dict has symmetrical entries
    """
    #create a dict with unique_contig names as keys and scaffold end names as values
    inverse_dict = {v[0]: k for k, v in scaff_end_dict.items()}

    pairs_dict = {}
    for end1 in scaff_end_dict:
        link_dict = {}
        scaff_no = int(end1.split("_")[1])
        scaffold = scaffold_list[scaff_no]
        if scaffold[1] != other_end(scaffold[-2]):
            link_dict = best_pairing(end1, scaff_end_dict, inverse_dict, hit_list, scaff_l_min_score, scaff_r_min_score)

        end_dict = {}
        for link in link_dict:
            end = inverse_dict[link]
            if end != end1: #Do not allow ends to join to themselves
                end_dict[end] = link_dict[link]
        pairs_dict[end1] = end_dict
    
    return pairs_dict


def make_links_dict(pairs_dict):
    """     
    Creates links_dict by pruning pairs_dict to a single best link for each scaffold end      
    """
    links_dict = {}
    for end1 in pairs_dict:
    
        if (end1 in pairs_dict) and (len(pairs_dict[end1])) > 0:
            best_pair = max(pairs_dict[end1], key = pairs_dict[end1].get)
            
            if best_pair in pairs_dict and len(pairs_dict[best_pair]) > 0:
                
                if max(pairs_dict[best_pair], key = pairs_dict[best_pair].get) == end1:
                    links_dict[end1] = best_pair
                    links_dict[best_pair] = end1
    return links_dict


def join_scaffolds(first_end, new_scaffold, new_end, links_dict, scaffold_list, used_scaffs):
    """
    Extends a scaffold by adding other scaffolds that can be linked by using links_dict
    to find matching end sequences 
    Calls other_side(), other_end()
    """    
    while new_end not in used_scaffs:

        if new_end in links_dict and len(links_dict[new_end]) > 0:
            next_scaff_start = links_dict[new_end]
            if next_scaff_start != other_side(first_end):
                ns = next_scaff_start.split("_")
                next_scaff_number = int(ns[1])
                next_scaff_dir = ns[0]
                next_scaffold = scaffold_list[next_scaff_number]
                if next_scaff_dir == "right":
                    next_scaffold.reverse()
                if new_scaffold[-2] != other_end(next_scaffold[1]):
                    new_scaffold = new_scaffold[:-1] + [other_end(new_scaffold[-2])] + [other_end(next_scaffold[1])] + next_scaffold[1:]
                else:
                    new_scaffold = new_scaffold[:-1] + next_scaffold[1:]
                used_scaffs.append(new_end)
                used_scaffs.append(next_scaff_start)
                new_end = other_side(next_scaff_start)
            else: 
                new_scaffold[-1] = "join_circle"
                used_scaffs.append(new_end)
            
        else:
            new_scaffold[-1] = "join_not_found"
            used_scaffs.append(new_end)

    return new_scaffold    


def remove_duplicates(my_list):
    """
    Creates a new list with just the first occurrence of any repeated item
    """
    result = []
    for item in my_list:
        if item not in result:
            result.append(item)
    return result
          

def merge_scaffolds(scaff1, scaff2):
    """
    creates a single scaffold from two
    """
    if scaff1[-2] == other_end(scaff2[1]):
        new_scaff = scaff1[:-1] + scaff2[1:]
    elif scaff1[-3] == scaff2[1] and scaff1[-2] == scaff2[2]:
        new_scaff = scaff1[:-1] + scaff2[3:]
    else:
        new_scaff = scaff1[:-1] + [other_end(scaff1[-2]), other_end(scaff2[1])] + scaff2[1:]
    return new_scaff
    
            
def send_to_make_seq(scaffold_list, overlaps_dict):
    '''
    Sends a set of scaffolds to make_seq
    Calls make_seq()
    '''
    assembly_tuples = []
    for scaffold in scaffold_list:
        if len(scaffold) > 2:
            scaffold_tuple = make_seq(scaffold, overlaps_dict)
            #Discard empty sequences
            if len(scaffold_tuple[0]) != 0:
                assembly_tuples.append(scaffold_tuple)
    return assembly_tuples
            

def make_seq(scaffold, o_dict):
    """
    Uses a scaffold description to create the corresponding DNA sequence and a 'nice'
    compact description that is used in the scaffold_descriptions output file
    Called by send_to_make_seq()
    Calls initiate_seq(), extend_seq()
    """
    scaff_name = scaffold[0]
    sequence = []
    
    nice_scaff = "contigs__"
    
    scaff_string = str(scaffold)
    while scaffold:
    
        if len(scaffold) == 1:
            #This should never happen!
            paf("\nWARNING: odd number of elements in scaffold!")
            paf("scaffold is: " + scaff_string)
            nice_scaff += "WARNING:_odd_number_of_elements_in_scaffold!"
            sequence.description = scaff_name
            return sequence, nice_scaff

        end1 = scaffold.pop(0)
        end2 = scaffold.pop(0)
        
        if end1[0:4] != "five" and end1[0:5] != "three":
            if end2 in repeat_contigs and end2[0:10] == "threeprime":
                #Only attach a repeat if connected by fiveprime end,
                # to avoid creating duplicate copies
                ''' this condition has been removed!
                end1 = scaffold.pop(0)
                end2 = scaffold.pop(0)
                #threeprime ends of repeats are not attached
                if end2[0:4] != "five" and end2[0:5] != "three": end2 = other_end(end1)
                '''
                
            if "dummy" in end2:
                end1 = scaffold.pop(0)
                end2 = scaffold.pop(0)

            if end2[0:4] != "five" and end2[0:5] != "three":
                #This should never happen! 
                paf("\nWARNING: scaffold not included in assembly!")
                paf("scaffold is: " + scaff_string)
                paf("end1 is: " + str(end1))
                paf("end2 is: " + str(end2)+ "\n")
                nice_scaff += "scaffold.not.included.in.assembly!" + str(end1) + "." + str(end2)
                sequence.description = scaff_name
                return sequence, nice_scaff
            else:
                sequence, nice_scaff = initiate_seq(end2, nice_scaff)
        elif (end2 != "link_circular") and ("dummy" not in end1):
            sequence, nice_scaff = extend_seq(sequence, end0, end1, o_dict, nice_scaff)
        end0 = end2
        
    sequence.description = scaff_name
    
    return sequence, nice_scaff
    
    
def initiate_seq(end, nice_scaff):
    """
    Starts a new scaffold sequence and corresponding 'nice' description
    Called by make_seq()
    """
    end_name = end.split("prime_")
    
    contig_number = end.split("_")[contig_num_pos + 1]
    
    contig = contigs_dict[end_name[1]]
    if end_name[0] == "five":
        contig = contig.reverse_complement()
        nice_scaff += contig_number + "r"
    else:
    	nice_scaff += contig_number + "f"
    seq = contig
    return seq, nice_scaff
    
         
def extend_seq(seq, end0, end, o_dict, nice_scaff):
    """
    Extends a scaffold sequence and the corresponding 'nice' description by 1 contig
    Called by make_seq()
    """
    end_name = end.split("prime_")
    contig_number = end.split("_")[contig_num_pos + 1]
    contig = contigs_dict[end_name[1]]

    if end_name[0] == "three":
        contig = contig.reverse_complement()
        contig_number += "r"
    else:
    	contig_number += "f"
    	
    if end in o_dict[end0]:
        overlap = o_dict[end0][end]        
        seq += contig[overlap:]
        nice_scaff += ":" + contig_number
    else: 
    	seq += "NNNNNNNNNNNNNNNNNNNN" + contig
    	nice_scaff += "." + contig_number
    
    return seq, nice_scaff
         
########################################################
# MAIN                                                 #
########################################################

paf("Running " + script_file)
paf("version " + version)


data_headers = "strain \tmedian_coverage \tall_contigs \tdiscarded_contigs \tclean_contigs \tunique_cont_ends \trepeat_cont_ends \tlinkless_cont_ends \tunplaced_cont_ends \trep_count \trep_contigs \tchromosomal_contigs \tchromosomal_scaffolds \tplasmid_scaffolds \tfragment_scaffolds \ttotal_scaffolds \tchromosomal_length \tplasmid_length \tfragment_length \ttotal_length"
with open(data_table_file, "a") as outfile:
    outfile.write(data_headers + "\n")

for strain in strain_list:

    seq_file = source_folder + strain + file_extension
    discards_file = working_folder + "discarded_contigs_" + version + "/" + strain + "_" + version + "_discards.fas"
    clean_contigs_file = working_folder + "clean_contigs_" + version + "/" + strain + "_" + version + "_clean_contigs.fas"
    overlaps_txt = working_folder + "contig_overlaps_" + version + "/"  + strain +"_" + version + "_overlaps.txt"
   
    tag_file = blast_files_folder + strain + "_" + version + "_tags.fas"
    tag_blast_output = blast_files_folder + strain + "_" + version + "_tag_blast_ref8genomes.tab"
    #scaffolded_contigs = working_folder + "scaffolded_contigs_lists_" + version + "/" + strain + "_" + version + "_scaffolded_contigs_list.txt"
    scaffold_seq_file = working_folder + "scaffolds_" + version + "/" + strain + "_" + version + "_scaffolds.fas"
    scaffold_description_file = working_folder + "scaffold_descriptions_" + version + "/" + strain + "_" + version + "_scaffold_contigs.txt"
    contig_count_file = working_folder + "contig_counts_" + version + "/" + strain + "_" + version + "_contig_counts.txt"
    rep_contig_file = working_folder + "marker_contigs_" + version  + "/" + strain + "_" + version + "_rep_contigs.txt"
    chr_contig_file = working_folder + "marker_contigs_" + version  + "/" + strain + "_" + version + "_chr_contigs.txt"

    rep_blast_output = blast_files_folder + strain + "_RepA_blast.tab"
    chr_blast_output = blast_files_folder + strain + "_chr_blast.tab"

    paf("\n\n*** STARTING STRAIN " + strain + " at " + str(datetime.datetime.now()) + "\n")
    data_string = strain + "\t"
    
    contigs_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
    
    #Clean the data, discarding contaminant contigs with low coverage or homopolymers
    #################################################################################

    median_coverage = get_median_coverage(contigs_dict)
    paf("median coverage = " + str(median_coverage))
    data_string += str(median_coverage) + "\t"
    coverage_cutoff = median_coverage * coverage_threshold
    data_string += str(len(contigs_dict)) + "\t"   

    discards = {}
    for contig in list(contigs_dict):
        this_seq = contigs_dict[contig].seq
        #discard homopolymers
        if this_seq.count(this_seq[0]) == len(this_seq):
            discards[contig] = contigs_dict.pop(contig)
        #discard low-coverage contigs
        elif coverage_threshold:
            coverage = float(contig.split("_")[5])
            if coverage < coverage_cutoff:
                discards[contig] = contigs_dict.pop(contig)
    with open(discards_file, "w") as outfile:
        SeqIO.write(discards.values(), outfile, "fasta")
    with open(clean_contigs_file, "w") as outfile:
        SeqIO.write(contigs_dict.values(), outfile, "fasta")
    data_string += str(len(discards)) + "\t" + str(len(contigs_dict)) + "\t"
    
    #Make an overlaps dict (assembly graph) and classify contigs as unique or repeat
    ################################################################################
            
    contigs_list = SeqIO.parse(clean_contigs_file, "fasta")        
    contig_ends = end_list(contigs_list, 0, overlap)
    overlaps_dict = matching_ends(contig_ends, max_mismatches)

    with open(overlaps_txt, "w") as ofile:
        for key, value in overlaps_dict.items():
            ofile.write(key) 
            ofile.write(str(value))
            ofile.write("\n")
    
    #Save a version of overlaps_dict that includes head-to-tail self matches
    complete_overlaps_dict = copy.deepcopy(overlaps_dict)
        
    #remove 'self' hits from overlaps_dict
    for contig in overlaps_dict:
        if other_end(contig) in overlaps_dict[contig]:
            overlaps_dict[contig].pop(other_end(contig))

    unique_contigs = []
    repeat_contigs = []
    unlinked_contigs = []

    for key, value in overlaps_dict.items():
        if len(value) == 1:
            unique_contigs.append(key)
        if len(value) > 1:
            repeat_contigs.append(key)
        if len(value) == 0:
            unlinked_contigs.append(key) 
    all_contigs = unique_contigs + repeat_contigs + unlinked_contigs
        
    contig_tags = []
    for end_name in all_contigs:
        name_parts = end_name.split("prime_")
        end = name_parts[0]
        seq_id = name_parts[1]
        if seq_id in contigs_dict:
            contig = contigs_dict[seq_id]
            if end == "five":
                tag = contig[overlap:overlap + tag_length]
            if end == "three":
                tag = contig.reverse_complement()[overlap : (overlap + tag_length)]
            tag.id = end + "prime_" + contig.id
            contig_tags.append(tag)
    
    outfile = open(tag_file, "w")
    SeqIO.write(contig_tags, outfile, "fasta")
    outfile.close()

    #Run a blast search of the reference genomes with the tag_file as query
    subprocess.call("blastn -db " + reference_db + " -query " + tag_file +" -out " + 
    tag_blast_output +" -outfmt 6 -max_hsps 1 -evalue 1E-40 -num_threads 6", shell=True)

    #Make a table of hits from the blast output
    hit_list = []
    with open(tag_blast_output) as blast_matches:
        for hit in blast_matches:
            hit_data = hit.rstrip("\n").split("\t")
            hit_list.append(hit_data)
     
    #contigs with one link on one side but more than one on the other are repeats
    for contig_end in unique_contigs:
        if other_end(contig_end) in repeat_contigs:
            repeat_contigs.append(contig_end)
    for contig_end in repeat_contigs:
        if contig_end in unique_contigs:
            unique_contigs.remove(contig_end) 
    
            
    #contigs longer than 10kb are treated as unique even if they have >1 link,
    #which can arise if they are next to a polymorphic contig set
    long_repeats = []
    for contig_end in repeat_contigs:
        if int(contig_end.split("_")[4]) >= 10000:
            unique_contigs.append(contig_end)
            long_repeats.append(contig_end)
    for contig_end in unique_contigs:
        if contig_end in repeat_contigs:
            repeat_contigs.remove (contig_end)
    
    #Make sure both ends of unique contigs are in the list        
    for contig_end in unique_contigs:
        if other_end(contig_end) not in unique_contigs:
            unique_contigs.append(other_end(contig_end))
            
    #contigs that are unlinked but have blast hits need to be included in assembly
    linkless_contigs = []
    #contigs that are unlinked and have no hits in reference genomes will be unplaced
    unplaced_contigs = []
    for contig_end in unlinked_contigs:
        if other_end(contig_end) in unlinked_contigs:
            found_hit = False
    
            for hit in hit_list:
                if hit[0] == contig_end and float(hit[2]) >= 90:
                    if contig_end not in linkless_contigs:
                        linkless_contigs.append(contig_end)
                    if other_end(contig_end) not in linkless_contigs:
                        linkless_contigs.append(other_end(contig_end))
                    found_hit = True
                    break
            if found_hit == False:
                if contig_end not in unplaced_contigs:
                    unplaced_contigs.append(contig_end)
                if other_end(contig_end) not in unplaced_contigs:
                    unplaced_contigs.append(other_end(contig_end))
            
    if verbose:
        with open(working_folder + strain + "_unique_" + version + ".txt", "w") as ofile:
            for rec in unique_contigs:
                ofile.write(rec + "\n")  
        with open(working_folder + strain  + "_repeat_" + version + ".txt", "w") as ofile:
            for rec in repeat_contigs:
                ofile.write(rec + "\n")  
        with open(working_folder + strain  + "_unlinked_" + version + ".txt", "w") as ofile:
            for rec in unlinked_contigs:
                ofile.write(rec + "\n")  
        with open(working_folder + strain  + "_linkless_" + version + ".txt", "w") as ofile:
            for rec in linkless_contigs:
                ofile.write(rec + "\n")  
        with open(working_folder + strain  + "_unplaced_" + version + ".txt", "w") as ofile:
            for rec in unplaced_contigs:
                ofile.write(rec + "\n")  
    
    data_string += str(len(unique_contigs)) + "\t" +   str(len(repeat_contigs)) + "\t" + str(len(linkless_contigs)) + "\t" + str(len(unplaced_contigs)) + "\t"   
    unique_contigs_master = copy.deepcopy(unique_contigs)  
    

    #Identify contigs with plasmid rep sequences or core chromosomal genes
    ######################################################################     
      
    paf("finding plasmid and chromosomal genes at " + str(datetime.datetime.now()))
    subprocess.call("makeblastdb -in " + clean_contigs_file + " -input_type fasta -dbtype nucl -logfile logfile.txt", shell=True)
    
    #Run a blast search of the contigs with a file of RepA standards as query
    subprocess.call("tblastn -db " + clean_contigs_file + " -query " + rep_file +" -out " + 
    rep_blast_output +" -outfmt 6 -qcov_hsp_perc 90.0  -num_threads 6", shell=True)
    
    #Make a dict of hits from the RepA blast output
    #The values are the names of the rep type(s) in that contig
    rep_contig = {}
    rep_count = 0
    with open(rep_blast_output) as blast_matches:
        for hit in blast_matches:
            hit_data = hit.rstrip("\n").split("\t")
            if float(hit_data[2]) >= 95.0:
                rep_count += 1
                rep_type = hit_data[0].split("_")[1]
                contig = hit_data[1]
                if rep_type in str(rep_contig):
                    paf("DUPLICATED RepA HITS! " + rep_type) #Should not happen!
                if contig in rep_contig:
                    rep_contig[contig] += "+" + rep_type
                else: 
                    rep_contig[contig] = rep_type
     
    with open(rep_contig_file, "w") as ofile:
        for contig in rep_contig:
            ofile.write(strain + "\t" + rep_contig[contig] + " \t" + contig)
            ofile.write("\n")
               
    data_string += str(rep_count) + "\t" + str(len(rep_contig)) + "\t"
    
    #Run a blast search of the contigs with a file of chromosomal core genes as query
    subprocess.call("blastn -db " + clean_contigs_file + " -query " + chr_file +" -out " + 
    chr_blast_output +" -outfmt 6 -qcov_hsp_perc 90.0  -num_threads 6", shell=True)
    
    #Make a dict of hits from the chr blast output
    #The values are the number of chromosomal genes identified in that contig
    #NB: this earlier code is partly duplicated by order_chromosomal_contigs().
    chr_contig = {}
    with open(chr_blast_output) as blast_matches:
        for hit in blast_matches:
            hit_data = hit.rstrip("\n").split("\t")
            #Only use hits with >=90% identity 
            if float(hit_data[2]) >= 90.0: 
                contig = hit_data[1]
                if contig in chr_contig:
                    chr_contig[contig] += 1
                else: 
                    chr_contig[contig] = 1
     
    with open(chr_contig_file, "w") as ofile:
        for contig in chr_contig:
            ofile.write(strain + "\t" + contig + "\t" + str(chr_contig[contig]))
            ofile.write("\n")
    data_string += str(len(chr_contig)) + "\t"
    
    #Assemble chromosomal scaffolds using conserved contig order and overlaps
    #########################################################################    
            
    paf("starting scaffolding at " + str(datetime.datetime.now()))
                                               
    #Create an ordered list of chromosomal contigs based on hits in Blast output         
    ordered_chr_contigs = order_chromosomal_contigs(chr_blast_output)
    
    for item in ordered_chr_contigs:
        if item[0] in str(rep_contig):
            paf("WARNING: chromosomal contig also has rep: removing from ordered_chr_contigs")
            paf(str(item)) 
            ordered_chr_contigs.remove(item)
            
    for item in ordered_chr_contigs:
        if item[0] in str(repeat_contigs): 
            paf("WARNING: chromosomal contig is a repeat: ")
            paf(str(item))
     
    chr_links = {}
    for c1 in ordered_chr_contigs:
        c0 = ordered_chr_contigs[ordered_chr_contigs.index(c1) - 1]
        link_l = ("threeprime_" if c0[1] == 1 else "fiveprime_") + c0[0]
        link_r = ("fiveprime_" if c1[1] == 1 else "threeprime_") + c1[0]
        chr_links[link_l] = link_r
        chr_links[link_r] = link_l
    
    #Create versions of unique_contigs and overlaps_dict that have the plasmid rep contigs
    #removed so that plasmid contigs cannot be 'cointegrated' into the chromosome
    
    x_rep_o_dict = copy.deepcopy(overlaps_dict)
    x_rep_u_contigs = copy.deepcopy(unique_contigs)
    
    for rep_c in rep_contig.keys():
        rep_c_end = "fiveprime_" + rep_c
        censor_contig(rep_c_end, x_rep_u_contigs, x_rep_o_dict)
                    
    #Go through the contigs in order, extending each rightwards into a scaffold
    remaining_chr_contigs = copy.deepcopy(ordered_chr_contigs)
    chr_scaffold_list = []
    n = 1    
    current_contig = get_next_contig(remaining_chr_contigs)
    header = "chromosome-" + '{:02d}'.format(n)
    scaffold = [header, current_contig]
    
    while remaining_chr_contigs:

        next_contig = get_next_contig(remaining_chr_contigs)
        
        #allow scaffold to grow, but only as far as the next chr contig                
        new_scaffold = grow_scaffold(scaffold, x_rep_o_dict, other_end(next_contig), repeat_contigs)                
        if new_scaffold[-1][0:4] == "link": #scaffold stopped before reaching next chr contig        
            if len(new_scaffold) == 2:  #no extension beyond the initial contig
                new_scaffold = [header, "fiveprime_dummy", other_end(current_contig), current_contig, "threeprime_dummy", "linkless"]
                current_contig = next_contig
            else:
                #skip over the other chr contigs that were already added to the scaffold
                while other_end(next_contig) in new_scaffold:
                    if remaining_chr_contigs:                 
                        next_contig = get_next_contig(remaining_chr_contigs)
                    else: next_contig = ""
                current_contig = next_contig
                                
            chr_scaffold_list.append(new_scaffold)          
            n += 1
            header = "chromosome-" + '{:02d}'.format(n)
            scaffold = [header, current_contig]
            
            if current_contig and not remaining_chr_contigs: #add the last chr_contig if not in previous scaffold 
                new_scaffold = [header, "fiveprime_dummy", other_end(current_contig), current_contig, "threeprime_dummy", "linkless"]
                chr_scaffold_list.append(new_scaffold)          
                

        else: scaffold = new_scaffold   #try to extend the scaffold further

    #Extend each scaffold leftwards
    new_chr_scaffold_list = []        
    for scaffold in chr_scaffold_list:
        header = scaffold[0]
        
        if scaffold[1] == "fiveprime_dummy":
            scaffold.pop(1)
            scaffold.pop(1)
    
        scaffold.reverse()
        scaffold[-1] = other_end(scaffold[-2])
        
        scaffold = grow_scaffold(scaffold, x_rep_o_dict, x_rep_u_contigs, repeat_contigs)
        scaffold.reverse()
        scaffold[0] = header
        
        if scaffold[-2] == "threeprime_dummy":
            if len(scaffold) == 4:
                scaffold.insert(1, other_end(scaffold[1]))
                scaffold.insert(1, "fiveprime_dummy")
            else:
                scaffold.pop(-2)
                scaffold.pop(-2)
        new_chr_scaffold_list.append(scaffold)

    chr_scaffold_list = new_chr_scaffold_list 
    
    for scaff in chr_scaffold_list:
        scaff[-1] = "chr_end" 

    #Remove scaffolds that duplicate part of their neighbour
    chr_scaffold_list = purge_redundancy(chr_scaffold_list)
    
    #Split scaffolds that have conflicting placement of unique contigs
    chr_scaffold_list = unique_contigs_are_unique(chr_scaffold_list, unique_contigs)

    #Get the scaffolds back in order    
    new_chr_scaff_list = []
                    
    for chr_cont in ordered_chr_contigs:
        if chr_cont[0] not in str(chr_scaffold_list):
            paf("WARNING: chr contig missing: " + chr_cont[0])
        else:
            for scaff in chr_scaffold_list:
                if chr_cont[0] in str(scaff):
                    fp = "fiveprime_" + chr_cont[0]
                    tp = "threeprime_" + chr_cont[0]
                    if (scaff[1] == fp or scaff[-2] == tp 
                    or (fp in scaff and tp in scaff and scaff.index(fp) > scaff.index(tp))):
                        cont_dir = -1
                    else: cont_dir = 1
                    if cont_dir * chr_cont[1] == -1:
                        scaff.reverse()
                    if scaff not in new_chr_scaff_list:
                        new_chr_scaff_list.append(scaff)             
    
    #Get scaffolds in correct orientation and then renumber
    chr_scaffold_list = new_chr_scaff_list
    for n in range(len(chr_scaffold_list)):
        if "chromosome" in chr_scaffold_list[n][-1] or "chr_end" in chr_scaffold_list[n][0]:
            chr_scaffold_list[n].reverse()
        header = "chromosome-" + '{:02d}'.format(n+1)
        chr_scaffold_list[n][0] = header 
        
    #Remove unique chromosomal contigs from unique list and overlaps dict
    edited_overlaps_dict = copy.deepcopy(overlaps_dict)
    for scaffold in chr_scaffold_list:
        for contig_end in scaffold:
            if contig_end in unique_contigs:
                censor_contig(contig_end, unique_contigs, edited_overlaps_dict)
                                       
    #Assemble plasmid scaffolds starting with each rep contig   
    #########################################################
    
    plasmid_scaffold_list = []

    for rep_c in rep_contig:
        temp_u_contigs = copy.deepcopy(unique_contigs)
        temp_o_dict = copy.deepcopy(edited_overlaps_dict)
     
        #Censor the other rep contigs so that cointegrates cannot form
        for r_c in rep_contig:
            if r_c != rep_c:
                r_c_end = "fiveprime_" + r_c
                censor_contig(r_c_end, temp_u_contigs, temp_o_dict)
    
        rep_type = rep_contig[rep_c]
        new_start = "threeprime_" + rep_c
       
        seed = ["start", new_start]
        scaffold = grow_scaffold(seed, temp_o_dict, temp_u_contigs, repeat_contigs)
                
        if scaffold[0] == "start":
            scaffold.reverse()
            scaffold[-1] = other_end(scaffold[-2])
        
            if scaffold[1] == scaffold[-1]:
                scaffold[-1] = "link_circular"
            else:
                scaffold = grow_scaffold(scaffold, temp_o_dict, temp_u_contigs, repeat_contigs)
        
        scaffold[0] = "plasmid-" + rep_type
        if len(scaffold) == 2:
            scaffold = [scaffold[0], "fiveprime_dummy", new_start, other_end(new_start), "threeprime_dummy", "linkless"] 
        plasmid_scaffold_list.append(scaffold)

    ucau_plasmid_scaffold_list = unique_contigs_are_unique(plasmid_scaffold_list, unique_contigs_master)
    
    #find scaffolds that still contain plasmid rep genes and move the rest to frags list
    new_plasmid_scaffold_list = []
    frag_scaffs_from_plasmids = []
    for scaff in ucau_plasmid_scaffold_list:
        for contig in rep_contig:
            if contig in str(scaff):
                new_plasmid_scaffold_list.append(scaff)
        if scaff not in new_plasmid_scaffold_list:
            frag_scaffs_from_plasmids.append(scaff)
    plasmid_scaffold_list = new_plasmid_scaffold_list            
                
    #Remove unique plasmid contigs from unique list only
    for scaffold in plasmid_scaffold_list:
        for contig_end in scaffold:
            if contig_end in unique_contigs:
                unique_contigs.remove(contig_end)

    #Take remaining unique contigs and extend them into 'fragment' scaffolds
    ########################################################################
    
    frag_scaffold_list = []
    while len(unique_contigs) > 0:
        new_start = unique_contigs.pop()          
        seed = ["start", new_start]
        header = "*** Seed contig: " + new_start + " ***"
        scaffold = grow_scaffold(seed, edited_overlaps_dict, unique_contigs, repeat_contigs)
                
        #Extend each scaffold leftwards                
        if scaffold[1] == "fiveprime_dummy":
            scaffold.pop(1)
            scaffold.pop(1)
    
        scaffold.reverse()
        scaffold[-1] = other_end(scaffold[-2])
        
        if scaffold[1] == scaffold[-1]:
            scaffold[-1] = "link_circular"
        else:
            scaffold = grow_scaffold(scaffold, edited_overlaps_dict, unique_contigs, repeat_contigs)
        scaffold.reverse()
        scaffold[0] = header
        
        if scaffold[-2] == "threeprime_dummy":
            if len(scaffold) == 4:
                scaffold.insert(1, other_end(scaffold[1]))
                scaffold.insert(1, "fiveprime_dummy")
            else:
                scaffold.pop(-2)
                scaffold.pop(-2)
        
        if len(scaffold) == 2:
            scaffold = [scaffold[0], "fiveprime_dummy", new_start, other_end(new_start), "threeprime_dummy", "linkless"] 

        frag_scaffold_list.append(scaffold)
        for contig_end in scaffold:
            if contig_end in unique_contigs:
                unique_contigs.remove(contig_end)
                
    #add fragments trimmed from plasmid scaffolds by ucau
    frag_scaffold_list.extend(frag_scaffs_from_plasmids)             

    #add the remaining linkless contigs to frag_scaffold_list
    scaffold_list = chr_scaffold_list + plasmid_scaffold_list + frag_scaffold_list            
    for contig_end in linkless_contigs:
        if contig_end not in (str(scaffold_list)) and contig_end[0:4] == "five":
            ll_scaffold = ["linkless", "fiveprime_dummy", contig_end, other_end(contig_end), "threeprime_dummy", "linkless"]
            frag_scaffold_list.append(ll_scaffold)
            
    #use ucau on combined plasmid and fragment scaffolds, then sort out plasmids again                                
    fp_scaffold_list = unique_contigs_are_unique(plasmid_scaffold_list + frag_scaffold_list, unique_contigs_master)
    
    new_plasmid_scaffold_list = []
    new_frag_scaffold_list = []
    for scaff in fp_scaffold_list:
        for contig in rep_contig:
            if contig in str(scaff):
                new_plasmid_scaffold_list.append(scaff)
        if scaff not in new_plasmid_scaffold_list:
            new_frag_scaffold_list.append(scaff)
    plasmid_scaffold_list = new_plasmid_scaffold_list
    frag_scaffold_list = new_frag_scaffold_list            
     
    scaffold_list = chr_scaffold_list + plasmid_scaffold_list + frag_scaffold_list            
       
    if verbose:
        with open(working_folder + strain + "_stage1_scaffolds_" + version + ".txt", "w") as ofile:

            for n in range(len(scaffold_list)):
                ofile.write("stage1_" + str(n) + str(scaffold_list[n]))
                ofile.write("\n\n")
    
    #Join the prelim chromosomal scaffolds based on best adjacent matches in reference genomes
    ##########################################################################################

    paf("joining preliminary scaffolds at " + str(datetime.datetime.now()))
        
    #create a dict with scaffold left and right ends as keys and (unique_contig, flange_length) tuples as values
    scaff_end_dict = unique_scaffold_ends(chr_scaffold_list, repeat_contigs, contigs_dict, overlap)

    #list connections that are supported by homology in reference genomes
    pairs_dict = make_pairs_dict(chr_scaffold_list, scaff_end_dict, hit_list)
    
    #merge adjacent chromosomal scaffolds if they are linked in pairs_dict
    # (even if they have higher-scoring links elsewhere)
    lcsl = len(chr_scaffold_list)
    to_delete = ["deleted"] #place-marker for deleted scaffolds
    for n in range(lcsl-1):
        m = n+1
        end_a = "right_" + str(n)
        end_b = "left_" + str(m)
        if end_a in pairs_dict and end_b in pairs_dict[end_a]:
            new_scaffold = merge_scaffolds(chr_scaffold_list[n], chr_scaffold_list[m])
            chr_scaffold_list[n] = to_delete
            chr_scaffold_list[m] = new_scaffold
            
    while to_delete in chr_scaffold_list:
        chr_scaffold_list.remove(to_delete)

    #Check for join between last scaffold and first
    end_a = "right_" + str(lcsl-1)
    if end_a in pairs_dict and "left_0" in pairs_dict[end_a]:
        if len(chr_scaffold_list) == 1: 
            chr_scaffold_list[0][0] = "chromosome-circular"
        else: 
            new_scaffold = merge_scaffolds(chr_scaffold_list[-1], chr_scaffold_list[0])
            chr_scaffold_list[0] = new_scaffold
            chr_scaffold_list.pop()
    
    #Find fragments that patch gaps in the chromosome
    
    #create a links_dict for chr + plasmid + frag     
    cpf_scaff_list = chr_scaffold_list + plasmid_scaffold_list + frag_scaffold_list
    scaff_end_dict = unique_scaffold_ends(cpf_scaff_list, repeat_contigs, contigs_dict, overlap)
    pairs_dict = make_pairs_dict(cpf_scaff_list, scaff_end_dict, hit_list)

    links_dict = make_links_dict(pairs_dict)
    
    #Check whether any fragment scaffold is the best link for the left or right of 
    #each chromosomal scaffold.  If so, merge them.
    lcsl = len(chr_scaffold_list)
    lpsl = len(plasmid_scaffold_list)
    for n in range(lcsl):
        for c_side in ["left", "right"]:        
            c_end = c_side + "_" + str(n)
            if c_end in links_dict: 
                top_link = links_dict[c_end]
                f_side, f_scaff_no = top_link.split("_")
                f_scaff_no = int(f_scaff_no)
                if f_scaff_no in range(lcsl+lpsl+1, len(cpf_scaff_list)):
                    c_scaff = cpf_scaff_list[n]
                    f_scaff = cpf_scaff_list[f_scaff_no]
                    if "deleted" not in f_scaff:
                        if c_side == "left":
                            if f_side == "left": f_scaff.reverse()
                            new_scaffold = merge_scaffolds(f_scaff, c_scaff)
                        if c_side == "right":
                            if f_side == "right": f_scaff.reverse()
                            new_scaffold = merge_scaffolds(c_scaff, f_scaff)
                        chr_scaffold_list[n] = new_scaffold
                        cpf_scaff_list[f_scaff_no] = to_delete

    while to_delete in cpf_scaff_list:
        cpf_scaff_list.remove(to_delete)
    frag_scaffold_list = cpf_scaff_list[lcsl+lpsl+1:]
    
    #Join adjacent chr scaffolds that share end contig or pair of end contigs    
    lcsl = len(chr_scaffold_list)
    if lcsl > 1:
        to_delete = ["deleted"] #place-marker for deleted scaffolds
        for n in range(-1,lcsl-1):
            m = n+1
            if ("dummy" not in str(chr_scaffold_list[n] + chr_scaffold_list[m])
            and "deleted" not in str(chr_scaffold_list[n] + chr_scaffold_list[m])):
                if (chr_scaffold_list[n][-2] == other_end(chr_scaffold_list[m][1]) 
                or (chr_scaffold_list[n][-3] == chr_scaffold_list[m][1] 
                and chr_scaffold_list[n][-2] == chr_scaffold_list[m][2])):
                    new_scaffold = merge_scaffolds(chr_scaffold_list[n], chr_scaffold_list[m])
                    chr_scaffold_list[n] = to_delete
                    chr_scaffold_list[m] = new_scaffold
            
        while to_delete in chr_scaffold_list:
            chr_scaffold_list.remove(to_delete)

    #Renumber the scaffolds
    for n in range(len(chr_scaffold_list)):
        header = "chromosome-" + '{:02d}'.format(n+1)
        chr_scaffold_list[n][0] = header 

    #Join preliminary plasmid and fragment scaffolds
    ################################################
            
    other_scaffold_list =  plasmid_scaffold_list + frag_scaffold_list 

    #create a dict with scaffold left and right ends as keys and (unique_contig, flange_length) tuples as values
    scaff_end_dict = unique_scaffold_ends(other_scaffold_list, repeat_contigs, contigs_dict, overlap)

    #list connections that are supported by homology in reference genomes
    pairs_dict = make_pairs_dict(other_scaffold_list, scaff_end_dict, hit_list)        
    links_dict = make_links_dict(pairs_dict)            

    #use links_dict to add scaffolds to the right, then flip and add to the other end
    new_scaffold_list = []    
    used_scaffs = []
    for n in range(len(other_scaffold_list)):
        first_end = "right_" + str(n)
    
        if first_end not in used_scaffs:
            new_scaffold = other_scaffold_list[n]
            title = new_scaffold[0]
            new_end = first_end
        
            if new_end in links_dict:
                new_scaffold = join_scaffolds(first_end, new_scaffold, new_end, links_dict, other_scaffold_list, used_scaffs)
        
                new_end = other_side(first_end)                
                new_beginning = other_side(used_scaffs[-1])
                if new_end in links_dict:
                    new_scaffold.reverse()
                    new_scaffold[-1] = other_end(new_scaffold[-2])               
                    new_scaffold = join_scaffolds(new_beginning, new_scaffold, new_end, links_dict, other_scaffold_list, used_scaffs)
                    new_scaffold.reverse()
                    new_scaffold[0] = title

            new_scaffold_list.append(new_scaffold)

    #find unique contigs used more than once and merge the scaffolds containing them
    new_scaffold_list = unique_contigs_are_unique(new_scaffold_list, unique_contigs_master)

    #find any contigs that somehow got left out
    contigs_list = SeqIO.parse(clean_contigs_file, "fasta")        
    for contig in contigs_list:
    
        if contig.id not in str(new_scaffold_list + chr_scaffold_list):
            lost_scaffold = ["fragment", "fiveprime_dummy", "fiveprime_"+contig.id, "threeprime_"+contig.id, "threeprime_dummy", "lost"]
            new_scaffold_list.append(lost_scaffold)
    
    if verbose:      
        with open(working_folder + strain + "_stage2_scaffolds_" + version + ".txt", "w") as ofile:
            for n in range(len(new_scaffold_list)):
                ofile.write("stage2_" + str(n) + str(new_scaffold_list[n]))
                ofile.write("\n\n")
         
    #Find the rep contigs again in the new scaffolds and name scaffolds as plasmids    
    for scaff in new_scaffold_list:
        scaff[0] = "fragment"
    rep_contig_location = {}
    for contig in rep_contig:
        contig_end = "fiveprime_" + contig
        finds = find_unique_contig(contig_end, new_scaffold_list)
        if not finds:
            contig_end = "threeprime_" + contig
            finds = find_unique_contig(contig_end, new_scaffold_list)
    
        if finds:
            rep_type = rep_contig[contig]
            rep_contig_location[rep_type] = finds[0]
        else: paf(str(contig) + " NOT FOUND!")

    for rep_type in rep_contig_location:
        rep_scaff = rep_contig_location[rep_type]
        if rep_scaff[0][:7] == "plasmid":
            rep_scaff[0] += "+" + rep_type
        else: 
            rep_scaff[0] = "plasmid-" + rep_type
            
    #Create a new plasmid scaffold list and remove plasmids from the fragment list
    plasmid_scaffold_list = []
    for scaff in new_scaffold_list:
        if "plasmid" in scaff[0]:
            plasmid_scaffold_list.append(scaff)              
    for rep_scaff in plasmid_scaffold_list:
        if rep_scaff in new_scaffold_list: new_scaffold_list.remove(rep_scaff)     
        if rep_scaff[1] == other_end(rep_scaff[-2]):
            if rep_scaff[1] != 'fiveprime_dummy':
                rep_scaff[-1] = "circular"
                rep_scaff[0] += "-circular"        

    plasmid_scaffold_list.sort()
    frag_scaffold_list = new_scaffold_list
    
    #remove any fragment scaffolds that include unique contigs already used in chr or plasmids
    temp_fsl = copy.deepcopy(frag_scaffold_list)
    
    for scaff in frag_scaffold_list:
        for contig_end in scaff:
            if contig_end in unique_contigs_master:
                if contig_end in str(chr_scaffold_list + plasmid_scaffold_list):
                    if scaff in temp_fsl: 
                        temp_fsl.remove(scaff)
                    
    frag_scaffold_list = temp_fsl
    
    #bring back any contigs that were lost entirely from the assembly
    contigs_list = SeqIO.parse(clean_contigs_file, "fasta")        
    for contig in contigs_list:    
        if contig.id not in str(chr_scaffold_list + plasmid_scaffold_list + frag_scaffold_list):
            lost_scaffold = ["fragment", "fiveprime_dummy", "fiveprime_"+contig.id, "threeprime_"+contig.id, "threeprime_dummy", "lost"]
            frag_scaffold_list.append(lost_scaffold)          
    
    if verbose:
        genome_scaffold_list = chr_scaffold_list + plasmid_scaffold_list + frag_scaffold_list
        with open(working_folder + strain + "_final_scaffolds_" + version + ".txt", "w") as ofile:

            for n in range(len(genome_scaffold_list)):
                ofile.write(str(n) + str(genome_scaffold_list[n]))
                ofile.write("\n\n")
                
    

    #Use the scaffold lists to assemble DNA sequences and create various output files
    #################################################################################    
    
    paf("starting assembly at " + str(datetime.datetime.now()))
            
    chromosome_tuples = send_to_make_seq(chr_scaffold_list, overlaps_dict)
    plasmid_tuples = send_to_make_seq(plasmid_scaffold_list, overlaps_dict)
    fragment_tuples = send_to_make_seq(frag_scaffold_list, overlaps_dict)
    fragment_tuples = sorted(fragment_tuples, key=lambda x: len(x[0]), reverse = True)
    
    chr01_file = working_folder + "chromosome-01.fas"
    with open(chr01_file , "w") as outfile:
        SeqIO.write(chromosome_tuples[0][0], outfile, "fasta")   
        
    #split the first scaffold so that chr-00 ends just before the ATG start of dnaA and
    #chr-01 starts 127 bases upstream of the ATG (i.e. 127 base overlap)
        
    #run blast-2-sequences to find location of dnaA in chromosome-01.
    subprocess.call("blastn -subject " + chr01_file + " -query " + dnaA_file +" -out " + 
    dnaA_blast_out +" -outfmt 6 -max_hsps 1 -evalue 1E-40 ", shell=True)
    with open(dnaA_blast_out) as blast_matches:
        hit = list(blast_matches)[0]
        hit_data = hit.rstrip("\n").split("\t")
        dnaA_start_pos = int(hit_data[8]) - int(hit_data[6])

    #create the new scaffold sequences        
    chr01 = chromosome_tuples[0][0]    
    new_chr00 = chr01[:dnaA_start_pos]
    new_chr01 = chr01[dnaA_start_pos -127:]
    new_chr00.description = "chromosome-00"
    new_chr01.description = "chromosome-01"

    #create the new nice scaffold descriptions for the scaffold_contig file    
    first_contig = ordered_chr_contigs[0]
    first_contig_no = first_contig[0].split("_")[contig_num_pos]
    fc_dir = first_contig[1]
    if fc_dir == 1: d = "f"
    elif fc_dir == -1: d = "r"
    fc = str(first_contig_no) + d
    chr01_nice = chromosome_tuples[0][1]
    #need to ensure we find, e.g. 1f and not 11f or 21f, so include the preceding spacer
    if ":" + fc in chr01_nice:
        cut_pos = chr01_nice.find(":" + fc)
    elif "." + fc in chr01_nice:
        cut_pos = chr01_nice.find("." + fc)
    elif "_" + fc in chr01_nice:
        cut_pos = chr01_nice.find("_" + fc)
    else:
        paf("WARNING! first contig not found: " + fc)
        paf(chr01_nice)
        cut_pos = 9       
    new_chr00_nice = "contigs__" + chr01_nice[9:cut_pos + 1] + "="
    new_chr01_nice = "contigs__=" + chr01_nice[cut_pos + 1:]

    #put the new scaffold tuples into the list    
    new_chr00_tuple = (new_chr00, new_chr00_nice)
    new_chr01_tuple = (new_chr01, new_chr01_nice)
    chromosome_tuples[0] = new_chr01_tuple
    chromosome_tuples.append(new_chr00_tuple)
    
    genome_tuples = chromosome_tuples + plasmid_tuples + fragment_tuples

    assembly = []
    scaff_description_list = []
    n = 0
    for scaffold in genome_tuples:
        n += 1
        scaffold[0].id = strain + "_scaf_" + str(n) 
        scaffold[0].description += "_len_" + str(len(scaffold[0]))
        scaffold_desc = scaffold[0].id + "_" + scaffold[0].description + "  " + scaffold[1]
        assembly.append(scaffold[0])
        scaff_description_list.append(scaffold_desc)

    with open(scaffold_seq_file, "w") as outfile:
        SeqIO.write(assembly, outfile, "fasta")

    with open(scaffold_description_file, "w") as outfile:
        for scaffold in scaff_description_list:
            outfile.write(scaffold + "\n")

    contig_totals = {}
    for scaff in scaff_description_list:        
        cont_contents = scaff.split('__')[1]
        contigs = re.split('[=:\.fr]+', cont_contents)
        for contig in contigs:
            if contig:
                if contig in contig_totals:
                    contig_totals[contig] +=1
                else:
                    contig_totals[contig] = 1
                        
    for contig_end in unique_contigs_master:
        parts = contig_end.split("_")
        if parts[0] == "fiveprime":
            if parts[2] not in contig_totals:
                paf("WARNING: unique contig missing: ")
                paf(contig_end)
            elif contig_totals[parts[2]] > 1 and contig_end not in long_repeats:
                paf("WARNING: unique contig included more than once:")
                paf(contig_end + " found " +  str(contig_totals[parts[2]]) + " times")
                    
    scaff_data = {'chromosome':[0,0], 'plasmid':[0,0], 'fragment':[0,0], 'total':[0,0]}
    for scaff in scaff_description_list:
        
        scaff_title = scaff.split()[0]
        title_parts = scaff_title.split('_scaf_')
        scaff_info = title_parts[1].split('_')
        scaff_type = scaff_info[1].split('-')[0]
        scaff_type = scaff_type.split('+')[0]
        scaff_length = scaff_info[3]
        scaff_data[scaff_type][0] += 1
        scaff_data[scaff_type][1] += int(scaff_length)
        scaff_data['total'][0] += 1
        scaff_data['total'][1] += int(scaff_length)
    for n in range(0,2):    
        data_string += str(scaff_data['chromosome'][n]) + "\t" + str(scaff_data['plasmid'][n]) + "\t" + str(scaff_data['fragment'][n]) + "\t" + str(scaff_data['total'][n]) + "\t"
                
    all_contigs_list = SeqIO.parse(seq_file, "fasta")
    contig_count_list = []
    for contig in all_contigs_list:
    
            if contig.id in str(unique_contigs_master): contig_class = "unique  "
            elif contig.id in str(repeat_contigs): contig_class = "repeat  "
            elif contig.id in str(unlinked_contigs): contig_class = "unlinked"
            elif contig.id in str(linkless_contigs): contig_class = "linkless"
            elif contig.id in str(unplaced_contigs): contig_class = "unplaced"
            else: contig_class = "unused  "
    
            contig_number = contig.id.split("_")[1]    
            if contig.id in discards:
                contig_count_list.append((contig_number, contig.id, contig_class, "discarded"))
            elif contig_number in contig_totals:
                contig_count_list.append((contig_number, contig.id, contig_class, str(contig_totals[contig_number])))
            else:
                contig_count_list.append((contig_number, contig.id, contig_class, "missing"))
                        
    contig_count_list.sort(key = lambda x: int(x[0]))
        
    with open(contig_count_file, "w") as outfile:
        for contig in contig_count_list:
            outfile.write(contig[0] + "\t" + contig[1] + "\t" + contig[2] + "\t" + contig[3] + "\n")
            
    with open(data_table_file, "a") as outfile:
        outfile.write(data_string + "\n")
       
    paf("finished strain " + strain + " at " + str(datetime.datetime.now()))
     

        

