import argparse
import glob
import gzip
from collections import namedtuple
import bottleneck as bn
import numpy as np
import pandas as pd
from Bio import AlignIO
from scipy import linalg

Assignment = namedtuple("Assignment", ["a", "score", "a_excl_none", "score_excl_none"])

def parse_pop_map(file_name = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, sara_id, origin, country in zip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Origin2']):
        pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country}
    return pop_map

def normalize(matrix, axis=0):
    """Normalizes a matrix' columns or rows to mean 0 and std 1."""
    mean = np.mean(matrix, axis=axis)
    std = np.std(matrix, axis=axis)
    if axis == 1:
        mean = mean[:, None]
        std = std[:, None]

    matrix = (matrix - mean) / std
    return np.nan_to_num(matrix)

def replace_column_nans_by_mean(matrix):
    # Set the value of gaps/dashes in each column to be the average of the other values in the column.
    nan_indices = np.where(np.isnan(matrix))
    # Note: bn.nanmean() instead of np.nanmean() because it is a lot(!) faster.
    column_nanmeans = np.nanmean(matrix, axis=0)

    # For each column, assign the NaNs in that column the column's mean.
    # See http://stackoverflow.com/a/18689440 for an explanation of the following line.
    matrix[nan_indices] = np.take(column_nanmeans, nan_indices[1])


def account_for_structure(in_glob, out_dir, maf=0.01):
    #with open("/home/agb/data/rhizo/full_dataset/orthologs/strain_list.txt") as strain_list_in:
    #    strain_list = np.array([l.strip() for l in strain_list_in])

    # used a core gene as the base of my strains:
    data = np.load('/Users/PM/Desktop/New_data/group_snps/group1254.fna.snps.npz')
    strain_list = data["strains"]

    #print(strain_list)
    print('number of strains')
    print(len(strain_list))

    snp_matrices = []
    matrix_lengths = []
    matrix_file_paths = []
    strain_list_masks = []
    snps_to_remove = []

    print("Reading data and calculating covariance matrix...")
    for i, in_file_name in enumerate(glob.iglob(in_glob)):
        #print(in_file_name)
        with np.load(in_file_name) as data:
            matrix, strains = data["matrix"], data["strains"]
            minor_frequencies = data["minor_frequencies"]

        # Skip SNPs with low MAF.
        columns_to_remove = minor_frequencies < maf
        matrix = matrix[:, ~columns_to_remove]

   
        # If a strain is in the group multiple times, throw away one of them (we know they are both similar in sequence,
        # otherwise they wouldn't be grouped).
        unique_strains, unique_strain_idx = np.unique(strains, return_index=True)

        # We are discarding smaller groups to prevent having to impute values for the missing strains in the small
        # groups.
        if len(unique_strain_idx) >= 100:
            continue

        # The SNP matrices are assumed to be sorted by strain. Create a NxM matrix (N = # strains, M = # SNPs) with the
        # correct rows filled in by the data from the SNP file.
        strains_list_mask = np.in1d(strain_list, unique_strains, assume_unique=True)

        full_matrix = np.empty((len(strain_list), matrix.shape[1]))
        full_matrix[:] = np.NAN
        full_matrix[strains_list_mask, :] = matrix[unique_strain_idx, :]

        snp_matrices.append(full_matrix)
        strain_list_masks.append(strains_list_mask)
        matrix_lengths.append(full_matrix.shape[1])
        matrix_file_paths.append(in_file_name)
        snps_to_remove.append(columns_to_remove)

    snp_boundaries = np.cumsum(matrix_lengths).tolist()

    print("Normalizing genotype matrix...")
    full_genotype_matrix = np.hstack(snp_matrices)
    del snp_matrices
    replace_column_nans_by_mean(full_genotype_matrix)
    full_genotype_matrix = normalize(full_genotype_matrix, axis=1)
    print("Calculating genotype matrix covariance...")
    cov = np.cov(full_genotype_matrix)
    print("Finding inverse and sqrt of covariance matrix...")
    # Genetically identical individuals results in singular (i.e. non-invertible). This can happen for a subset of the
    # data but should not happen for the entire dataset. Use the pseudo inverse instead. When a matrix is invertible,
    # its pseudo inverse is its inverse anyway. Similarly, pseudo inverses might not be positive definite, meaning that
    # we can't use Cholesky decomposition. If that happens, we can use SciPy's linalg.sqrtm() method (I don't actually
    # know if that is equivalent). Anyway, use linalg.sqrtm(linalg.pinv(cov)) when dealing with small sample sizes.
    inv_cov_sqrt = linalg.cholesky(linalg.inv(cov))
    print("Calculating pseudo SNPS...")
    pseudo_snps = np.column_stack(np.dot(inv_cov_sqrt, col) for col in full_genotype_matrix.T)
    del full_genotype_matrix

    #import matplotlib.pyplot as plt
    #plt.matshow(cov)
    #plt.show()
    #plt.matshow(np.cov(pseudo_snps))
    #plt.show()

    print("Updating files...")
    # Extract the original genes from the large pseudo SNP matrix.
    for i, (start, end) in enumerate(zip([0] + snp_boundaries, snp_boundaries)):
        strains_list_mask = strain_list_masks[i]
        snps = pseudo_snps[strains_list_mask, start:end]
        strains = strain_list[strains_list_mask]
        snps_to_discard = snps_to_remove[i]

        file_path = matrix_file_paths[i]
        file_name = file_path[file_path.rindex("/") + 1:]

        with np.load(file_path) as data:
            np.savez_compressed("{}/{}".format(out_dir, file_name), matrix=snps, strains=strains,
                                positions=data["positions"][~snps_to_discard],
                                minor_frequencies=data["minor_frequencies"][~snps_to_discard],
                                alignment_length=data["alignment_length"])

account_for_structure(in_glob= '/Users/PM/Desktop/New_data/group_snps/group*.snps.npz', out_dir = '/Users/PM/Desktop/New_data/group_snps_corrected_structure')