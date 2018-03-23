import argparse
import glob
import gzip
from collections import namedtuple

#import bottleneck as bn
import numpy as np
from Bio import AlignIO
#from scipy import linalg


def get_encoded_non_bi_allelic_columns(strain_to_seq_map, locus):
    nrows = len(strain_to_seq_map.keys())
    strains = sorted(strain_to_seq_map.keys())

    columns = []
    for fixed_ref in ['A', 'C', 'G', 'T']:
        fixed_ref_seen = False
        column = np.zeros(nrows)
        for s_idx, s in enumerate(strains):
            allele = strain_to_seq_map[s][locus]
            value = 0
            if allele != "-":
                if allele != fixed_ref:
                    value = 1
                else:
                    fixed_ref_seen = True
            else:
                # Gap/dash
                value = np.nan

            column[s_idx] = value

        if fixed_ref_seen:
            columns.append(column)

    return columns


def unique_with_counts(ar):
    """A simplified version of numpy.unique() focusing on minimizing the conditionals to maximize performance."""
    ar.sort()
    flag = np.concatenate(([True], ar[1:] != ar[:-1]))
    idx = np.concatenate(np.nonzero(flag) + ([ar.size],))
    return ar[flag], np.diff(idx)


def build_complete_matrix(aln_matrix, skip_multi_allelic_variants=False, min_minor_frequency=0.):
    nrows, ncols = aln_matrix.shape
    encoded_matrix = np.zeros((nrows, ncols), dtype=float)
    gap_positions = aln_matrix == b'-'
    gap_counts = gap_positions.sum(axis=0)
    encoded_matrix[gap_positions] = np.nan

    minor_frequencies = []
    # An array mapping the entries in the columns array to 1-based positions in the sequences.
    column_idx_to_keep = []
    for locus_idx, locus in enumerate(aln_matrix.T):
        column = encoded_matrix[:, locus_idx]
        unique, counts = unique_with_counts(locus[~gap_positions[:, locus_idx]])

        if len(unique) < 2:
            # Not a variant.
            continue
        if skip_multi_allelic_variants and len(unique) > 2:
            continue

        # Multi-allelic sites are encoded treated as bi-allelic with the major allele(s) as 1 and the rest as 0.
        freqs_excl_gaps = counts / nrows
        max_freq_idx = np.argmax(freqs_excl_gaps)
        ref, ref_freq = unique[max_freq_idx], freqs_excl_gaps[max_freq_idx]
        ref_is_major = ref_freq >= sum(freqs_excl_gaps) / 2

        if ref_is_major:
            gap_freq = gap_counts[locus_idx] / nrows
            minor_frequency = 1 - ref_freq - gap_freq
        else:
            minor_frequency = ref_freq
        if minor_frequency < min_minor_frequency:
            # The minor allele(s) should have a frequency higher than min_minor_frequency.
            continue

        if ref_is_major:
            column[locus == ref] = 1
        else:
            column[np.logical_and(~gap_positions[:, locus_idx], locus != ref)] = 1
        # Minor alleles are 0 because of the initialization to 0.

        column_idx_to_keep.append(locus_idx)
        minor_frequencies.append(minor_frequency)

    if len(column_idx_to_keep) == 0:
        return None, None, None

    encoded_matrix = encoded_matrix[:, column_idx_to_keep]
    positions = np.array(column_idx_to_keep) + 1
    return encoded_matrix, positions, np.array(minor_frequencies)


def columns_all_equal(data):
    # Assume all NaNs have been removed/replaced by the column means.
    # If a column is all 0s or all 1s, the sum of the column is 0 or equal to the number of rows in the column.
    colsum = np.sum(data, axis=0)
    nrows = data.shape[0]
    # noinspection PyTypeChecker
    return np.isclose(colsum, 0) | np.isclose(colsum, nrows)


def replace_column_nans_by_mean(matrix):
    # Set the value of gaps/dashes in each column to be the average of the other values in the column.
    nan_indices = np.where(np.isnan(matrix))
    # Note: bn.nanmean() instead of np.nanmean() because it is a lot(!) faster.
    column_nanmeans = np.nanmean(matrix, axis=0)

    # For each column, assign the NaNs in that column the column's mean.
    # See http://stackoverflow.com/a/18689440 for an explanation of the following line.
    matrix[nan_indices] = np.take(column_nanmeans, nan_indices[1])


def normalize_columns(matrix):
    """Normalizes a matrix' columns to mean 0 and std 1."""
    mean = np.mean(matrix, axis=0)
    std = np.std(matrix, axis=0)
    matrix = (matrix - mean) / std
    return np.nan_to_num(matrix)


def normalize(matrix, axis=0):
    """Normalizes a matrix' columns or rows to mean 0 and std 1."""
    mean = np.mean(matrix, axis=axis)
    std = np.std(matrix, axis=axis)
    if axis == 1:
        mean = mean[:, None]
        std = std[:, None]

    matrix = (matrix - mean) / std
    return np.nan_to_num(matrix)


def specialize_complete_matrix(matrix, strain_indices, positions):
    subset_rows = np.copy(matrix[strain_indices, :])
    replace_column_nans_by_mean(subset_rows)

    # Now that we have removed some rows, we might end up with non-SNP positions (the columns are all 0s or 1s).
    # Filter out those columns.
    columns_to_keep = ~columns_all_equal(subset_rows)
    subset_matrix = subset_rows[:, columns_to_keep]
    if subset_matrix.shape[1] == 0:
        return None, None

    return normalize(subset_matrix), positions[columns_to_keep]


def get_intersection_indices(l1, l2):
    n1 = len(l1)
    n2 = len(l2)

    indices1 = []
    indices2 = []

    i = 0
    j = 0
    while i < n1 and j < n2:
        if l1[i] > l2[j]:
            j += 1
        elif l1[i] < l2[j]:
            i += 1
        else:
            indices1.append(i)
            indices2.append(j)
            i += 1
            j += 1

    return indices1, indices2


def calculate_variant_matrices(in_dir, in_file_name, matrices_out_file, strains_out_file, in_file_format="clustal"):
    variant_matrices_out_file = "{}/{}".format(in_dir, matrices_out_file)
    variant_matrices_strains_out_file = "{}/{}".format(in_dir, strains_out_file)

    with open("{}/{}".format(in_dir, in_file_name)) as in_file_h, np.errstate(divide='ignore', invalid='ignore'):
        alns = tuple(AlignIO.parse(in_file_h, in_file_format))
        print("Creating variant matrices for each gene family...")
        matrices = []
        strain_lists = []
        for i, aln in enumerate(alns):
            print("\tCreating variant matrix for family {}...".format(i))
            strain_to_seq_map = {seq.id.split("|")[1]: str(seq.seq) for seq in aln}
            strains = strain_to_seq_map.keys()
            matrix, _ = build_complete_matrix(len(aln[0]), strain_to_seq_map)
            matrices.append(matrix)
            strain_lists.append(strains)

    print("Saving variant matrices to disk...")
    np.save(variant_matrices_out_file, matrices)
    np.save(variant_matrices_strains_out_file, strain_lists)


def corrcoef_mat_vec(m, v):
    """Calculates correlation coefficients between the columns of a matrix and fixed vector (both with mean 0)."""
    r_num = m.T.dot(v)
    return r_num / len(v)


def corrcoef_mat_mat(m1, m2):
    """Calculates the correlation coefficients between all column pairs of two matrices.
    Assumes the columns have mean 0 and std 1."""
    r = np.zeros((m1.shape[1], m2.shape[1]))

    for idx, m2_col in enumerate(m2.T):
        r[:, idx] = corrcoef_mat_vec(m1, m2_col)

    return r

def build_snp_matrices(in_glob, in_format, out_dir, verbose=False, sort_by_strain=True, replace_nans_by_mean=True):
    for i, in_file_name in enumerate(glob.iglob(in_glob)):
        print(in_file_name)
        aln = AlignIO.parse(in_file_name, in_format).__next__()
        print("Processing family {} ({})...".format(i, len(aln)))

        aln_matrix = np.array([list(str(seq.seq)) for seq in aln], dtype=np.character, order="F")
        #print(aln_matrix)
        alignment_length = aln_matrix.shape[1]

        if sort_by_strain:
            # In case a strain has multiple members in the group, use both strain and gene id as the sort key.
            strains = ["".join(seq.id.split("|")[1:]) for seq in aln]
            print(strains)
            aln_matrix = aln_matrix[np.argsort(strains), :]
            # For the strain array that is saved with the SNP matrix, don't include the gene id in the strain list.
            strains = sorted([seq.id.split("|")[1] for seq in aln])
        else:
            strains = [seq.id.split("|")[1] for seq in aln]

        matrix, positions, minor_frequencies = build_complete_matrix(aln_matrix)

        if matrix is None or matrix.shape[1] == 1:
            if verbose:
                print("\tFamily skipped because there is zero or 1 SNPs.")
            continue

        if replace_nans_by_mean:
            # Replace the gaps/dashes (NaNs) in each column by the column average.
            replace_column_nans_by_mean(matrix)
        #print(matrix)

        out_file_name = "{}{}.snps.npz".format(out_dir, in_file_name[in_file_name.rfind("/"):])
        #print(out_file_name)
        np.savez_compressed(out_file_name, matrix=matrix, positions=positions, minor_frequencies=minor_frequencies,
                            alignment_length=alignment_length, strains=strains)
        if verbose:
            print("\tSNPs saved as {}.".format(out_file_name))


print(build_snp_matrices(in_glob = "/Users/PM/Desktop/New_data/group_alns/*.fna", in_format = "fasta", out_dir = "/Users/PM/Desktop/New_data/", verbose=True, sort_by_strain=True,
                       replace_nans_by_mean=True))

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument("--in_glob", dest="in_glob", type=str,
#                         default="/Users/PM/Desktop/group_alns/*.fna")
#     parser.add_argument("--snp_output_dir", dest="snp_output_dir", type=str,
#                         default="/Users/PM/Desktop/group_snps/")
#     args = parser.parse_args()

#     main(args)
