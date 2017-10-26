import argparse
import glob
import gzip
from collections import namedtuple

import bottleneck as bn
import numpy as np
from Bio import AlignIO
from scipy import linalg

Assignment = namedtuple("Assignment", ["a", "score", "a_excl_none", "score_excl_none"])


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
    column_nanmeans = bn.nanmean(matrix, axis=0)

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


def ld(in_dir, matrices_file, strains_file, part=1, num_parts=1):
    with np.errstate(divide='ignore', invalid='ignore'):
        print("Loading variant matrices and strain lists...")
        matrices = np.load("{}/{}".format(in_dir, matrices_file))
        strain_lists = np.load("{}/{}".format(in_dir, strains_file))

        out_file = "{}/ld_part{}_{}.txt.gz".format(in_dir, part, num_parts)
        idx_range = range(part - 1, len(matrices), num_parts)

        with gzip.open(out_file, "w") as partial_out:
            for aln1_idx in idx_range:
                aln1_full_matrix = matrices[aln1_idx]
                if aln1_full_matrix is None:
                    continue

                # Holds specialized matrices for aln1 indexed by strain indices.
                specialized_matrix_cache = {}

                for j, aln2_full_matrix in enumerate(matrices[aln1_idx:]):
                    aln2_idx = aln1_idx + j

                    if aln2_full_matrix is None:
                        continue

                    common_strains_indices1, common_strains_indices2 = get_intersection_indices(strain_lists[aln1_idx],
                                                                                                strain_lists[aln2_idx])

                    if len(common_strains_indices1) == 0:
                        continue

                    common_strains_indices1_set = frozenset(common_strains_indices1)
                    if common_strains_indices1_set in specialized_matrix_cache:
                        aln1_matrix = specialized_matrix_cache[common_strains_indices1_set]
                    else:
                        aln1_matrix = specialize_complete_matrix(aln1_full_matrix, common_strains_indices1)
                        specialized_matrix_cache[common_strains_indices1_set] = aln1_matrix
                    if aln1_matrix is None:
                        continue

                    if aln1_idx == aln2_idx:
                        aln2_matrix = aln1_matrix
                    else:
                        aln2_matrix = specialize_complete_matrix(aln2_full_matrix, common_strains_indices2)
                        if aln2_matrix is None:
                            continue

                    correlation_matrix = corrcoef_mat_mat(aln1_matrix, aln2_matrix)
                    correlation_matrix = np.power(correlation_matrix, 2)
                    correlation = np.mean(correlation_matrix)
                    num_strains_in_common = len(common_strains_indices1)
                    num_aln1_snps, num_aln2_snps = correlation_matrix.shape
                    partial_out.write("{}\t{}\t{:.5f}\t{}\t{}\t{}\n".format(aln1_idx, aln2_idx, correlation,
                                                                            num_strains_in_common, num_aln1_snps,
                                                                            num_aln2_snps).encode())


def combine_parts(in_dir, strains_in_file, ld_matrix_out_file, num_strains_out_file, num_parts=1):
    ld_matrix_out_file = "{}/{}".format(in_dir, ld_matrix_out_file)
    num_strains_out_file = "{}/{}".format(in_dir, num_strains_out_file)

    out_files = ["{}/ld_part{}_{}.txt.gz".format(in_dir, part, num_parts) for part in range(1, num_parts + 1)]
    # Use the length of the array f strain lists to figure out the dimensions of the resulting LD matrix.
    strain_lists = np.load("{}/{}".format(in_dir, strains_in_file))
    n = len(strain_lists)

    ld_matrix = np.zeros((n, n))
    ld_matrix[:] = np.nan
    num_strains_matrix = np.zeros((n, n), dtype=int)

    for out_file in out_files:
        print("Opening {}...".format(out_file))
        with gzip.open(out_file) as out_file_h:
            for raw_g1, raw_g2, raw_correlation, raw_num_strains_in_common, raw_num_snps1, raw_num_snps2 in map(
                    lambda l: l.decode().strip().split("\t"), out_file_h):
                g1 = int(raw_g1)
                g2 = int(raw_g2)
                correlation = float(raw_correlation)
                ld_matrix[g1, g2] = correlation
                ld_matrix[g2, g1] = correlation
                num_strains_in_common = int(raw_num_strains_in_common)
                num_strains_matrix[g1, g2] = num_strains_in_common
                num_strains_matrix[g2, g1] = num_strains_in_common

    np.save(ld_matrix_out_file, ld_matrix)
    np.save(num_strains_out_file, num_strains_matrix)


def d_prime(matrix):
    d_prime_matrix = np.zeros((matrix.shape[1], matrix.shape[1]))
    for i, site1 in enumerate(matrix.T):
        for j, site2 in enumerate(matrix[:, i:].T):
            site1_2 = list(zip(site1, site2))
            site1_2_freq = {p: site1_2.count(p) / len(site1_2) for p in set(site1_2)}
            haplotype_freq = np.zeros((2, 2))
            haplotype_freq[0, 0] = site1_2_freq[(0, 0)] if (0, 0) in site1_2_freq else 0
            haplotype_freq[1, 0] = site1_2_freq[(1, 0)] if (1, 0) in site1_2_freq else 0
            haplotype_freq[0, 1] = site1_2_freq[(0, 1)] if (0, 1) in site1_2_freq else 0
            haplotype_freq[1, 1] = site1_2_freq[(1, 1)] if (1, 1) in site1_2_freq else 0
            # Normalize to sum 1 to remove the gap frequencies.
            haplotype_freq_sum = np.sum(haplotype_freq)
            # The sum of the haplotypes might be 0 if site 1 has gaps in all positions that site 2 has non-gaps or vice
            # versa.
            if haplotype_freq_sum != 0:
                haplotype_freq /= haplotype_freq_sum

            f_site1_0 = haplotype_freq[0, 0] + haplotype_freq[0, 1]
            # f_site1_1 = haplotype_freq[1, 0] + haplotype_freq[1, 1]
            f_site2_0 = haplotype_freq[0, 0] + haplotype_freq[1, 0]
            # f_site2_1 = haplotype_freq[0, 1] + haplotype_freq[1, 1]
            d = haplotype_freq[0, 0] - f_site1_0 * f_site2_0
            if d < 0:
                d_max = min(f_site1_0 * f_site2_0, (1 - f_site1_0) * (1 - f_site2_0))
            else:
                d_max = min(f_site1_0 * (1 - f_site2_0), (1 - f_site1_0) * f_site2_0)
            if d_max == 0:
                d = 0
            else:
                d = d / d_max
            d_prime_matrix[i, j + i] = d
            d_prime_matrix[j + i, i] = d

    return d_prime_matrix


def build_snp_matrices(in_glob, in_format, out_dir, verbose=False, sort_by_strain=True, replace_nans_by_mean=True):
    for i, in_file_name in enumerate(glob.iglob(in_glob)):
        aln = AlignIO.parse(in_file_name, in_format).__next__()
        print("Processing family {} ({})...".format(i, len(aln)))

        aln_matrix = np.array([list(str(seq.seq)) for seq in aln], dtype=np.character, order="F")

        alignment_length = aln_matrix.shape[1]

        if sort_by_strain:
            # In case a strain has multiple members in the group, use both strain and gene id as the sort key.
            strains = ["".join(seq.id.split("|")[1:]) for seq in aln]
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

        out_file_name = "{}{}.snps.npz".format(out_dir, in_file_name[in_file_name.rfind("/"):])
        np.savez_compressed(out_file_name, matrix=matrix, positions=positions, minor_frequencies=minor_frequencies,
                            alignment_length=alignment_length, strains=strains)
        if verbose:
            print("\tSNPs saved as {}.".format(out_file_name))


def intragenic_ld(in_glob, out_file_name, assignments, min_assignment_score=0., min_minor_frequency=0.,
                  min_alignment_members=0, min_alignment_length=0, max_alignment_length=None,
                  minor_frequency_bins=np.linspace(0.0, 1.0, num=21), verbose=False):
    distance_to_ld = {}

    for i, in_file_name in enumerate(glob.iglob(in_glob)):
        group_name = in_file_name[in_file_name.rfind("/") + 1:].replace(".fna.snps.npz", "")
        print("Processing family {} ({})...".format(i, group_name))

        assignment = assignments[group_name]

        if assignment.score_excl_none < min_assignment_score:
            if verbose:
                print("\tSkipping family because its assignment score is low ({}).".format(assignment.score_excl_none))
                continue

        with np.load(in_file_name) as data:
            matrix, strains = data["matrix"], data["strains"]
            positions, minor_frequencies = data["positions"], data["minor_frequencies"]
            alignment_length = data["alignment_length"]

        if alignment_length < min_alignment_length or (
                        max_alignment_length is not None and alignment_length > max_alignment_length):
            if verbose:
                print("\tSkipping family as the alignment length {} is not in [{}, {}].".format(alignment_length,
                                                                                                min_alignment_length,
                                                                                                max_alignment_length))
            continue
        elif matrix.shape[0] < min_alignment_members:
            if verbose:
                print("\tSkipping family as it has fewer than {} members.".format(min_alignment_members))
            continue

        if verbose:
            print("\tDiscarding low-frequency minor alleles...")
        columns_to_discard = minor_frequencies < min_minor_frequency

        if np.all(columns_to_discard):
            print("\t\tAll SNPs discarded. Skipping family.")
            continue

        matrix = matrix[:, ~columns_to_discard]
        positions = positions[~columns_to_discard]
        minor_frequencies = minor_frequencies[~columns_to_discard]

        if verbose:
            print("\tCalculating the correlation matrix ({})...".format(matrix.shape))
        # Use either np.corrcoef, d_prime or corrcoef_mat_mat (with normalize()).
        # d_prime is less sensitive to minor allele frequencies but also more expensive to calculate.
        # correlation_matrix = np.abs(d_prime(matrix))
        # np.corrcoef only uses numpy functions, however, it is also a bit slower than our own implementation.
        # correlation_matrix = np.corrcoef(matrix, rowvar=0) ** 2
        # corrcoef_mat_mat calculates the correlation between the columns of two matrices. It's quite fast.
        matrix = normalize(matrix)
        correlation_matrix = corrcoef_mat_mat(matrix, matrix) ** 2

        if verbose:
            print("\tCalculating positional information...")

        # Calculate the pairwise distance between entries in the position list.
        pos_vector = np.matrix(positions)
        # distance_matrix[i, j] is the difference between the positions of the SNPs giving rise to
        # correlation_matrix[i, j]. Note that distance_matrix is a symmetric matrix with diagonal 0.
        distance_matrix = np.asarray(np.abs(pos_vector - pos_vector.T))

        # Calculate whether or not some positions j and k are both 3rd bases of codons.
        end_positions = np.mod(positions, 3) == 0
        codon_end_pairs = np.logical_and(end_positions, end_positions[:, None])

        if verbose:
            print("\tBinning SNP comparisons according to minor allele frequency...")

        # Bin each pairwise SNP comparison according to minor allele frequency.
        minor_freq_bin_idx = np.digitize(minor_frequencies, minor_frequency_bins, right=True)

        if verbose:
            print("\tGroup LD according to distance...")

        # For each pairwise distance, get the corresponding LD.
        current_distance_to_ld = {}
        for j, dist_row in enumerate(distance_matrix):
            minor_freq_bin_1 = minor_frequency_bins[minor_freq_bin_idx[j]]

            for k, dist in enumerate(dist_row[j + 1:]):
                minor_freq_bin_2 = minor_frequency_bins[minor_freq_bin_idx[k + j + 1]]

                codon_end_pair = codon_end_pairs[j, k + j + 1]
                key = (dist, codon_end_pair, minor_freq_bin_1, minor_freq_bin_2, assignment.a_excl_none)
                existing_sum, existing_n = current_distance_to_ld[key] if key in current_distance_to_ld else (0, 0)
                current_distance_to_ld[key] = (existing_sum + correlation_matrix[j, k + j + 1], existing_n + 1)

        # Update the overall average. Do it here to minimize the number of multiplications necessary.
        for key, (current_sum, current_n) in current_distance_to_ld.items():
            old_n, old_avg = distance_to_ld[key] if key in distance_to_ld else (0, 0)
            new_n = old_n + current_n
            distance_to_ld[key] = (new_n, (old_n * old_avg + current_sum) / new_n)

    with open(out_file_name, "w") as out:
        # Rows of (dist, avg ld for dist, data points, freq bin 1, freq bin 2, codon_end_pair, assignment).
        out.write("\n".join("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(dist, avg_ld, n, freq_bin_1, freq_bin_2,
                                                                "T" if codon_end_pair else "F", assignment)
                            for (dist, codon_end_pair, freq_bin_1, freq_bin_2, assignment), (n, avg_ld) in
                            distance_to_ld.items()))


def get_group_assignments(in_file_name):
    with open(in_file_name) as in_file:
        return {name: Assignment(a, float(score), a_excl_none, float(score_excl_none))
                for name, a, score, a_excl_none, score_excl_none in (l.strip().split("\t") for l in in_file)}


def intergenic_ld_parts(in_glob, out_file_syn, out_file_non_syn, part, num_parts, verbose=False):
    assert part <= num_parts
    part -= 1

    if verbose:
        print("Caching group files...")

    # Caching uses approx. 2.5 GB for the 16,000 groups when caching "strains" and "matrix" and "positions".
    data_cache = {}
    file_names = glob.glob(in_glob)
    for in_file_name in file_names[part:]:
        with np.load(in_file_name) as data:
            data_cache[in_file_name] = {k: data[k] for k in ["strains", "matrix", "positions"]}

    if verbose:
        print("Caching done.")

    n = len(file_names)
    syn_p = np.zeros((n, n))
    np.fill_diagonal(syn_p, 1)
    non_syn_p = np.zeros((n, n))
    np.fill_diagonal(non_syn_p, 1)
    strains_in_common = np.zeros((n, n), dtype=int)
    num_pos_syn_1 = np.zeros((n, n), dtype=int)
    num_pos_syn_2 = np.zeros_like(num_pos_syn_1)
    num_pos_non_syn_1 = np.zeros((n, n), dtype=int)
    num_pos_non_syn_2 = np.zeros_like(num_pos_non_syn_1)

    jobs = range(part, n, num_parts)

    for idx_1 in jobs:
        print("Processing group id {}...".format(idx_1))

        in_file_name_1 = file_names[idx_1]
        data_1 = data_cache[in_file_name_1]
        # Fill out the diagonal of the strains_in_common matrix.
        strains_in_common[idx_1, idx_1] = len(data_1["strains"])

        for idx_2 in range(idx_1 + 1, len(file_names)):
            in_file_name_2 = file_names[idx_2]
            data_2 = data_cache[in_file_name_2]

            if verbose:
                print("Calculating the strain intersection between the two groups...")
            common_strains_indices1, common_strains_indices2 = get_intersection_indices(data_1["strains"],
                                                                                        data_2["strains"])

            if len(common_strains_indices1) == 0:
                if verbose:
                    print("The groups have no strains in common. Skipping pairwise comparison.")
                continue

            # Consider caching specialized matrices by their strain indexes.
            matrix_1, matrix_1_positions = specialize_complete_matrix(data_1["matrix"], common_strains_indices1,
                                                                      data_1["positions"])
            if matrix_1 is None:
                if verbose:
                    print("Specializing group 1 resulted in a non-SNP matrix. Skipping pairwise comparison.")
                continue
            matrix_2, matrix_2_positions = specialize_complete_matrix(data_2["matrix"], common_strains_indices2,
                                                                      data_2["positions"])
            if matrix_2 is None:
                if verbose:
                    print("Specializing group 2 resulted in a non-SNP matrix. Skipping pairwise comparison.")
                continue

            if verbose:
                print("Calculating correlation matrix...")

            matrix_1_syn_pos = np.mod(matrix_1_positions, 3) == 0
            matrix_2_syn_pos = np.mod(matrix_2_positions, 3) == 0

            num_matrix_1_syn_pos = np.sum(matrix_1_syn_pos)
            num_matrix_2_syn_pos = np.sum(matrix_2_syn_pos)
            if num_matrix_1_syn_pos > 0 and num_matrix_2_syn_pos > 0:
                matrix_1_syn = matrix_1[:, matrix_1_syn_pos]
                matrix_2_syn = matrix_2[:, matrix_2_syn_pos]
                correlation_matrix_syn = np.mean(np.power(corrcoef_mat_mat(matrix_1_syn, matrix_2_syn), 2))
                syn_p[idx_1, idx_2] = correlation_matrix_syn
                syn_p[idx_2, idx_1] = correlation_matrix_syn

                num_pos_syn_1[idx_1, idx_2] = num_matrix_1_syn_pos
                num_pos_syn_1[idx_2, idx_1] = num_matrix_1_syn_pos
                num_pos_syn_2[idx_2, idx_1] = num_matrix_2_syn_pos
                num_pos_syn_2[idx_1, idx_2] = num_matrix_2_syn_pos

            num_matrix_1_non_syn_pos = len(matrix_1_positions) - num_matrix_1_syn_pos
            num_matrix_2_non_syn_pos = len(matrix_2_positions) - num_matrix_2_syn_pos
            if num_matrix_1_non_syn_pos > 0 and num_matrix_2_non_syn_pos > 0:
                matrix_1_non_syn = matrix_1[:, ~matrix_1_syn_pos]
                matrix_2_non_syn = matrix_2[:, ~matrix_2_syn_pos]
                correlation_matrix_non_syn = np.mean(np.power(corrcoef_mat_mat(matrix_1_non_syn, matrix_2_non_syn), 2))
                non_syn_p[idx_1, idx_2] = correlation_matrix_non_syn
                non_syn_p[idx_2, idx_1] = correlation_matrix_non_syn

                num_pos_non_syn_1[idx_1, idx_2] = num_matrix_1_non_syn_pos
                num_pos_non_syn_1[idx_2, idx_1] = num_matrix_1_non_syn_pos
                num_pos_non_syn_2[idx_2, idx_1] = num_matrix_2_non_syn_pos
                num_pos_non_syn_2[idx_1, idx_2] = num_matrix_2_non_syn_pos

            strains_in_common[idx_1, idx_2] = len(common_strains_indices1)
            strains_in_common[idx_2, idx_1] = len(common_strains_indices1)

    np.savez_compressed("{}_part{}_{}".format(out_file_syn, part + 1, num_parts), p=syn_p,
                        strains_in_common=strains_in_common, num_pos_1=num_pos_syn_1, num_pos_2=num_pos_syn_2)
    np.savez_compressed("{}_part{}_{}".format(out_file_non_syn, part + 1, num_parts), p=non_syn_p,
                        strains_in_common=strains_in_common, num_pos_1=num_pos_non_syn_1, num_pos_2=num_pos_non_syn_2)


def account_for_structure(in_glob, out_dir, maf=0.0):
    with open("/home/agb/data/rhizo/full_dataset/orthologs/strain_list.txt") as strain_list_in:
        strain_list = np.array([l.strip() for l in strain_list_in])

    snp_matrices = []
    matrix_lengths = []
    matrix_file_paths = []
    strain_list_masks = []
    snps_to_remove = []

    print("Reading data and calculating covariance matrix...")
    for i, in_file_name in enumerate(glob.iglob(in_glob)):
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
        if len(unique_strain_idx) < 100:
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

    import matplotlib.pyplot as plt
    plt.matshow(cov)
    plt.show()
    plt.matshow(np.cov(pseudo_snps))
    plt.show()

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


def combine_intergenic_ld_parts(in_glob, out_file_name):
    combined_p = None
    combined_strains_in_common = None
    combined_num_pos_1 = None
    combined_num_pos_2 = None
    for file_name in glob.iglob(in_glob):
        print("\tReading {}".format(file_name))
        with np.load(file_name) as data:
            p = data["p"]
            strains_in_common = data["strains_in_common"]
            num_pos_1 = data["num_pos_1"]
            num_pos_2 = data["num_pos_2"]

            if combined_p is None:
                combined_p = p
                combined_strains_in_common = strains_in_common
                combined_num_pos_1 = num_pos_1
                combined_num_pos_2 = num_pos_2
            else:
                entries_to_replace = p > combined_p
                combined_p[entries_to_replace] = p[entries_to_replace]
                combined_strains_in_common[entries_to_replace] = strains_in_common[entries_to_replace]
                combined_num_pos_1[entries_to_replace] = num_pos_1[entries_to_replace]
                combined_num_pos_2[entries_to_replace] = num_pos_2[entries_to_replace]

    np.savez_compressed(out_file_name, p=combined_p, strains_in_common=combined_strains_in_common,
                        num_pos_1=combined_num_pos_1, num_pos_2=combined_num_pos_2)


def main(args):
    build_snp_matrices("/home/agb/data/rhizo/full_dataset/orthologs/group_alns/*.fna", "fasta",
                       "/home/agb/data/rhizo/full_dataset/orthologs/group_snps", verbose=True, sort_by_strain=True,
                       replace_nans_by_mean=False)
    # with np.errstate(divide='ignore', invalid='ignore'):
    #     account_for_structure("/home/agb/data/rhizo/full_dataset/orthologs/group_intergenic_snps/group*.snps.npz",
    #                           "/home/agb/data/rhizo/full_dataset/orthologs/group_snps_no_structure", maf=0.1)
    # assignments = get_group_assignments("/home/agb/data/rhizo/full_dataset/orthologs/200_rhizobium.poff.assigned.txt")
    # intragenic_ld("/home/agb/data/rhizo/full_dataset/orthologs/group_snps_no_structure/group*.snps.npz",
    #               "/home/agb/data/rhizo/full_dataset/ld/test/intragenic_ld_full_100members.tsv", assignments,
    #               min_alignment_members=100, min_alignment_length=0, max_alignment_length=None, verbose=False)

    # build_snp_matrices(args.in_glob, "fasta", args.snp_output_dir, verbose=False, sort_by_strain=True,
    #                    replace_nans_by_mean=False)
    # intergenic_ld_parts(args.in_glob, "./intergenic_ld_parts/syn_ld.npz", "./intergenic_ld_parts/non_syn_ld.npz",
    #                     part=args.part, num_parts=args.num_parts, verbose=False)
    # combine_intergenic_ld_parts("./intergenic_ld_parts/syn_ld.npz_part*_*", "./syn_ld.npz")
    # combine_intergenic_ld_parts("./intergenic_ld_parts/non_syn_ld.npz_part*_*", "./non_syn_ld.npz")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_glob", dest="in_glob", type=str,
                        default="/home/agb/data/rhizo/full_dataset/orthologs/group_intergenic_snps/group*.snps.npz")
    parser.add_argument("--snp_output_dir", dest="snp_output_dir", type=str,
                        default="/home/agb/data/rhizo/full_dataset/orthologs/group_intergenic_snps")
    parser.add_argument("--part", dest="part", default=1, type=int)
    parser.add_argument("--num_parts", dest="num_parts", default=1, type=int)
    args = parser.parse_args()

    main(args)
