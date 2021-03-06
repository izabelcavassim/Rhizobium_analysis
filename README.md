Introduction
============

This document was created as a repository of the scripts used in the analysis of the article: Symbiosis genes show a unique pattern of introgression and selection within a Rhizobium leguminosarum species complex (doi: https://doi.org/10.1099/mgen.0.000351).

The data used in this analysis is available in the following link: 

   https://www.dropbox.com/sh/6fceqmwfa3p3fm6/AAAkFIRCf7ZxgO1a4fHv3FeOa?dl=0
   
   The folder contains the following data:
   * **Group_fnas**: Fasta files of orthologous genes (no-alignment) 
   * **Group_alns**: Aligned fasta files of orthologous genes 
   * **corrected_snps**: snp files of polymorphic genes corrected by population structure
   * **ani_sorted_by_genospecies_282_sorted**: Average nucleotide identity (ANI) across 282 conserved core genes
   * **README** file explaining the rest of the data provided

It can also be accessed in the open repository:
https://doi.org/10.6084/m9.figshare.11568894.v1

Genome assemblies are found at NCBI under the project: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA510726/

Genome assembly: Spades
-----------------------

``` bash
#!/bin/bash
for f in ../reads/*_1_trimmed.fastq.gz; do
    FILE=${f#../reads/*}
    STRAIN=${FILE%*_1_trimmed.fastq.gz}
    if [ -d $STRAIN"/"$STRAIN"_spades_careful" ]; then
        echo $STRAIN": _spades_careful directory exists. Skipping strain."
        continue
    fi

    mkdir -p $STRAIN
    qx --no-scratch -c 8 -m 32g -t 05:00:00 "/project/clover/faststorage/tools/SPAdes-3.6.2-Linux/bin/spades.py --careful -t 8 -m 32 -o $STRAIN"/"$STRAIN"_spades_careful" -1 $f -2 ${f%*1_trimmed.fastq.gz}"2_trimmed.fastq.gz" -s ${f%*1_trimmed.fastq.gz}"U1_trimmed.fastq.gz" -s ${f%*1_trimmed.fastq.gz}"U2_trimmed.fastq.gz"";
done
```

Assembly stats: QUAST
---------------------

``` bash
#!/bin/bash
# For each strain, create a quast report over SPAdes, a5, Discovar and CISA assemblies.
for f in contigs/*.fasta; do
    FILE_NAME=${f#*/}
    STRAIN_NAME=${FILE_NAME%*.fasta}
    qx --no-scratch -c 1 -t 00:10:00 "../../tools/quast-3.2/quast.py -t 1 -o assemblies/"$STRAIN_NAME"/quast contigs/"$STRAIN_NAME".fasta assemblies/"$STRAIN_NAME"/spades/contigs.fasta assemblies/"$STRAIN_NAME"/a5/"$STRAIN_NAME".contigs.fasta assemblies/"$STRAIN_NAME"/discovar/a.final/a.lines.fasta assemblies/"$STRAIN_NAME"/cisa/CISA_Out.fa"
done
```

Gene identification: Prokka
---------------------------

``` bash
#!/bin/bash
# For each contig run prokka, many output files are generated: .gff, .gbk, .fna, .faa and so on.
for f in */*.contigs.fna; do
    STRAIN=$(IFS="/"; set -- $f; echo $1)
    if [ -f $STRAIN"/"$STRAIN".gff" ]; then
    echo $STRAIN": gff file exists. Skipping strain."
    continue
    fi
    cd $STRAIN
    qx --no-scratch -t 10:00:00 -c 4 "/project/clover/faststorage/tools/prokka-1.11/bin/prokka -cpus=4 *.contigs.fna"
    cd ..
done
```

``` python
#Same analysis using gwf:

from gwf import Workflow
import glob
import os

gwf = Workflow()

def prokka_run(directory, ID):
    outdir = f'{directory}/{ID}'
    infile = f'{directory}/{ID}.contigs.fna'
    inputs = [f'{directory}/{ID}.contigs.fna']
    outputs = [f'{outdir}/{ID}{x}' for x in ['.gff', '.gbk', '.fna', '.faa']]
    #print('\n'.join([outdir, infile, "".join(outfiles)]))
    options = {
        #'inputs': [infile],
        #'outputs': outfiles,
        'memory': '4g',
        'cores': '4',
        'walltime': '10:00:00',
        'account': 'NChain'
    }

    spec = f"/project/NChain/faststorage/tools/prokka-1.12/bin/prokka --outdir {outdir} --prefix {ID}  {infile}"
    print(spec)

    return inputs, outputs, options, spec

# For each fasta file run prokka, many output files are generated: .gff, .gbk, .fna, .faa and so on.
# Find all the subdirectories existent in that directory
d = '.'
folders = [os.path.join(d, o) for o in os.listdir(d) 
                    if os.path.isdir(os.path.join(d,o))]

for strain in folders: 
    print(strain)
    if strain[2] == '3':
        workflow_name = 'prokka{}'.format(strain[2:6])
        gwf.target_from_template(workflow_name, prokka_run(directory = strain, ID = strain[2::])) 
```

Gene alignment: Clustalo
------------------------

``` bash
#!/bin/bash                                                                                                          
for f in data_test_clustaol/group*.fna; do
    OUT="group_alns"/${f#*/*}
    echo "$f - $OUT"
    qx -A nchain --no-scratch -t 00:10:00 -m 1g "/home/mica16/anaconda2/envs/py36/bin/python codon_aware_clustal.py \
$f /project/NChain/faststorage/tools/clustalo > $OUT"
done
```

Codon aware alignment
---------------------

``` python

# Group fasta files are aligned in a codon aware matter. The clustalO runs with default parameters. 

import argparse
import io
import subprocess
from itertools import groupby

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def fasta_iter(fasta_file):
    """
    Given a FASTA file, yield tuples of header, sequence.
    Author: Brent Pedersen (https://www.biostars.org/p/710/#1412)
    """
    fh = open(fasta_file)
    # Ditch the boolean (x[0]) and just keep the header or sequence since we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # Drop the ">".
        header = header.__next__()[1:].strip()
        # Join all sequence lines to one in case it spans multiple lines.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield header, seq

def seq_to_str(header, seq):
    return ">{}\n{}\n".format(header, seq)

def back_translate(fasta_aln, original_sequences):
    string_handle = io.StringIO()
    string_handle.write(fasta_aln)
    string_handle.seek(0)
    outer_alignments = AlignIO.parse(string_handle, "fasta")
    back_translated_records = []
    for alignments in outer_alignments:
        for aln in alignments:
            original_seq = original_sequences[aln.id]

            back_translated = []
            gaps = 0
            for i in range(len(aln.seq)):
                if aln.seq[i] == "-":
                    gaps += 1
                    back_translated.append("---")
                else:
                    idx = (i - gaps) * 3
                    back_translated.append(original_seq[idx:(idx + 3)])

            back_translated_records.append(SeqRecord(Seq("".join(back_translated), generic_dna), aln.id))

    return MultipleSeqAlignment(back_translated_records, alphabet=generic_dna)


def msa_to_str(msa):
    # msa.format("fasta") prints unecessary things such as description. This is a simpler version.
    return "\n".join(seq_to_str(aln.id, aln.seq) for aln in sorted(msa, key=lambda aln: aln.id))


def main():
    parser = argparse.ArgumentParser(description="Read a fasta file and align the sequences in it using ClustalO in a "
                                                 "codon-aware style.")
    parser.add_argument("in_file", type=str, help="the input fasta file")
    parser.add_argument("clustal_o_exec", type=str, help="the path to the ClustalO executeable", default="clustalo")
    parser.add_argument('--output_aa', dest='output_aa', action='store_true')
    parser.add_argument("--threads", type=int, help="the number of threads to use", default=1)
    args = parser.parse_args()

    original_sequences = {}
    for h, s in fasta_iter(args.in_file):
        orig_h = h
        i = 2
        while h in original_sequences:
            h = "{}_{}".format(orig_h, i)
            i += 1
        original_sequences[h] = s

    translated_sequences = [(header, Seq(seq).translate()) for header, seq in original_sequences.items()]
    p = subprocess.Popen([args.clustal_o_exec, "--outfmt=fasta", "--seqtype=Protein",
                          "--threads={}".format(str(args.threads)), "--infile=-"],
                         stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    aln_out = p.communicate(input="\n".join(seq_to_str(h, s) for h, s in translated_sequences).encode())
    if args.output_aa:
        print(aln_out[0].decode())
    else:
        back_translated_alignment = back_translate(aln_out[0].decode(), original_sequences)
        print(msa_to_str(back_translated_alignment))


if __name__ == '__main__':
    main()
```

Orthologous identification: ProteinOrtho
----------------------------------------

``` bash
#!/bin/bash
# For each protein fasta file, go through each step of proteinortho, with synteny flag enabled

qx --no-scratch -t 00:30:00 -c 16 -w "/home/agb/tools/proteinortho_v5.13/proteinortho5.pl -cpus=16 -synteny -step=1 data/*.faa"

for i in {1..199}; do
    qx --no-scratch -t 24:00:00 -c 4 "/home/agb/tools/proteinortho_v5.13/proteinortho5.pl -cpus=4 -synteny -step=2 -startat=$i -stopat=$i data/*.faa"
done

qx --no-scratch -t 4:00:00 -c 8 "/home/agb/tools/proteinortho_v5.13/proteinortho5.pl -cpus=8 -synteny -step=3 data/*.faa"
```

SNP calling
-----------

``` python
# Script first made by Asger (ld.py), and then improved by Bjarni (rhizob_ld.py ). Files are in scripts_methods.

import argparse
import glob
import gzip
from collections import namedtuple

#import bottleneck as bn
import numpy as np
from Bio import AlignIO
from scipy import linalg


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
        print(aln)
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
        print(out_file_name)
        np.savez_compressed(out_file_name, matrix=matrix, positions=positions, minor_frequencies=minor_frequencies,
                            alignment_length=alignment_length, strains=strains)
        if verbose:
            print("\tSNPs saved as {}.".format(out_file_name))


print(build_snp_matrices(in_glob = "/home/mica16/NChain/faststorage/rhizobium/orthologs_proteinortho/group_alns/*.fna", in_format = "fasta", out_dir = "./group_snps", verbose=True, sort_by_strain=True, replace_nans_by_mean=False))
                       
```
