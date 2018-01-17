Introduction
============

This document was created as a repository of the scripts used in the analysis of the Rhizobium paper.

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

Gene alignment: ClustalO
------------------------

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

SNP calling
-----------

``` python
# Script first made by Asger (ld.py), and then improved by Bjarni (rhizob_ld.py ). Files are in scripts_methods.

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

Orthologous identification: OrthoMCL
------------------------------------

OrthoMCL analysis were executed with PosgresSQL instead of mySQL. More details on how to run OrthoMCL using PostgreSQL can be found in <https://bitbucket.org/asger/orthomcl-postgresql>. The output of OrthoMCL it was then parsed into fasta files. The script `disambiguate_orthomcl_groups.py` is then used to filter for paralagous genes.

``` python

import argparse
import glob
from itertools import groupby

from GFFReader import gff_iter
from disambiguate_orthomcl_groups import get_disambiguated_groups, Gene


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


def reverse_complement(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])


def get_genes_from_gff(gff_glob):
    genes = {}
    for gff_file_name in glob.iglob(gff_glob):
        gff_glob_file_name = gff_file_name[gff_file_name.rfind("/") + 1:]
        strain_name = gff_glob_file_name[:gff_glob_file_name.find(".")][0:4]
        print("Reading genes for strain {} (file {})...".format(strain_name, gff_file_name))

        genes[strain_name] = {}

        previous_gene = None
        for gff_entry in gff_iter(gff_file_name, type_filter="CDS"):
            gene_id = gff_entry.attributes["ID"].replace("PROKKA_", "")
            gene = Gene(gene_id, gff_entry.sequence, gff_entry.start, gff_entry.end, gff_entry.strand == "-")
            if previous_gene is not None and previous_gene.seq_name == gff_entry.sequence:
                gene.neighbor_l = previous_gene
                previous_gene.neighbor_r = gene

            genes[strain_name][gff_entry.attributes["ID"]] = gene
            previous_gene = gene

    return genes


def read_orthomcl_groups(file_name, ortho_name_to_strain_name=None):
    groups = []
    with open(file_name) as file_h:
        for line in file_h:
            if line.isspace():
                continue

            line_members = line.strip().split(" ")
            if ortho_name_to_strain_name is not None:
                group_members = [(ortho_name_to_strain_name[m.split("|")[0]], m.split("|")[1])
                                 for m in line_members[1:]]
            else:
                group_members = [(m.split("|")[0], m.split("|")[1]) for m in line_members[1:]]
            group_id = line_members[0][:-1]
            groups.append((group_id, group_members))

    return groups


def main():
    parser = argparse.ArgumentParser(description="Read an OrthoMCL output file and output each group as a fasta file.")
    parser.add_argument("in_file", type=str, help="the input (groups.txt) file")
    parser.add_argument("fna_glob", type=str, help="a glob pattern pointing to all fna files to load."
                                                   " Remember to quote your argument to stop wildcard expansion.")
    parser.add_argument("gff_glob", type=str, help="a glob pattern pointing to all the gff files to load."
                                                   " Remember to quote your argument to stop wildcard expansion.")
    parser.add_argument("u_out_files_dir", type=str, help="the output directory of the unambiguous group fasta files")
    parser.add_argument("a_out_files_dir", type=str, help="the output directory of the ambiguous group fasta files")

    args = parser.parse_args()
    assert args.fna_glob.startswith("~"), "You cannot use '~' in the fna glob pattern."
    assert args.gff_glob.startswith("~"), "You cannot use '~' in the gff glob pattern."

    # Extract the sequences from each strain's fasta file. Assume we can keep it all in memory!
    sequences = {}
    print("Using fna glob '{}'.".format(args.fna_glob))
    for strain_file_name in glob.iglob(args.fna_glob):
        fna_glob_file_name = strain_file_name[strain_file_name.rfind("/") + 1:]
        strain_name = fna_glob_file_name[:fna_glob_file_name.find(".")][0:4]
        print("Reading sequences for strain {} (file {})...".format(strain_name, strain_file_name))
        sequences[strain_name] = {}
        for seq_record in fasta_iter(strain_file_name):
            seq_id = seq_record[0].split(" ")[0]
            sequences[strain_name][seq_id] = seq_record[1]

    print("Using gff glob '{}'.".format(args.gff_glob))
    genes = get_genes_from_gff(args.gff_glob)

    assert set(sequences.keys()) == set(genes.keys())

    groups = read_orthomcl_groups(args.in_file)
    groups, ambiguous_groups, unambiguous_groups = get_disambiguated_groups(groups, genes, verbose=True)

    one_based = True
    start_pos_offset = 1 if one_based else 0

    for group_id, group_members in unambiguous_groups:
        write_group_fasta_file(args.u_out_files_dir, genes, group_id, group_members, sequences, start_pos_offset)
    for group_id, group_members in ambiguous_groups:
        write_group_fasta_file(args.a_out_files_dir, genes, group_id, group_members, sequences, start_pos_offset)


def write_group_fasta_file(out_files_dir, genes, group_id, group_members, sequences, start_pos_offset):
    fasta_strs = []
    for strain, gs in group_members.items():
        for gene_name in gs:
            gene = genes[strain][gene_name]
            start, end = min(gene.start, gene.end), max(gene.start, gene.end)
            gene_subseq = sequences[strain][gene.seq_name][(start - start_pos_offset):end]
            if gene.reverse:
                gene_subseq = reverse_complement(gene_subseq)

            header = "{}|{}|{}".format(group_id, strain, gene.name)
            fasta_strs.append(">{}\n{}\n".format(header, gene_subseq))
    with open("{}/{}.fna".format(out_files_dir, group_id), "w") as out_file:
        out_file.write("\n".join(fasta_strs))


if __name__ == '__main__':
    main()

# for f in group_fasta_files/u/group1000.fna; do
#     OUT="group_fasta_files/u/alns/"${f#group_fasta_files/u/*}
#     echo "$f - $OUT"
#     qx --no-scratch -t 00:15:00 "/home/agb/tools/prank -d=$f -o=$OUT -codon"
# done
```

Alignments of orthologous groups:Clustal Omega
----------------------------------------------

``` bash

#!/bin/bash
for f in u/*.fna; do
    OUT="u/normal_alns/"${f#u/*}
    echo "$f - $OUT"
    qx -A nchain --no-scratch -t 00:10:00 -m 1g "/home/agb/tools/clustalo-1.2.0 --outfmt=fasta --threads=1 --infile=$f > $OUT"
done
for f in a/*.fna; do
    OUT="a/normal_alns/"${f#a/*}
    echo "$f - $OUT"
    qx -A nchain --no-scratch -t 01:00:00 -m 1g "/home/agb/tools/clustalo-1.2.0 --outfmt=fasta --threads=1 --infile=$f > $OUT"
done
```

Useful bash commands
--------------------

``` bash

# Changing file names:

for file in *.fas; do mv "$file" "${file/_0_1p_scaffolds.fas/.contigs.fna}"; done

# Creating a folder for each file:

for file in *.contigs.fna; do
  suffix=".contigs.fna"
  dir=${file%$suffix}
  echo $dir
  mkdir -p "./$dir" &&
  cp -iv "$file" "./$dir"
done

# Activate the gwf envinroment 
gwf config set backend slurm

# Running prokka with careful parameters 
--metagenome      Improve gene predictions for highly fragmented genomes (default OFF)
--proteins [X]    Fasta file of trusted proteins to first annotate from (default '') 
The symbiotic fasta file was addedls


# Remove all folders of a directory
rm -R -- */

# Change the permission of files
chmod 777 filename

# Running proteinortho:
#!/bin/bash                                                                                                                                            
# For each protein fasta file, go through each step of proteinortho, with synteny flag enabled                                                         

qx --no-scratch -t 00:30:00 -c 16 -w "/project/NChain/faststorage/tools/proteinortho_v5.16_final/proteinortho5.pl -cpus=16 -synteny -step=1 genomes/da\
ta/*.faa"

for i in {1..199}; do
    qx --no-scratch -t 24:00:00 -c 4 "/project/NChain/faststorage/tools/proteinortho_v5.16_final/proteinortho5.pl  -cpus=4 -synteny -step=2 -jobs=$i/2\
00 genomes/data/*.faa"
done

qx --no-scratch -t 4:00:00 -c 8 "/project/NChain/faststorage/tools/proteinortho_v5.16_final/proteinortho5.pl  -cpus=8 -synteny -step=3 genomes/data/*.\
faa"

## Parsing proteinortho output


## Aligning using clustalo (fix code, i need to source python)
#!/bin/bash                                                                                                                                            
for f in data_test_clustaol/group*.fna; do
    OUT="group_alns"/${f#*/*}
    echo "$f - $OUT"
    qx -A nchain --no-scratch -t 00:10:00 -m 1g "/home/mica16/anaconda3/bin/python codon_aware_clustal.py $f ./project/clover/faststorage/tools/clustalo > $OUT"
done

## Detecting SNPs
```
