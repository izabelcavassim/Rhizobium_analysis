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
