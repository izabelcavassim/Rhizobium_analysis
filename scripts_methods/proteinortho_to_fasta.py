import argparse
import glob
from itertools import groupby

from GFFReader import gff_iter


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
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])


def get_genes_from_gff(gff_glob):
    genes = {}
    for gff_file_name in glob.iglob(gff_glob):
        gff_glob_file_name = gff_file_name[gff_file_name.rfind("/") + 1:]
        strain_name = gff_glob_file_name[:gff_glob_file_name.find(".")]

        genes[strain_name] = {}

        for gff_entry in gff_iter(gff_file_name, type_filter="CDS"):
            genes[strain_name][gff_entry.attributes["ID"]] = gff_entry

    return genes


def read_proteinortho_groups(file_name):
    groups = []
    with open(file_name) as file_h:
        _ = file_h.__next__()

        for group_idx, l in enumerate(file_h):
            members_fields = l.strip().split("\t")[3:]
            members = []
            for member_field in members_fields:
                members.extend([(s, g) for s, g in (tuple(m.split("_")) for m in member_field.split(",") if m != "*")])

            groups.append((group_idx, members))

    return groups


def write_group_fasta_file(out_files_dir, genes, group_id, group_members, sequences):
    group_name = "group{}".format(group_id + 1)
    fasta_strs = []
    for strain, gene_name in group_members:
        gene = genes[strain]["{}_{}".format(strain, gene_name)]
        start, end = min(gene.start, gene.end), max(gene.start, gene.end)
        gene_subseq = sequences[strain][gene.sequence][(start - 1):end]
        if gene.strand == "-":
            gene_subseq = reverse_complement(gene_subseq)

        header = "{}|{}|{}".format(group_name, strain, gene.attributes["ID"])
        fasta_strs.append(">{}\n{}\n".format(header, gene_subseq))
    with open("{}/{}.fna".format(out_files_dir, group_name), "w") as out_file:
        out_file.write("\n".join(fasta_strs))


def main(group_file_name, out_dir, gff_glob, fna_glob):
    sequences = {}
    print("Reading fasta files...")
    for strain_file_name in glob.iglob(fna_glob):
        fna_glob_file_name = strain_file_name[strain_file_name.rfind("/") + 1:]
        strain_name = fna_glob_file_name[:fna_glob_file_name.find(".")]
        sequences[strain_name] = {}
        for seq_record in fasta_iter(strain_file_name):
            seq_id = seq_record[0].split(" ")[0]
            sequences[strain_name][seq_id] = seq_record[1]

    print("Reading genes...")
    genes = get_genes_from_gff(gff_glob)
    print("Reading groups...")
    groups = read_proteinortho_groups(group_file_name)

    print("Writing groups...")
    for group_id, group_members in groups:
        write_group_fasta_file(out_dir, genes, group_id, group_members, sequences)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read a Proteinortho output file and output each group as a fasta"
                                                 " file.")
    parser.add_argument("in_file", type=str, help="the input (groups.txt) file")
    parser.add_argument("fna_glob", type=str, help="a glob pattern pointing to all fna files to load."
                                                   " Remember to quote your argument to stop wildcard expansion.")
    parser.add_argument("gff_glob", type=str, help="a glob pattern pointing to all the gff files to load."
                                                   " Remember to quote your argument to stop wildcard expansion.")
    parser.add_argument("out_dir", type=str, help="the output directory")

    args = parser.parse_args()
    main(args.in_file, args.out_dir, args.gff_glob, args.fna_glob)
