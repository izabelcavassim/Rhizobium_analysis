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
