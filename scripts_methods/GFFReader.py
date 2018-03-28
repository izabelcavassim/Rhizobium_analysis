from collections import namedtuple

__author__ = 'Asger Bachmann (agb@birc.au.dk)'

GFFEntry = namedtuple('GFFEntry', ['sequence', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                   'attributes'])


def gff_iter(gff_file_name, type_filter=None):
    """
    Yields GFFEntry named tuples from a GFF file. Skips comments.
    Does not support escaped values in any columns, except for tab (%09) in the attributes.
    Stops iteration at the end of file or when "##FASTA" is read.
    :param gff_file_name: the file name of the GFF file to read.
    :param type_filter: if set, return only entries of this type (e.g. gene, CDS). List of string values or single
                        string value. Case-sensitive.
    :return: GFFEntry named tuples ('sequence', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
             'attributes'). attributes is a key-value dictionary of type <string, string>.

    Example:
    >>> for entry in gff_iter("file.gff", type_filter=["CDS", "gene"]:
    >>>     print "({}, {}, {}): {}".format(entry.sequence, entry.start, entry.end, entry.attributes["ID"])
    sequence id, start and end positions, and ID of CDS' and genes in file.gff
    """
    if type_filter is not None and type(type_filter) is not list:
        type_filter = [type_filter]

    with open(gff_file_name) as gff:
        for gff_line in gff:
            if gff_line.startswith("#"):
                if gff_line.strip() == "##FASTA":
                    break
                else:
                    continue

            line_information = gff_line.split("\t")
            if type_filter is not None:
                if line_information[2] not in type_filter:
                    continue

            start = int(line_information[3])
            end = int(line_information[4])
            score = None if line_information[5] == "." else float(line_information[5])
            phase = None if line_information[7] == "." else int(line_information[7])

            attrs = {}
            for attr in line_information[8].split(";"):
                if attr.isspace():
                    continue

                attr_split = attr.split("=")
                #print attr_split
                attrs[attr_split[0]] = attr_split[1].strip().replace("%09", "\t")   

            yield GFFEntry(line_information[0], line_information[1], line_information[2], start, end, score,
                           line_information[6], phase, attrs)
