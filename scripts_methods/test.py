import dendropy

seqstr = """\
    #NEXUS

    BEGIN TAXA;
        DIMENSIONS NTAX=13;
        TAXLABELS a1 a2 a3 b1 b2 b3 c1 c2 c3 c4 c5 d1 d2;
    END;
    BEGIN CHARACTERS;
        DIMENSIONS NCHAR=7;
        FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=.;
        MATRIX
            a1 ACCTTTG
            a2 ACCTTTG
            a3 ACCTTTG
            b1 ATCTTTG
            b2 ATCTTTG
            b3 ACCTTTG
            c1 ACCCTTG
            c2 ACCCTTG
            c3 ACCCTTG
            c4 ACCCTTG
            c5 ACCCTTG
            d1 ACAAAAG
            d2 ACCAAAG
        ;
    END
    """
seqs = dendropy.DnaCharacterMatrix.get_from_string(seqstr, 'nexus')
taxon_namespace = seqs.taxon_namespace
print(dir(taxon_namespace))


def mf(t):
    return t.label[0]

tax_parts = taxon_namespace.partition(membership_func=mf)

tax_pop_label_map = {}
for t in taxon_namespace:
    tax_pop_label_map[t] = t.label[0]

tax_parts = taxon_namespace.partition(membership_dict=tax_pop_label_map)


dict_map = parse_pop_map()
seqstr = """\
    #NEXUS

     BEGIN TAXA;
         DIMENSIONS NTAX=13;
         TAXLABELS 3404 3233 3311 3329;
     END;
     BEGIN CHARACTERS;
        DIMENSIONS NCHAR=7;
        FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=.;
        MATRIX
            3404 ACCTTTG
            3233 ACCTTTG
            3311 ACCTTTG
            3329 ATCTTTG
        ;
    END
    """

seqstr = """\
    #NEXUS

     BEGIN TAXA;
         DIMENSIONS NTAX=13;
         TAXLABELS a1 a2 b1 b2;
     END;
     BEGIN CHARACTERS;
        DIMENSIONS NCHAR=7;
        FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=.;
        MATRIX
            a1 ACCTTTG
            a2 ACCTTTG
            b1 ACCTTTG
            b2 ATCTTTG
        ;
    END
    """

seqs = dendropy.DnaCharacterMatrix.get_from_string(seqstr, 'nexus')
taxon_namespace = seqs.taxon_namespace
tax_pop_label_map = {}
for t in taxon_namespace:
    tax_pop_label_map[t] = t.label[0]

tax_parts = taxon_namespace.partition(membership_dict=tax_pop_label_map)
