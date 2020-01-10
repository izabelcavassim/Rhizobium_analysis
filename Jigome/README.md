# Author
This script was developed by J. Peter W. Young and applied in the article: Symbiosis genes show a unique pattern of introgression and selection within a Rhizobium leguminosarum species complex (doi: https://doi.org/10.1101/526707). 
Please contact the authors for details.

# genomics
Genome analysis of rhizobia and other bacteria

# Jigome
This script takes sets of contigs assembled by SPAdes and carries out further assembly to 
create larger scaffolds.  Assembly is guided by the connections in the contig graph and by
homologous links found a set of reference genomes.  The chromosomal contigs are oriented
and ordered using a set of core genes that have the same order in all members of a set of
fully-assembled genomes. Plasmid replication genes are identified and used to label the
corresponding scaffolds.
The current version is specific for genomes of Rhizobium leguminosarum, but it could be modified for use with other bacteria.

Additional files provided are the reference sequences that Jigome uses: the 20 RepA types, the conserved chromosomal core genes in their conserved order, the dnaA sequence, and the reference genomes (for which only the contig headers are given - the actual genome sequences are too large, but can mostly be downloaded from GenBank).
