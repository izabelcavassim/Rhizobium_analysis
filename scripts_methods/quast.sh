#!/bin/bash
# For each strain, create a quast report over SPAdes, a5, Discovar and CISA assemblies.
for f in contigs/*.fasta; do
    FILE_NAME=${f#*/}
    STRAIN_NAME=${FILE_NAME%*.fasta}
    qx --no-scratch -c 1 -t 00:10:00 "../../tools/quast-3.2/quast.py -t 1 -o assemblies/"$STRAIN_NAME"/quast contigs/"$STRAIN_NAME".fasta assemblies/"$STRAIN_NAME"/spades/contigs.fasta assemblies/"$STRAIN_NAME"/a5/"$STRAIN_NAME".contigs.fasta assemblies/"$STRAIN_NAME"/discovar/a.final/a.lines.fasta assemblies/"$STRAIN_NAME"/cisa/CISA_Out.fa"
done
