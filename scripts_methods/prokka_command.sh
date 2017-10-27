#!/bin/bash
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
