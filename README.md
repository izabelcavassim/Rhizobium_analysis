README
================
Maria Izabel cavassim Alves
10/26/2017

-   [Introduction](#introduction)
    -   [Genome assembly: Spades](#genome-assembly-spades)
    -   [Assembly stats: QUAST](#assembly-stats-quast)
    -   [Gene identification: Prokka](#gene-identification-prokka)
    -   [Gene alignment: ClustalO](#gene-alignment-clustalo)
    -   [SNP calling](#snp-calling)
    -   [Orthologous identification: ProteinOrtho](#orthologous-identification-proteinortho)
    -   [Orthologous identification: OrthoMCL](#orthologous-identification-orthomcl)

Introduction
============

This document was created as a resource of the scripts used in the analysis of Rhizobium paper.

Genome assembly: Spades
-----------------------

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

Assembly stats: QUAST
---------------------

Gene identification: Prokka
---------------------------

Gene alignment: ClustalO
------------------------

SNP calling
-----------

Orthologous identification: ProteinOrtho
----------------------------------------

Orthologous identification: OrthoMCL
------------------------------------
