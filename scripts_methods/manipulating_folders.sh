#!/bin/bash

# Making a directory containing all the fasta files
mkdir data2
for i in */*/*.faa; do cp -iv $i "./data2"; done
