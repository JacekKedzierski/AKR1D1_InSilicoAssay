#!/bin/bash
mkdir good
mkdir good/input
for filename in input/*.pdb; do
    plip -f "$filename" -x
    var=$(grep -20 -n "<longname>UNK</longname>" report.xml | grep -e "<num_hbd>0</num_hbd>" -e "<num_unpaired_hba>0</num_unpaired_hba>" | wc -l)
    if [ $var = "2" ]; then
        cp "$filename" good/"$filename"
    fi
done