#!/bin/bash
mkdir output/HappyHB
for filename in output/MaeGz2Mae/*.pdb; do
    report="$filename.xml"
    plip -f "$filename" -x --name "$filename"
    var=$(grep -20 -n "<longname>UNK</longname>" "$report" | grep -e "<num_hbd>0</num_hbd>" -e "<num_unpaired_hba>0</num_unpaired_hba>" | wc -l)
    echo $var
    if [ $var = "2" ]; then
        cp "$filename" output/HappyHB/.
    fi
done