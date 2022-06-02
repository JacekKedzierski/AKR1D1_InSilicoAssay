mkdir -p output/ligands/smi

while read line; do
  line=($line)
  echo "${line[0]}" > "output/CreateSmi/${line[1]}.smi"
done < $1
