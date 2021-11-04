#!/bin/bash

python ../scripts/create_fasta_from_ncbi_datasets.py 50 archaea,bacteria,viruses,eukaryotes > fasta.log 2>&1
cd ../db
for a in ../meta_data_20210924/contamination_*.fsa
do
  n=$(basename $a)
  n=${n%.fsa}
  echo "makeblastdb -dbtype nucl -in ${a} -title ${n} -out ${n} -taxid_map ../meta_data_20210924/${n}_taxid -parse_seqids"
done | parallel -j 6 > parallel.log 2>&1

exit 0
