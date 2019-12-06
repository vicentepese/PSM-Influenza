#!bin/bash

# Get allele 
allele=$1

cd netMHCIIpan-3.2/tmp
echo "$allele"
../netMHCIIpan -f ./tmpSeq.fasta -s -a  $allele > ./tmpOut.out
