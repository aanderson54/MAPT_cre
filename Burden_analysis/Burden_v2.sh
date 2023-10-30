#!/bin/bash

# Read input VCF
VCFID=$1

# Extract genotype information from the VCF using SnpSift
for i in ${VCFID}.vcf; do 
	cat ${i} | java -Xmx64G -jar /cluster/home/ncochran/bin/snpEff-4.3s/snpEff/SnpSift.jar extractFields - GEN[*].GT > ${VCFID}_Gen.txt 
	sed -i 's/|/\//g;s/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s/1\/0/1/g; s/\.\/\./0/g; s/0\/\./0/g; s/\.\/0/0/g' ${VCFID}_Gen.txt 
done

# Transpose, replace tabs with +, and find the sum 
awk '(NR>1)' ${VCFID}_Gen.txt | datamash transpose | sed 's/\t/+/g' | bc > ${VCFID}_wGen_xPose-Sums.txt

# Check file generated from last step. If 0, print 0; otherwise print 1. This ensures that individuals with more than one qualifying variant are only counted once.
awk '{if ($1=="0") print "0"; else print "1";}' ${VCFID}_wGen_xPose-Sums.txt > ${VCFID}_wGen_xPose-Sums-Collapse.txt

# Clean up temp files
rm ${VCFID}_Gen.txt
rm ${VCFID}_wGen_xPose-Sums.txt