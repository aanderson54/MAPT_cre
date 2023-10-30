#!/bin/bash

module load cluster/bedtools/2.28.0
module load cluster/R/4.2.1

### Change paths to local path for users files
snp_sift="/cluster/home/ncochran/bin/snpEff_5.0/SnpSift.jar"
burden_script="Burden_v2.sh"
skat_script="2023.10.18-SKAT-ADSP_GenoTools.R"


regions_lists=("MAPT-Overall", "MAPT-ManualCuration", "MAPT-AutomatedSelection")

# Change to working directory
cd /cluster/home/ncochran/MAPT-Working/9-28-23_FreshRun


### Ensure there are VCFs, FAM, and covariate files for each region list in this working directory (FAM and cov files should be formatted for SKAT script)
### The naming scheme for the vcfs should follow this: {region list}_AF-0.001_CADD-10.vcf


# Filter VCFs based on CADD score for each region list
for j in 20; do 
	for i in "${regions_lists[@]}"; do 
		java -Xmx64G -jar ${snp_sift} filter "(CADD_1.6_phred>${j})" ${i}_AF-0.001_CADD-10.vcf > ${i}_AF-0.001_CADD-${j}.vcf 
	done
done

# Filter VCFs based on AF values for each region list
for k in 0.0001; do 
	for j in 10 20; do 
		for i in "${regions_lists[@]}"; do 
			java -Xmx64G -jar ${snp_sift} filter "!(Bravo_AF>${k}) & (AF<${k})" ${i}_AF-0.001_CADD-${j}.vcf > ${i}_AF-${k}_CADD-${j}.vcf 
		done 
	done 
done

# Filter VCFs based on AC values for each region list: note AC=1 is hard-coded at this allele frequency level given the number of individuals in either AMP-PD or ADSP.
for k in 0.00001; do 
	for j in 10 20; do 
		for i in "${regions_lists[@]}"; do 
			java -Xmx64G -jar ${snp_sift} filter "!(Bravo_AF>${k}) & (AC=1)" ${i}_AF-0.0001_CADD-${j}.vcf > ${i}_AF-${k}_CADD-${j}.vcf 
		done 
	done 
done

# Filter VCFs based on Impact for each region list - this step provides non-coding variant only files.
for k in 0.001 0.0001 0.00001; do 
	for j in 10 20; do 
		for i in "${regions_lists[@]}"; do 
			java -Xmx64G -jar ${snp_sift} filter "!((ANN[0].IMPACT = 'HIGH') | (ANN[0].IMPACT = 'MODERATE'))" ${i}_AF-${k}_CADD-${j}.vcf > ${i}_AF-${k}_CADD-${j}_NC.vcf 
		done 
	done 
done

# Create a list of vcf files
ls MAPT*AF*CADD*.vcf | sed -e "s/\.vcf//" > InputList.txt

job_ids=() # Array to hold job IDs

# Submit jobs for each input in the list and store job IDs
while read i; do 
	job_id=$(sbatch --nodes 1 --mem=20G --wrap "bash $burden_script ${i}" | awk '{print $NF}')
	job_ids+=("$job_id")
done < InputList.txt

# Wait for jobs to finish
for id in "${job_ids[@]}"; do
	while squeue -j "$id" 2>&1 >/dev/null; do 
	sleep 60
    done
done

# Create lists of FAM files and output files from burden script
ls *.fam > FamInputs.txt
ls *_xPose-Sums-Collapse.txt > MAPTinputs.txt

# Run the SKAT R script for each combination of fam file and output file
while read j; do 
	while read i; do 
		Rscript ${skat_script} ${i} Covars.cov ${j} 
	done < FamInputs.txt 
done < MAPTinputs.txt

# Add all SKAT results to the combined output
cat *SKAT-output.txt > CombinedResults.txt