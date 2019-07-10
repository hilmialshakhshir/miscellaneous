# This pipeline is specifically made to investigate Allele Specific Expression in RNA-seq data
This pipeline is made so that a BAM file is provided as input and a CSV file is output containing the base counts for heterozygous sites within the BAM file.

## What is provided within this repository
1. deepTools (v 2.0) docker image
2. GATK4 (latest version) docker image

## Steps in this pipeline
1. GATK4 AddOrReplaceReadGroups - This step is essential for GATK4's ASEReadCounter step. The BAM file is modified with ReadGroups, and a new index file is provided.
2. deeptTools computeGCBias - GC content within the BAM file is calculated. The parameters of this algorithm must be changed from run to run depending on the parameters used for sequencing (e.g. effective genome size, fragment length)
3. deeptTools correctGCBias - This step utilizes the GCBias matrix output from Step 2 and corrects for GC bias in the provided bam file.
4. GATK ASEReadCounter - Base pair counts are extracted from heterozygous sites from the GC corrected BAM file. The output is in csv format and will be present in the output directory.

## Structure of the pipeline
 - Input directory : BAM files should be placed here
 - Reference directory : Reference genome and VCF files should be placed here
 - Output directory : This is where the output should be found at the end of the analysis
 - Temporary directory : The intermediate steps produce large files that are required to be saved nor directly manipulated. This directory is emptied at the end of the analysis.
 
 ![alt text](https://github.com/hilmialshakhshir/DCM/blob/master/pipeline.png)
 
 ## Future work
 - Automate the pipeline : The pipeline can be set up as to automatically analyze BAM files as they are placed in the Input directory.
 - Designate the reference files and VCF files fixed variables : The user can then select the variable that suits their dataset.
 - Add additional functions to the pipeline depending on the scope of the analysis.
