# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)
library(readr)

## Location of the package with the data files.
#base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "/Users/hilmi/DCMUM/snp_maf.csv";
snps_location_file_name = "/Users/hilmi/DCMUM/snpspos.csv";

# Gene expression file name
expression_file_name = "/Users/hilmi/Documents/gene_vst.csv";
gene_location_file_name = "/Users/hilmi/DCMUM/genepos_final.csv";

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "/Users/hilmi/DCMUM/cvrtv2.csv";

# Output file name
output_file_name_cis = tempfile("ciseQTL", fileext = ".csv");
output_file_name_tra = tempfile("transeQTL", fileext = ".csv");

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 5e-2; 
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = ",";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = ",";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name,
                       header = TRUE, sep = ",", stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name,
                        header = TRUE, sep = ",", stringsAsFactors = FALSE);

me.cvrt = Matrix_eQTL_main(
        snps = snps, 
        gene = gene, 
        cvrt = cvrt,
        output_file_name  = output_file_name_tra,
        pvOutputThreshold = pvOutputThreshold_tra,
        useModel = useModel, 
        errorCovariance = errorCovariance, 
        verbose = TRUE, 
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos, 
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

message('Analysis done in: ', me.cvrt$time.in.sec, ' seconds');
message('Detected local eQTLs:');
show(me.cvrt$cis$eqtls);
message('Detected distant eQTLs:');
show(me.cvrt$trans$eqtls);

## Plot the histogram of local and distant p-values

plot(me.cvrt);
