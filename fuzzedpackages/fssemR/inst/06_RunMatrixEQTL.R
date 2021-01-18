# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R.
#
setwd("./data")
library(methods)
library(MatrixEQTL)

## Settings

# Linear model to use, modelANOVA or modelLINEAR
useModel = modelLINEAR; # modelANOVA or modelLINEAR

# Genotype file name
SNP_file_name = 'SNP.txt';

# Gene expression file name
expression_file_name = 'GE.txt';

# Covariates file name
# Set to character() for no covariates
# covariates_file_name = character();
covariates_file_name = 'Covariates.txt';
## covariates_file_name = character();

# Output file name
output_file_name = 'eQTL_results_R.txt';
output_file_name.cis = 'cis_eQTL_results_R.txt'

# Only associations significant at this level will be output
pvOutputThreshold = 1e-4;
pvOutputThreshold.cis = 1e-4;

# Error covariance matrix
# Set to character() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data
tic_load = proc.time()[3];

snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = 'Nh' ;# denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 10000; # read file in pieces of 10,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = '\t'; # the TAB character
gene$fileOmitCharacters = 'NA'; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 10000; # read file in pieces of 10,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = '\t'; # the TAB character
cvrt$fileOmitCharacters = 'NA'; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
cvrt$fileSliceSize = snps$nCols()+1; # read file in one piece
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

snp.pos  = readRDS("SNP.rds")
snp.pos$pos = as.numeric(snp.pos$pos)
gene.pos = readRDS("GE.rds")

toc_load = proc.time()[3];
#cat('eQTL time: ', toc_load-tic_load, ' sec\n');

## Run the analysis
{
    tic_eqtl = proc.time()[3];
    me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = output_file_name,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      pvOutputThreshold.cis = pvOutputThreshold.cis,
      output_file_name.cis = output_file_name.cis,
      snpspos = snp.pos,
      genepos = gene.pos,
      cisDist = 1e6
    )

    toc_eqtl = proc.time()[3];
}
# cat('eQTL time: ', toc_eqtl-tic_eqtl, ' sec\n');
cat('\n\n');
show(data.frame(load = toc_load-tic_load, eQTL = toc_eqtl-tic_eqtl))
