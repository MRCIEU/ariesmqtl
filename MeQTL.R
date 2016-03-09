#!/usr/bin/Rscript

library(MatrixEQTL)

args <- commandArgs(trailingOnly = TRUE)

# see MeQTL for settings
useModel             = modelLINEAR; # one of: modelANOVA, modelLINEAR, or modelLINEAR_CROSS

SNP_file_name        = args[1];
expression_file_name = args[2];
covariates_file_name = character()  # set to character() for no covariates
output_file_name     = args[3];

errorCovariance      = numeric();   # set to numeric() for no covariance
pvOutputThreshold    = as.numeric(args[4]);


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA";
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 2000;
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";
cvrt$fileOmitCharacters = "NA";
cvrt$fileSkipRows = 1;
cvrt$fileSkipColumns = 1;
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);


## added r2 to output
if (nrow(me$all$eqtls) > 0) {
  assoc <- read.csv(output_file_name, sep="\t", head=T)
  assoc$r2 <- assoc$t.stat / sqrt(me$param$dfFull + assoc$t.stat^2)
  
  write.table(assoc, file=output_file_name, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
} else{
  ## remove if no associations ...
  unlink(output_file_name)
}


cat('MeQTL Done: ', me$time.in.sec, ' Seconds', '\n');
