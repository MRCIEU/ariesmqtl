#!/usr/bin/Rscript

library(argparse)

library(foreach)
library(doMC)

registerDoMC(detectCores())

parser <- ArgumentParser()

parser$add_argument("--meth", type="character")
parser$add_argument("--keep-probes", type="character")
parser$add_argument("--covars", type="character")

Args <- parser$parse_args()


# read data
print(paste("reading meth: ", Args$meth, sep=""))
load(Args$meth) # M
dim(M)

print(paste("reading keep-probes: ", Args$meth, sep=""))
if (!is.null(Args$keep_probes)) {
  Args$keep_probes <- scan(Args$keep_probes, what=character())
} else {
  Args$keep_probes <- rownames(M)
}
dim(Args$keep_probes)

print("reading geno: ./parsed_data/ALSPAC.1.tab")
geno <- read.csv("./parsed_data/ALSPAC.1.tab", head=T, row.names=1, sep="\t", check.names=F) # pick one at random (others have same header ...)
dim(geno)

print(paste("reading covars: ", Args$covars, sep=""))
covar <- read.csv(Args$covars, head=T, row.names=1, sep="\t")
dim(covar)


# FETCH COMMON ALNs ...
ALNs <- colnames(M)

print(paste("ALNs in ARIES: ", length(ALNs), sep=""))

diff <- setdiff(ALNs, rownames(geno))
print(paste("Not Genotyped: ", length(diff), sep=""))
print(diff)
ALNs <- intersect(ALNs, rownames(geno))

covar <- na.omit(covar) # remove those with no PCs ...
diff  <- setdiff(ALNs, rownames(covar))
print(paste(" Missing Covar (PCs?): ", length(diff), sep=""))
print(diff)
ALNs  <- intersect(ALNs, rownames(covar))

print(paste("ALNs in FINAL DATASET: ", length(ALNs), sep=""))


# REGRESS COVARIATES & WRITE DATA
M     <- M[Args$keep_probes, ALNs, drop=F]

print("FINAL METH:")
dim(M)

covar <- covar[ALNs,]
covar$Batch <- as.factor(covar$Batch)

if (length(which(names(covar)=="Sex")) == 1) {
  covar$Sex <- as.factor(covar$Sex)
}

print("Running Regression ...")

M <- M[rowSums(is.na(M)) != ncol(M),,drop=F] # remove rows where all NA

M <- apply(M, 1, function(xv) { xv[is.na(xv)] <- mean(xv, na.rm=TRUE); return(xv) }) # NAs -> row mean

fit <- lm(M ~ ., data=covar)

r <- data.frame(fit$residuals)
colnames(r) <- colnames(M)
r <- t(r) # transpose for MeQTL

write.table(r, file="./parsed_data/probes.tab", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

p <- list.files(path="./parsed_data/", pattern="ALSPAC.*.tab") # loop imputed files
foreach(i=1:length(p)) %dopar% {
  q  <- p[i]
  print(q)

  df <- read.csv(paste("./parsed_data/", q, sep=""), head=T, row.names=1, sep="\t", check.names=F)
  df <- t(df[ALNs,,drop=F])

  write.table(df, file=paste("./parsed_data/", q, sep=""), sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
}
