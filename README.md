# ALSPAC meQTL

This package provides the framework used within ALSPAC/ARIES for running meQTL analysis

## Dependencies

Python 2.7+ (do not use Python 3)

PLINK v1.9 needs to be on your path

The following R Packages are also required - these can be installed using the ```install.packages``` command during an interactive R session:

- argparse
- foreach
- doMC
- parallel
- MatrixEQTL

## Using ALSPAC_meQTL

To use the package, just type

```
python run.py
```

## Arguments

You can see the list of accepted arguments using the ```--help``` parameter, i.e.

```
python run.py --help
```

### Required Arguments

```--geno```: path to the binary ped file, i.e. PLINK "bfile"

```--meth```: path to the mathylation data **_in Rdata_** format

```--covars```: path to the covariates used to residualize the methylation data.  

**Note:** covariate data must be a tab-delimited file with samples as rows and covariates as columns.  Also, "Sex" and "Batch" **_MUST_** be covariates (_proprocessing.R_ expects them and will fail if they are omitted)

```--tp```: a file listing the sample IDs in the ARIES time point to analyse

```--fo```: where the output is written

### Optional Arguments

```--keep-snps```: path to a file containing the list of SNPs to analyze (defaults to all SNPs)

```--keep-probes```: path to a file containing a list of CpG probes to analyze (defaults to all probes)

```--pv```: the p-value threshold used when reporting results (defaults to 1e-5)

```--distance```: the distance used when flagging associations as _trans_ (defaults to 1Mb)

### Debugging Arguments

These arguments do not affect the results, but do provide developers with some debugging options:

```--batch```: the batch size to use (defaults to 10,000)

```--skip-annotation```: skip the final annotation step to save some time

```--keep-intermediate-data```: do not delete intermediate data files when finished

## Examples

- run meQTL analysis in 15up for **_all SNPs/probe combination_**:

```
python run.py --geno <path_to_geno> --meth <path_to_meth> --covars <path_to_covars> --tp 15up --fo <path_to_output>
```

- run meQTL analysis in FOM for **_all SNPs, but a handful of probes_**:

```
python run.py --geno <path_to_geno> --meth <path_to_meth> --covars <path_to_covars> --tp FOM --fo <path_to_output> --keep-probes <path_to_probes>
```

- run meQTL analysis in F7 for **_all probes, but a handful of SNPs_**:

```
python run.py --geno <path_to_geno> --meth <path_to_meth> --covars <path_to_covars> --tp F7 --fo <path_to_output> --keep-snps <path_to_snps>
```
