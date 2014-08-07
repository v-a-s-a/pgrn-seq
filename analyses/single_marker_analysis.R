# simple and multiple regression of the phenotype

##############################################################################
## PACKAGES ##################################################################
##############################################################################
## Downgrading to R 3.0.2 version of R to maintain package compatibility
.libPaths('~/R/x86_64-pc-linux-gnu-library/3.0/')
library(GenABEL)
library(ggplot2)

##############################################################################
## CONSTANTS #################################################################
##############################################################################
PHENO.FN <- 'data/rdata/pheno.Rdata'
EXOME.GENO.FN <- 'data/rdata/exome.geno.Rdata'
SEQ.GENO.FN <- 'data/rdata/seq.geno.Rdata'

##############################################################################
## FUNCTIONS #################################################################
##############################################################################

##############################################################################
## ANALYSIS ##################################################################
##############################################################################

## load genotype data
if(!file.access(EXOME.GENO.FN)==0 | !file.access(SEQ.GENO.FN)==0) {
  ## genotype data has not yet been processed
  source('analyses/plink_genotype_preprocessing.R')
} else {
  ## load processed genotype data
  load(EXOME.GENO.FN)
  load(SEQ.GENO.FN)
}

## run scans
basic.model <- ANC ~ 1
full.model <- ANC ~ sex + site + dose
#score.results <- qtscore(scale(ANC) ~ sex + site + dose, seq.geno, trait="gaussian")
seq.reg.results <- mlreg(basic.model, seq.geno, trait="gaussian")
exome.reg.results <- mlreg(basic.model, exome.geno, trait="gaussian")

top.exome <- head(exome.reg.results[order(exome.reg.results[,"P1df"]),])

## summarize
source('scripts/qqunif.r')
qqunif(seq.reg.results[,"P1df"])

