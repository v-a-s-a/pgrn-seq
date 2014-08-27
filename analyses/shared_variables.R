### Keep shared variables in one place

## Dataframe of clean phenotype data
PHENO.FN <- '../data/rdata/pheno.Rdata'

## GenABEL genotype data objects
### FULL: all samples included
### EURO: only samples in european cluster of PC plot
EURO.EXOME.GENO.FN <- '../data/rdata/euro.exome.geno.Rdata'
FULL.EXOME.GENO.FN <- '../data/rdata/full.exome.geno.Rdata'
FULL.SEQ.GENO.FN <- '../data/rdata/full.seq.geno.Rdata'
EURO.SEQ.GENO.FN <- '../data/rdata/euro.seq.geno.Rdata'

## Models
basic.model <- trans.fun(ANC) ~ 1
full.model <- trans.fun(ANC) ~ sex + site + as.numeric(dose)
pca.model <- trans.fun(ANC) ~ sex + site + as.numeric(dose) + eigenvect.1 + eigenvect.2 + eigenvect.3 + eigenvect.4 + eigenvect.5

