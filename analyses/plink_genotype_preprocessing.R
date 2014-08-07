## process PLINK files from exome chip and sequence data for single marker analyses

##############################################################################
## PACKAGES ##################################################################
##############################################################################
.libPaths('~/R/x86_64-pc-linux-gnu-library/3.0/')
library(GenABEL)
library(ggplot2)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

##############################################################################
## CONSTANTS #################################################################
##############################################################################
PHENO.FN <- 'data/rdata/pheno.Rdata'
EXOME.GENO.FN <- 'data/rdata/exome.geno.Rdata'
SEQ.GENO.FN <- 'data/rdata/seq.geno.Rdata'

##############################################################################
## FUNCTIONS #################################################################
##############################################################################

## keep a subset of samples in a PLINK dataset; external system command
plink.subset <- function(originalBase, subsetBase, samples.fn) {
  subsetCmd <- paste('plink',
               '--bfile', originalBase,
               '--keep', samples.fn,
               '--recode',
               '--out', subsetBase)
  system(subsetCmd)
}

##############################################################################
## GENOTYPE ##################################################################
##############################################################################

## load and process the phenotype information
if(!file.access(PHENO.FN)==0) {
  ## phenotype data has not yet been generated -- do that now
  source('analyses/phenotype_preprocessing.R')
} else {
  ## load alrady processed phenotype file
  load(PHENO.FN)
}

## subset original PLINK files to only samples with non-missing phenotype data
nonmissingSamples.fn <- 'data/annotation/nonmissing_phenotype_samples.txt'
write.table( data.frame(fid=pheno.df$id, iid=pheno.df$id),
              nonmissingSamples.fn,
              row.names=F,
              col.names=F,
              sep='\t')

## convert sequence PLINK data to GenABEL
seqPlink <- 'data/qc_report/consensus/consensus'
plink.subset(seqPlink, paste0(seqPlink, '.nomiss'), nonmissingSamples.fn)
genabel.geno <- 'data/qc_report/consensus/consensus.genabel'
convert.snp.ped(pedfile = paste0(seqPlink, '.nomiss.ped'),
                mapfile = paste0(seqPlink, '.nomiss.map'),
                outfile = genabel.geno)

## convert exome PLINK data to GenABEL
exomePlink <- 'data/exomechip_data/innocenti_082613'
plink.subset(exomePlink, paste0(exomePlink, '.nomiss'), nonmissingSamples.fn)
exome.geno <- 'data/exomechip_data/innocenti_082613.genabel'
convert.snp.ped(pedfile = paste0(exomePlink, '.nomiss.ped'),
                mapfile = paste0(exomePlink, '.nomiss.map'),
                outfile = exome.geno)

## load and clean data
## subset markers by: minor allele frequency and callrate
## subset samples by: PCA cluster (i.e. genetic ethnicity)
genabel.pheno <- "data/annotation/genable.pheno"
write.table(pheno.df, file=genabel.pheno, row.names=F)
## Seqeunce data
consensus <- load.gwaa.data(phenofile = genabel.pheno, genofile = genabel.geno)
seq.snpsubset <- check.marker(data = consensus,
                  maf = 0.05,
                  callrate = 0.95,
                  idsubset = as.character(pheno.df$id[which(pheno.df$cluster==1 | pheno.df$cluster==3)]))
seq.geno <- consensus[seq.snpsubset$idok, seq.snpsubset$snpok]
## Exome data
exome <- load.gwaa.data(phenofile = genabel.pheno, genofile = exome.geno)
exome.snpsubset <- check.marker(data = exome,
                    maf = 0.05,
                    callrate = 0.95,
                    idsubset = as.character(pheno.df$id[which(pheno.df$cluster==1 | pheno.df$cluster==3)]))
exome.geno <- exome[exome.snpsubset$idok, exome.snpsubset$snpok]

## dump data to disk
save(exome.geno, file=EXOME.GENO.FN)
save(seq.geno, file=EXOME.GENO.FN)

