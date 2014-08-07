## Downgrading to R 3.0.2 version of R to maintain package compatibility
.libPaths('~/R/x86_64-pc-linux-gnu-library/3.0/')
library(GenABEL)
library(ggplot2)
library(mclust)
library(SNPRelate)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

## weird sample ID format for italian samples
removeACO <- function(x) {
  return(gsub('\\+AC0', '', x))
}
## formatting of dosages
removemgm2 <- function(x) {
  return(gsub(' mg\\/m2' ,'', x))
}
## correct inconsistency in Erasmus sample ID encoding
homogenizeSampleIds <- function(x) {
  return(gsub('M0', 'MO', x))
}

## read in phenotype (ANC) and covariate data
anc <- read.table('annotation/ANC_Nadir_FULL_061614.csv', sep='\t', header=T, stringsAsFactors=F)
## clean up the italian sample IDs
anc$Sample <- unlist(lapply(anc$Subject_ID, removeACO) )
anc$Sample <- unlist(lapply(anc$Sample, homogenizeSampleIds) )
covars <- read.table('annotation/Covariates_Vasa_newIDs.csv', head=T, sep='\t', stringsAsFactors=F)
covars$SM_ID <- as.character(covars$SM_ID)
clean.pheno <- merge(anc, covars, by.x="Sample", by.y="Sample.ID")
nomiss.pheno <- clean.pheno[which(!is.na(clean.pheno$ANC_nadir_1000cells.per.mcL)  & !is.na(clean.pheno$SM_ID) ),]

## make factor variables out of the covariates
nomiss.pheno$Study.site <- as.factor(nomiss.pheno$Study.site)
## homogenize the sex encodings
nomiss.pheno$Sex[which(nomiss.pheno$Sex=='Male')] <- '1'
nomiss.pheno$Sex[which(nomiss.pheno$Sex=='M')] <- '1'
nomiss.pheno$Sex[which(nomiss.pheno$Sex=='Female')] <- '0'
nomiss.pheno$Sex[which(nomiss.pheno$Sex=='F')] <- '0'
nomiss.pheno$Sex <- as.factor(nomiss.pheno$Sex)
nomiss.pheno$Race[which(nomiss.pheno$Race == "")] <- NA
## clean up the dose variables
nomiss.pheno$Dose..mg. <- as.numeric(unlist(lapply(list(nomiss.pheno$Dose..mg.), removemgm2)))

## write a genabel phenotype file
pheno.df <- data.frame(id=nomiss.pheno$SM_ID,
                       sex=nomiss.pheno$Sex,
                       ANC=nomiss.pheno$ANC_nadir_1000cells.per.mcL,
                       site=nomiss.pheno$Study.site,
                       race=nomiss.pheno$Race,
                       dose=nomiss.pheno$Dose..mg.)

## run PCA and cluster samples
snpgdsBED2GDS(bed.fn='exome_data/innocenti_082613.bed',
              bim.fn='exome_data/innocenti_082613.bim',
              fam.fn='exome_data/innocenti_082613.fam',
              out.gdsfn='exome_data/innocenti_082613.gds')
exome.geno <- openfn.gds(filename='exome_data/innocenti_082613.gds')

## prune for LD
ld.snpset <- snpgdsLDpruning(exome.geno)

## run PCA
pca <- snpgdsPCA(exome.geno, snp.id = unlist(ld.snpset), maf = 0.05)
pca.df <- data.frame(sample.id=pca$sample.id, eigenvect=pca$eigenvect)

## pull sample covariates to identify sources of variation
pheno.df <- merge(pca.df, pheno.df, by.x = "sample.id", by.y = "id")
names(pheno.df)[which(names(pheno.df) == "sample.id")] <- "id"

## run clustering
cluster <- Mclust(pheno.df[,c("eigenvect.1", "eigenvect.2")])
pheno.df$cluster <- cluster$classification

## keep a subset of samples in a PLINK dataset; external system command
plink.subset <- function(originalBase, subsetBase, samples.fn) {
  subsetCmd <- paste('plink',
               '--bfile', originalBase,
               '--keep', samples.fn,
               '--recode',
               '--out', subsetBase)
  system(subsetCmd) 
}

## subset original PLINK files to only samples with non-missing phenotype data
nonmissingSamples.fn <- 'annotation/nonmissing_phenotype_samples.txt'
write.table( data.frame(fid=pheno.df$id, iid=pheno.df$id),
              nonmissingSamples.fn,
              row.names=F,
              col.names=F,
              sep='\t')
exomePlink <- 'exome_data/innocenti_082613'
plink.subset(exomePlink, paste0(exomePlink, '.nomiss'), nonmissingSamples.fn)
seqPlink <- 'qc_report/consensus/consensus'
plink.subset(seqPlink, paste0(seqPlink, '.nomiss'), nonmissingSamples.fn)

## convert PLINK data to GenABEL
genabel.geno <- 'qc_report/consensus/consensus.genabel'
convert.snp.ped(pedfile = paste0(seqPlink, '.nomiss.ped'),
                mapfile = paste0(seqPlink, '.nomiss.map'),
                outfile = genabel.geno)
exome.geno <- 'exome_data/innocenti_082613.genabel'
convert.snp.ped(pedfile = paste0(exomePlink, '.nomiss.ped'),
                mapfile = paste0(exomePlink, '.nomiss.map'),
                outfile = exome.geno)

## load and clean data
genabel.pheno <- "annotation/genable.pheno"
write.table(pheno.df, file=genabel.pheno, row.names=F)
consensus <- load.gwaa.data(phenofile = genabel.pheno, genofile = genabel.geno)
exome <- load.gwaa.data(phenofile = genabel.pheno, genofile = exome.geno)
## subset markers by: minor allele frequency and callrate
## subset samples by: PCA cluster (i.e. genetic ethnicity)
seq.snpsubset <- check.marker(data = consensus,
                  maf = 0.05,
                  callrate = 0.95,
                  idsubset = as.character(pheno.df$id[which(pheno.df$cluster==1 | pheno.df$cluster==3)]))
seq.geno <- consensus[seq.snpsubset$idok, seq.snpsubset$snpok]
exome.snpsubset <- check.marker(data = exome,
                    maf = 0.05,
                    callrate = 0.95,
                    idsubset = as.character(pheno.df$id[which(pheno.df$cluster==1 | pheno.df$cluster==3)]))
exome.geno <- exome[exome.snpsubset$idok, exome.snpsubset$snpok]

## run a few scans
basic.model <- ANC ~ 1
full.model <- ANC ~ sex + site + dose
#score.results <- qtscore(scale(ANC) ~ sex + site + dose, seq.geno, trait="gaussian")
seq.reg.results <- mlreg(basic.model, seq.geno, trait="gaussian")
exome.reg.results <- mlreg(basic.model, exome.geno, trait="gaussian")

top.exome <- head(exome.reg.results[order(exome.reg.results[,"P1df"]),])

## summarize
source('scripts/qqunif.r')
qqunif(reg.results[,"P1df"])

