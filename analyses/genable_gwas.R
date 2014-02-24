library(GenABEL)
library(ggplot2)
library(mclust)
library(SNPRelate)

## weird sample ID format for italian samples
removeACO <- function(x) {
  return(gsub('\\+AC0', '', x))
}
## formatting of dosages
removemgm2 <- function(x) {
  return(gsub(' mg\\/m2' ,'', x))
}

## read in phenotype and covariate data
pheno <- read.table('annotation/ANC_Nadir_vasyaModified.csv', sep='\t', header=T, stringsAsFactors=F)
## clean up the italian sample IDs
pheno$Sample <- unlist(lapply(pheno$Sample, removeACO) )
pheno2 <- read.table('annotation/Covariates_Vasa_newIDs.csv', head=T, sep='\t', stringsAsFactors=F)
pheno2$SM_ID <- as.character(pheno2$SM_ID)
clean.pheno <- merge(pheno, pheno2, by.x="Sample", by.y="Sample.ID")
nomiss.pheno <- clean.pheno[which(!is.na(clean.pheno$ANC.Nadir..ul.AC0.1.)  & !is.na(clean.pheno$SM_ID) ),]

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
pheno.df <- data.frame(id=nomiss.pheno$SM_ID, sex=nomiss.pheno$Sex, ANC=nomiss.pheno$ANC.Nadir..ul.AC0.1., site=nomiss.pheno$Study.site, race=nomiss.pheno$Race, dose=nomiss.pheno$Dose..mg. )

## run PCA and cluster samples
snpgdsBED2GDS(bed.fn='exome_data/innocenti_082613.bed',
              bim.fn='exome_data/innocenti_082613.bim',
              fam.fn='exome_data/innocenti_082613.fam',
              'exome_data/innocenti_082613.gds')
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

## convert PLINK data to GenABEL
genabel.geno <- 'qc_report/consensus/consensus.genabel'
convert.snp.ped(pedfile='qc_report/consensus/consensus.nomiss.ped',
                mapfile='qc_report/consensus/consensus.nomiss.map',
                outfile = genabel.geno)
exome.geno <- 'exome_data/innocenti_082613.genabel'
convert.snp.ped(pedfile='exome_data/innocenti_082613.nomisspheno.ped',
                mapfile='exome_data/innocenti_082613.nomisspheno.map',
                outfile = exome.geno)

## load and clean data
genabel.pheno <- "annotation/genable.pheno"
write.table(pheno.df, file=genabel.pheno, row.names=F)
consensus <- load.gwaa.data(phenofile = genabel.pheno, genofile = genabel.geno)
exome <- load.gwaa.data(phenofile = genabel.pheno, genofile = exome.geno)
## subset markers by: minor allele frequency and callrate
## subset samples by: PCA cluster (i.e. genetic ethnicity)
seq.snpsubset <- check.marker(data = consensus, maf = 0.05, callrate = 0.95, idsubset = as.character(pheno.df$id[which(pheno.df$cluster==1 & pheno.df$ANC>1)]))
seq.geno <- consensus[seq.snpsubset$idok, seq.snpsubset$snpok]
exome.snpsubset <- check.marker(data = exome, maf = 0.05, callrate = 0.95, idsubset = as.character(pheno.df$id[which(pheno.df$cluster==1 & pheno.df$ANC>1)]))
exome.geno <- exome[exome.snpsubset$idok, exome.snpsubset$snpok]

## run a few scans
basic.model <- scale(ANC) ~ 1
full.model <- scale(ANC) ~ sex + site + dose
#score.results <- qtscore(scale(ANC) ~ sex + site + dose, seq.geno, trait="gaussian")
seq.reg.results <- mlreg(full.model, seq.geno, trait="gaussian")
exome.reg.results <- mlreg(full.model, exome.geno, trait="gaussian")

## summarize
source('scripts/qqunif.r')
qqunif(reg.results[,"P1df"])

