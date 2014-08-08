## tidying genotype and phenotype data

##############################################################################
## PACKAGES ##################################################################
##############################################################################
.libPaths('~/R/x86_64-pc-linux-gnu-library/3.0/')
library(ggplot2)
library(mclust)
library(SNPRelate)

##############################################################################
## CONSTANTS #################################################################
##############################################################################
PHENO.FN <- 'data/rdata/pheno.Rdata' # dataframe of processed phenotype data
exome.fbase <- 'data/exomechip_data/innocenti_082613'

##############################################################################
## FUNCTIONS #################################################################
##############################################################################

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

##############################################################################
## PHENOTYPE #################################################################
##############################################################################

## read in phenotype (ANC) and covariate data
anc <- read.table('data/annotation/ANC_Nadir_FULL_061614.csv', sep='\t', header=T, stringsAsFactors=F)
## clean up the italian sample IDs
anc$Sample <- unlist(lapply(anc$Subject_ID, removeACO) )
anc$Sample <- unlist(lapply(anc$Sample, homogenizeSampleIds) )
covars <- read.table('data/annotation/Covariates_Vasa_newIDs.csv', head=T, sep='\t', stringsAsFactors=F)
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

## cleanup SNPRelate files if they exist
if(!file.access(paste0(exome.fbase, '.gds'))) unlink(paste0(exome.fbase, '.gds'))

## run PCA and cluster samples
snpgdsBED2GDS(bed.fn=paste0(exome.fbase, '.bed'),
              bim.fn=paste0(exome.fbase, '.bim'),
              fam.fn=paste0(exome.fbase, '.fam'),
              out.gdsfn=paste0(exome.fbase, '.gds'))
exome.geno <- openfn.gds(filename=paste0(exome.fbase, '.gds'))

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

## dump phenotype data to disk
save(pheno.df, file=PHENO.FN)
