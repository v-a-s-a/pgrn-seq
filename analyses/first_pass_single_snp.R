library(SNPRelate)
library(ggplot2)


## weird sample ID format for italian samples
removeACO <- function(x) {
  return(gsub('\\+AC0', '', x))
}


## read in phenotype and covariate data
pheno <- read.table('annotation/ANC_Nadir_vasyaModified.csv', sep='\t', header=T, stringsAsFactors=F)
## clean up the italian sample IDs
pheno$Sample <- unlist(lapply(pheno$Sample, removeACO) )
pheno2 <- read.table('annotation/Covariates_Vasa_newIDs.csv', head=T, sep='\t', stringsAsFactors=F)
pheno2$SM_ID <- as.character(pheno2$SM_ID)
clean.pheno <- merge(pheno, pheno2, by.x="Sample", by.y="Sample.ID")
nomiss.pheno <- clean.pheno[which(!is.na(clean.pheno$ANC.Nadir..ul.AC0.1.) & !is.na(clean.pheno$SM_ID) ),]
nomiss.pheno <- nomiss.pheno[which(nomiss.pheno$ANC.Nadir..ul.AC0.1. != 0),]

## make factor variables out of the covariates
nomiss.pheno$Study.site <- as.factor(nomiss.pheno$Study.site)
## homogenize the sex encodings
nomiss.pheno$Sex[which(nomiss.pheno$Sex=='Male')] <- 'M'
nomiss.pheno$Sex[which(nomiss.pheno$Sex=='Female')] <- 'F'
nomiss.pheno$Sex <- as.factor(nomiss.pheno$Sex)
nomiss.pheno$Race[which(nomiss.pheno$Race == "")] <- NA


## read in genotype data
#snpgdsVCF2GDS(vcf.fn='../variant_data/consensus.vcf',outfn.gds='../variant_data/consensus.gds')
snpgdsBED2GDS(bed.fn='qc_report/consensus/consensus.bed', bim.fn='qc_report/consensus/consensus.bim', fam.fn='qc_report/consensus/consensus.fam', 'variant_data/consensus.gds')
genofile <- openfn.gds(filename='variant_data/consensus.gds')

## run PCA for use as covariates
ld.snpset <- snpgdsLDpruning(genofile)
pca <- snpgdsPCA(genofile, snp.id = unlist(ld.snpset), num.thread = 12, maf = 0.05, eigen.cnt=3)
pca.df <- data.frame(sample.id=pca$sample.id, eigenvect=pca$eigenvect)

tmp <- merge(pca.df, nomiss.pheno, by.x = "sample.id", by.y = "SM_ID")
nomiss.pheno <- tmp


## subset to reported ethnicity
nomiss.pheno <- nomiss.pheno[which(nomiss.pheno$Race == "Caucasian"),]

## remove outlier phenotypes
nomiss.pheno <- nomiss.pheno[which(log10(nomiss.pheno$ANC.Nadir..ul.AC0.1.) > 1.5),]

## standardize phenotype
nomiss.pheno$ANC.Nadir..ul.AC0.1. <- scale(nomiss.pheno$ANC.Nadir..ul.AC0.1., center=TRUE, scale=TRUE)


## remove extremely rare variants
snpset <- snpgdsSelectSNP(genofile, maf = 0.01, missing.rate=0.05)
X <- snpgdsGetGeno(genofile, snp.id = snpset)
colnames(X) <- snpset 
rownames(X) <- read.gdsn(index.gdsn(genofile, "sample.id"))
## scale to colMean=0 and colSd=1
X <- scale(X, center = TRUE, scale = TRUE)
## set missing values to new column mean, i.e. 0.0
X[is.na(X)] <- 0.0
nomiss.X <- X[which(rownames(X) %in% nomiss.pheno$sample.id),]


## order the phenotypes
ordered.pheno <- nomiss.pheno[match(nomiss.pheno$sample.id, rownames(nomiss.X)),]

refcovar <- function(index) {
  return(paste0('covar[[', index, ']]'))
}

extendmodel <- function(model, value) {
  return( paste(model, '+', value) )
}



## incorporate covariate information
minimal.lm <- function(X, Y, covar=NULL) {

  ## construct full model expression
  model <- paste('Y', '~', 'X')
  
  if (!is.null(covar)) {
    model.covars <- lapply(X=list(1:length(covar)), FUN=refcovar)
    model <- Reduce(extendmodel, unlist(model.covars), init=model)
  }

  ## fit model and extraxt pvalue
  lm.fit <- lm(as.formula(model))
  return(lm.fit)
}


## return the pcalue of a model fit using anova
lmp.anova <- function(modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  afit <- anova(modelobject)
  p <- afit[["Pr(>F)"]][1]
  return(p)
}

## return a pvalue from a fit
##   This function may depend on how the model was specified
lmp.fit <- function(lm.fit) {
  return(summary(lm.fit)$coefficients[2,4])
}


## specify the model and apply it over the SNP genotypes
#model.covar <- data.frame(sex=ordered.pheno$Sex, site=ordered.pheno$Study.site, ev1=ordered.pheno$eigenvect.1, ev2=ordered.pheno$eigenvect.2)
model.covar2 <- data.frame(sex=ordered.pheno$Sex, site=ordered.pheno$Study.site)
fits <- apply(nomiss.X, 2, minimal.lm, Y = ordered.pheno$ANC.Nadir..ul.AC0.1., covar=model.covar2 )
#raw.fits <- apply(nomiss.X, 2, minimal.lm, Y = ordered.pheno$ANC.Nadir..ul.AC0.1.)
pvals <- unlist(lapply(fits, lmp.fit))
resids <- unlist(lapply(fits, resid))

source('scripts/qqunif.r')
names(pvals) <- colnames(nomiss.X)
qqunif(pvals)

