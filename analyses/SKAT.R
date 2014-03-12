## first pass at a gene-based test using the SKAT package

library(VariantAnnotation)
library(SNPRelate)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(SKAT)
library(doMC)
library(foreach)

source('scripts/qqunif.r')

## load in genotype data
vcf <- readVcf('variant_data/consensus.vcf.gz', 'hg19')
vcfRows <- rowData(vcf)

## annotate variants with Entrez gene IDs
loc <- locateVariants(vcfRows, TxDb.Hsapiens.UCSC.hg19.knownGene, AllVariants())
names(loc) <- NULL
annot <- as.data.frame(loc)
annot$names <- names(vcfRows)[ annot$QUERYID ]
annot <- annot[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID")]

## pull the genotype matrix
XsnpMat <- genotypeToSnpMatrix(vcf)
Xraw <- as(object=XsnpMat$genotypes, Class="numeric")

## remove snps with no variation and high missingness


## load up the phenotype data
## NOTE: This phenotype file was created for the GenABEL GWAS. It is not raw!
##        It has removed some missing data, and standardized some sample names.
pheno.df <- read.table('annotation/genable.pheno', header=T)
pheno.df$dose <- as.numeric(pheno.df$dose)
pheno.df$sex <- as.factor(pheno.df$sex)
## subset to european-like individuals
eur.pheno <- pheno.df[which(pheno.df$cluster == 1),]

## match phenotype and genotype matrices
X <- Xraw[as.character(eur.pheno$id),]

## Convert to a sparse matrix by reporting the dosage of the minor allele.
## The snpMatrix class defaults to reporting dosage of the alternate allele.
flip.allele <- function(x) {
  if (all(is.na(x))) {
    return(x)
  } else if (mean(x, na.rm=TRUE) > 1) {
  ## alleles need to be flipped
  ## 0 -> 2
  ## 1 -> 1
  ## 2 -> 0
  x[x == 0] <- 3
  x[x == 2] <- 0
  x[x == 3] <- 2
  }
  return(x)
}
sparseX <- apply(X, 2, flip.allele)


## run SKAT test for a gene
genes <- unique(annot$GENEID[which(!is.na(annot$GENEID))])
attach(eur.pheno)
null.model <- SKAT_Null_Model(ANC ~ sex + site + dose, out_type="C", Adjustment=FALSE)
#res <- data.frame(gene = c(),
#                  burden.asymp.p = c(),
#                  burden.resample.p = c(),
#                  rarecommon.asymp.p = c(),
#                  rarecommon.resample.p = c())
#for (i in seq_along(genes)) {
#  print(i)
#  gene <- genes[i]
#  gene.snps <- annot$names[which(annot$GENEID == gene)]
#  Xgene <- sparseX[,gene.snps]
#  try({
#    burden <- SKAT(Xgene, null.model)
#    rarecommon <- SKAT_CommonRare(Xgene, null.model)
#    res <- rbind(res, data.frame(gene=gene,
#                                 burden.asymp.p = burden$p.value,
#                                 burden.resample.p = Get_Resampling_Pvalue(burden)$p.value
#                                 rarecommon.asymp.p = rarecommon$p.value,
#                                 rarecommon.resample.p = Get_Resampling_Pvalue(rarecommon)$p.value)
#    })
#}

registerDoMC(cores=8)
res <- foreach(i = 1:length(genes), .combine = rbind) %dopar% {
  print(i)
  gene <- genes[i]
  gene.snps <- annot$names[which(annot$GENEID == gene)]
  Xgene <- sparseX[,gene.snps]
  try({
    burden <- SKAT(Xgene, null.model)
    data.frame(gene=gene,
               burden.asymp.p = burden$p.value)
    })
}
#res[,2] <- lapply(res[,2], as.numeric)

## calculate some permutation based p-values for top hits
permute.p <- function() {


}

## add gene symbols for interpretation
Symbol2id <- as.list( org.Hs.egSYMBOL2EG )
id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
names(id2Symbol) <- unlist(Symbol2id)
res$symbol <- id2Symbol[ as.character(res$gene) ]
