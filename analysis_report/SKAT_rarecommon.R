## first pass at a gene-based test using the SKAT package

library(VariantAnnotation)
library(SNPRelate)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## load in genotype data
vcf <- readVcf('./consensus.vcf.gz', 'hg19')
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

## load up the phenotype data
## NOTE: This phenotype file was created for the GenABEL GWAS. It is not raw!
##        It has removed some missing data, and standardized some sample names.
pheno.df <- read.table('genable.pheno', header=T)
## subset to european-like individuals
eur.pheno <- pheno.df[which(pheno.df$cluster == 1),]

## match phenotype and genotype matrices
X <- Xraw[as.character(eur.pheno$id),]

## run SKAT test for a gene
pvals <- c()
genes <- unique(annot$GENEID[which(!is.na(annot$GENEID))])
attach(eur.pheno)
null.skat <- SKAT_Null_Model(ANC ~ sex + site + dose, out_type="C", data=eur.pheno)
for (i in seq_along(genes)) {
  gene <- genes[i]
  gene.snps <- annot$names[which(annot$GENEID == gene)]
  Xgene <- X[,gene.snps]
  if (!is.matrix(Xgene)) next
  #pvals <- c(pvals, SKAT(Xgene, obj)$p.value)
  pvals <- c(pvals, SKAT_CommonRare(Xgene, obj)$p.value)
}



## load in phenotype data
## we need a way to annotate genotypes by gene or unit
