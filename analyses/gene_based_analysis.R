## gene-based scan using the SKAT package

##############################################################################
## PACKAGES ##################################################################
##############################################################################
.libPaths('~/R/x86_64-pc-linux-gnu-library/3.0/')
library(VariantAnnotation, quietly = T)
library(SNPRelate, quietly = T)
library(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = T)
library(org.Hs.eg.db, quietly = T)
library(SKAT, quietly = T)
library(doMC, quietly = T)
library(foreach, quietly = T)

source('scripts/qqunif.r')

##############################################################################
## CONSTANTS #################################################################
##############################################################################
PHENO.FN <- 'data/rdata/pheno.Rdata'

##############################################################################
## FUNCTIONS #################################################################
##############################################################################

## Convert to a sparse matrix by reporting the dosage of the minor allele.
## The snpMatrix class defaults to reporting dosage of the alternate allele.
flip.allele <- function(x) {
  if (all(is.na(x))) {
    return(x)
  } else if (mean(x, na.rm=TRUE) > 1) {
  ## alleles need to be flipped
  x[x == 0] <- 3
  x[x == 2] <- 0
  x[x == 3] <- 2
  }
  return(x)
}

## parallel SKAT scan
parallel.skat.scan <- function(genes, pheno, model, snpMat, annot, n.core = 1) {
  ## compute the SKAT null model
  null.model <- with(eur.pheno, SKAT_Null_Model(as.formula(model), out_type="C", Adjustment=TRUE))
  registerDoMC(cores=n.core)
  res <- foreach(i = 1:length(genes), .combine = rbind) %dopar% {
    gene <- genes[i]
    gene.snps <- annot$names[which(annot$GENEID == gene)]
    Xgene <- sparseX[,gene.snps]
    ## record mean MAF
    try({skat <- SKAT(Xgene, null.model)
         data.frame(gene=gene,
                  p = as.numeric(skat$p.value),
                  Q = as.numeric(skat$Q),
                  nsnps = as.numeric(skat$param$n.marker.test))
      })
  }
  ## cast types into numerical results
  res$p <- as.numeric(res$p)
  res$Q <- as.numeric(res$Q)
  res$nsnps <- as.numeric(res$nsnps)
  clean.res <- res[which(!is.na(res$p)),]

  ## add gene symbols for interpretation
  Symbol2id <- as.list( org.Hs.egSYMBOL2EG )
  id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
  names(id2Symbol) <- unlist(Symbol2id)
  clean.res$symbol <- id2Symbol[ as.character(clean.res$gene) ]

  return(clean.res)
}

## calculate some permutation based p-values for top hits
## should also try generating one large null distribution, and testing genes against that
permute.p <- function(pheno.df, X, n, Q) {
  null.q <- rep(0, n)
  #null.skat <- SKAT_Null_Model(ANC ~ sex + site + dose, out_type="C", Adjustment=F)
  null.skat <- SKAT_Null_Model(ANC ~ 1, out_type="C", Adjustment=F)
  for(i in 1:n) {
    print(i)
    ## shuffle phenotype and covariate values
    idx <- sample(1:nrow(pheno.df))
    ## fit null model
    skat <- SKAT(X[idx,], null.skat)
    null.q[i] <- skat$Q
  }
  return( 1 - (length(which(Q > null.q)) / n) )
}

## annotate a snp matrix
annotate.genotypes <- function(vcfFile) {
  ## load in genotype data
  vcf <- readVcf(vcfFile, 'hg19')
  vcfRows <- rowData(vcf)
  ## annotate variants with Entrez gene IDs
  loc <- locateVariants(vcfRows, TxDb.Hsapiens.UCSC.hg19.knownGene, AllVariants())
  names(loc) <- NULL
  annot <- as.data.frame(loc)
  annot$names <- names(vcfRows)[ annot$QUERYID ]
  return(annot[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID")])
}

## make a sparse genotype matrix from a vcf file
load.snp.mat <- function(vcf, pheno) {
  ## pull the genotype matrix
  vcf <- readVcf(vcfFile, 'hg19')
  vcfRows <- rowData(vcf)
  XsnpMat <- genotypeToSnpMatrix(vcf)
  Xraw <- as(object=XsnpMat$genotypes, Class="numeric")
  ## match phenotype and genotype matrices
  X <- Xraw[as.character(pheno$id),]
  sparseX <- apply(X, 2, flip.allele)
  return(sparseX)
}

##############################################################################
## ANALYSIS ##################################################################
##############################################################################

## load phenotype data
if(!file.access(PHENO.FN)==0) {
  ## phenotype data has not yet been generated -- do that now
  source('analyses/phenotype_preprocessing.R')
} else {
  ## load alrady processed phenotype file
  load(PHENO.FN)
}
## subset to european-like individuals
eur.pheno <- pheno.df[which(pheno.df$cluster == 1 | pheno.df$cluster == 3),]


######## PGRNseq Analysis ###################
## load in genotype data
vcfFile <- 'consensus_seq_variants/consensus.geno.vcf.gz'
annot <- annotate.genotypes(vcfFile)
sparseX <- load.snp.mat(vcfFile, eur.pheno)

## identify the unique genes we have captured here
genes <- unique(annot$GENEID[which(!is.na(annot$GENEID))])

## run the scan for the PGRNseq data
clean.res <- parallel.skat.scan(genes = genes,
              pheno = eur.pheno,
              model = "ANC ~ 1",
              snpMat = sparseX,
              annot = annot)
top.hits <- clean.res[order(clean.res$p),]
gene.X <- sparseX[,annot$names[which(annot$GENEID == top.hits[1, "gene"])]]

## Some initial permutation tests
#permute.p(pheno.df = eur.pheno, X = gene.X, n = 10000, Q = top.hits[1, "Q"])


########## Exome chip analysis ###############
#exomeFile <- 'exome_data/innocenti_082613_clean_variant_chr.vcf.gz'
#exomeannot <- annotate.genotypes(exomeFile)
#exomegenes <- unique(exomeannot$GENEID[which(!is.na(exomeannot$GENEID))])
#sparseX <- load.snp.mat(exomeFile, eur.pheno)
#exome.res <-  parallel.skat.scan(genes = exomegenes, pheno = eur.pheno, model = "ANC ~ 1", snpMat = sparseX, annot = exomeannot)
