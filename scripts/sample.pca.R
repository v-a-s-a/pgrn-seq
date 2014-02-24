## sample genotype PCA

library(gdsfmt)
library(SNPRelate)
library(ggplot2)

## load data
bed <- 'qc_report/consensus/consensus.bed'
bim <- 'qc_report/consensus/consensus.bim' 
fam <- 'qc_report/consensus/consensus.fam'
gds <- 'qc_report/consensus/consensus.gds'
snpgdsBED2GDS(bed.fn = bed, bim.fn = bim, fam.fn= fam, out.gdsfn=gds )

## prune for LD
geno <-  openfn.gds(gds)
ld.snpset <- snpgdsLDpruning(geno)

## run PCA
pca <- snpgdsPCA(geno, snp.id = unlist(ld.snpset), maf = 0.05)
pca.df <- data.frame(sample.id=pca$sample.id, eigenvect=pca$eigenvect)

## pull sample covariates to identify sources of variation
covar <- read.table('annotation/Covariates_Vasa_newIDs.csv', header=TRUE, sep='\t')

pca <- merge(pca.df, covar, by.x = "sample.id", by.y = "SM_ID")


## plot results
plt <- ggplot(pca, aes(x=eigenvect.1, y=eigenvect.2, colour=Race)) +
        geom_point() +
        scale_color_brewer(palette="Set2") +
        labs(x="PC1", y="PC2")

show(plt)
