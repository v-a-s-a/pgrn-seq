## Check out the exome chip data PCA

library(ggplot2)
library(SNPRelate)
library(mclust)

## import data into GDS
snpgdsBED2GDS(bed.fn='../exome_data/innocenti_082613.bed', bim.fn='../exome_data/innocenti_082613.bim', fam.fn='../exome_data/innocenti_082613.fam', '../exome_data/innocenti_082613.gds')
genofile <- openfn.gds(filename='../exome_data/innocenti_082613.gds')

## perform PCA
ld.snpset <- snpgdsLDpruning(genofile)
pca <- snpgdsPCA(genofile, snp.id = unlist(ld.snpset), num.thread = 12, maf = 0.05, eigen.cnt=3)
pca.df <- data.frame(sample.id=pca$sample.id, eigenvect=pca$eigenvect)

## cluster samples in PCA space to pull europeans
cluster <- Mclust(data=pca.df[c("eigenvect.1", "eigenvect.2")])
pca.df$cluster <- as.factor(cluster$classification)

## annotate with self reported ancestry
covar <- read.table('../annotation/Covariates_Vasa_newIDs.csv', header=TRUE, sep='\t')
pca <- merge(pca.df, covar, by.x = "sample.id", by.y = "SM_ID")


print(ggplot(pca, aes(x=eigenvect.1, eigenvect.2, color=Race)) + geom_point())
