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
ld.snpset <- snpgdsLDpruning(geno, ld.threshold=0.2)

## run IBD analysis
ibd <- snpgdsIBDMLE(geno, snp.id=unlist(ld.snpset), maf=0.05, missing.rate=0.05, num.thread=12)

## pull out IBD coefficients
ibd.coeff <- snpgdsIBDSelection(ibd)

## plot
plt <- ggplot(ibd.coeff, aes(x=k0, y=k1)) + geom_point()
show(plt)



