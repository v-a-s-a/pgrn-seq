## graph site and sample concordance 

library(ggplot2)

site.dat <- read.table('../qc_report/consensus/consensus.site.concord.txt')
sample.dat <- read.table('../qc_report/consensus/consensus.sample.concord.txt')


sample.plt <- ggplot(sample.dat, aes(x=as.numeric(V2))) +
                geom_histogram(binwidth=0.0035) +
                theme_bw() +
                xlab("# of concordant sites / total # of sites") +
                ylab("# of samples") +
                labs(title="Variant Site Concordance Rate per Sample")
show(sample.plt)

x11()
site.plt <- ggplot(site.dat, aes(x=as.numeric(V2))) +
                geom_histogram(binwidth=0.015) +
                theme_bw() +
                xlab("# of concordant genotypes / total possible genotypes") +
                ylab("# of sites") +
                labs(title="Genotype Concordance Rate per Variant Site")
show(site.plt)


