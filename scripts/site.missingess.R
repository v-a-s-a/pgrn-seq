## plot mendelian error per sample

library(plyr)
library(ggplot2)

parse_imen <- function(x) {
  read.table(x, header=T)
}


sets <- list('atlas', 'gatk', 'freebayes', 'consensus')
files <- list('../qc_report/atlas/atlas.pgrn.proc.ontarget.minQ10.imiss',
              '../qc_report/gatk/gatk.pgrn.ontarget.minQ10.imiss',
              '../qc_report/freebayes/freebayes.pgrn.ontarget.minQ10.imiss',
              '../qc_report/consensus/consensus.imiss')

## read in tables
dat <- lapply(files, parse_imen)

## drop unnecessary columns
dat <- Map( function(x) x[,!(names(x) %in% c('CHR','N_MISS','N_GENO'))], dat )

## uniquify errors column names
#dat <- Map( function(dataf, newname) rename(dataf, c("F_MISS"=newname)), dat, sets )


## Merge them into a single data frame properly ordered by sample
#dat <- Reduce( function(...) merge(..., by.x=c('SNP'), by.y=c('SNP'), suffixes), dat )

## add caller names to data frame
dat <- Map( function(dataf, name) dataf[,'caller'] <- rep(name, length(dataf[,1])), dat, sets )

## make into one big dataframe for plotting
plt.dat <- ldply(dat, data.frame)

plt <- ggplot(plt.dat, aes(x=values)) + facet_grid( ind ~ . ,scale= "free") 
plt <- plt + geom_histogram(binwidth=0.003)
plt <- plt + xlab('# of sites with this rate of missingness') + ylab('Rate of missing genotypes')
plt <- plt + labs('Variant Site Genotype Missingness Spectra')
plt <- plt + theme_bw()

show(plt)
