library(RColorBrewer)
library(plyr)
library(ggplot2)

parse_imen <- function(x) {
  read.table(x, header=T)
}

## Make a data frame of locus missingness data
make_fpath <- function(base, basedir) {
  fpath <- paste( basedir, base, '.imiss', sep = '' )
  return(fpath)
}

## MAIN
proj.dir <- '/nas40t0/vasya/pgrn_seq/qc_report/'
call.sets = list( 'atlas/atlas.pgrn.proc.ontarget.minQ10',
                  'freebayes/freebayes.pgrn.ontarget.minQ10',
                  'consensus/consensus',
                  'gatk/gatk.pgrn.ontarget.minQ10',
                  'original_calls/innocenti_pgrn_irinotecan_1_multisample_filtered' )

callers <- list('Atlas', 'Freebayes', 'CGES', 'GATK', 'Original')
imiss.fpaths <- mapply(make_fpath, call.sets, rep(proj.dir, length(call.sets)) )
imiss.dat <- lapply( imiss.fpaths, read.table, header=TRUE )


## drop unnecessary columns
dat <- Map( function(x) x[,!(names(x) %in% c('MISS_PHENO','N_MISS','N_GENO'))], imiss.dat )

## uniquify errors column names
dat <- Map( function(dataf, newname) rename(dataf, c("F_MISS"=newname)), dat, callers )

## Merge them into a single data frame properly ordered by sample
dat <- Reduce( function(...) merge(..., by.x=c('FID','IID'), by.y=c('FID', 'IID'), suffixes), dat )

## friendly form for plotting
plt.dat <- stack( dat, select = unlist(callers) )

## keep coloring consistent between plots
myColors <- brewer.pal(5,"Set1")
names(myColors) <- unlist(callers)
colScale <- scale_fill_manual(values = myColors)


plt <- ggplot(plt.dat, aes(x=values, fill=ind)) + facet_grid( ind ~ . ,scale= "free") 
plt <- plt + geom_histogram(binwidth=0.003)
plt <- plt + xlab('# of samples with this rate of missingness') + ylab('rate of missing genotypes')
plt <- plt + labs('Sample Genotype Missingness Spectra')
plt <- plt + theme_bw()
plt <- plt + colScale


show(plt)
