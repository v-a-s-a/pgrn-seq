## plot out variant counts

library(ggplot2)

## Make a data frame of locus missingness data
make_fpath <- function(base, basedir) {
  fpath <- paste( basedir, base, '.lmiss', sep = '' )
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
lmiss.fpaths <- mapply(make_fpath, call.sets, rep(proj.dir, length(call.sets)) )
lmiss.dat <- lapply( lmiss.fpaths, read.table, header=TRUE )

m.var <- data.frame( mvar = unlist( lapply(lmiss.dat, nrow) ),
                    caller = unlist(callers))
m.var$order_callers <- reorder(m.var$caller, m.var$mvar)

plt <- ggplot( m.var, aes(x=order_callers, y=mvar, fill=caller) ) +
        geom_bar(stat="identity") +
        labs(x="Caller", y="Number of Variants") +
        theme_bw() +
        scale_fill_brewer(palette="Set1")

show(plt)






