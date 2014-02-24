library(ggplot2)

#' Make a data frame of TsTV summary data
make_fpath <- function(base, basedir) {
  fpath <- paste( basedir, base, '.tstv.TsTv.summary', sep = '' )
  return(fpath)
}

#' Calculate the TsTv ratio
extract_tstv <- function(x) {
  return( x$COUNT[7] / x$COUNT[8] )
}

## MAIN

proj.dir <- '/nas40t0/vasya/pgrn_seq/qc_report/'
call.sets = list( 'atlas/atlas.pgrn.proc.ontarget.minQ10',
                  'freebayes/freebayes.pgrn.ontarget.minQ10',
                  'consensus/consensus',
                  'gatk/gatk.pgrn.ontarget.minQ10',
                  'original_calls/innocenti_pgrn_irinotecan_1_multisample_filtered' )



callers <- list('Atlas', 'Freebayes', 'CGES', 'GATK', 'Original')
tstv.fpaths <- mapply(make_fpath, call.sets, rep(proj.dir, length(call.sets)) )
tstv.dat <- lapply( tstv.fpaths, read.table, header=TRUE )


tstv <- data.frame( tstv = unlist( lapply(t(tstv.dat), extract_tstv) ),
                    caller = unlist(callers))
tstv$order_callers <- reorder(tstv$caller, tstv$tstv)

plt <- ggplot( tstv, aes(x=order_callers, y=tstv, fill=caller) ) +
        geom_bar(stat="identity") +
        labs(x="Caller", y="Ts/Tv") +
        theme_bw() + 
        scale_fill_brewer(palette="Set1")

show(plt)

  

