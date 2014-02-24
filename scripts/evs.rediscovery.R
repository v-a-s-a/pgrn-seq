## plot Exome Variant Server rediscovery
library(ggplot2)

#' Make a data frame of TsTV summary data
make_fpath <- function(base, basedir) {
  fpath <- paste( basedir, base, '.evs', sep = '' )
  return(fpath)
}

proj.dir <- '/nas40t0/vasya/pgrn_seq/qc_report/'
call.sets = list( 'atlas/atlas.pgrn.proc.ontarget.minQ10',
                  'freebayes/freebayes.pgrn.ontarget.minQ10',
                  'consensus/consensus',
                  'gatk/gatk.pgrn.ontarget.minQ10',
                  'original_calls/innocenti_pgrn_irinotecan_1_multisample_filtered' )


set.names <- list( 'Atlas', 'Freebayes', 'CGES', 'GATK', 'Original' )
evs.fpaths <-  mapply(make_fpath, call.sets, rep(proj.dir, length(call.sets)) )
evs.dat <- lapply( evs.fpaths, read.table )

evs <- data.frame( evs = unlist( lapply(evs.dat, subset, select=V2) ),
                   set = unlist(set.names) )
evs$order_sets <- reorder(evs$set, evs$evs)


plt <- ggplot( evs, aes(x=order_sets, y=evs, fill=set) ) +
        geom_bar(stat="identity") +
        labs(x="Caller", y="EVS rediscovery rate") +
        theme_bw() +
        scale_fill_brewer(palette="Set1")

show(plt)
