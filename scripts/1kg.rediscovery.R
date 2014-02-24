## plot Exome Variant Server rediscovery
library(ggplot2)

#' Make a data frame of 1000 genomes rediscovery data
make_fpath <- function(base, basedir) {
  fpath <- paste( basedir, base, '.kg', sep = '' )
  return(fpath)
}

proj.dir <- '/nas40t0/vasya/pgrn_seq/qc_report/'
call.sets = list( 'atlas/atlas.pgrn.proc.ontarget.minQ10',
                  'freebayes/freebayes.pgrn.ontarget.minQ10',
                  'consensus/consensus',
                  'gatk/gatk.pgrn.ontarget.minQ10',
                  'original_calls/innocenti_pgrn_irinotecan_1_multisample_filtered' )


set.names <- list( 'Atlas', 'Freebayes', 'CGES', 'GATK', 'Original' )
kg.fpaths <-  mapply(make_fpath, call.sets, rep(proj.dir, length(call.sets)) )
kg.dat <- lapply( kg.fpaths, read.table )

kg <- data.frame( kg = unlist( lapply(kg.dat, subset, select=V2) ),
                   set = unlist(set.names) )
kg$order_sets <- reorder(kg$set, kg$kg)


plt <- ggplot( kg, aes(x=order_sets, y=kg, fill=set) ) +
        geom_bar(stat="identity") +
        labs(x="Caller", y="1000 Genomes rediscovery rate") +
        theme_bw() +
        scale_fill_brewer(palette="Set1")

show(plt)
