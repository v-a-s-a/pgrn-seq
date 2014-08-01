#!/usr/bin/python

import vcf
import matplotlib.pyplot as plt
import numpy as np


def __main__():
  '''
  Check genotype concordance between PGRNseq and Illumina exome chip data.
  '''

  pgrnReader = vcf.Reader(open("../variant_data/consensus.geno.vcf.gz"))
  exomeReader = vcf.Reader(open("../exome_data/innocenti_082613_clean.vcf.gz"))

  commonSamples =  set(pgrnReader.samples).union(set(exomeReader.samples))
  numSamples = len(commonSamples)  

  overlap = set()
  complement = set()

  siteConcordance = list()

  totalGeno = 0
  totalConcordant = 0

  ## iterate through PGRN variants, and query exome file
  for pgrnRec in pgrnReader:
    ## try to fetch corresponding variant in the exome file
    try:
      exomeRec = exomeReader.fetch(pgrnRec.CHROM.strip('chr'), pgrnRec.POS)
    except KeyError:    
      pass

    ## keep variants found in overlap and complement -- move on if record does not overlap
    if exomeRec:
      overlap.add(pgrnRec)
      ## look at genotype concordance
      concordant = 0
      for sample in commonSamples:
        if pgrnRec.genotype(sample)['GT'] == exomeRec.genotype(sample)['GT']:
          concordant = concordant + 1
      totalConcordant = totalConcordant + concordant
      totalGeno = totalGeno + numSamples
      siteConcordance.append( float(concordant) / numSamples )
    else:
      complement.add(pgrnRec)

  print "Number of overlapping variant sites:", len(overlap)
  print "Number of non-overlapping variant sites:", len(complement)
  print "Genotype concordance:", str(float(totalConcordant) / float(totalGeno))

  plt.hist(np.array(siteConcordance), 100) 
  plt.show()

if __name__ == "__main__":
  __main__()

