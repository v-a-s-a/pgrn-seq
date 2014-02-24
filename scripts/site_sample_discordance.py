'''
Calculate site and sample discordance.
'''

import vcf as pyvcf
import optparse as opt


def __main__():

  ## parse command line arguments
  parser = opt.OptionParser()
  parser.add_option('--vcf', dest='vcfFile')
  parser.add_option('--out', dest='outBase')
  (opts, args) = parser.parse_args()

  reader = pyvcf.Reader(open(opts.vcfFile))
  
  ## init data structures
  sampleConcordance = dict()
  for sample in reader.samples:
    sampleConcordance[sample] = 0
  siteConcordance = dict()
  
  ## init file outputs
  sampleConcordFile = open(opts.outBase + '.sample.concord.txt', 'w')
  siteConcordFile = open(opts.outBase + '.site.concord.txt', 'w')

  siteCount = 0
  total = 0
  for rec in reader:
    siteCount += 1    

    ## update sample concordance/discordance
    concordant = 0
    for sample in reader.samples:
      status = rec.genotype(sample)['CN']
      genotype = rec.genotype(sample)['GT']
      if status == 'T':
        sampleConcordance[sample] += 1
        concordant += 1
      else:
        total += 1

    ## compute site concordance and write out to file
    siteConcordRate = float(concordant) / float(len(reader.samples))
    print >> siteConcordFile, '\t'.join( [rec.ID, str(siteConcordRate)] )
  

  print '# concordant:', total

  ## calculate sample concordance rate and write to file
  for sample in reader.samples:
    concordRate = float(sampleConcordance[sample]) / float(siteCount)
    print >> sampleConcordFile, '\t'.join( [sample, str(concordRate)] )

  ## close file handles
  sampleConcordFile.close()
  siteConcordFile.close()


if __name__ == '__main__':
  __main__()






