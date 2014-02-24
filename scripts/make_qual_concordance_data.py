'''
Make a table of data which shows GATK QUAL versus concordance
'''

import vcf as pyvcf
import math

consensusVcf = '/nas40t0/vasya/autism/157_sample_publication_data/uchicago_autism_exomes_12-02-13/consensus_clean.vcf.gz'
gatkVcf = '/nas40t0/vasya/autism/157_sample_publication_data/uchicago_autism_exomes_12-02-13/gatk_clean.vcf.gz'

consensusReader = pyvcf.Reader(open(consensusVcf))
gatkReader = pyvcf.Reader(open(gatkVcf))


for rec in consensusReader:
  ## calculcate concordance
  concord =  [ 1.0 for sample in rec.samples if sample['CN']=='T' ]
  perc = sum(concord) / len(rec.samples)
  if perc == 1.0: continue
  concord = -math.log10( 1.0 - perc)  

  ## pull corresponding record from the GATK vcf
  print concord, gatkReader.fetch(rec.CHROM.strip('chr'), rec.POS).QUAL
