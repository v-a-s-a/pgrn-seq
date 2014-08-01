'''
Check which samples we are using for the IRB update.
'''

irbFile = 'innocenti_pgrn_QC_for_Vasa_021114.csv'
studyFile = 'current_consensus_samples.txt'

samplesUsed = set()
for sample in open(studyFile):
  sample = sample.strip()
  samplesUsed.add(sample)

irbConn = open(irbFile)
allowed = set()
for line in irbConn:
  line = line.strip().split()
  sample = line[3]
  allowed.add(sample)

for sample in samplesUsed:
  if not sample in allowed: print sample

