import os

lookup = dict()
for (root, dir, files) in os.walk('/nas40t2/innocentiData84/bamFiles'):
  bams = [ f.split('.') for f in files if f.endswith('bam') ]
  for bam in bams:
    id = bam[0]
    sm = bam[1]
    lookup[id] = sm
    lookup[sm] = id

covarFile = open('../annotation/Covariates_Vasa.csv')
newcovarFile = open('../annotation/Covariates_Vasa_newIDs.csv', 'w')
header = covarFile.next().strip().split('\t')
header.append('SM_ID')
print >> newcovarFile, '\t'.join(header)
for line in covarFile:
  line = line.strip().split('\t')
  id = line[0]
  sm = lookup.get(id)

  if sm:
    line.append(sm)
  else:
    line.append('NA')
    print id, sm

  print >> newcovarFile, '\t'.join(line)

