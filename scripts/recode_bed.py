'''
Recode some of the custom contigs into canonical form.
'''

oldbed = open('bed/PGRN-seq_ALL_target_Final.bed')
newbed = open('bed/PGRN-seq_ALL_target_Final_cononical.bed', 'w')


for line in oldbed:
  line = line.lstrip('chr')
  line = line.split()
  line[0] = line[0].split('_')[0]
  print >> newbed, '\t'.join(line)






