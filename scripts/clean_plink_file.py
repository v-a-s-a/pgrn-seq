'''
Our exome chip results contain some funny allele codes that are really messing up down stream analyses.

Essentially, we are interested in iterating concurrently over the .ped and .map file. Any time a "D" or "I" allele is encountered, we want to remove the entire site from the data set.
'''

def __main__():

  ## open ped and map files
  pedFile = '../exome_data/innocenti_082613.ped'

  ## iterate through ped file excluding the sites slated for removal
  newPed = open('../exome_data/../exome_data/innocenti_082613_fixedalleles.ped','w')
  for line in open(pedFile):
    line = line.strip().split()
    newGeno = [ "0" if x=="I" or x=="D" else x for x in line[7:] ]
    newLine = line[:7] + newGeno
    print >> newPed, '\t'.join(newLine) 


if __name__ == "__main__":
  __main__()

