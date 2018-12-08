import os
import io
import ast
import os.path

chrom = 22

with open('/home/raul/Documents/genomes/data/chr'+str(chrom)+'_PUR_inds.txt', 'r') as ff:
   PUR_inds = [line.rstrip('\n') for line in ff]
