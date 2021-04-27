#################################################################################
# 
# NAME: build_short_reads.py
#
# DESC: Generates short reads from genome file.
# 
#################################################################################

import sys
from random import randint

if len(sys.argv) == 5:
   fname = sys.argv[1]
   N = int(sys.argv[2])
   mini = int(sys.argv[3])
   maxi = int(sys.argv[4])
else:
   print("Usage: <genome_filename> <number_of_reads> <min_length> <max_length>")
   sys.exit(0)


with open(fname) as fp:
   genome = ""
   line = fp.readline()
   while (line):
      genome += line
      line = fp.readline()


for i in range(N):
   s = randint(mini, maxi)
   n = randint(0, len(genome) - s)
   print(genome[n:n+s])