#################################################################################
# 
# NAME: build_genome.py
#
# DESC: Generates a random string of ACGTs.
# 
#################################################################################

import sys
from random import randint

if len(sys.argv) == 2:
   N = int(sys.argv[1])
else:
   print("Usage: <Length_of_Genome>")
   sys.exit(0)

nucl_dict = {
   0: "A",
   1: "C",
   2: "G",
   3: "T",
   #"A": 0,
   #"C": 1,
   #"G": 2,
   #"T": 3
}

for i in range(N):
   print(nucl_dict[randint(0,3)], end="")
print("")