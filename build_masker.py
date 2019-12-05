#################################################################################
# 
# NAME: build_genome.py
#
# DESC: Generates a random string of ACGTs.
# 
#################################################################################

datasize = 8 * 8

print("DATA_TYPE masker = { ", end="")

for i in range(0, datasize+1, 2):

   b = "0b"

   for j in range(i):
      b+= "1"

   for j in range(i, datasize):
      b += "0"

   print(b, end="")

   if (i < datasize ):
      print(", ", end="")

print(" };")
