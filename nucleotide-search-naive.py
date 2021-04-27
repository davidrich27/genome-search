#!/usr/bin/env python
##################################################################################
# NAME: nucleotide-search.py
# DESC: Naive search of short reads in larger genome sequence.
##################################################################################

import sys

def print_all_indices_where_short_read_occurs(genome, r):
  results_for_read = []
  index = 0
  while True:
    index = genome.find(r,index)
    if index == -1:
      break
    index += 1
    results_for_read.append(index)
  print(''.join([str(x)+' ' for x in results_for_read]))

def main(argv):
  if len(argv) != 2:
    print('usage: <genome_fname> <reads_fname>')
  else:
    genome_fname, reads_fname = argv
    genome = open(genome_fname).read()
    genome = genome.replace('\n', '')
    for read in open(reads_fname).readlines():
      read = read.replace('\n','')
      print_all_indices_where_short_read_occurs(genome, read)

if __name__=='__main__':
  main(sys.argv[1:])
