/***********************************************************************
* 
*  NAME:  nucleotide_search_unpacked.cpp
*
*  DESC:  Search for short reads in genome file.
*
************************************************************************/

// imports
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// local imports 
#include "Clock.hpp"

using namespace std;

int main(int argc, char *argv[]) 
{
   string genome_fname, sreads_fname;
   string genome, sread; 
   vector<string> sreads;
   string line;
   char c_i, c_v;
   bool id;
   int n,i,j;
   vector<vector<int>> results;

   if (argc == 3) 
   {
      genome_fname = argv[1];
      sreads_fname = argv[2];

      // cout << "GENOME FILE: " << genome_fname << endl;
      // cout << "SHORT READS FILE: " << sreads_fname << endl;
   } 
   else 
   {
      cout << "Usage: <genome_filename> <short_reads_filename>" << endl;
      exit(0);
   }

   Clock c;
   c.tick();

   // read in genome
   {
      ifstream inf(genome_fname);
      while (inf.get(c_i))
      {
         genome.push_back(c_i);
      }
      inf.close();
   }

   // read in short reads
   {
      ifstream inf(sreads_fname);
      while (inf)
      {
         getline(inf, line);
         sreads.push_back(line);

         // for each short read, add a result vector
         vector<int> result;
         results.push_back(result);
      }
      inf.close();
   }

   // print genome and short reads
   // cout << "=== GENOME ===" << endl << genome << endl;
   // cout << "=== SHORT READS ===" << endl;
   // for (int i = 0; i < sreads.size(); ++i)
   // {
   //    cout << sreads[i] << endl;
   // }

   // search for each shortread in genome
   for (n = 0; n < sreads.size(); ++n)
   {
      sread = sreads[n];
      for (i = 0; i < genome.size(); ++i)
      {
         if (genome[i] == sread[0]) 
         {
            id = true;
            for (j = 1; j < sread.size(); ++j)
            {
               if (genome[i+j] != sread[j])
               {
                  id = false;
                  break;
               }
            }
            if (id == true)
            {
               results[n].push_back(i+1);
            }
         }
      }
   }

   c.ptock();

   // print results
   // cout << "=== RESULTS ===" << endl;
   for (n = 0; n < results.size(); ++n)
   {
      for (i = 0; i < results[n].size(); ++i)
      {
         cout << results[n][i] << " ";
      }
      cout << endl;
   }

   exit(0);
}