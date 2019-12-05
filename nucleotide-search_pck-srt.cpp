/***********************************************************************
* 
*  NAME:  nucleotide_search.cpp
*
*  DESC:  Pack nucleotide into 
*
************************************************************************/

// imports
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

// local imports 
#include "../Clock.hpp"

// data type of packed blocks
#define DATA_TYPE unsigned long
// data size of packed blocks (in bits)
#define DATA_WIDTH (sizeof(DATA_TYPE) * CHAR_BIT)
// number of nucleotides that fit in the data type (divide by 2)
#define NUCL_WIDTH (DATA_WIDTH >> 1)

using namespace std;

// used so I can find the original index in string after sorting suffixes
typedef struct {
   DATA_TYPE data;
   unsigned int index;
} SUFFIX;

// encoder for finding the bit-packed values of nucleotides
int nucl_encoder[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
DATA_TYPE masker[] = { 0b0000000000000000000000000000000000000000000000000000000000000000, 0b1100000000000000000000000000000000000000000000000000000000000000, 0b1111000000000000000000000000000000000000000000000000000000000000, 0b1111110000000000000000000000000000000000000000000000000000000000, 0b1111111100000000000000000000000000000000000000000000000000000000, 0b1111111111000000000000000000000000000000000000000000000000000000, 0b1111111111110000000000000000000000000000000000000000000000000000, 0b1111111111111100000000000000000000000000000000000000000000000000, 0b1111111111111111000000000000000000000000000000000000000000000000, 0b1111111111111111110000000000000000000000000000000000000000000000, 0b1111111111111111111100000000000000000000000000000000000000000000, 0b1111111111111111111111000000000000000000000000000000000000000000, 0b1111111111111111111111110000000000000000000000000000000000000000, 0b1111111111111111111111111100000000000000000000000000000000000000, 0b1111111111111111111111111111000000000000000000000000000000000000, 0b1111111111111111111111111111110000000000000000000000000000000000, 0b1111111111111111111111111111111100000000000000000000000000000000, 0b1111111111111111111111111111111111000000000000000000000000000000, 0b1111111111111111111111111111111111110000000000000000000000000000, 0b1111111111111111111111111111111111111100000000000000000000000000, 0b1111111111111111111111111111111111111111000000000000000000000000, 0b1111111111111111111111111111111111111111110000000000000000000000, 0b1111111111111111111111111111111111111111111100000000000000000000, 0b1111111111111111111111111111111111111111111111000000000000000000, 0b1111111111111111111111111111111111111111111111110000000000000000, 0b1111111111111111111111111111111111111111111111111100000000000000, 0b1111111111111111111111111111111111111111111111111111000000000000, 0b1111111111111111111111111111111111111111111111111111110000000000, 0b1111111111111111111111111111111111111111111111111111111100000000, 0b1111111111111111111111111111111111111111111111111111111111000000, 0b1111111111111111111111111111111111111111111111111111111111110000, 0b1111111111111111111111111111111111111111111111111111111111111100, 0b1111111111111111111111111111111111111111111111111111111111111111 };

/* print bit-packed data in binary */
void print_bin(DATA_TYPE data)
{
   cout << "0b";
   DATA_TYPE bit = 0;
   for (int i = 0; i < DATA_WIDTH; ++i)
   {
      bit = data >> (DATA_WIDTH - 1);
      cout << bit;
      data <<= 1;
   }
   cout << endl;
}

// print bit-packed data in {A,C,G,T}
void print_sread(DATA_TYPE sread, int len)
{
   DATA_TYPE mask = masker[1];
   DATA_TYPE ch = 0;
   for (int i = 0; i < len; ++i)
   {
      // print_bin(sread);
      ch = (mask & sread) >> (DATA_WIDTH - 2);
      sread <<= 2;
      // print_bin(ch);
      switch(ch) 
      {
         case 0:
            cout << 'A';
            break;
         case 1:
            cout << 'C';
            break;
         case 2:
            cout << 'G';
            break;
         case 3:
            cout << 'T';
            break;
         default:
            cout << 'X';
      }
   }
   cout << endl;
}

// compare two suffixes
bool compareSuffix(const SUFFIX &a, const SUFFIX &b)
{
   return a.data < b.data;
}

int main(int argc, char *argv[]) 
{
   string genome_fname, sreads_fname;
   vector<SUFFIX> genome, genome_edge;
   vector<DATA_TYPE> sreads; 
   vector<unsigned int> sreads_len;
   DATA_TYPE sread, sread_top, sread_btm, sgenome, mask, gen;
   unsigned int sread_len, idx, top_idx, btm_idx;
   char c_i, c_v;
   int n,i,j;
   vector< vector<unsigned int> > results;

   if (argc == 3) 
   {
      genome_fname = argv[1];
      sreads_fname = argv[2];

      // cout << "DATA_WIDTH: " << DATA_WIDTH << endl;
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
      DATA_TYPE block = 0;
      unsigned int b = 0;
      ifstream inf(genome_fname);

      // fill block with data before adding it to the genome (TODO: check if genome length < 32)?
      for (i = 0; i < (NUCL_WIDTH-1); ++i)
      {
         inf.get(c_i);
         block <<= 2;
         block += nucl_encoder[c_i];
      }

      // interior blocks
      while (inf.get(c_i))
      {
         if (c_i == '\n')
            break;

         block <<= 2;
         block += nucl_encoder[c_i];
         // print_bin(block);
         ++b;
         genome.push_back({block, b});
      }
      inf.close();

      // reached end of sequence (ending data is null)
      for (i=0; i<(DATA_WIDTH-1); i += 2)
      {
         block <<= 2;
         // print_bin(block);
         ++b;
         genome_edge.push_back({block, b});
      }
   }

   // read in short reads
   {
      DATA_TYPE block = 0;
      unsigned int b = 0;
      ifstream inf(sreads_fname);
      while (inf.get(c_i))
      {
         // short read is <32 bits, so don't worry about overflowing block
         if (c_i != '\n')
         {
            block <<= 2;
            block += nucl_encoder[c_i];
            b++;
         }
         else 
         {
            // when we reach end of each short read
            block <<= (DATA_WIDTH - b*2); // shift data to left side of register
            sreads.push_back(block);
            sreads_len.push_back(b);
            b = 0;
            block = 0;

            // add a result vector for that short read
            vector<unsigned int> result;
            results.push_back(result);
         }
      }
      inf.close();
   }

   // sort genome 
   sort(genome.begin(), genome.end(), [](const SUFFIX &a, const SUFFIX &b) { return a.data < b.data; });

   // // verify sorted
   // for (i = 0; i < genome.size()-1; ++i)
   // {
   //    assert(genome[i].data < genome[i+1].data);
   // }

   // cout << "SORTED GENOME:" << endl;
   // for (i = 0; i < genome.size(); ++i)
   // {
   //    print_bin(genome[i].data);
   // }

   // search for each shortread in genome
   for (n = 0; n < sreads.size(); ++n)
   {
      // cout << "n=" << n << "..." << endl;
      sread = sreads[n];
      sread_len = sreads_len[n];
      mask = masker[sread_len];

      // smallest value containing sread
      sread_top = sread;
      // largest value containing sread
      sread_btm = sread | ~mask;

      // cout << "LEN: " << sread_len << endl;
      // cout << "TOP: ";
      // print_bin(sread_top);
      // cout << "BTM: ";
      // print_bin(sread_btm);
      // cout << "MSK: ";
      // print_bin(mask);

      // perform binary search for matching internal genome suffixes
      top_idx = genome.size() >> 1;
      btm_idx = genome.size() >> 1;
      for (i = (genome.size() >> 2); i >= 1; i >>= 1)
      {
         // find last element smaller than match
         if (sread_top <= genome[top_idx].data) 
         {
            top_idx -= i;
         }
         else
         {
            top_idx += i;
         }

         // find first element larger than match
         if (sread_btm >= genome[btm_idx].data)
         {
            btm_idx += i;
         }
         else
         {
            btm_idx -= i;
         }
      }

      // cout << "(" << top_idx << "," << btm_idx << ")" << endl;

      // final correction
      while (sread_top > genome[top_idx].data)
      {
         if (top_idx >= genome.size()-1)
            break;
         top_idx += 1;
      }
      while (sread_top <= genome[top_idx].data)
      {
         if (top_idx <= 0)
            break;
         top_idx -= 1;
      }
      while (sread_btm < genome[btm_idx].data)
      {
         if (btm_idx <= 0)
            break;
         btm_idx -= 1;
      }
      while (sread_btm >= genome[btm_idx].data) 
      {
         if (btm_idx >= genome.size()-1)
            break;
         btm_idx += 1;
      }

      // cout << "(" << top_idx << "," << btm_idx << ")" << endl;
      // cout << "T-E: ";
      // print_bin(genome[top_idx].data);
      // cout << "B-E: ";
      // print_bin(genome[btm_idx].data);

      // check top edge
      if (sread == (genome[top_idx].data & mask) )
      {
         results[n].push_back(genome[top_idx].index);
      }

      // add all elements suffixes in search range to results
      for (i = top_idx+1; i < btm_idx; i++)
      {
         results[n].push_back(genome[i].index);
      }

      // check bottom edge
      if (sread == (genome[btm_idx].data & mask) && top_idx != btm_idx)
      {
         results[n].push_back(genome[btm_idx].index);
      }

      // genome substrings with null chars at end
      unsigned int gen_len = NUCL_WIDTH;
      for (i = 0; i < NUCL_WIDTH; ++i)
      {
         --gen_len;

         // if read length is longer than remain genome seq, we can't have a match
         if ( sread_len > gen_len ) 
         {
            break;
         } 

         gen = genome_edge[i].data;
         idx = genome_edge[i].index;
         gen &= mask;

         if (gen == sread) 
         {
            results[n].push_back(idx);
         }
      }
   }
   // exit(0);

   c.ptock();

   // print results
   // cout << "=== RESULTS ===" << endl;
   for (n = 0; n < results.size(); ++n)
   {
      vector<unsigned int> result = results[n];
      sort(result.begin(), result.end());

      for (i = 0; i < result.size(); ++i)
      {
         cout << result[i] << " ";
      }
      cout << endl;
   }

   exit(0);
}