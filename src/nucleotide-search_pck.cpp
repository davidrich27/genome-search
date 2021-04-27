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

// local imports 
#include "../Clock.hpp"

// data type of packed blocks
#define DATA_TYPE unsigned long
// data size of packed blocks (in bits)
#define DATA_WIDTH (sizeof(DATA_TYPE) * CHAR_BIT)
// number of nucleotides that fit in the data type (divide by 2)
#define NUCL_WIDTH (DATA_WIDTH >> 1)

using namespace std;

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

int main(int argc, char *argv[]) 
{
   string genome_fname, sreads_fname;
   vector<DATA_TYPE> sreads, genome, genome_edge; 
   vector<int> sreads_len;
   DATA_TYPE sread, sgenome, mask, gen;
   int sread_len;
   char c_i, c_v;
   int n,i,j;
   vector<vector<int>> results;

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
      int b = 0;
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
         // // for faster encoding, remove character checking
         // if (nucl_encoder[c_i] == -1) 
         // {
         //    cerr << "ERROR: Bad Nucleotide Character" << endl;
         //    exit(1);
         // }

         if (c_i == '\n')
            break;

         block <<= 2;
         block += nucl_encoder[c_i];
         // print_bin(block);
         genome.push_back(block);
      }
      inf.close();

      // reached end of sequence (ending data is null)
      for (i=0; i<(DATA_WIDTH-1); i += 2)
      {
         block <<= 2;
         // print_bin(block);
         genome_edge.push_back(block);
      }
   }

   // read in short reads
   {
      DATA_TYPE block = 0;
      int b = 0;
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
            vector<int> result;
            results.push_back(result);
         }
      }
      inf.close();
   }

   // // print genome and short reads
   // cout << "=== GENOME === (" << genome.size() << ")" << endl; 
   // for (i = 0; i < genome.size(); ++i)
   // {
   //    print_sread(genome[i], NUCL_WIDTH);
   // }
   // cout << "=== SHORT READS ===" << endl;
   // for (i = 0; i < sreads.size(); ++i)
   // {
   //    print_sread(sreads[i], sreads_len[i]);
   // }

   // search for each shortread in genome
   for (n = 0; n < sreads.size(); ++n)
   {
      sread = sreads[n];
      sread_len = sreads_len[n];
      mask = masker[sread_len];

      // genome substrings that are > DATA_WIDTH from the end
      for (i = 0; i < genome.size(); ++i)
      {
         gen = genome[i];
         gen &= mask; 

         if (gen == sread) 
         {
            results[n].push_back(i+1);
         }
      }

      // genome substrings with null chars at end
      int gen_len = NUCL_WIDTH;
      for (i = 0; i < NUCL_WIDTH; ++i)
      {
         --gen_len;

         // if read length is longer than remain genome seq, we can't have a match
         if ( sread_len > gen_len ) 
         {
            break;
         } 

         gen = genome_edge[i];
         gen &= mask;

         if (gen == sread) 
         {
            results[n].push_back(i+genome.size()+1);
         }
      }
   }

   c.ptock();

   // print results
   // cout << "=== RESULTS ===" << endl;
   for (n = 0; n < results.size(); ++n)
   {
      vector<int> result = results[n];
      for (i = 0; i < result.size(); ++i)
      {
         cout << result[i] << " ";
      }
      cout << endl;
   }

   exit(0);
}