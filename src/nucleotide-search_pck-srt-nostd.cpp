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
#include <stdlib.h>
#include <time.h>

// local imports
#include "../Clock.hpp"

// data type of packed blocks
#define DATA_TYPE unsigned long
#define SORT_DATA SUFFIX
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

typedef struct {
   DATA_TYPE data;
   unsigned int len;
} SREAD;

typedef struct {
   unsigned int len;
   unsigned int size;
   SUFFIX *data;
} VECTOR_SUFFIX;

typedef struct {
   unsigned int len;
   unsigned int size;
   SREAD *data;
} VECTOR_SREAD;

typedef struct {
   unsigned int len;
   unsigned int size;
   unsigned int *data;
} RESULT;

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
      switch (ch)
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

void suffix_selection_sort(SUFFIX *src, unsigned int N)
{
   unsigned int i, j;
   unsigned int min_idx;
   SUFFIX min, tmp;

   for (i = 0; i < N; i++)
   {
      min = src[i];
      min_idx = i;

      for (j = i; j < N; j++)
      {
         if (min.data > src[j].data)
         {
            min = src[j];
            min_idx = j;
         }
      }

      tmp = src[min_idx];
      src[min_idx] = src[i];
      src[i] = tmp;
   }
}

void suffix_mergesort(SUFFIX *source, SUFFIX *buffer, const unsigned long N, const unsigned long K )
{
   if (N <= K)
   {
      suffix_selection_sort(source, N);
      return;
   }

   SUFFIX *source_2 = source + N / 2;
   suffix_mergesort(source, buffer, N / 2, K);
   suffix_mergesort(source_2, buffer, N - N / 2, K);
   // Merge sorted halves into buffer:
   unsigned long i = 0, j = 0, buffer_ind = 0;

   while (i < N / 2 && j < (N - N / 2))
   {
      if (source[i].data < source_2[j].data)
      {
         buffer[buffer_ind] = source[i];
         ++i;
      }
      else
      {
         // In case of equality, order doesn't matter, so use this case:
         buffer[buffer_ind] = source_2[j];
         ++j;
      }
      ++buffer_ind;
   }
   // Copy remaining values:
   for (; i < N / 2; ++i, ++buffer_ind)
      buffer[buffer_ind] = source[i];
   for (; j < N - N / 2; ++j, ++buffer_ind)
      buffer[buffer_ind] = source_2[j];
   // Copy back sorted list from buffer:
   for (i = 0; i < N; ++i)
      source[i] = buffer[i];
}

void suffix_quicksort(SUFFIX *arr, int lo, int hi)
{
   /* if we have recursed to single element, return */
   if (hi - lo <= 1)
   { 
      return; 
   }

   // printf("...quicksort, depth=%d, lo=%d, hi=%d\n", depth, lo, hi);
   int depth_1, depth_2;
   int k = lo;   /* k is pivot element */
   int i = lo;   /* i is left index */
   int j = hi-1; /* j is right index */
   
   SUFFIX pivot = arr[k]; 
   // printf("pivot: %.0f\n", pivot);
   // array_Print(arr, lo, hi);
   SUFFIX swp;

   while (i < j)
   {
      /* decrement j until we find a j <= pivot */
      if (arr[j].data > pivot.data)
      {
         j--;
      }
      else
      {
         while (i < j)
         {
            /* then, increment until we find an i > pivot */
            if (arr[i].data <= pivot.data)
            {
               i++;            
            }
            /* then, swap values at i and j */
            else 
            {
               // printf("swp: %d <-> %d\n", i, j);
               std::swap( arr[i], arr[j] );
               // array_Print(arr, lo, hi);
               j--;
               break;
            }
         }
      }
   }

   // printf("end quicksort, i=%d, j=%d\n", i, j);

   /* swap value at i with pivot element */
   arr[lo] = arr[i];
   arr[i] = pivot;
   // array_Print(arr, lo, hi);

   /* recurse on sub-problems */
   suffix_quicksort( arr, lo, i );
   suffix_quicksort( arr, i+1, hi );

   return;
}

void results_quicksort(unsigned int *arr, int lo, int hi)
{
   /* if we have recursed to single element, return */
   if (hi - lo <= 1)
   { 
      return; 
   }

   // printf("...quicksort, depth=%d, lo=%d, hi=%d\n", depth, lo, hi);
   int depth_1, depth_2;
   int k = lo;   /* k is pivot element */
   int i = lo;   /* i is left index */
   int j = hi-1; /* j is right index */
   
   unsigned int pivot = arr[k]; 
   // printf("pivot: %.0f\n", pivot);
   // array_Print(arr, lo, hi);
   unsigned int swp;

   while (i < j)
   {
      /* decrement j until we find a j <= pivot */
      if (arr[j] > pivot)
      {
         j--;
      }
      else
      {
         while (i < j)
         {
            /* then, increment until we find an i > pivot */
            if (arr[i] <= pivot)
            {
               i++;            }
            /* then, swap values at i and j */
            else 
            {
               // printf("swp: %d <-> %d\n", i, j);
               std::swap( arr[i], arr[j] );
               // array_Print(arr, lo, hi);
               j--;
               break;
            }
         }
      }
   }

   // printf("end quicksort, i=%d, j=%d\n", i, j);

   /* swap value at i with pivot element */
   arr[lo] = arr[i];
   arr[i] = pivot;
   // array_Print(arr, lo, hi);

   /* recurse on sub-problems */
   results_quicksort( arr, lo, i );
   results_quicksort( arr, i+1, hi );

   return;
}

void results_sort(unsigned int *src, unsigned int N)
{
   unsigned int i, j;
   unsigned int min_idx;
   unsigned int min, tmp;

   for (i = 0; i < N; i++)
   {
      min = src[i];
      min_idx = i;

      for (j = i; j < N; j++)
      {
         if (min > src[j])
         {
            min = src[j];
            min_idx = j;
         }
      }

      tmp = src[min_idx];
      src[min_idx] = src[i];
      src[i] = tmp;
   }
}

int main(int argc, char *argv[])
{
   string genome_fname, sreads_fname;
   VECTOR_SUFFIX genome, genome_edge, buffer;
   VECTOR_SREAD sreads;
   DATA_TYPE sread, sread_top, sread_btm, sgenome, mask, gen;
   unsigned int idx, top_idx, btm_idx, sread_len;  // indices for binary searches
   char c_i;
   int n, i, j;
   RESULT *results;

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
      genome.size = 256;
      genome.data = (SUFFIX *) malloc(genome.size * sizeof(SUFFIX));

      DATA_TYPE block = 0;
      unsigned int b = 0;
      ifstream inf(genome_fname);

      // fill block with data before adding it to the genome (TODO: check if genome->length < 32)?
      for (i = 0; i < (NUCL_WIDTH - 1); ++i)
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

         genome.data[b] = {block, b + 1};
         ++b;

         // resize genome if needed
         if (b >= genome.size)
         {
            genome.size *= 2;
            genome.data = (SUFFIX *) realloc(genome.data, genome.size * sizeof(SUFFIX));
         }
      }
      genome.len = b;
      inf.close();

      // reached end of sequence (ending data is null)
      genome_edge.len = DATA_WIDTH;
      genome_edge.size = DATA_WIDTH;
      genome_edge.data = (SUFFIX *) malloc(genome_edge.len * sizeof(SUFFIX));
      for (i = 0; i < ((DATA_WIDTH - 1) >> 1); i += 1)
      {
         block <<= 2;
         genome_edge.data[i] = {block, b + 1};
         ++b;
      }
   }

   // read in short reads
   {
      sreads.size = 128;
      sreads.len = 0;
      sreads.data = (SREAD *) malloc(sreads.size * sizeof(SREAD));

      DATA_TYPE block = 0;
      unsigned int b = 0;
      ifstream inf(sreads_fname);

      while (inf.get(c_i))
      {
         // all short read are <32 bits, so don't worry about overflowing block
         if (c_i != '\n')
         {
            block <<= 2;
            block += nucl_encoder[c_i];
            b++;
         }
         else
         {
            // when we reach end of each short read
            block <<= (DATA_WIDTH - b * 2); // shift data to left side of register

            // add short read to array
            sreads.data[sreads.len] = {block, b};
            sreads.len++;

            // resize sreads and results if necessary
            if (sreads.len >= sreads.size)
            {
               sreads.size *= 2;
               sreads.data = (SREAD *) realloc(sreads.data, sreads.size * sizeof(SREAD));
            }

            b = 0;
            block = 0;
         }
      }
      inf.close();

      // results
      results = (RESULT *) malloc(sreads.size * sizeof(RESULT));
      for (i = 0; i < sreads.len; ++i)
      {
         // add result array for short read
         results[i].size = 128;
         results[i].data = (unsigned int *) malloc(128 * sizeof(unsigned int));
         results[i].len = 0;
      }
   }

   // // print genome
   // cout << "GENOME SUFFIXES (UNSORTED): " << endl;
   // for (int i=0; i<genome.len; i++)
   // {
   //    cout << genome.data[i].index << "\t";
   //    print_bin(genome.data[i].data);
   // }

   // // print sread
   // cout << endl << "SHORT READS: " << endl;
   // for (int i=0; i<sreads.len; i++)
   // {
   //    cout << sreads.data[i].len << "\t";
   //    // print_sread(sreads.data[i].data, sreads.data[i].len);
   //    print_bin(sreads.data[i].data);
   // }

   // sort genome
   // buffer.size = genome.size;
   // buffer.data = (SUFFIX *) malloc(buffer.size * sizeof(SUFFIX));
   // suffix_mergesort(genome.data, buffer.data, genome.len, 16);
   suffix_quicksort(genome.data, 0, genome.len);
   // Quicksort(genome.data, genome.len, 0);

   // cout << "SORTED GENOME:" << endl;
   // for (i = 0; i < genome.len; ++i)
   // {
   //    cout << genome.data[i].index << "\t";
   //    print_bin(genome.data[i].data);
   // }
   // cout << "EDGE GENOME:" << endl;
   // for (i = 0; i < NUCL_WIDTH-1; ++i)
   // {
   //    cout << genome_edge.data[i].index << "\t";
   //    print_bin(genome_edge.data[i].data);
   // }

   // verify sorted
   // for (i = 0; i < genome.len - 1; ++i)
   // {
   //    assert(genome.data[i].data < genome.data[i + 1].data);
   // }

   // search for each shortread in genome
   for (n = 0; n < sreads.len; ++n)
   {
      // cout << "n=" << n << "..." << endl;
      sread = sreads.data[n].data;
      sread_len = sreads.data[n].len;
      mask = masker[sread_len];

      // smallest value containing sread
      sread_top = sread;
      // largest value containing sread
      sread_btm = sread | ~mask;

      // cout << "\nLEN:\t" << sreads.data[n].len << endl;
      // cout << "TOP:\t";
      // print_bin(sread_top);
      // cout << "BTM:\t";
      // print_bin(sread_btm);
      // cout << "MSK:\t";
      // print_bin(mask);

      // perform binary search for matching internal genome suffixes
      top_idx = genome.len >> 1;
      btm_idx = genome.len >> 1;
      for (i = (genome.len >> 2); i >= 1; i >>= 1)
      {
         // find last element smaller than match
         if (sread_top <= genome.data[top_idx].data)
         {
            top_idx -= i;
         }
         else
         {
            top_idx += i;
         }

         // find first element larger than match
         if (sread_btm >= genome.data[btm_idx].data)
         {
            btm_idx += i;
         }
         else
         {
            btm_idx -= i;
         }
      }

      // cout << "\t(" << top_idx << "," << btm_idx << ")" << endl;

      // final correction
      while (sread_top > genome.data[top_idx].data)
      {
         if (top_idx >= genome.len - 1)
            break;
         top_idx += 1;
      }
      while (sread_top <= genome.data[top_idx].data)
      {
         if (top_idx <= 0)
            break;
         top_idx -= 1;
      }
      while (sread_btm < genome.data[btm_idx].data)
      {
         if (btm_idx <= 0)
            break;
         btm_idx -= 1;
      }
      while (sread_btm >= genome.data[btm_idx].data)
      {
         if (btm_idx >= genome.len - 1)
            break;
         btm_idx += 1;
      }

      // cout << "\t(" << top_idx << "," << btm_idx << ")" << endl;
      // cout << "T-E:\t";
      // print_bin(genome.data[top_idx].data);
      // cout << "B-E:\t";
      // print_bin(genome.data[btm_idx].data);

      // check top edge
      if (sread == (genome.data[top_idx].data & mask) )
      {
         results[n].data[results[n].len] = genome.data[top_idx].index;
         results[n].len++;
      }

      // add all elements suffixes in search range to results
      for (i = top_idx + 1; i < btm_idx; i++)
      {
         results[n].data[results[n].len] = genome.data[i].index;
         results[n].len++;
         // resize results data if necessary
         if (results[n].len >= results[n].size) {
            results[n].size *= 2;
            results[n].data = (unsigned int *) realloc(results[n].data, results[n].size * sizeof(unsigned int));
         }
      }

      // check bottom edge
      if (sread == (genome.data[btm_idx].data & mask) )
      {
         results[n].data[results[n].len] = genome.data[btm_idx].index;
         results[n].len++;
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

         gen = genome_edge.data[i].data;
         idx = genome_edge.data[i].index;
         gen &= mask;

         if (gen == sread)
         {
            results[n].data[results[n].len] = genome_edge.data[i].index;
            results[n].len++;
            // resize results data if necessary
            if (results[n].len >= results[n].size) {
               results[n].size *= 2;
               results[n].data = (unsigned int *) realloc(results[n].data, results[n].size * sizeof(unsigned int));
            }
         }
      }
   }

   c.ptock();

   // print results
   // cout << "=== RESULTS ===" << endl;
   for (n = 0; n < sreads.len; ++n)
   {
      results_quicksort(results[n].data, 0, results[n].len);
      for (i = 0; i < results[n].len; ++i)
      {
         cout << results[n].data[i] << " ";
      }
      cout << endl;
   }

   exit(0);
}