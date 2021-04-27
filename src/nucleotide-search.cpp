/***********************************************************************
*  NAME:       nucleotide_search.cpp
*  PURPOSE:    Search for >32 character short reads in larger genome file.
*  METHOD:     - Encodes 4-alphabet {A,C,G,T} genome sequence into bitpacked 2-bits {00,01,10,11} per character.
*              - Builds suffix array of genome file by generating a suffix starting at every position in genome and sorting array.
*              - For each short read, creates minimum and maximum prefix containing short read by padding with 0s or 1s, resp.
*              - Perform binary search for minimum prefix and maximum prefix.  
*              - If results are found, the inclusive range of maximum and minimum search within suffix array all match short read.
*
************************************************************************/

// imports
#include <iostream>
#include <fstream>
#include <string>
#include <climits>

// local imports
#include "Clock.hpp"

// data type of bitpacked blocks
#define SUFFIX_DATATYPE unsigned long
// data size of packed blocks (in bits)
#define DATA_WIDTH (sizeof(SUFFIX_DATATYPE) * CHAR_BIT)
// number of nucleotides that fit in the data type (divide by 2)
#define NUCL_WIDTH (DATA_WIDTH >> 1)

using namespace std;

/* bitpacked suffix and index into genome */
typedef struct {
   SUFFIX_DATATYPE data;
   unsigned int index;
} SUFFIX;

/* bitpacked short reads */
typedef struct {
   SUFFIX_DATATYPE data;
   unsigned int len;
} SREAD;

/* resizable suffix array of genome */
typedef struct {
   unsigned int len;
   unsigned int size;
   SUFFIX *data;
} VECTOR_SUFFIX;

/* short read vector */
typedef struct {
   unsigned int len;
   unsigned int size;
   SREAD *data;
} VECTOR_SREAD;

/* result vector of matching data references */
typedef struct {
   unsigned int len;
   unsigned int size;
   unsigned int *data;
} RESULT;

/* encoder for quickly translating from 8-bit ascii values to the bit-packed values of nucleotides  */
int nucl_encoder[] = { 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 
};

/* bitmasks for selecting a extracting prefix of bit-packed data */
SUFFIX_DATATYPE masker[] = { 
    0b0000000000000000000000000000000000000000000000000000000000000000, 0b1100000000000000000000000000000000000000000000000000000000000000, 
    0b1111000000000000000000000000000000000000000000000000000000000000, 0b1111110000000000000000000000000000000000000000000000000000000000, 
    0b1111111100000000000000000000000000000000000000000000000000000000, 0b1111111111000000000000000000000000000000000000000000000000000000, 
    0b1111111111110000000000000000000000000000000000000000000000000000, 0b1111111111111100000000000000000000000000000000000000000000000000, 
    0b1111111111111111000000000000000000000000000000000000000000000000, 0b1111111111111111110000000000000000000000000000000000000000000000, 
    0b1111111111111111111100000000000000000000000000000000000000000000, 0b1111111111111111111111000000000000000000000000000000000000000000, 
    0b1111111111111111111111110000000000000000000000000000000000000000, 0b1111111111111111111111111100000000000000000000000000000000000000, 
    0b1111111111111111111111111111000000000000000000000000000000000000, 0b1111111111111111111111111111110000000000000000000000000000000000, 
    0b1111111111111111111111111111111100000000000000000000000000000000, 0b1111111111111111111111111111111111000000000000000000000000000000, 
    0b1111111111111111111111111111111111110000000000000000000000000000, 0b1111111111111111111111111111111111111100000000000000000000000000, 
    0b1111111111111111111111111111111111111111000000000000000000000000, 0b1111111111111111111111111111111111111111110000000000000000000000, 
    0b1111111111111111111111111111111111111111111100000000000000000000, 0b1111111111111111111111111111111111111111111111000000000000000000, 
    0b1111111111111111111111111111111111111111111111110000000000000000, 0b1111111111111111111111111111111111111111111111111100000000000000, 
    0b1111111111111111111111111111111111111111111111111111000000000000, 0b1111111111111111111111111111111111111111111111111111110000000000, 
    0b1111111111111111111111111111111111111111111111111111111100000000, 0b1111111111111111111111111111111111111111111111111111111111000000, 
    0b1111111111111111111111111111111111111111111111111111111111110000, 0b1111111111111111111111111111111111111111111111111111111111111100, 
    0b1111111111111111111111111111111111111111111111111111111111111111 
};

/* print bit-packed data in binary */
void print_bin(SUFFIX_DATATYPE data)
{
   cout << "0b";
   SUFFIX_DATATYPE bit = 0;
   for (int i = 0; i < DATA_WIDTH; ++i)
   {
      bit = data >> (DATA_WIDTH - 1);
      cout << bit;
      data <<= 1;
   }
   cout << endl;
}

/* print bit-packed data in ascii {A,C,G,T} */
void print_sread(SUFFIX_DATATYPE sread, int len)
{
   SUFFIX_DATATYPE mask = masker[1];
   SUFFIX_DATATYPE ch = 0;
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

/* selection sort array <arr> on range (<lo> inclusive, <hi> exclusive) */ 
void suffix_selectsort(SUFFIX*__restrict arr, unsigned int lo, unsigned int hi)
{
   SUFFIX tmp;
   SUFFIX_DATATYPE min;
   unsigned int min_idx;

   for (unsigned int i=lo; i<hi; ++i) {
      min = arr[i].data;
      min_idx = i;
      for (unsigned int j=i; j<hi; ++j) {
         if (min > arr[j].data) {
            min = arr[j].data;
            min_idx = j;
         }
      }
      tmp = arr[i];
      arr[i] = arr[min_idx];
      arr[min_idx] = tmp;
   }

   return;
}

/* quiksort array <arr> on range (<lo> inclusive, <hi> exclusive) */ 
void suffix_quicksort(SUFFIX*__restrict arr, unsigned int lo, unsigned int hi)
{
   /* if we have recursed to single element, return */
   if (hi - lo <= 1)
   { 
      // suffix_selectsort(arr, lo, hi);
      return; 
   }

   // printf("...quicksort, depth=%d, lo=%d, hi=%d\n", depth, lo, hi);
   unsigned int depth_1, depth_2;
   unsigned int k = lo;   /* k is pivot element */
   unsigned int i = lo;   /* i is left index */
   unsigned int j = hi-1; /* j is right index */
   
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

   /* swap value at i with pivot element */
   arr[lo] = arr[i];
   arr[i] = pivot;
   // array_Print(arr, lo, hi);

   /* recurse on sub-problems */
   suffix_quicksort( arr, lo, i );
   suffix_quicksort( arr, i+1, hi );

   return;
}

/* selection sort of results array <arr> on range (<lo> inclusive, <hi> exclusive) */ 
void results_selectsort(unsigned int*__restrict arr, unsigned int lo, unsigned int hi)
{
   unsigned int tmp;
   unsigned int min, min_idx;

   for (unsigned int i=lo; i<hi; ++i) {
      min = arr[i];
      min_idx = i;
      for (unsigned int j=i; j<hi; ++j) {
         if (min > arr[j]) {
            min = arr[j];
            min_idx = j;
         }
      }
      tmp = arr[i];
      arr[i] = arr[min_idx];
      arr[min_idx] = tmp;
   }

   return;
}

/* quicksort of results array <arr> on range (<lo> inclusive, <hi> exclusive) */
void results_quicksort(unsigned int*__restrict arr, unsigned int lo, unsigned int hi)
{
   /* if we have recursed to single element, return */
   if (hi - lo <= 1)
   { 
      // results_selectsort(arr, lo, hi);
      return; 
   }

   // printf("...quicksort, depth=%d, lo=%d, hi=%d\n", depth, lo, hi);
   unsigned int depth_1, depth_2;
   unsigned int k = lo;   /* k is pivot element */
   unsigned int i = lo;   /* i is left index */
   unsigned int j = hi-1; /* j is right index */
   
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
            /* then, increment i until we find an i > pivot */
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

   /* swap value at i with pivot element */
   arr[lo] = arr[i];
   arr[i] = pivot;
   // array_Print(arr, lo, hi);

   /* recurse on sub-problems */
   results_quicksort( arr, lo, i );
   results_quicksort( arr, i+1, hi );

   return;
}

/*    
 *
 */
int main(int argc, char *argv[])
{
   string genome_fname;          /* input filename: full genome */
   string sreads_fname;          /* input filename: short reads */
   VECTOR_SUFFIX genome;         /* target: suffix array of genome */
   VECTOR_SUFFIX genome_edge;    /* target: suffix array of tail genome, which does not fill vector */
   VECTOR_SUFFIX buffer;         /* buffer for reading in suffix arrays */
   VECTOR_SREAD sreads;          /* queries: vector of short genome reads */
   SUFFIX_DATATYPE sread, sread_top, sread_btm, sgenome, mask, gen;
   unsigned int idx, top_idx, btm_idx, sread_len;  /* indices for binary searches */
   char char_at_i;               /* buffer character for reading in files */
   int n, i, j;                  /* indices for sorting */
   RESULT*__restrict results;

   /* parse commandline arguments */
   if (argc == 3)
   {
      genome_fname = argv[1];
      sreads_fname = argv[2];
   }
   else
   {
      /* output help in_fileo if incorrect number of args given */
      cout << "Usage: <genome_filename> <short_reads_filename>" << endl;
      exit(0);
   }

   /* start trial timer */
   Clock c;
   c.tick();

   /* read in target genome */
   {
      genome.size = 256;
      genome.data = (SUFFIX *) malloc(genome.size * sizeof(SUFFIX));

      SUFFIX_DATATYPE block = 0;
      unsigned int b = 0;

      /* open genome file */
      ifstream in_file(genome_fname);

      /* fill block with data before adding it to the genome */
      /* (TODO: check if genome->length < 32)? */
      for (i = 0; i < (NUCL_WIDTH - 1); ++i)
      {
         in_file.get(char_at_i);
         block <<= 2;
         block += nucl_encoder[char_at_i];
      }

      /* interior blocks of data */
      while (in_file.get(char_at_i))
      {
         /* check for end of line */
         if (char_at_i == '\n')
            break;

         /* shift bitpacked sequence left by two (one bit-packed character) */
         block <<= 2;
         block += nucl_encoder[char_at_i];

         /* add suffix to array with position in sequence */
         genome.data[b] = {block, b + 1};
         ++b;

         /* resize genome suffix array if needed */
         if (b >= genome.size)
         {
            genome.size *= 2;
            genome.data = (SUFFIX *) realloc(genome.data, genome.size * sizeof(SUFFIX));
         }
      }
      genome.len = b;
      in_file.close();

      /* reached end of sequence (ending data is null) */
      genome_edge.len = DATA_WIDTH;
      genome_edge.size = DATA_WIDTH;
      genome_edge.data = (SUFFIX *) malloc(genome_edge.len * sizeof(SUFFIX));

      /* store all partial suffixes (not full block) in separate list */
      for (i = 0; i < ((DATA_WIDTH - 1) >> 1); i += 1)
      {
         block <<= 2;
         genome_edge.data[i] = {block, b + 1};
         ++b;
      }
   }

   /* read in query short reads */
   {
      sreads.size = 128;
      sreads.len = 0;
      sreads.data = (SREAD *) malloc(sreads.size * sizeof(SREAD));

      SUFFIX_DATATYPE block = 0;
      unsigned int b = 0;

      /* open short reads file */
      ifstream in_file(sreads_fname);

      while (in_file.get(char_at_i))
      {
         // all short read are <32 bits, so don't worry about overflowing block
         if (char_at_i != '\n')
         {
            block <<= 2;
            block += nucl_encoder[char_at_i];
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
      in_file.close();

      /* create results array */
      results = (RESULT *) malloc(sreads.size * sizeof(RESULT));
      for (i = 0; i < sreads.len; ++i)
      {
         // add result array for short read
         results[i].size = 128;
         results[i].data = (unsigned int *) malloc(128 * sizeof(unsigned int));
         results[i].len = 0;
      }
   }

   /* sort genome suffix array */
   suffix_quicksort(genome.data, 0, genome.len);

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

      /* smallest value containing sread => query, with ending padded by 0s */
      sread_top = sread;
      /* largest value containing sread => query, with ending padded by 1s */
      sread_btm = sread | ~mask;

      /* perform binary search for matching internal genome suffixes */
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

      /* final correction for off-by-one errors */
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

      /* check if top edge is a match */
      if (sread == (genome.data[top_idx].data & mask) )
      {
         results[n].data[results[n].len] = genome.data[top_idx].index;
         results[n].len++;
      }

      /* add all suffixes in internal search range to results */
      for (i = top_idx + 1; i < btm_idx; i++)
      {
         results[n].data[results[n].len] = genome.data[i].index;
         results[n].len++;
         /* resize results data if necessary */
         if (results[n].len >= results[n].size) {
            results[n].size *= 2;
            results[n].data = (unsigned int *) realloc(results[n].data, results[n].size * sizeof(unsigned int));
         }
      }

      /* check bottom edge is a match */
      if (sread == (genome.data[btm_idx].data & mask) )
      {
         results[n].data[results[n].len] = genome.data[btm_idx].index;
         results[n].len++;
      }

      /* genome substrings with null chars at end */
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

   /* end trial timer */
   c.ptock();

   /* print results */
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
