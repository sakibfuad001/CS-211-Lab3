/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1, p, n) - 1)
#define size_blk 1048576

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   int   local_first;
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   char  *local_prime_marked;
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;
   unsigned long int  local_prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */
   unsigned long int  local_prime_size;
	unsigned long long int low_value_blk;
   unsigned long long int high_value_blk;


   MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);

   /* Add you code here  */

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

	local_prime_size  = (int) sqrt((double)(n)) - 1;
	low_value = (2 + BLOCK_LOW(id, p, n - 1)) + ((2 + BLOCK_LOW(id, p, n - 1)) + 1) % 2;
   high_value = (2 + BLOCK_HIGH(id, p, n - 1)) - ((2 + BLOCK_HIGH(id, p, n - 1)) + 1) % 2;
   size = 1 + ((high_value - low_value) / 2);
   proc0_size = ((n/2) - 1) / p;
	unsigned long int check = (int) sqrt((double) n/2);

   if ((proc0_size + 2) < check){
   	if (!id) printf("Too many processes.\n");
      MPI_Finalize();
      exit(1);
   }

	/* Allocate this process's share of the array. */

	marked = (char*) malloc(size);
  	local_prime_marked = (char *) malloc (local_prime_size);

  	if (local_prime_marked == NULL || marked == NULL) {   
   	printf("Cannot allocate enough memory.\n");
      MPI_Finalize();
      exit(1);
   }

	local_prime = 2;


   for (i = 0; i < local_prime_size; i++) local_prime_marked[i] = 0;
	for (i = 0; i < size; i++) marked[i] = 0;

   index = 0;
   do {
   	local_first = (local_prime * local_prime) - 2;
      for (i = local_first; i < local_prime_size; i += local_prime){
      	local_prime_marked[i] = 1;
		}

      while (local_prime_marked[++index] == 1);

      local_prime = 2 + index;
   } while (local_prime * local_prime <= n);

   low_value_blk = low_value;
   high_value_blk = 2 * (size_blk - 1) + low_value_blk;

   do {
   	prime = 3;
		index = 0;

      while (prime * prime <= high_value_blk) {
      	if (prime * prime > low_value_blk) {
         	first = (prime * prime - low_value_blk) / 2;
			} else {
				first = (low_value_blk % prime) == 0 ? 0 : ((prime - (low_value_blk % prime) + low_value_blk / prime % 2 * prime) / 2);
         }

         for (i = first + (low_value_blk - low_value) / 2; i < ((high_value_blk - low_value) / 2) + 1; i = i + prime) marked[i] = 1;
         while(local_prime_marked[++index] == 1);
         prime = index + 2;
      }

      low_value_blk = 2 + high_value_blk;
      high_value_blk =  2 * (size_blk - 1) + low_value_blk;

      if (high_value < high_value_blk) high_value_blk = high_value;
    } while (low_value_blk <= high_value);

    count = 0;
    for (i = 0; i < size; i++)
        if (!(marked[i])) count++;
    if (!(id)) count++; 
    if (p > 1)
        MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

   /* Stop the timer */

   elapsed_time += MPI_Wtime();

   /* Print the results */

   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);

   }
   MPI_Finalize ();
   return 0;
}

