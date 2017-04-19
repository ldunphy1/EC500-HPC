/* RCB and ESW, 2016-02-05, Comparing different search algorithms.
   Various codes are borrowed from the web.
   Updated 2017-01-30. */

#include <iostream>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include "search.h"
#include "sort.h"

#define GIG 1000000000

using std::cout; 
using std::endl;

timespec diff(timespec start, timespec end)
{
  timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0)
  {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

// A function to print arrays.
void printArray(int list[], int arraySize)
{
  cout << " --------" << endl;
  for (int i = 0; i < arraySize; i++)
  {
    cout << list[i] <<  " " << endl;
  }
  cout << " --------" << endl;
}


int main()
{
	 // Iterators.
	int i, j;

	// Timing structures.
  	timespec time1, time2, timediff; 

	// Counters 
	int BinCount; // counter for binary search. 
	int LinCount; // counter for linear search. 
	int DctCount; // counter for dictionary search.
	int time_seconds = time(0);
	srand(time_seconds%100); // seed the random number generator.
	int N = 10000000;

	// Timing storage.
  	std::vector<double> timing_list_bin;
  	std::vector<double> timing_list_lin;
  	std::vector<double> timing_list_dct;

	for(j=10;j<N;j*=10)
	{
	  	int Nsize = j; 
	  	int * list = new  int[Nsize];

	  	// Fill the array with random numbers.
	  	for(i = 0; i < Nsize; i++)
	  	{
			list[i] =  rand()%Nsize;
	  	}

	  	// Randomly choose an element to search for.
	  	int find = list[rand()%Nsize];

	  	// Start timer.
    	clock_gettime(CLOCK_REALTIME, &time1);

	  	// First, search with a linear search. A linear search doesn't depend
	  	// on a list being sorted.
	  	i = linearSearch(list,find ,Nsize, LinCount);

	  	// End timer.
    	clock_gettime(CLOCK_REALTIME, &time2);

    	// Get difference.
    	timediff = diff(time1, time2);

    	// Save time in seconds to vector.
    	timing_list_lin.push_back(((double)GIG * (double)timediff.tv_sec + (double)timediff.tv_nsec)/((double)GIG));

	  	// Print the number of iterations compared to the average.
	  	//printf("Size %d is done. The iteration count was %d\n", j, LinCount);

	 	// Next, let's look at a binary search and a dictionary search.
	  	// This requires a sorted list. Feel free to play with these if
	  	// you want. 
	  	// WARNING: the selection sort takes a long time if Nsize gets
	  	// into the millions... you probably have better things to do.
	  
	  	mergeSort(list, Nsize);  

	  	// Start timer.
    	clock_gettime(CLOCK_REALTIME, &time1);

	  	//Perform a binary search. 
	  	i = binarySearch(list,find ,Nsize, BinCount); 

	  	// End timer.
    	clock_gettime(CLOCK_REALTIME, &time2);

    	// Get difference.
    	timediff = diff(time1, time2);

    	// Save time in seconds to vector.
    	timing_list_bin.push_back(((double)GIG * (double)timediff.tv_sec + (double)timediff.tv_nsec)/((double)GIG));

	  	// Print the number of iterations compared to the average.
	  	//printf("Size %d is done. The interation count was %d\n", j, BinCount);

	  	// Start timer.
    	clock_gettime(CLOCK_REALTIME, &time1);

	  	// Perform a dictionary search. 
	  	i = dictionarySearch(list,find ,Nsize, DctCount); 

	  	// End timer.
    	clock_gettime(CLOCK_REALTIME, &time2);

    	// Get difference.
    	timediff = diff(time1, time2);

    	// Save time in seconds to vector.
    	timing_list_dct.push_back(((double)GIG * (double)timediff.tv_sec + (double)timediff.tv_nsec)/((double)GIG)); 

	  	// Print the number of iterations compared to the average.
	  	//printf("Size %d is done. The iteration count was %d\n", j, DctCount);

	  	// Free the allocated memory.
	  	delete [] list;
	}

	// Print timings.
	  int counter = 0;
	  for (j = 10; j <N; j*=10)
	  {
	    printf("%d %.8e %.8e %.8e\n", j, timing_list_lin[counter], timing_list_bin[counter], timing_list_dct[counter]);
	    counter++;
	  }

	return 0;
}



