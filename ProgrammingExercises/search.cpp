/* RCB and ESW, 2016-02-05, Comparing different search algorithms.
   Various codes are borrowed from the web.
   Updated 2017-01-30. */

#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>

#include "search.h"
#include "sort.h"

using std::cout; 
using std::endl;

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

	// Counters 
	int BinCount; // counter for binary search. 
	int LinCount; // counter for linear search. 
	int DctCount; // counter for dictionary search.
	int time_seconds = time(0);
	srand(time_seconds%100); // seed the random number generator.
	int N = 10000000;
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
	  int find = list[rand()%Nsize] ;

	  // First, search with a linear search. A linear search doesn't depend
	  // on a list being sorted.
	  cout << Nsize<<" ";
	  i = linearSearch(list,find ,Nsize, LinCount);

	  // Print the number of iterations compared to the average.
	  cout <<LinCount <<" ";

	  // Next, let's look at a binary search and a dictionary search.
	  // This requires a sorted list. Feel free to play with these if
	  // you want. 
	  // WARNING: the selection sort takes a long time if Nsize gets
	  // into the millions... you probably have better things to do.
	  
	  mergeSort(list, Nsize);  

	  //Perform a binary search. 
	  i = binarySearch(list,find ,Nsize, BinCount); 

	  // Print the number of iterations compared to the average.
	  cout <<BinCount <<" ";

	  // Perform a dictionary search. 
	  i = dictionarySearch(list,find ,Nsize, DctCount);  

	  // Print the number of iterations compared to the average.
	  cout<<DctCount <<endl; 

	  // Free the allocated memory.
	  delete [] list;
	}
  return 0;
}



