#include <iostream>
#include <iomanip> 
#include <cmath>
#include <float.h> // For long doubles. 

using std::cout;
using std::cin;
using std::endl;
using std::ios;

using std::abs;

// With finite precision, we can't really "exactly" find
// the zero. This define sets the tolerance: how small we 
// want the relative error between the true zero and our
// approximate zero to be. 
//#define TOL 1.0E-15

// This function performs one iteration of Newton's method
// and returns a new guess (x - f(x)/f'(x) -> x_new).
// For now, you need to hard-code the numerical function and
// its derivative. 
long double newton(long double x, long double A, int n);

// This function performs one iteration of bisection
// and updates either "min" or "max" (note how they are both
// passed by reference), and returns the current "midpoint".
// Again, you need to hard-code the numerical function. Bisection
// does not require a derivative. 
long double bisection(long double A, long double & min, long double & max, int n);

int main()
{
  // Declare variables to hold the current guess
  // and relative error. 
  long double x = 0.0, fractional_error = 0.5;
  
  // Declare a variable to hold "A". 
  long double A = 2;
  
  // Declare a counter.
  int count,n;
  long double N,TOL;
  
  // Fix the output format.
  for(int j=1;j<=15;j++)
  {
	  N = pow(10,j);
	  TOL = 1/N;
	  cout<<N<<" ";
	  for(int i=2;i<=3;i++)
	  {
		  n=i;
		  long double min, max;
		  min = 0.0;
		  max = A;
		  count = 0;

		  do 
		  {
			count++;
			x =  bisection( A,  min, max, n);
			fractional_error = 0.5*abs((pow(x,n))/A-1);
		  }
		  while(fractional_error > TOL);
		  
		  cout<<count<<" ";
	  }
	  cout<<endl;
  }
  cout<<endl;
  return  0;
}

long double bisection(long double A, long double & min, long double & max, int n)
{
  long double x  = (min + max)/2.0;
  if((pow(x,n)-A) < 0.0)
    min = x;
  else
    max = x;

  return x;
}