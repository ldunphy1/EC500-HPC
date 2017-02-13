#include <iostream>
#include <iomanip> 
#include <cmath>
#include <float.h>
#include <stdlib.h>
#define TOL 1.0E-15

using std::cout;
using std::cin;
using std::endl;
using std::ios;
using std::abs;

long double newton(long double x, long double A, int n)
{
  return x*(1 - 1.0/n) + A/(n*pow(x,n-1));
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

int main()
{
  // Declare variables to hold the current guess and relative error. 
  long double x = 0.0, fractional_error = 0.5;
  
  // Declare a variable to hold "A". 
  long double A;
  
  // Declare a counter and size
  int count, n; 
  
  // Fix the output format.
  cout.setf(ios::fixed,ios::floatfield); 
  cout.precision(20);
  
  // Describe the problem and prompt the user:
  cout << " Compute the nth root by Newton's Method and the bisection method to a tolerance " << TOL << "." << endl;
  cout <<"Give a number A: ";
  cin >> A;
  cout <<"Give a root n: ";
  cin >> n;
  // Choose an initial guess for Newton's method: in this case, A/2. Set the output precision as well.
  x = A/2;
  count = 0;
  cout.precision(20);
  do 
  {
	  count++;
	  x = newton(x, A, n);
	  fractional_error = 0.5*abs((pow(x,n))/A-1);
  }
  while(fractional_error > TOL && count < 100000);
  cout.precision(40);
  cout << "After " << count << " iterations, Newton's method " << endl;
  cout << "gave = "  << x  <<  " vs cmath = "  << pow(A,(1.0/n)) << endl; 
  
  /*   Compare with bisection */
  cout<< "Bisection Method starting with min = 0 and max = A\n";
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
  cout.precision(40);
  cout << "Bisection's value  in " << count << " iterations " << endl;
  cout << "gave = "  << x  <<  " vs cmath = "  << pow(A,(1.0/n)) << endl;
  return  0;
}
