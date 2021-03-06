#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

double func(double x)
{ return sin(x);}
double derivative(double x)
{ return cos(x); }
double forward_diff(double x, double h)
{ return (func(x+h)-func(x))/h; }
double backward_diff(double x, double h)
{ return (func(x)-func(x-h))/h;}
double central_diff(double x, double h)
{ return (func(x+h/2)-func(x-h/2))/h;}
int main(int argc, char** argv)
{
	double h=1.0;
	const double x = 1.0;
	for (int i = 0; i < 17; i++)
	{ 
		h*=0.1;
		// Print h, the forward, backward, and central difference of sin(x) at 1,
		// as well as the exact derivative cos(x) at 1.
		// Should you use spaces or tabs as delimiters? Does it matter?
		
		cout << setprecision(15) << h << " " << forward_diff(x,h) << " " << backward_diff(x,h) << " " << central_diff(x,h) << " "<< cos(1.0) << endl;
	}
	return 0;
}
