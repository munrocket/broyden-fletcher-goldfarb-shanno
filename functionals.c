#include <stdio.h>
#include <math.h>
#include "functionals.h"

#define DIFF_STEP 1e-6

//numerical differantiation 
double differential(double_in_double f, double x)
{
    double h = DIFF_STEP;
    return (f(x - 2*h) - 8*f(x - h) + 8*f(x + h) - f(x + 2*h)) / 12 / h;    //O(h^4) approximation
    //return (f(x + h) - f(x - h)) / 2 / h;                                   //O(h^2) approximation
}

//numerical integration with O(h^4)
double simpson(double_in_double f, double left, double right, int k)
{
    int i;
    double h = (right - left) / k;
    double Sf = f(left) + f(right); 
    for(i = 1; i < 2*k; i++)
        Sf += (i % 2 + 1) * 2 * f(i * h / 2);
    return Sf * h / 6;
}

//golden section search
double golden_section(double_in_double q, double left, double right, double epsilon)
{
	double a, b;		        //left < a < b < right
	double qa, qb;
	double psi;			        //psi = 1 / phi
	double x;			        
	int i;
	psi = 2/(1 + sqrt(5));
	b = left + psi * (right - left);
	a = right - psi * (right - left);
	qa = q(a);
	qb = q(b);
	for (i = 0; b - a > epsilon; i++)
	{
		if (qa > qb)
		{
			left = a;
			a = b;
			b = left + psi * (right - left);
			qa = qb;
			qb = q(b);
		}
		else
		{
			right = b;
			b = a;
			a = right - psi * (right - left);
			qb = qa;
			qa = q(a);
		}
	}
	x = (b + a) / 2;
	return x;
}

double max(double a, double b)
{
    return (a > b) ? a : b;
}

double min(double a, double b)
{
    return (a < b) ? a : b;
}
#undef DEBUG
