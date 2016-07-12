#include <stdio.h>
#include <math.h>
#include "functionals.h"

#define DIFF_STEP 1e-6

double differential(double_in_double f, double x)
{
    double h = DIFF_STEP;
    return (f(x - 2*h) - 8*f(x - h) + 8*f(x + h) - f(x + 2*h)) / 12 / h;    //4-й порядок точности
    //return (f(x + h) - f(x - h)) / 2 / h;                                   //2-й порядок точности
}

double simpson(double_in_double f, double left, double right, int k)
{
    int i;
    double h = (right - left) / k;
    double Sf = f(left) + f(right);     //Sf = f(x1) + 4*f(x2) + 2*f(x3) + ... + 4*f(x<2k-1>) + f(x<2k>)
    for(i = 1; i < 2*k; i++)
        Sf += (i % 2 + 1) * 2 * f(i * h / 2);
    return Sf * h / 6;
}

double golden_section(double_in_double q, double left, double right, double epsilon)
{
	double a, b;		        //золотые сечения left < a < b < right
	double qa, qb;	            //значение функции q(a) q(b)
	double psi;			        //psi=1/phi, где phi золотая пропорция
	double x;			        //точка минимума
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
