#include <iostream>
#include <cmath>
#include "matrix.h"
#include "functionals.h"

#define GOLD_EPS 1e-9           ///accuracy in one dimension optimization (optimal when >= 2e-9)
#define DIFF_STEP 1e-5          ///numerical differantiation step (>= 1e-6, otherwise observational error increases)
#define EPSILON 5e-9            ///target accuracy
#define UNIM_STEP 0.1           ///first step in s_unimodal funciton
#define TAX 1e4                 ///penalty coefficient (if 1/TAX > 2*DIFF_STEP numerical derivative in boundary become wrong)

#define _X_0_ -2.0, -2.0        ///fist point (be sure that it is double, double)
#define DEBUG (0)               ///debug \ production mode

/// global variables
vector x_k(2);
matrix D(2, 2);

/// exterior penalty function with constraints (g[i] <= 0)
double penalty(vector x)
{
    double pnt = 0;
    pnt += pow(max(0, x[0]*x[0] + 2*(x[1]-2) - 8), 2);                  //first constraint
    pnt += pow(max(0, 1 - x[0]*x[0] - (x[1]-2)*(x[1]-2)), 2);           //second coinstraint
    return TAX * pnt;
}

/// target function
double f(vector x)
{
    return x[0]*x[0] + x[1]*x[1] + penalty(x);                          //paraboloid with penalty funciton
    //return pow(x[0] + 1, 2) + pow(x[1] - 1, 2) + penalty(x);          //shifted paraboloid with penalty
    //return 100 * pow(x[1] - x[0]*x[0], 2) + pow(1 - x[0], 2);         //Rosenbrock function
    //return - (pow(x[0]*x[0]+x[1]-11, 2) + pow(x[0]+x[1]*x[1]-7, 2));  //Himmelblau function
    //return exp(x[0]*x[0] +  2*x[1]*x[1]) - x[0] - x[1] + penalty(x);  //task 1
}

/// gradient of target function with numerical differantiation
vector grad_f(vector x)
{
    double h = DIFF_STEP;
    vector g(2);
    vector dx0(2);
    vector dx1(2);
    init(dx0, h, 0.0);
    init(dx1, 0.0, h);
    g[0] = (f(x - 2*dx0) - 8*f(x - dx0) + 8*f(x + dx0) - f(x + 2*dx0)) / 12 / h;
    g[1] = (f(x - 2*dx1) - 8*f(x - dx1) + 8*f(x + dx1) - f(x + 2*dx1)) / 12 / h;
    //g[0] = (f(x + dx0) - f(x - dx0)) / 2 / h;
    //g[1] = (f(x + dx1) - f(x - dx1)) / 2 / h;
    return g;
}

/// function that used by golden section search and return next f(x<k+1>)
inline double f_x_kk(double s)
{
    extern vector x_k;
    extern matrix D;
    return f(x_k - s*D*grad_f(x_k));
}

/// function that simply increases the interval as long as f(x) decreases and return unimodal interval
void s_unimodal(double& s_min, double& s_max)
{
    extern vector x_k;
    double step = UNIM_STEP / norm(grad_f(x_k))
    double s0 = 0;
    double s1 = step;
    double s2 = 2 * step;
    if (f_x_kk(s0) > f_x_kk(s1)) {
        for(int i = 0; f_x_kk(s1) > f_x_kk(s2) && i < 1000; i++) {
            s0 = s1;
            s1 = s2;
            s2 *= 2;
        }
    }
    s_min = s0;
    s_max = s2;
}

int main()
{
    extern vector x_k;                                                  ///current point
    extern matrix D;                                                    ///oproximated Hesse matrix in BFGS (D<0> == E)
    int k;
    matrix A(2, 2);
    matrix B(2, 2);
    vector x_kk(2);
    vector u(2);
    vector v(2);
    double s;                                                           ///variable in one dimension optimization
    double s_min, s_max;

    k = 1;
    x_k = vector(2);
    init(x_k, _X_0_);
    D = identity(2);

/// 1 step
    s_unimodal(s_min, s_max);
    s = golden_section(f_x_kk, s_min, s_max, GOLD_EPS);
	x_kk = x_k - s*D*grad_f(x_k);
	printf("#       X               Y               Z        n\n");     ///print format
    printf("%15.12lf %15.12lf %15.12lf  %d\n", x_k[0], x_k[1], f(x_k), k-1);
    printf("%15.12lf %15.12lf %15.12lf  %d\n", x_kk[0], x_kk[1], f(x_kk), k);

/// 2..n step
    do
    {
        u = x_kk - x_k;                                                 ///u<k> = x<k+1> - x<k>
        v = grad_f(x_kk) - grad_f(x_k);                                 ///v<k>
        A = u * transpose(u) / (transpose(u) * v);                      ///A<k> = u * u^t / (u^t * v)
        B = (D * v * transpose(v) * D) / (- (transpose(v) * D * v));    ///B<k> = -(D * v * v^t * D) / (v^t * D * v)
        D += A + B;                                                     ///D<k+1>
        k++;
        x_k = x_kk;                                                     
        s_unimodal(s_min, s_max);
        s = golden_section(f_x_kk, s_min, s_max, GOLD_EPS);
        x_kk = x_k - s*D*grad_f(x_k);                                   ///x<k+1>
        printf("%15.12lf %15.12lf %15.12lf  %d\n", x_kk[0], x_kk[1], f(x_kk), k);
    } while (norm(x_kk - x_k) > EPSILON);

    putchar('\n');
    return 0;
}
