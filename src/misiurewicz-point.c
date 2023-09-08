/*


====================================
A preperiodic Misiurewicz point c of preperiod q  and period p  = M_{q,p} satisfies:

F^{q+p}(0,c) = F^q(0,c)

----------------------------------

gcc misiurewicz-point.c -lm -Wall -Wextra


*/
#include <stdio.h> // fprintf
#include <math.h>
#include <complex.h>


/*
input
* aproximation c0
* preperiod q  of critical point
* period p 
* maximal number of iterations n

output c ( Misiurewicz point  M_{q,p})

A naive implementation of Newton’s root finding method iterations for a Misiurewicz point takes the form:
https://mathr.co.uk/web/m-misiurewicz-naive.html

*/
double _Complex m_misiurewicz_naive
    (double _Complex c0, int q, int p, int n)
{
    double _Complex c = c0;
    for (int m = 0; m < n; ++m)
    {
        double _Complex z = 0;
        double _Complex dc = 0;
        double _Complex zp = 0;
        double _Complex dcp = 0;
        for (int i = 0; i < q + p; ++i)
        {
            if (i == q)
            {
                zp = z;
                dcp = dc;
            }
            dc = 2 * z * dc + 1;
            z = z * z + c;
        }
        c = c - (z - zp) / (dc - dcp);
    }
    return c;
}


/*
input
* aproximation c0
* preperiod q  of critical point
* period p 
* maximal number of iterations n

output c ( Misiurewicz point  M_{q,p})

A full implementation of Newton’s root finding method iterations for a Misiurewicz point takes the form:
https://mathr.co.uk/web/m-misiurewicz-full.html

Applying Newton’s root finding method iterations to G(c)=0  finds the preperiodic Misiurewicz point, 
with a larger basin of attraction than the “naive” method. 
However it may not converge as precisely, one approach is to first use this “full” method and then use the “naive” method on its output.


*/

double _Complex m_misiurewicz_full
    (double _Complex c0, int q, int p, int n)
{
    double _Complex c = c0;
    for (int m = 0; m < n; ++m)
    {
        double _Complex z = 0;
        double _Complex dc = 0;
        double _Complex zp = 0;
        double _Complex dcp = 0;
        double _Complex h = 1;
        double _Complex dh = 0;
        for (int i = 0; i < p; ++i)
        {
            dc = 2 * z * dc + 1;
            z = z * z + c;
        }
        for (int i = 0; i < q; ++i)
        {
            double _Complex k = z - zp;
            h = h * k;
            dh = dh + (dc - dcp) / k;
            dc = 2 * z * dc + 1;
            z = z * z + c;
            dcp = 2 * zp * dcp + 1;
            zp = zp * zp + c;
        }
        dh = dh * h;
        double _Complex g = z - zp;
        double _Complex dg = dc - dcp;
        double _Complex f = g / h;
        double _Complex df = (dg * h - g * dh) / (h * h);
        c = c - f / df;
    }
    return c;
}

int main(void){

	double _Complex c;
	double _Complex c0 = 1.1*I;
	int q = 2; // preperiod of critical point
	int p = 2; 
	int n = 100000;
	
	c = m_misiurewicz_naive(c0, q, p, n);
	fprintf(stdout, "Misiurewicz point M_{%d,%d} = c = %.16f%+.16f (naive)\n", q, p, creal(c), cimag(c));
	
	c = m_misiurewicz_full(c0, q, p, n);
	fprintf(stdout, "Misiurewicz point M_{%d,%d} = c = %.16f%+.16f (full)\n", q, p, creal(c), cimag(c));
	

	return 0;
	
}
