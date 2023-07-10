/*

gcc nucleus.c -lm -Wall -Wextra 

Computing c just by iterating fc  could take many 1000s of iterations, especially when multiplier λ  is close to 1 = point is near boundary of componennt

You may need very small epsilon and very large n, otherwise for example c = -3/4+10^{-10}  will probably give an incorrect period of 2  instead of the correct period of 1
, which error will compound to an incorrect interior distance estimate (for example distance 3.8e-8 with your method (epsilon 1e-12, n 79,573,343
) instead of 2e-10 with my method).



*/

#include <stdio.h>
#include <math.h>
#include <complex.h>


/* 

input
c0 = a reasonable starting guess  for Newton's method
period 
mMax = maximal nuimbetr of Newton iterations ( steps)

output: c = nucleus

*/
double complex m_nucleus(const double complex c0, const int period, const int mMax)
{
	double complex c = c0;
	for (int m = 0; m < mMax; ++m)
	{
		double complex z = 0;
		double complex dc = 0;
		for (int i = 0; i < period; ++i)
		{
			dc = 2 * z * dc + 1;
			z = z * z + c;
		}
		c = c - z / dc;
	}
	return c;
}


int main(void){

	int period = 1;
	double complex nucleus = m_nucleus(0.0, period, 100);
	fprintf(stdout, "nucleus = %.16f %+.16f period = %d\n ", creal(nucleus), cimag(nucleus), period);




	return 0;
}
