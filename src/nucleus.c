/*

gcc nucleus.c -lm -Wall -Wextra 
*/

#include <stdio.h>
#include <math.h>
#include <complex.h>

double complex m_nucleus(double complex c0, int period, int mMax)
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
