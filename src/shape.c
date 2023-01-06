/*
   1 // mandelbrot-numerics -- numerical algorithms related to the Mandelbrot set
   2 // Copyright (C) 2015-2018 Claude Heiland-Allen
   3 // License GPL3+ http://www.gnu.org/licenses/gpl.html


*/

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>


#define kMax 17 // number of examples, see plane_examples



// plane_center_x	plane_center_y	plane_radius	componenet_period
double plane_examples[kMax][4] = {
	{0.0,		+0.0,		0.8,		1}, 
	{+0.2925755,	-0.0149977, 	0.00025,	32}, 
	{-1.763,  	+0.0,		0.016,		3}, 
	{-0.15842, 	+1.03335, 	0.01,		4},  
	{+0.358431,	+ 0.643507,	0.006,		5},  
	{+0.442990,	+0.373727,	0.005,		6}, 
	{+0.432259,	+0.227315,	0.003,		7}, 
	{+0.404879,	+0.146216,	0.002,		8}, 
	{+0.378631,	+0.098841,	0.001,		9}, 
	{+0.356854, 	+0.069659,	0.001,		10},
	{+0.339454,	+0.050823,	0.001,		11},
	{+0.325631,	+0.038164,	0.001,		12},
	{-1.0,		0.0,		1.0,		2},
	{-0.122561166876654, +0.744861766619744, 0.5,    3},
	{0.282271390766914, +0.530060617578525, 0.5, 4},
	{-1.772892903381624, +0.0,  0.2,  6},
	{-1.757783060083098,  +0.013796143369485, 0.2, 9}
	
	 
};

enum m_shape { m_cardioid, m_circle };
  typedef enum m_shape m_shape;
  
  
double _Complex m_d_shape_estimate(double _Complex c, int p)
{
	double _Complex z = c;
	double _Complex dc = 1;
	double _Complex dz = 1;
	double _Complex dcdc = 0;
	double _Complex dcdz = 0;
	for (int i = 1; i < p; ++i)
	{
		dcdc = 2 * (z * dcdc + dc * dc);
		dcdz = 2 * (z * dcdz + dc * dz);
		dc = 2 * z * dc + 1;
		dz = 2 * z * dz;
		z = z * z + c;
	}
	double _Complex s = - (dcdc / (2 * dc) + dcdz / dz) / (dc * dz);
	return s;
}

m_shape m_d_shape_discriminant(double _Complex shape) {
  if (cabs(shape) < cabs(shape - 1)) {
    return m_cardioid;
  } else {
    return m_circle;
  }
}

m_shape m_d_shape(double _Complex nucleus, int period) {
  return m_d_shape_discriminant(m_d_shape_estimate(nucleus, period));
}


int CheckShape(int k){

	complex double nucleus = plane_examples[k][0] + I*plane_examples[k][1];
	int period = (int) plane_examples[k][3];
	
	m_shape shape;
	shape = m_d_shape(nucleus, period);
	
	printf("nucleus = %.16e %+.16e\t period = %d\t", creal(nucleus), cimag(nucleus), period);
	switch (shape) {
      		case m_cardioid: printf("pseudocardioid %.1f %.1f\n", creal(shape), cimag(shape)); return 0;
      		case m_circle:   printf("pseudocircle %.1f %.1f\n",   creal(shape), cimag(shape)); return 0;
    	}
	return 0;

}



int main(void){



int k =0; // 
  
	
	
	for (k=0; k<kMax; ++k) {
		CheckShape(k);
	}




}
