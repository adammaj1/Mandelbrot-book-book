/* 

mandelbrot-numerics -- numerical algorithms related to the Mandelbrot set
Copyright (C) 2015-2018 Claude Heiland-Allen
License GPL3+ http://www.gnu.org/licenses/gpl.html

gcc m_d_nucleus.c -lm -Wall -Wextra 
./a.out



*/
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

 
 
const double pi = 3.141592653589793;
const double twopi = 6.283185307179586;

// 
const double epsilon = 4.440892098500626e-16;

// epsilon^2
static const double epsilon2 = 1.9721522630525295e-31;
 
 
 enum m_newton { m_failed, m_stepped, m_converged };
typedef enum m_newton m_newton;
 
int sgn(double z) {
  if (z > 0) { return  1; }
  if (z < 0) { return -1; }
  return 0;
}

bool odd(int a) {
  return a & 1;
}

double cabs2(double _Complex z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

bool cisfinite(double _Complex z) {
  return isfinite(creal(z)) && isfinite(cimag(z));
}

 // ******************************************************************************************
// **************************** Naive ****************************************************************
// *****************************************************************************************


m_newton m_d_nucleus_naive_step(double _Complex *c_out, double _Complex c_guess, int period) {
  double _Complex z = 0;
  double _Complex dc = 0;
  for (int i = 0; i < period; ++i) {
    dc = 2 * z * dc + 1;
    z = z * z + c_guess;
  }
#if 0
  if (cabs2(dc) <= epsilon2) {
    *c_out = c_guess;
    return m_converged;
  }
#endif
  double _Complex c_new = c_guess - z / dc;
  double _Complex d = c_new - c_guess;
  if (cabs2(d) <= epsilon2) {
    *c_out = c_new;
    return m_converged;
  }
  if (cisfinite(d)) {
    *c_out = c_new;
    return m_stepped;
  } else {
    *c_out = c_guess;
    return m_failed;
  }
}

m_newton m_d_nucleus_naive(double _Complex *c_out, double _Complex c_guess, int period, int maxsteps) {
  m_newton result = m_failed;
  double _Complex c = c_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_nucleus_naive_step(&c, c, period))) {
      break;
    }
  }
  *c_out = c;
  return result;
}

// ******************************************************************************************
// ********************************************************************************************
// *****************************************************************************************


m_newton m_d_nucleus_step(double _Complex *c_out, double _Complex c_guess, int period) {
  double _Complex z = 0;
  double _Complex dc = 0;
  double _Complex h = 1;
  double _Complex dh = 0;
  for (int i = 1; i <= period; ++i) {
    dc = 2 * z * dc + 1;
    z = z * z + c_guess;
    // reject lower periods
    if (i < period && period % i == 0)
    {
      h = h * z;
      dh = dh + dc / z;
    }
  }
  // build function
  dh = dh * h;
  double _Complex g = z;
  double _Complex dg = dc;
  double _Complex f = g / h;
  double _Complex df = (dg * h - g * dh) / (h * h);
  // newton step
  double _Complex c_new = c_guess - f / df;
  // check convergence
  double _Complex d = c_new - c_guess;
  if (cabs2(d) <= epsilon2) {
    *c_out = c_new;
    return m_converged;
  }
  if (cisfinite(d)) {
    *c_out = c_new;
    return m_stepped;
  } else {
    *c_out = c_guess;
    return m_failed;
  }
}

m_newton m_d_nucleus(double _Complex *c_out, double _Complex c_guess, int period, int maxsteps) {
  m_newton result = m_failed;
  double _Complex c = c_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_nucleus_step(&c, c, period))) {
      break;
    }
  }
  *c_out = c;
  return result;
}

// ******************************************************************************************
// ************************   Main  ********************************************************************
// *****************************************************************************************


int main(void){

	int period = 1;
	double _Complex c0 = 0.0; // A reasonable starting guess for Newton's method is within the atom domain of the component.
	double complex nucleus; // output
	
	m_newton result = m_d_nucleus( &nucleus, c0, period, 100);
	if (result == m_converged )
		{ fprintf(stdout, "nucleus = %.16f %+.16f period = %d\n ", creal(nucleus), cimag(nucleus), period); }
		else {fprintf(stdout, "Newton method failed = %d\n ", result);}
		
return 0;
}
