// mandelbrot-numerics -- numerical algorithms related to the Mandelbrot set
// Copyright (C) 2015-2018 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

/*

https://www.mrob.com/pub/muency/period.html

Finding the Period of a mu-Atom
If the approximate location and size of a mu-atom are known, its period can be found using the following method, which is similar to the Synchronous-Orbit Algorithm.
1. Identify a polygon with at least three vertices, such that the nucleus is known to be inside the polygon. (In practice, a square is the simplest polygon to use.)
2. Define values of Ci for each vertex, and set initial values of Zi=0. Then, iterate Zi2+Ci for all values of Zi, one step at a time.
3. At each iteration step, count how many edges of the polygon cross the positive half of the real axis. (Actually, either half of either axis will work). If an odd number of edges cross the positive real axis, then we know the origin is inside the polygon. (This method works because of the Jordan Curve Theorem.)
4. Look and see how many iterations you had to do to get to the set of vertices that contains the origin. This is the period.
5. This method will find the lowest period within the starting polygon. Therefore, if the location of the nucleus is known only approximately, the method will still work.




gcc b.c -std=c99 -Wall -Wextra -pedantic -lm


*/



#include <stdio.h> // fprintf
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


static inline int sgn(double z) {
  if (z > 0) { return  1; }
  if (z < 0) { return -1; }
  return 0;
}

static inline bool odd(int a) {
  return a & 1;
}

static inline double cabs2(double _Complex z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

static inline bool cisfinite(double _Complex z) {
  return isfinite(creal(z)) && isfinite(cimag(z));
}

//****************************************************************
//**************** Box period *************************************
//*****************************************************************



static double cross(double _Complex a, double _Complex b) {
  return cimag(a) * creal(b) - creal(a) * cimag(b);
}

static bool crosses_positive_real_axis(double _Complex a, double _Complex b) {
  if (sgn(cimag(a)) != sgn(cimag(b))) {
    double _Complex d = b - a;
    int s = sgn(cimag(d));
    int t = sgn(cross(d, a));
    return s == t;
  }
  return false;
}

static bool surrounds_origin(double _Complex a, double _Complex b, double _Complex c, double _Complex d) {
  return odd
    ( crosses_positive_real_axis(a, b)
    + crosses_positive_real_axis(b, c)
    + crosses_positive_real_axis(c, d)
    + crosses_positive_real_axis(d, a)
    );
}

typedef struct  {
  double _Complex c[4];
  double _Complex z[4];
  int p;
} m_d_box_period ;

m_d_box_period *m_d_box_period_new(double _Complex center, double radius) {
  m_d_box_period *box = (m_d_box_period *) malloc(sizeof(*box));
  if (! box) {
    return 0;
  }
  box->z[0] = box->c[0] = center + ((-radius) + I * (-radius));
  box->z[1] = box->c[1] = center + (( radius) + I * (-radius));
  box->z[2] = box->c[2] = center + (( radius) + I * ( radius));
  box->z[3] = box->c[3] = center + ((-radius) + I * ( radius));
  box->p = 1;
  return box;
}

void m_d_box_period_delete(m_d_box_period *box) {
  if (box) {
    free(box);
  }
}

bool m_d_box_period_step(m_d_box_period *box) {
  if (! box) {
    return false;
  }
  bool ok = true;
  for (int i = 0; i < 4; ++i) {
    box->z[i] = box->z[i] * box->z[i] + box->c[i];
    ok = ok && cisfinite(box->z[i]);
  }
  box->p = box->p + 1;
  return ok;
}

bool m_d_box_period_have_period(const m_d_box_period *box) {
  if (! box) {
    return true;
  }
  return surrounds_origin(box->z[0], box->z[1], box->z[2], box->z[3]);
}

int m_d_box_period_get_period(const m_d_box_period *box) {
  if (! box) {
    return 0;
  }
  return box->p;
}

int m_d_box_period_do(double _Complex center, double radius, int maxperiod) {
  m_d_box_period *box = m_d_box_period_new(center, radius);
  if (! box) {
    return 0;
  }
  int period = 0;
  for (int i = 0; i < maxperiod; ++i) {
    if (m_d_box_period_have_period(box)) {
      period = m_d_box_period_get_period(box);
      break;
    }
    if (! m_d_box_period_step(box)) {
      break;
    }
  }
  m_d_box_period_delete(box);
  return period;
}



//****************************************************************
//**************** Main  *************************************
//*****************************************************************


int main(void){
	
	complex double c = -0.120972062945854+0.643951893407125*I; //
	complex double dc = 1.0;
	int maxiters = 1000;
	
	/* find the period of a nucleus within a large box : radius = 4.0 * cabs(dc)
	uses Robert P. Munafo's Jordan curve method  */
	int p = m_d_box_period_do(c, 4.0 * cabs(dc), maxiters);
	fprintf(stdout, "c = %.16f %+.16f \t period = %d\n", creal(c), cimag(c),p);
	
	 // refine the nucleus location (uses Newton's method)
        if (m_converged == m_d_nucleus(&c0, c0, p, 16))
        
      

	return 0;

}
