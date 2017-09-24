#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenMP
/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/blob/HEAD:/book/


gcc exterior-coordinates.c -lm -Wall -fopenmp

./a.out >ec.pgm

*/


const double pi = 3.141592653589793;

double cnorm(double _Complex z)
{
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

double cabs0(double _Complex z)
{
  return fabs(creal(z)) + fabs(cimag(z));
}

int main()
{
  int aa = 4;
  int w = 800 * aa;
  int h = 800 * aa;
  int n = 1024;
  double r = 2;
  double r2 = 25 * 25 * 25 * 25;
  unsigned char *img = malloc(w * h);
  #pragma omp parallel for
  for (int j = 0; j < h; ++j)
  {
    double y = (h/2 - (j + 0.5)) / (h/2) * r;
    for (int i = 0; i < w; ++i)
    {
      double x = (i + 0.5 - w/2) / (h/2) * r;
      double _Complex c = x + I * y;
      double _Complex z = 0;
      int k;
      for (k = 0; k < n; ++k)
      {
        z = z * z + c;
        if (cnorm(z) > r2)
          break;
      }
      double g = 1;
      if (k < n)
      {
        double _Complex e = fmod(1 + carg(z) / (2 * pi), 1) + I * (2 - log(cnorm(z))/log(r2));
        if (cabs(e - (0.25 + 0.25 * I)) < 0.2) g = 0;
        if (cabs0(e - (0.66 + 0.66 * I)) < 0.3) g = 0.67;
        if (cabs(e - (0.17 + 0.75 * I)) < 0.1) g = 0.33;
      }
      img[j * w + i] = 255 * g;
    }
  }
  printf("P5\n%d %d\n255\n", w, h);
  fwrite(img, w * h, 1, stdout);
  free(img);
  return 0;
}
