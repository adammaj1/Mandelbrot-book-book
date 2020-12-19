#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenM

/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/blob/HEAD:/book/


gcc interior-distance.c -lm -Wall -fopenmp

./a.out > id.pgm

*/

const double pi = 3.141592653589793;

 double cnorm(double _Complex z)
{
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

double m_interior_distance(double _Complex z0, double _Complex c, int p) {
  double _Complex z = z0;
  double _Complex dz = 1;
  double _Complex dzdz = 0;
  double _Complex dc = 0;
  double _Complex dcdz = 0;
  for (int m = 0; m < p; ++m)
  {
    dcdz = 2 * (z * dcdz + dz * dc);
    dc = 2 * z * dc + 1;
    dzdz = 2 * (dz * dz + z * dzdz);
    dz = 2 * z * dz;
    z = z * z + c;
  }
  return (1 - cnorm(dz)) / cabs(dcdz + dzdz * dc / (1 - dz));
}

double _Complex m_attractor(double _Complex w0, double _Complex c, int p, int n)
{
  double _Complex w = w0;
  for (int m = 0; m < n; ++m)
  {
    double _Complex z = w;
    double _Complex dz = 1;
    for (int i = 0; i < p; ++i)
    {
      dz = 2 * z * dz;
      z = z * z + c;
    }
    w = w - (z - w) / (dz - 1);
  }
  return w;
}

double m_distance(int N, double R, double _Complex c)
{
  double _Complex dc = 0;
  double _Complex z = 0;
  double m = 1.0 / 0.0;
  int p = 0;
  for (int n = 1; n <= N; ++n)
  {
    dc = 2 * z * dc + 1;
    z = z * z + c;
    if (cabs(z) > R)
      return 2 * cabs(z) * log(cabs(z)) / cabs(dc);
    if (cabs(z) < m)
    {
      m = cabs(z);
      p = n;
      double _Complex z0 = m_attractor(z, c, p, 64);
      double _Complex w = z0;
      double _Complex dw = 1;
      for (int k = 0; k < p; ++k)
      {
        dw = 2 * w * dw;
        w = w * w + c;
      }
      if (cabs(dw) <= 1)
        return m_interior_distance(z0, c, p);
    }
  }
  return 0;
}

int main()
{
  int aa = 4;
  int w = 800 * aa;
  int h = 800 * aa;
  int n = 1024;
  double r = 2;
  double px = r / (h/2);
  double r2 = 25 * 25;
  unsigned char *img = malloc(w * h);
  #pragma omp parallel for schedule(static, 1)
  for (int j = 0; j < h; ++j)
  {
    double y = (h/2 - (j + 0.5)) / (h/2) * r;
    for (int i = 0; i < w; ++i)
    {
      double x = (i + 0.5 - w/2) / (h/2) * r;
      double _Complex c = x + I * y;
      double de = m_distance(n, r2, c);
      double g = tanh(de / px / aa);
      img[j * w + i] = 255 * g;
    }
  }
  printf("P5\n%d %d\n255\n", w, h);
  fwrite(img, w * h, 1, stdout);
  free(img);
  return 0;
}
