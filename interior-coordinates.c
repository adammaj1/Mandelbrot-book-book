#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenMP
/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/blob/HEAD:/book/binary-decomposition.c


gcc interior-coordinates.c -lm -Wall -fopenmp

./a.out >ic.pgm

*/




const double pi = 3.141592653589793;

double cnorm(double _Complex z)
{
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

void hsv2rgb(double h, double s, double v, int *red, int *grn, int *blu) {
  double i, f, p, q, t, r, g, b;
  int ii;
  if (s == 0.0) { r = g = b = v; } else {
    h = 6 * (h - floor(h));
    ii = i = floor(h);
    f = h - i;
    p = v * (1 - s);
    q = v * (1 - (s * f));
    t = v * (1 - (s * (1 - f)));
    switch(ii) {
      case 0: r = v; g = t; b = p; break;
      case 1: r = q; g = v; b = p; break;
      case 2: r = p; g = v; b = t; break;
      case 3: r = p; g = q; b = v; break;
      case 4: r = t; g = p; b = v; break;
      default:r = v; g = p; b = q; break;
    }
  }
  *red = fmin(fmax(255 * r + 0.5, 0), 255);
  *grn = fmin(fmax(255 * g + 0.5, 0), 255);
  *blu = fmin(fmax(255 * b + 0.5, 0), 255);
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







int main()
{
  int aa = 4;
  int w = 800 * aa;
  int h = 800 * aa;
  int n = 1024;
  double r = 2;
  double r2 = 25 * 25 * 25 * 25;
  double dc0 = r / (h/2);
  unsigned char *img = malloc(3 * w * h);
  #pragma omp parallel for schedule(dynamic, 1)
  for (int j = 0; j < h; ++j)
  {
    double y = (h/2 - (j + 0.5)) / (h/2) * r;
    for (int i = 0; i < w; ++i)
    {
      double x = (i + 1.0 - w/2) / (h/2) * r;
      double _Complex c = x + I * y;
      double _Complex z = 0;
      double _Complex dc = 0;
      double mz = 1.0 / 0.0;
      double _Complex b = 0;
      double de = -1;
      int k;
      for (k = 1; k < n; ++k)
      {
        dc = 2 * z * dc + 1;
        z = z * z + c;
        double zp = cabs(z);
        if (zp < mz)
        {
          mz = zp;
          double _Complex w = m_attractor(z, c, k, 16);
          double _Complex dw = 1;
          for (int m = 0; m < k; ++m)
          {
            dw = 2 * w * dw;
            w = w * w + c;
          }
          if (cabs(dw) <= 1)
          {
            b = dw;
            break;
          }
        }
        if (cnorm(z) > r2)
        {
          de = 2 * cabs(z) * log(cabs(z)) / cabs(dc * dc0);
          break;
        }
      }
      double hue = 0;
      double sat = 0;
      double val = 1;
      if (k < n && de < 0)
      {
        hue = carg(b) / (2 * pi);
        sat = cnorm(b);
        val = 1;
      }
      if (de >= 0)
      {
        hue = 0;
        sat = 0;
        val = tanh(de / aa);
      }
      int red, grn, blu;
      hsv2rgb(hue, sat, val, &red, &grn, &blu);
      img[3 * (j * w + i) + 0] = red;
      img[3 * (j * w + i) + 1] = grn;
      img[3 * (j * w + i) + 2] = blu;
    }
  }
  printf("P6\n%d %d\n255\n", w, h);
  fwrite(img, 3 * w * h, 1, stdout);
  free(img);
  return 0;
}
