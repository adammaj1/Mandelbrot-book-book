#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenM

/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/blob/HEAD:/book/


gcc atom-domains.c -lm -Wall -fopenmp

./a.out >ad.pgm

*/


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

int main()
{
  double phi = (sqrt(5) + 1) / 2;
  double gold = 1 / (phi * phi);
  int aa = 4;
  int w = 800 * aa;
  int h = 800 * aa;
  int n = 1024;
  double r = 2;
  double r2 = r * r;
  double dc0 = r / (h/2);
  unsigned char *img = malloc(3 * w * h);
  #pragma omp parallel for schedule(dynamic, 1)
  for (int j = 0; j < h; ++j)
  {
    double y = (h/2 - (j + 0.5)) / (h/2) * r;
    for (int i = 0; i < w; ++i)
    {
      double x = (i + 0.5 - w/2) / (h/2) * r;
      double _Complex c = x + I * y;
      double _Complex z = 0;
      double _Complex dc = 0;
      double mz = 1.0 / 0.0;
      int m = 0;
      double de = -1;
      int k;
      for (k = 1; k < n; ++k)
      {
        dc = 2 * z * dc + 1;
        z = z * z + c;
        double z2 = cnorm(z);
        if (z2 < mz)
        {
          mz = z2;
          m = k;
        }
        if (z2 > r2)
        {
          de = 2 * cabs(z) * log(cabs(z)) / cabs(dc * dc0);
          break;
        }
      }
      double hue = 1 - fmod(m * gold, 1);
      double sat = 0.25;
      double val = de >= 0 ? tanh(de / aa) : 1;
      int red, grn, blu;
      hsv2rgb(hue, sat, val, &red, &grn, &blu);
      img[3*(j * w + i)+0] = red;
      img[3*(j * w + i)+1] = grn;
      img[3*(j * w + i)+2] = blu;
    }
  }
  printf("P6\n%d %d\n255\n", w, h);
  fwrite(img, 3 * w * h, 1, stdout);
  free(img);
  return 0;
}
