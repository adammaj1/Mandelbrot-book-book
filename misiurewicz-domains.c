#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenM

/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/tree/HEAD:/book


gcc misiurewicz-domains.c -lm -Wall -fopenmp

./a.out > md.pgm


cd existing_folder
git init
git remote add origin git@gitlab.com:adammajewski/mandelbrot-book_book.git
git add misiurewicz-domainsc
git commit -m "misiurewicz-domains"
git push -u origin master

*/


static double cabs2(complex double z) {
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
  int maxiters = 4096;
  double er2 = 65536;
  double _Complex c0 = -1.7548776662466929;
  double r = 1.9035515913132437e-02 * 2;
  int period = 1;
  unsigned char *img = malloc(3 * w * h);
  #pragma omp parallel for schedule(static, 1)
  for (int j = 0; j < h; ++j)
  {
    double y = (h/2 - (j + 0.5)) / (h/2);
    for (int i = 0; i < w; ++i)
    {
      double x = ((i + 0.5) - w/2) / (h/2);
      double _Complex c = c0 + r * (x + I * y);
      double dc0 = r / (h/2);

      double _Complex z = c;
      double _Complex dc = 1;
      double _Complex zp = c;
      int pp = 0;
      double mp2 = 1.0 / 0.0;
      double mz2 = 1.0 / 0.0;
      double de = -1;
      for (int n = 0; n < period; ++n)
      {
            double z2 = cabs2(z);
            if (z2 > er2)
            {
              de = sqrt(z2 / cabs2(dc * dc0)) * log(z2);
              break;
            }
            dc = 2 * z * dc + 1;
            z = z * z + c;
      }
      if (de < 0)
      {
            for (int n = 0; n < maxiters - period; ++n) {
              double z2 = cabs2(z);
              if (z2 < mz2) {
                mz2 = z2;
              }
              if (z2 > er2) {
                de = sqrt(z2 / cabs2(dc * dc0)) * log(z2);
                break;
              }
              double p2 = cabs2(z - zp);
              if (p2 < mp2) {
                mp2 = p2;
                pp = n;
              }
              dc = 2 * z * dc + 1;
              z = z * z + c;
              zp = zp * zp + c;
            }
      }
      double hue = 0, sat = 0, val = 1;
      if (de >= 0)
      {
        hue = 1 - fmod(gold * pp, 1);
        sat = 0.25;
        val = tanh(de / aa);
      }
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
