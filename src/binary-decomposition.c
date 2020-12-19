#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenM

/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/blob/HEAD:/book/


gcc binary-decomposition.c -lm -Wall -fopenmp

./a.out >bd.pgm

*/


 double cnorm(double _Complex z)
{
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

int main()
{
  //int aa = 4;
  int w = 2000 ; // width in piels
  int h = 2000 ; // height in pixels
  int kMax = 1024; // iteration max
  
  double r = 2; // radius: https://en.wikibooks.org/wiki/Fractals/Computer_graphic_techniques/2D/plane#radius
  
  double r2 = 25 * 25; // escape_radius^2
  
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
      for (k = 0; k < kMax; ++k)
      {
        z = z * z + c;
        if (cnorm(z) > r2)
          break;
      }
      img[j * w + i] = (k < kMax && cimag(z) < 0) ? 0 : 255;
    }
  }
  printf("P5\n%d %d\n255\n", w, h);
  fwrite(img, w * h, 1, stdout);
  free(img);
  return 0;
}
