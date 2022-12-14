#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenM

/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/blob/HEAD:/book/

-----------------------------------------------------
hobold: Let's remove the iteration bands and change the relative sizes of black and white checkerboard areas. Now the structural similarity is more emphasized.
https://fractalforums.org/fractal-mathematics-and-new-theories/28/plotting-field-lines-during-iteration/4233

The white curve in my first image is a polyline connecting the (numerically approximated) corners of the checkers of some specific iteration band. I believe the limit was 27 iterations before the smallest "square" became too tiny to distinguish adjacent vertices with double precision.

The pattern in the 2nd image was done by looking at the final value of the iterated Z_n, just after it escapes. The usual checkerboarded "binary decomposition" looks at just the sign of the imaginary part. But you really can choose any axis through the origin and color based on what side of the axis you end up on. Or you can choose a sector smaller than 180 degrees like I did:
Code: [Select]
return (fabs(z.x)*0.1 < fabs(z.y));

The bailout radius needs to be reasonably large for those field lines to align nicely between iteration bands. Because this analogy with field lines is only strictly true for an infinite bailout. So with a minimal bailout radius of 2.0, the binary decomposition ends up being very visibly distorted and misaligned.
------------------------------------------------------------------------





gcc m.c -lm -Wall -Wextra -fopenmp

./a.out >f.pgm

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
  
  double r = 2; // plane radius: https://en.wikibooks.org/wiki/Fractals/Computer_graphic_techniques/2D/plane#radius
  // The bailout radius needs to be reasonably large for those field lines to align nicely between iteration bands
  double ER = 2000; // Escape Radius = b
  double ER2 = ER * ER; // escape_radius^2
  
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
        if (cnorm(z) > ER2)
          break;
      }
      
      //  fabs(z.x)*0.1 < fabs(z.y)   
      img[j * w + i] = (k < kMax && fabs(creal(z))*0.1 < fabs(cimag(z))) ? 255 : 0;
    }
  }
  printf("P5\n%d %d\n255\n", w, h);
  fwrite(img, w * h, 1, stdout);
  free(img);
  return 0;
}

