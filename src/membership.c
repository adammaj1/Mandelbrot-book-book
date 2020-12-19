#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenM

/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/tree/HEAD:/book


gcc membership.c -lm -Wall -fopenmp

./a.out > m.pgm


cd existing_folder
git init
git remote add origin git@gitlab.com:adammajewski/mandelbrot-book_book.git
git add membership.c
git commit -m "Initial commit"
git push -u origin master

*/



 double cnorm(double _Complex z)
{
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

int main()
{
  int aa = 4;
  int w = 800 * aa;
  int h = 800 * aa;
  int n = 1024;
  double r = 2;
  double r2 = r * r;
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
      img[j * w + i] = (k == n) ? 0 : 255;
    }
  }
  printf("P5\n%d %d\n255\n", w, h);
  fwrite(img, w * h, 1, stdout);
  free(img);
  return 0;
}
