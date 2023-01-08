/*

https://mathr.co.uk/mandelbrot/mu-atom/

mu-atom mapping
about
period p hyperbolic components of the Mandelbrot set can each be mapped conformally to the unit disc, by the derivative d/dz of the periodic limit cycle where f_c^p(z_0) = z_0.



gcc mu-atom.c -std=c99 -Wall -Wextra -pedantic -O3  -lm
./a.out 1 > mu-atom-1.ppm 



printf "run the compiled program\n"
for k in {1..10}
do
  ./a.out "$k" > "$k".ppm 
done


*/


#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



double cnorm(double _Complex z)
{
  double x = creal(z);
  double y = cimag(z);
  return x * x + y * y;
}

#define W 640
#define H 360
#define N 1024

int main(int argc, char **argv)
{
  if (argc != 2) {
        fprintf(stderr, "Usage: ./a.out period > mu-atom-period.ppm\n");
        exit(EXIT_FAILURE);
    }

  
  
  double _Complex c0 = -0.75;
  double r0 = 3;
 
  int period = (int) strtol(argv[1], NULL, 10);
  double ER2 = 4;
  
  
  
  printf("P6\n%d %d\n255\n", W, H);
  for (int j = 0; j < H; ++j)
  for (int i = 0; i < W; ++i)
  {
    double x = ((i + 0.5) / W - 0.5) * W / H;
    double y = ((j + 0.5) / H - 0.5);
    double _Complex c = c0 + r0 * (x + I * y);
    double _Complex z = 0;
    double _Complex dc = 0;
    double M = 1.0 / 0.0;
    double DE = -1;
    int highlight = 0;
    for (int k = 1; k < N; ++k)
    {
        dc = 2 * z * dc + 1;
	z = z * z + c;
        if (cnorm(z) < M)
        {
          M = cnorm(z);
          double _Complex w = z;
          double _Complex du;
          for (int l = 0; l < 30; ++l)
          {
            double _Complex u = w;
            du = 1;
            for (int m = 0; m < k; ++m)
            {
              du = 2 * u * du;
              u = u * u + c;
            }
            w -= (u - w) / (du - 1);
          }
          if (cnorm(du) < 1)
          {
            DE = 0;
            highlight = k == period;
            break;
          }
        }
        if (cnorm(z) > ER2)
        {
          DE = sqrt(cnorm(z)/ER2) * log2(2*cnorm(z)/ER2) / (cabs(dc) * r0 / H);
          break;
        }
    }
 
    double MU = tanh(fmax(DE, 0));
    double red[3] = { 1, 0, 0 };
    // color of the pixel 
    for (int c = 0; c < 3; ++c)
    {
      putchar(255 * (highlight ? red[c] : MU));
    }
    
  } //  for (int i = 0; i < W; ++i)
  
  
  
  return 0;
}
