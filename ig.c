#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenM

/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/blob/HEAD:/book/

generating an improved grid ( in the exterior) 



HSV 
* Hue is circular, from red at 0, with increasing values running through the rainbow yellow
green blue violet, and then back through purple and pink to red again at 1, where the cycle repeats.
* Saturation blends colour intensity from greyscale at 0 to full blast at 1, 
* value controls brightness: black at 0 and full strength at 1


We use single precision float for the calculations: 24bits is more than enough, given each colour channel has only 8bits. 
We need to scale the output values from hsv2rgb() from [0..1] to [0..255], and we clamp them to the output
range to make sure no overflow glitches occur.





gcc ig.c -lm -Wall -fopenmp

./a.out > ig.pgm

*/



static double TwoPi=2.0*M_PI;

double c_arg(complex double z)
{
 double arg;
 arg = carg(z);
 if (arg<0.0) arg+= TwoPi ; 
 return arg; 
}

double c_turn(complex double z)
{
 double arg;
 arg = c_arg(z);
 return arg/TwoPi; 
}


const double pi = 3.141592653589793;

// 
double cnorm(double _Complex z) // https://stackoverflow.com/questions/6363247/what-is-a-complex-data-type-and-an-imaginary-data-type-in-c
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
  int aa = 4; // 
  int w = 800 * aa;
  int h = 800 * aa;
  int nMax = 1024;
  double r = 2;// radius of the plane 
  //double px = r / (h/2); // pixel size 
  
  double er = 512;
  double er2 = er * er; // escape_radius *escape_radius
  
  
  unsigned char *img = malloc(3 * w * h);
 
 
  #pragma omp parallel for schedule(static, 1)
  for (int j = 0; j < h; ++j)
  {
    double y = (h/2 - (j + 0.5)) / (h/2) * r;
    
    for (int i = 0; i < w; ++i)
    {
      double x = (i + 0.5 - w/2) / (h/2) * r;
      double _Complex c = x + I * y;
      double _Complex z = 0;
    
      // iteration    
      int n;
      for (n = 0; n < nMax; ++n)
      {
        z = z * z + c;
        if (cnorm(z) > er2) // abs(z)>er
          break;
      }
      
      // colour
      double hue = 0, sat = 0, val = 0; // interior
      
      
      if (n < nMax) // exterior 
      {
       
        double final_z_abs = log(cabs(z)) / log(er) - 1.0;
        int final_n = n;
        double continuous_escape_time = final_n - log2(final_z_abs + 1.0);
        double  final_z_arg = fmod( carg(z)/TwoPi + 1.0, 1.0); 
        
        // improved grid
        double k = pow ( 0.5 , 0.5 - final_z_abs ) ;
	double grid_weight = 0.05 ;
        int grid = 
        	grid_weight < final_z_abs &&
        	final_z_abs < 1.0 - grid_weight &&
        	grid_weight * k < final_z_arg && 
        	final_z_arg < 1.0 - grid_weight * k;
        
        
        //        
        hue = continuous_escape_time/64.0;
        sat = grid * 0.7; 
        val = 1.0;
      }
      
      // hsv to rgb conversion
      int red, grn, blu;
      hsv2rgb(hue, sat, val, &red, &grn, &blu);
      // 
      int k = 3*(j * w + i);
      img[k+0] = red;
      img[k+1] = grn;
      img[k+2] = blu;
    }
  }
 
 
  printf("P6\n%d %d\n255\n", w, h);
  fwrite(img, 3 * w * h, 1, stdout);
  free(img);
 
 
  return 0;
}
