#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenM

/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/blob/HEAD:/book/

generating an improved grid ( in the exterior) and removed noise near boundary 



In the book it is : 1.31 Distance estimation.

Our images look noisy and grainy near the boundary of the Mandelbrot set. The escape time
bands get closer and closer, while the pixel spacing is fixed. The pixel grid samples isolated
points of a mathematically abstract image defined on the continuous plane. The Nyquist-Shannon
sampling theorem shows that sampling isolated points from a continuum is a valid approximation
only so long as the values don’t change too quickly between the points. Aliasing occurs when
the values do change too quickly compared to the sampling rate, with the grainy noisy visual
effects as we have seen. Because the escape time bands increase in number without bound as we
approach the boundary of the Mandelbrot set, no sampling rate can be high enough.











gcc n.c -lm -Wall -fopenmp

./a.out > n.ppm

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
  double pixel_spacing = r / ( h/2.0); // = radius / ( height / 2.0 ) ;
  
  double er = 512; // R
  double er2 = er * er; // escape_radius *escape_radius = R*R = R2
  
  
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
      double _Complex dc = 0; // derivative with the respect to c 
    
      // iteration    
      int n;
      for (n = 0; n < nMax; ++n)
      {
      
        dc = 2 * z * dc + 1;
        z = z * z + c;
        
        if (cnorm(z) > er2) // abs(z)>er
          break;
      }
      
      // colour
      double hue = 0, sat = 0, val = 1.0; // interior
      
      
      if (n < nMax) // exterior 
      {
       
       
       double de = 2.0 * cabs(z) * log( cabs(z) ) / ( cabs( dc ) * pixel_spacing ) ; // distance estimation
       double final_z_abs = log(cabs(z)) / log(er) - 1.0; // not only cabs(z) 
       int final_n =  n; //
       double  final_z_arg = fmod( carg(z)/TwoPi + 1.0, 1.0); 
        
       double continuous_escape_time = final_n - log2(final_z_abs + 1.0);
        
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
        /*  The de here in your code is scaled by pixel spacing, so should be less than around 1 for pixels near the boundary. */
        val = tanh ( fmin ( fmax (  de , 0.0 ) , 4.0 ) ) ;
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
 
 
  printf("P6\n%d %d\n255\n", w, h); // ppm
  fwrite(img, 3 * w * h, 1, stdout);
  free(img);
 
 
  return 0;
}
