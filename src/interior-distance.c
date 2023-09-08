/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/blob/HEAD:/book/


gcc i.c -lm -Wall -Wextra -fopenmp

./a.out > id.pgm

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h> //OpenM


const double pi = 3.141592653589793;




// fast abs(complex number)
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


// find periodic point w = z using Newtom method
// n = number of Newton method steps
// p = period
// c = point of parameter  plane, and also parameter of fc = 2^2 + c function
// w0 reasonable starting guess for Newton's methods
// output : w = attractor = periodic point 
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

// compute distance estimation 
// R = escape radius = bailout valued
// N = iMax = maximal number of iteration
// c = point of parameter  plane, and also parameter of fc = 2^2 + c function
double m_distance(int N, double R, double _Complex c)
{
  double _Complex dc = 0;
  double _Complex z = 0;
  double m = 1.0 / 0.0;
  int p = 0;
  for (int n = 1; n <= N; ++n)
  {
    dc = 2 * z * dc + 1; // derivative 
     z = z * z + c; // fc
    if (cabs(z) > R)
    	// point is exterior = escapes
      // compute exterior distance estimate from z and dc
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
      	// point is interior with period p and known z0
        // compute interior distance estimate
        return m_interior_distance(z0, c, p);
    }
  }
  return 0;
}


int main()
{	
	double aspect_ratio = 1.0; // https://en.wikipedia.org/wiki/Aspect_ratio_(image)
	
	// integer coordinate ( pixel or screen or image coordinate)
	int width =  1200;
	int height = 1200;
	
	// double coordinate ( real world)
	double plane_radius = 2;
	double plane_center = -0.5;
	
	double cxMin = creal(plane_center) - plane_radius*aspect_ratio;	
	double cxMax = creal(plane_center) + plane_radius*aspect_ratio;	//
	double cyMin = cimag(plane_center) - plane_radius;	// inv
	double cyMax = cimag(plane_center) + plane_radius;	//
	double PixelWidth = (cxMax - cxMin)/width;
	double PixelHeight = (cyMax - cyMin)/height; // pixel_size = PixelWidth = PixelHeight
	
		
	int kMax = 1024; // maximal number of iterations
	double escape_radius = 2;
	//double escape_radius_2 = escape_radius * escape_radius;
  
	unsigned char *img = malloc(width * height); // image = dynamic array of colors
	
	
	#pragma omp parallel for
	for (int j = 0; j < height; ++j)
	{
		//double y = (height/2 - (j + 0.5)) / (height/2) * plane_radius;
		double y =  cyMin + j*PixelHeight; /* mapping from screen to world; reverse Y  axis */
		for (int i = 0; i < width; ++i)
		{
			// double x = plane_center + (i + 0.5 - width/2) / (height/2) * plane_radius;
			double x = cxMin + i*PixelWidth;
			double _Complex c = x + I * y; // parameter c of  fc(z) = z^2 + c
			
			double de = m_distance(kMax, escape_radius, c);
      			double g = tanh(de / PixelWidth );
      			img[j * width + i] = 255 * g;   // compute and save color to array
		}
	}
	
	
	// create pgm file using command redirection https://en.wikipedia.org/wiki/Redirection_(computing)
	printf("P5\n%d %d\n255\n", width, height); // header of P5 = binary pgm 
	fwrite(img, width * height, 1, stdout);
	free(img);
	return 0;


}
