/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/tree/HEAD:/book


gcc membership.c -lm -Wall -fopenmp

time ./a.out > m.pgm

convert m.pgm -resize 600x600 m.png


real	0m0.414s
user	0m1.045s
sys	0m0.008s

*/

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> //OpenM

/*
fork of 
mandelbrot-book	how to write a book about the Mandelbrot set by Claude Heiland-Alle
https://code.mathr.co.uk/mandelbrot-book/tree/HEAD:/book


gcc membership.c -lm -Wall -fopenmp

time ./a.out > m.pgm


convert m.pgm -resize 600x600 m.png


cd existing_folder
git init
git remote add origin git@gitlab.com:adammajewski/mandelbrot-book_book.git
git add membership.c
git commit -m "Initial commit"
git push -u origin master

*/


// fast abs(complex number)
 double cnorm(double _Complex z)
{
  return creal(z) * creal(z) + cimag(z) * cimag(z);
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
	double escape_radius_2 = escape_radius * escape_radius;
  
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
			double _Complex z = 0; // z = z0 = critical point 
			int k; // number of iterations
			for (k = 0; k < kMax; ++k)
			{
				z = z * z + c; // forward iteration of complex quadratic polynomial
				if (cnorm(z) > escape_radius_2) // if abs(z) > ER
					{break;}
			}
			img[j * width + i] = (k == kMax) ? 0 : 255; // compute and save color to array
		}
	}
	// create pgm file using command redirection https://en.wikipedia.org/wiki/Redirection_(computing)
	printf("P5\n%d %d\n255\n", width, height); // header of P5 = binary pgm 
	fwrite(img, width * height, 1, stdout);
	free(img);
	return 0;


}
