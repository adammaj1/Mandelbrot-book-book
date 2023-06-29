// http://mathr.co.uk/blog/2014-11-22_adaptive_supersampling_using_distance_estimate.html
// gcc -std=c99 -Wall -Wextra -pedantic -O3 -fopenmp -o adaptive-dem adaptive-dem.c -lm
// ./a.out > a.pgm // P5 = binary Portable GrayMap 
//  %s width height creal cimag radius maxiters depth > out.pgm
// 

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int counter;

double cabs2(complex double z) { return creal(z) * creal(z) + cimag(z) * cimag(z); }

double colour(complex double c, double pixel_spacing, int maxiters, int depth) {
  {
    #pragma omp atomic
    counter++;
  }
  double er2 = 512 * 512;
  complex double z = 0;
  complex double dc = 0;
  for (int i = 0; i < maxiters; ++i) {
    dc = 2 * z * dc + 1;
    z = z * z + c;
    if (cabs2(z) > er2) {
      double de = 2 * cabs(z) * log(cabs(z)) / (cabs(dc) * pixel_spacing);
      if (de <= 4 && depth > 0) {
        double p2 = pixel_spacing / 2;
        double p4 = pixel_spacing / 4;
        return 0.25 *
          ( colour(c + p4 + I * p4, p2, maxiters, depth - 1)
          + colour(c - p4 + I * p4, p2, maxiters, depth - 1)
          + colour(c + p4 - I * p4, p2, maxiters, depth - 1)
          + colour(c - p4 - I * p4, p2, maxiters, depth - 1)
          );
      } else {
        return tanh(de);
      }
    }
  }
  return 0;
}

void render(unsigned char *image, int width, int height, complex double c0, double radius, int maxiters, int depth) {
  double pixel_spacing = radius / (height/2.0);
  #pragma omp parallel for schedule(dynamic, 1)
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      complex double c = c0 + pixel_spacing * ((i + 0.5 - width/2) + I * (height/2 - j - 0.5));
      double v = colour(c, pixel_spacing, maxiters, depth);
      int g = fmin(fmax(255 * v, 0), 255);
      image[width * j + i] = g;
    }
  }
}

int main(int argc, char **argv) {
  if (argc != 8) {
    fprintf(stderr, "usage: %s width height creal cimag radius maxiters depth > out.pgm\n", argv[0]);
    return 1;
  }
  int width = atoi(argv[1]);
  int height = atoi(argv[2]);
  complex double c = atof(argv[3]) + I * atof(argv[4]);
  double radius = atof(argv[5]);
  int maxiters = atoi(argv[6]);
  int depth = atoi(argv[7]);
  unsigned char *image = malloc(width * height);
  counter = 0;
  render(image, width, height, c, radius, maxiters, depth);
  fprintf(stderr, "width = %d\t height = %d\t c = %f%+f\t radius = %f\t maxiters = %d\n", width, height, creal(c), cimag(c), radius, maxiters);
  fprintf(stderr, "depth = %d\tsamples per pixel ss = %.2f%%\n", depth, counter * 100.0 / (double)(width * height));
  printf("P5\n%d %d\n255\n", width, height);
  fwrite(image, width * height, 1, stdout);
  
   
  return 0;
}
