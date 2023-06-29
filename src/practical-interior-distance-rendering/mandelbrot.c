// http://mathr.co.uk/blog/2014-11-02_practical_interior_distance_rendering.html

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double pi = 3.141592653589793;
const double infinity = 1.0 / 0.0;
const double phi = 1.618033988749895; // (sqrt(5.0) + 1.0) / 2.0;
const double colour_modulus = 5.7581917135421046e-2; // (1.0 + 1.0 / (phi * phi)) / 24.0;
const double escape_radius_2 = 512.0 * 512.0;

const int BIAS_UNKNOWN = 0;
const int BIAS_INTERIOR = 1;
const int BIAS_EXTERIOR = 2;

const int ALGORITHM_PLAIN = 0;
const int ALGORITHM_UNBIASED = 1;
const int ALGORITHM_BIASED = 2;
const int ALGORITHM_ANALYSE = 4;

static inline double cabs2(complex double z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

static inline unsigned char *image_new(int width, int height) {
  return malloc(width * height * 3);
}

static inline void image_delete(unsigned char *image) {
  free(image);
}

static inline void image_save_ppm(unsigned char *image, int width, int height, const char *filename) {
  FILE *f = fopen(filename, "wb");
  if (f) {
    fprintf(f, "P6\n%d %d\n255\n", width, height);
    fwrite(image, width * height * 3, 1, f);
    fclose(f);
  } else {
    fprintf(stderr, "ERROR saving `%s'\n", filename);
  }
}

static inline void image_poke(unsigned char *image, int width, int i, int j, int r, int g, int b) {
  int k = (width * j + i) * 3;
  image[k++] = r;
  image[k++] = g;
  image[k  ] = b;
}

static inline void colour_hsv_to_rgb(double h, double s, double v, double *r, double *g, double *b) {
  double i, f, p, q, t;
  if (s == 0) { *r = *g = *b = v; } else {
    h = 6 * (h - floor(h));
    int ii = i = floor(h);
    f = h - i;
    p = v * (1 - s);
    q = v * (1 - (s * f));
    t = v * (1 - (s * (1 - f)));
    switch(ii) {
      case 0: *r = v; *g = t; *b = p; break;
      case 1: *r = q; *g = v; *b = p; break;
      case 2: *r = p; *g = v; *b = t; break;
      case 3: *r = p; *g = q; *b = v; break;
      case 4: *r = t; *g = p; *b = v; break;
      default:*r = v; *g = p; *b = q; break;
    }
  }
}

static inline void colour_to_bytes(double r, double g, double b, int *r_out, int *g_out, int *b_out) {
  *r_out = fmin(fmax(255 * r, 0), 255);
  *g_out = fmin(fmax(255 * g, 0), 255);
  *b_out = fmin(fmax(255 * b, 0), 255);
}

static inline void colour_mandelbrot(unsigned char *image, int width, int i, int j, int period, double distance) {
  double r, g, b;
  colour_hsv_to_rgb(period * colour_modulus, 0.5, tanh(distance), &r, &g, &b);
  int ir, ig, ib;
  colour_to_bytes(r, g, b, &ir, &ig, &ib);
  image_poke(image, width, i, j, ir, ig, ib);
}

static inline void colour_analysis(unsigned char *image, int width, int i, int j, int bias, int outcome) {
  int ir = 0, ig = 0, ib = 0;
  if (bias == outcome) {
    ig = 255;
  } else if (bias == BIAS_INTERIOR && outcome == BIAS_EXTERIOR) {
    ib = 255;
  } else if (bias == BIAS_EXTERIOR && outcome == BIAS_INTERIOR) {
    ir = 255;
  }
  image_poke(image, width, i, j, ir, ig, ib);
}

static inline int attractor(complex double *z_out, complex double *dz_out, complex double z0, complex double c, int period) {
  double epsilon_2 = 1e-20;
  complex double zz = z0;
  for (int j = 0; j < 64; ++j) {
    complex double z = zz;
    complex double dz = 1;
    for (int i = 0; i < period; ++i) {
      dz = 2.0 * z * dz;
      z = z * z + c;
    }
    complex double zz1 = zz - (z  - zz) / (dz - 1.0);
    if (cabs2(zz1 - zz) < epsilon_2) {
      *z_out = z;
      *dz_out = dz;
      return 1;
    }
    zz = zz1;
  }
  return 0;
}

static inline double interior_distance(complex double z0, complex double c, int period) {
  complex double z = z0;
  complex double dz = 1;
  complex double dzdz = 0;
  complex double dc = 0;
  complex double dcdz = 0;
  for (int p = 0; p < period; ++p) {
    dcdz = 2 * (z * dcdz + dz * dc);
    dc = 2 * z * dc + 1;
    dzdz = 2 * (dz * dz + z * dzdz);
    dz = 2 * z * dz;
    z = z * z + c;
  }
  return (1 - cabs2(dz)) / cabs(dcdz + dzdz * dc / (1 - dz));
}

struct partial {
  complex double z;
  int p;
};

static inline void render(unsigned char *image, int algorithm, int maxiters, int width, int height, complex double center, double radius) {
  double pixel_spacing = radius / (height / 2.0);
  #pragma omp parallel for schedule(dynamic, 1)
  for (int j = 0; j < height; ++j) {
    struct partial *partials = 0;
    int bias = BIAS_EXTERIOR, new_bias, npartials;
    if (algorithm & ALGORITHM_BIASED) {
      partials = malloc(maxiters * sizeof(struct partial));
    }
    for (int i = 0; i < width; ++i) {
      new_bias = BIAS_UNKNOWN;
      npartials = 0;
      double x = i + 0.5 - width / 2.0;
      double y = height / 2.0 - j - 0.5;
      complex double c = center + pixel_spacing * (x + I * y);
      complex double z = 0;
      complex double dc = 0;
      double minimum_z2 = infinity;
      int period = 0;
      for (int n = 1; n <= maxiters; ++n) {
        dc = 2 * z * dc + 1;
        z = z * z + c;
        double z2 = cabs2(z);
        if (z2 < minimum_z2) {
          minimum_z2 = z2;
          period = n;
          if (algorithm & (ALGORITHM_UNBIASED | ALGORITHM_BIASED)) {
            if (algorithm & ALGORITHM_UNBIASED || bias == BIAS_INTERIOR) {
              complex double z0 = 0, dz0 = 0;
              if (attractor(&z0, &dz0, z, c, period)) {
                if (cabs2(dz0) <= 1.0) {
                  if (algorithm & ALGORITHM_ANALYSE) {
                    colour_analysis(image, width, i, j, bias, BIAS_INTERIOR);
                  } else {
                    double distance = interior_distance(z0, c, period) / pixel_spacing;
                    colour_mandelbrot(image, width, i, j, period, distance);
                  }
                  new_bias = BIAS_INTERIOR;
                  break;
                }
              }
            } else if (algorithm & ALGORITHM_BIASED) {
              partials[npartials].z = z;
              partials[npartials].p = period;
              npartials++;
            }
          }
        }
        if (z2 >= escape_radius_2) {
          if (algorithm & ALGORITHM_ANALYSE) {
            colour_analysis(image, width, i, j, bias, BIAS_EXTERIOR);
          } else {
            double distance = sqrt(z2) * log(z2) / (cabs(dc) * pixel_spacing);
            colour_mandelbrot(image, width, i, j, period, distance);
          }
          new_bias = BIAS_EXTERIOR;
          break;
        }
      }
      if (algorithm & ALGORITHM_BIASED) {
        if (bias == BIAS_EXTERIOR && new_bias == BIAS_UNKNOWN) {
          for (int n = 0; n < npartials; ++n) {
            complex double z = partials[n].z;
            int period = partials[n].p;
            complex double z0 = 0, dz0 = 0;
            if (attractor(&z0, &dz0, z, c, period)) {
              if (cabs2(dz0) <= 1.0) {
                if (algorithm & ALGORITHM_ANALYSE) {
                  colour_analysis(image, width, i, j, bias, BIAS_INTERIOR);
                } else {
                  double distance = interior_distance(z0, c, period) / pixel_spacing;
                  colour_mandelbrot(image, width, i, j, period, distance);
                }
                new_bias = BIAS_INTERIOR;
                break;
              }
            }
          }
        }
      }
      if (new_bias == BIAS_UNKNOWN) {
        if (algorithm & ALGORITHM_ANALYSE) {
          colour_analysis(image, width, i, j, bias, BIAS_UNKNOWN);
        } else {
          if (algorithm & (ALGORITHM_UNBIASED | ALGORITHM_BIASED)) {
            colour_mandelbrot(image, width, i, j, period, 0.0);
          } else {
            colour_mandelbrot(image, width, i, j, period, 10.0);
          }
        }
        new_bias = BIAS_EXTERIOR;
      }
      bias = new_bias;
    }
    if (algorithm & ALGORITHM_BIASED) {
      free(partials);
    }
  }
}

int main(int argc, char **argv) {
  if (argc != 9) {
    fprintf(stderr,
      "usage: %s algorithm maxiters width height creal cimag radius filename\n"
      "values for algorithm:\n"
      "    0  plain (exterior only)\n"
      "    1  unbiased (exterior and interior, unoptimized)\n"
      "    2  biased (exterior and interior optimized by local connectedness)\n"
      "    6  analysis map of biased rendering\n"
      , argv[0]);
    return 1;
  }
  int algorithm = atoi(argv[1]);
  int maxiters = atoi(argv[2]);
  int width = atoi(argv[3]);
  int height = atoi(argv[4]);
  complex double center = atof(argv[5]) + I * atof(argv[6]);
  double radius = atof(argv[7]);
  const char *filename = argv[8];
  unsigned char *image = image_new(width, height);
  render(image, algorithm, maxiters, width, height, center, radius);
  image_save_ppm(image, width, height, filename);
  image_delete(image);
  return 0;
}
