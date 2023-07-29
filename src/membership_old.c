#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

static inline double cnorm(double _Complex z)
{
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

int main()
{
  int aa = 1;
  int w = 1200 * aa;
  int h = 1200 * aa;
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
      double x = -0.5+(i + 0.5 - w/2) / (h/2) * r - 0.75;
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
