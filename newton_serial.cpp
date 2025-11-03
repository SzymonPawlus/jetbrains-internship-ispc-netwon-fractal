#include <iostream>
static inline float fast_pow(float base, int exp) {
  // Fast exponentiation by squaring for integer exponents
  float result = 1.0f;
  while (exp > 0) {
    if (exp & 0x01) {
      result *= base;
    }
    base *= base;
    exp >>= 1;
  }
  return result;
}

static inline void complex_pow(float &xr, float &xi, int n) {
  float resr = 1.0f, resi = 0.0f;
  float baser = xr, basei = xi;
  int e = n;

  while (e > 0) {
    if (e & 1) {
      float tr = resr * baser - resi * basei;
      float ti = resr * basei + resi * baser;
      resr = tr;
      resi = ti;
    }
    e >>= 1;
    if (e) {
      float tr = baser * baser - basei * basei;
      float ti = 2.0f * baser * basei;
      baser = tr;
      basei = ti;
    }
  }
  xr = resr;
  xi = resi;
}

static inline float modulus_sqr(float xr, float xi) {
  return xr * xr + xi * xi;
}

static inline void complex_mul(float &ar, float &ai, float br, float bi) {
  float tr = ar * br - ai * bi;
  float ti = ar * bi + ai * br;
  ar = tr;
  ai = ti;
}

static inline void complex_inverse(float &xr, float &xi) {
  float denom = modulus_sqr(xr, xi);
  xr = xr / denom;
  xi = -xi / denom;
}

int newton_one(float z_re, float z_im, int maxIterations, int n, float &root_re,
               float &root_im) {
  // Calculate Newton's method for f(z) = z^n - 1
  bool inverted = false;
  float nf = static_cast<float>(n);

  for (int i = 0; i < maxIterations; ++i) {
    float new_re, new_im;
    float f_prime_re = z_re, f_prime_im = z_im;
    complex_pow(f_prime_re, f_prime_im, n - 1);
    float f_re = z_re, f_im = z_im;
    complex_mul(f_re, f_im, f_prime_re, f_prime_im);
    f_re -= 1.0f; // f(z) = z^n - 1
    complex_mul(f_prime_re, f_prime_im, n, 0);

    // Newton's method: z = z - f(z) / f'(z)
    float denom = modulus_sqr(f_prime_re, f_prime_im);
    if (denom == 0.0f)
      break;

    new_re = z_re - (f_re * f_prime_re + f_im * f_prime_im) / denom;
    new_im = z_im - (f_im * f_prime_re - f_re * f_prime_im) / denom;
    // Check for convergence
    float diff_re = new_re - z_re;
    float diff_im = new_im - z_im;
    if (diff_re * diff_re + diff_im * diff_im < 1e-6f) {
      if (inverted) {

      } else {
        root_re = new_re;
        root_im = new_im;
      }
      return i;
    }

    z_re = new_re;
    z_im = new_im;
    std::cout << "Iteration " << i << ": new_re = " << new_re
              << ", new_im = " << new_im << std::endl;
  }
  std::cout << "Did not converge: new_re = " << z_re << ", new_im = " << z_im
            << std::endl;
  return -1;
}

void newton_serial(int width, int height, float x0, float y0, float x1,
                   float y1, int maxIterations, float n, float root_re[],
                   float root_im[], int colour[]) {
  float dx = (x1 - x0) / width;
  float dy = (y1 - y0) / height;

  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      // Write which root for each pixel here using complex_pow function
      float xr = x0 + i * dx;
      float xi = y0 + j * dy;
      int iteration = 0;
      for (iteration = 0; iteration < maxIterations; iteration++) {
        // Compute f(z) and f'(z)
        float fr = xr;
        float fi = xi;
        complex_pow(fr, fi, n); // f(z) = z^n

        float fpr = n;
        float fpi = 0.0f;
        complex_pow(xr, xi, n - 1); // z^(n-1)
        fpr *= xr;
        fpi *= xi;

        // Newton's method: z = z - f(z)/f'(z)
        float denom = fpr * fpr + fpi * fpi;
        if (denom == 0.0f) {
          break; // Avoid division by zero
        }
        float tr = (fr * fpr + fi * fpi) / denom;
        float ti = (fi * fpr - fr * fpi) / denom;

        xr -= tr;
        xi -= ti;

        // Check for convergence to any root
        bool converged = false;
        for (int k = 0; k < n; k++) {
          float dr = xr - root_re[k];
          float di = xi - root_im[k];
          if (dr * dr + di * di < 1e-6f) {
            converged = true;
            break;
          }
        }
        if (converged) {
          break;
        }
      }
    }
  }
}
