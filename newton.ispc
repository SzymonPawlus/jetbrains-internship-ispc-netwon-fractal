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

static inline int newton_standard(float z_re, float z_im, int maxIterations,
                                  int n, float &root_re, float &root_im) {
  // Calculate Newton's method for f(z) = z^n - 1
  float nf = (float)n;

  for (int i = 0; i < maxIterations; ++i) {
    float new_re, new_im;
    // Heuristic to avoid NaNs for bigger |z|
    if (modulus_sqr(z_re, z_im) > pow(2.0f, 76.f / nf)) {
      new_re = z_re * (nf - 1.f) / nf;
      new_im = z_im * (nf - 1.f) / nf;
    } else {
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
    }
    // Check for convergence
    float diff_re = new_re - z_re;
    float diff_im = new_im - z_im;
    if (diff_re * diff_re + diff_im * diff_im < 1e-10f) {
      root_re = new_re;
      root_im = new_im;
      return i;
    }

    unmasked {
      z_re = new_re;
      z_im = new_im;
    }
  }
  return -1;
}

static inline int newton_radial(float z_re, float z_im, int maxIterations,
                                float n, float &root_re, float &root_im) {
  // Calculate Newton's method for f(z) = z^n - 1
  for (int i = 0; i < maxIterations; ++i) {
    // Calculate z^(n-1)
    float r = sqrt(z_re * z_re + z_im * z_im);
    float theta = atan2(z_im, z_re);
    float r_n_minus_1 = fast_pow(r, n - 1);
    float theta_n_minus_1 = (n - 1) * theta;
    float r_n = fast_pow(r, n);
    float theta_n = n * theta;
    float f_re = r_n_minus_1 * cos(theta_n_minus_1);
    float f_im = r_n_minus_1 * sin(theta_n_minus_1);

    // f'(z) = n * z^(n-1)
    float f_prime_re = n * f_re;
    float f_prime_im = n * f_im;
    // f(z) = z^n - 1
    f_re = r_n * cos(theta_n) - 1.0f;
    f_im = r_n * sin(theta_n);

    // Newton's method: z = z - f(z) / f'(z)
    float denom = f_prime_re * f_prime_re + f_prime_im * f_prime_im;
    if (denom == 0.0f)
      break;

    float new_re = z_re - (f_re * f_prime_re + f_im * f_prime_im) / denom;
    float new_im = z_im - (f_im * f_prime_re - f_re * f_prime_im) / denom;

    // Check for convergence
    float diff_re = new_re - z_re;
    float diff_im = new_im - z_im;
    // Print old and new values for debugging
    if (diff_re * diff_re + diff_im * diff_im < 1e-6f) {
      root_re = new_re;
      root_im = new_im;
      return i;
    }

    unmasked {
      z_re = new_re;
      z_im = new_im;
    }
  }
  return -1;
}

export void newton_ispc_standard(uniform int width, uniform int height,
                                 uniform float x0, uniform float y0,
                                 uniform float x1, uniform float y1,
                                 uniform int maxIterations, uniform int n,
                                 uniform float root_re[],
                                 uniform float root_im[], uniform int colour[],
                                 uniform int iterations[]) {
  float dx = (x1 - x0) / width;
  float dy = (y1 - y0) / height;

  for (uniform int j = 0; j < height; j++) {
    foreach (i = 0 ... width) {
      float z_re = x0 + i * dx;
      float z_im = y0 + j * dy;
      float r_re, r_im;
      int iters = newton_standard(z_re, z_im, maxIterations, n, r_re, r_im);
      iterations[j * width + i] = iters;
      if (iters == -1) {
        colour[j * width + i] = 0;
      } else {
        for (int k = 0; k < n; k++) {
          float diff_re = root_re[k] - r_re;
          float diff_im = root_im[k] - r_im;
          if (diff_re * diff_re + diff_im * diff_im < 1e-3f) {
            colour[j * width + i] = k + 1;
            break;
          }
        }
      }
    }
  }
}

export void newton_ispc_radial(uniform int width, uniform int height,
                               uniform float x0, uniform float y0,
                               uniform float x1, uniform float y1,
                               uniform int maxIterations, uniform int n,
                               uniform float root_re[], uniform float root_im[],
                               uniform int colour[], uniform int iterations[]) {
  float dx = (x1 - x0) / width;
  float dy = (y1 - y0) / height;

  for (uniform int j = 0; j < height; j++) {
    foreach (i = 0 ... width) {
      float z_re = x0 + i * dx;
      float z_im = y0 + j * dy;
      float r_re, r_im;
      int iters = newton_radial(z_re, z_im, maxIterations, n, r_re, r_im);
      iterations[j * width + i] = iters;
      if (iters == -1) {
        colour[j * width + i] = 0;
      } else {
        for (int k = 0; k < n; k++) {
          float diff_re = root_re[k] - r_re;
          float diff_im = root_im[k] - r_im;
          if (diff_re * diff_re + diff_im * diff_im < 1e-3f) {
            colour[j * width + i] = k + 1;
            break;
          }
        }
      }
    }
  }
}
