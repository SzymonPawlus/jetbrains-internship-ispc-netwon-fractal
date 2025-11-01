/*
  Copyright (c) 2010-2023, Intel Corporation

  SPDX-License-Identifier: BSD-3-Clause
*/

static inline int mandel(float c_re, float c_im, int count) {
  float z_re = c_re, z_im = c_im;
  int i;
  for (i = 0; i < count; ++i) {
    if (z_re * z_re + z_im * z_im > 4.)
      break;

    float new_re = z_re * z_re - z_im * z_im;
    float new_im = 2.f * z_re * z_im;
    unmasked {
      z_re = c_re + new_re;
      z_im = c_im + new_im;
    }
  }

  return i;
}

static inline int newton_radial(float z_re, float z_im, int maxIterations,
                                float n, float &root_re, float &root_im) {
  // Calculate Newton's method for f(z) = z^n - 1
  for (int i = 0; i < maxIterations; ++i) {
    // Calculate z^(n-1)
    float r = sqrt(z_re * z_re + z_im * z_im);
    float theta = atan2(z_im, z_re);
    float r_n_minus_1 = pow(r, n - 1);
    float theta_n_minus_1 = (n - 1) * theta;
    float r_n = pow(r, n);
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

export void mandelbrot_ispc(uniform float x0, uniform float y0,
                            uniform float x1, uniform float y1,
                            uniform int width, uniform int height,
                            uniform int maxIterations, uniform int output[]) {
  float dx = (x1 - x0) / width;
  float dy = (y1 - y0) / height;

  for (uniform int j = 0; j < height; j++) {
    // Note that we'll be doing programCount computations in parallel,
    // so increment i by that much.  This assumes that width evenly
    // divides programCount.
    foreach (i = 0 ... width) {
      // Figure out the position on the complex plane to compute the
      // number of iterations at.  Note that the x values are
      // different across different program instances, since its
      // initializer incorporates the value of the programIndex
      // variable.
      float x = x0 + i * dx;
      float y = y0 + j * dy;

      int index = j * width + i;
      output[index] = mandel(x, y, maxIterations);
    }
  }
}

typedef struct {
  uniform int width;
  uniform int height;
  uniform float x0;
  uniform float y0;
  uniform float x1;
  uniform float y1;
} size_params_t;

export void newton_ispc(uniform int width, uniform int height, uniform float x0,
                        uniform float y0, uniform float x1, uniform float y1,
                        uniform int maxIterations, uniform float n,
                        uniform float root_re[], uniform float root_im[],
                        uniform int colour[]) {
  float dx = (x1 - x0) / width;
  float dy = (y1 - y0) / height;

  for (uniform int j = 0; j < height; j++) {
    foreach (i = 0 ... width) {
      float z_re = x0 + i * dx;
      float z_im = y0 + j * dy;
      float r_re, r_im;
      int iters = newton_radial(z_re, z_im, maxIterations, n, r_re, r_im);
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
