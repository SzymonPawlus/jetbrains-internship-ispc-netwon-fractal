/*
  Copyright (c) 2010-2023, Intel Corporation

  SPDX-License-Identifier: BSD-3-Clause
*/

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
