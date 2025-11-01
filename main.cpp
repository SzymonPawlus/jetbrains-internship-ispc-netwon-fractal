/*
  Copyright (c) 2010-2011, Intel Corporation
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.


   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cmath>
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define NOMINMAX
#pragma warning(disable : 4244)
#pragma warning(disable : 4305)
#endif

#include "newton.h"
#include "timing.h"
#include <algorithm>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
using namespace ispc;

extern void mandelbrot_serial(float x0, float y0, float x1, float y1, int width,
                              int height, int maxIterations, int output[]);

const int default_colours[16] = {0x000000, 0x0000FF, 0x00FF00, 0x00FFFF,
                                 0xFF0000, 0xFF00FF, 0xFFFF00, 0xFFFFFF,
                                 0xC6C6C6, 0x840000, 0x008400, 0x848400,
                                 0x000084, 0x840084, 0x008484, 0x848484};

/* Write a PPM image file with the image of the Mandelbrot set */
static void writePPM(int *buf, int width, int height, const char *fn) {
  FILE *fp = fopen(fn, "wb");
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", width, height);
  fprintf(fp, "255\n");
  for (int i = 0; i < width * height; ++i) {
    // Map the iteration count to colors by just alternating between
    // two greys.
    int colour = default_colours[buf[i] % 16];

    for (int j = 0; j < 3; ++j) {
      unsigned char c = (unsigned char)(colour & 0xFF);
      colour >>= 8;
      fputc(c, fp);
    }
  }
  fclose(fp);
  printf("Wrote image file %s\n", fn);
}

static void calculate_roots(float *roots_re, float *roots_im, int n) {
  for (int k = 0; k < n; ++k) {
    float angle = 2.0f * 3.14159265358979323846f * float(k) / float(n);
    roots_re[k] = cosf(angle);
    roots_im[k] = sinf(angle);
  }
}

int main(int argc, char *argv[]) {
  static unsigned int test_iterations[] = {3, 3};
  unsigned int width = 4096;
  unsigned int height = 4096;
  float x0 = -1;
  float x1 = 1;
  float y0 = -1;
  float y1 = 1;
  int n = 7;

  if (argc > 1) {
    if (strncmp(argv[1], "--scale=", 8) == 0) {
      float scale = atof(argv[1] + 8);
      width *= scale;
      height *= scale;
    }
  }
  if ((argc == 3) || (argc == 4)) {
    for (int i = 0; i < 2; i++) {
      test_iterations[i] = atoi(argv[argc - 2 + i]);
    }
  }

  int maxIterations = 256;
  int *buf = new int[width * height];
  float *roots_re = new float[n];
  float *roots_im = new float[n];
  calculate_roots(roots_re, roots_im, n);

  //
  // Compute the image using the ispc implementation; report the minimum
  // time of three runs.
  //
  double minISPC = 1e30;
  for (int i = 0; i < test_iterations[0]; ++i) {
    reset_and_start_timer();
    newton_ispc(width, height, x0, y0, x1, y1, maxIterations, n, roots_re,
                roots_im, buf);
    double dt = get_elapsed_mcycles();
    printf("@time of ISPC run:\t\t\t[%.3f] million cycles\n", dt);
    minISPC = std::min(minISPC, dt);
  }

  printf("[mandelbrot ispc]:\t\t[%.3f] million cycles\n", minISPC);
  writePPM(buf, width, height, "mandelbrot-ispc.ppm");

  // Clear out the buffer
  for (unsigned int i = 0; i < width * height; ++i)
    buf[i] = 0;

  return 0;
}
