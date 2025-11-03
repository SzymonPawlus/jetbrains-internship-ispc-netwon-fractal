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
#include <vector>
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

struct HSV {
  float h;
  float s;
  float v;
};

struct RGB {
  float r;
  float g;
  float b;

  RGB(float r, float g, float b) : r(r), g(g), b(b) {}

  RGB(HSV &hsv) {
    double c = hsv.v * hsv.s;
    double x = c * (1 - std::fabs(std::fmod(hsv.h * 6, 2)));
    double m = hsv.v - c;
    if (hsv.h < 1.0 / 6) {
      r = c;
      g = x;
      b = 0;
    } else if (hsv.h < 2.0 / 6) {
      r = x;
      g = c;
      b = 0;
    } else if (hsv.h < 3.0 / 6) {
      r = 0;
      g = c;
      b = x;
    } else if (hsv.h < 4.0 / 6) {
      r = 0;
      g = x;
      b = c;
    } else if (hsv.h < 5.0 / 6) {
      r = x;
      g = 0;
      b = c;
    } else {
      r = c;
      g = 0;
      b = x;
    }
    r += m;
    g += m;
    b += m;
  }

  RGB operator*(float s) { return {r * s, g * s, b * s}; }

  void print(FILE *fp) {
    unsigned char r_c = static_cast<unsigned char>(r * 255);
    unsigned char g_c = static_cast<unsigned char>(g * 255);
    unsigned char b_c = static_cast<unsigned char>(b * 255);

    fputc(r_c, fp);
    fputc(g_c, fp);
    fputc(b_c, fp);
  }
};

extern void newton_serial(int width, int height, float x0, float y0, float x1,
                          float y1, int maxIterations, float n, float root_re[],
                          float root_im[], int colour[]);

/* Write a PPM image file with the image of the Mandelbrot set */
static void writePPM(int *buf, int width, int height, const char *fn,
                     std::vector<RGB> &root_colours) {
  int maxIterations = *std::max_element(buf, buf + width * height);
  FILE *fp = fopen(fn, "wb");
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", width, height);
  fprintf(fp, "255\n");
  for (int i = 0; i < width * height; ++i) {
    root_colours[buf[i]].print(fp);
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

static std::vector<RGB> generate_root_colours(int n) {
  std::vector<RGB> colours;
  colours.emplace_back(0, 0, 0);
  int step = 1;

  for (int i = 0; i < n; ++i) {
    float hue = std::fmod(i * step, n) / float(n);
    HSV hsv = {hue, 1.0, 1.0};
    colours.emplace_back(hsv);
  }
  return colours;
}

int main(int argc, char *argv[]) {
  static unsigned int test_iterations[] = {3, 3};
  unsigned int width = 2048;
  unsigned int height = 2048;
  float x0 = -2;
  float x1 = 2;
  float y0 = -2;
  float y1 = 2;
  int n = 3;

  if (argc > 1) {
    n = atoi(argv[1]);
  }

  int maxIterations = 256;
  int *buf = new int[width * height];
  float *roots_re = new float[n];
  float *roots_im = new float[n];
  calculate_roots(roots_re, roots_im, n);
  std::vector<RGB> colours = generate_root_colours(n);

  //
  // Compute the image using the ispc implementation; report the minimum
  // time of three runs.
  //
  double minISPC = 1e30;
  for (int i = 0; i < test_iterations[0]; ++i) {
    reset_and_start_timer();
    newton_ispc_standard(width, height, x0, y0, x1, y1, maxIterations, n,
                         roots_re, roots_im, buf);
    double dt = get_elapsed_mcycles();
    printf("@time of ISPC run:\t\t\t[%.3f] million cycles\n", dt);
    minISPC = std::min(minISPC, dt);
  }

  printf("[newon standard ispc]:\t\t[%.3f] million cycles\n", minISPC);
  writePPM(buf, width, height, "newton-ispc.ppm", colours);

  // Clear out the buffer
  for (unsigned int i = 0; i < width * height; ++i)
    buf[i] = 0;

  //
  // Compute the image using the radial ispc implementation
  //
  double minISPC_radial = 1e30;
  for (int i = 0; i < test_iterations[0]; ++i) {
    reset_and_start_timer();
    newton_ispc_radial(width, height, x0, y0, x1, y1, maxIterations, n,
                       roots_re, roots_im, buf);
    double dt = get_elapsed_mcycles();
    printf("@time of ISPC radial run:\t[%.3f] million cycles\n", dt);
    minISPC_radial = std::min(minISPC_radial, dt);
  }
  printf("[newton radial ispc]:\t\t[%.3f] million cycles\n", minISPC_radial);

  //
  // Compute the image using the serial implementation; report the minimum
  // time of three runs.
  //
  double minSerial = 1e30;
  for (int i = 0; i < test_iterations[1]; ++i) {
    reset_and_start_timer();
    newton_serial(width, height, x0, y0, x1, y1, maxIterations, n, roots_re,
                  roots_im, buf);
    double dt = get_elapsed_mcycles();
    printf("@time of serial run:\t\t[%.3f] million cycles\n", dt);
    minSerial = std::min(minSerial, dt);
  }
  printf("[newton serial]:\t\t[%.3f] million cycles\n", minSerial);

  return 0;
}
