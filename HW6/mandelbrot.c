%%cu
/*******************************************************************************
To compile: gcc -O3 -o mandelbrot mandelbrot.c -lm
To create an image with 4096 x 4096 pixels: ./mandelbrot 4096 4096
*******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cuda.h"

int writeMandelbrot(const char *fileName, int width, int height, float *img, int minI, int maxI);

#define MXITER 1000

/*******************************************************************************/
// Define a complex number
typedef struct {
  double x;
  double y;
}complex_t;

__device__ int testpoint(complex_t c) {
    int iter;
    complex_t z = c;

    for (iter = 0; iter < MXITER; iter++) {
        // real part of z^2 + c
        double tmp = (z.x * z.x) - (z.y * z.y) + c.x;
        // update with imaginary part of z^2 + c
        z.y = z.x * z.y * 2.0 + c.y;
        // update real part
        z.x = tmp;
        // check bound
        if ((z.x * z.x + z.y * z.y) > 4.0) {
            return iter;
        }
    }
    return iter;
}

__global__ void mandelbrotKernel(int Nre, int Nim, complex_t cmin, complex_t dc, float *count){

  int m = threadIdx.x + blockIdx.x * blockDim.x;
  int n = threadIdx.y + blockIdx.y * blockDim.y;

    if (m < Nre && n < Nim) {
        complex_t c;
        c.x = cmin.x + dc.x * m;
        c.y = cmin.y + dc.y * n;
        count[m + n * Nre] = (float)testpoint(c);
    }
}

/*******************************************************************************/
int main(int argc, char **argv){

  // to create a 4096x4096 pixel image
  // usage: ./mandelbrot 4096 4096

  int Nre = (argc==3) ? atoi(argv[1]): 4096;
  int Nim = (argc==3) ? atoi(argv[2]): 4096;

  // storage for the iteration counts
  float *count;
  count = (float*) malloc(Nre*Nim*sizeof(float));

  // Parameters for a bounding box for "c" that generates an interesting image
  // const float centRe = -.759856, centIm= .125547;
  // const float diam  = 0.151579;
  const float centRe = -0.5, centIm= 0;
  const float diam  = 3.0;

  complex_t cmin;
  complex_t cmax;
  complex_t dc;

  cmin.x = centRe - 0.5*diam;
  cmax.x = centRe + 0.5*diam;
  cmin.y = centIm - 0.5*diam;
  cmax.y = centIm + 0.5*diam;

  //set step sizes
  dc.x = (cmax.x-cmin.x)/(Nre-1);
  dc.y = (cmax.y-cmin.y)/(Nim-1);

  //clock_t start = clock(); //start time in CPU cycles

  // 1. allocate DEVICE array
  float *d_count;

  // 2. allocate memory on device
  cudaMalloc(&d_count, Nre * Nim * sizeof(float));

  // 3. create events
  cudaEvent_t start, end;
  cudaEventCreate(&start);
  cudaEventCreate(&end);

  // 4. calculate number of thread-blocks and threads per thread-block to use
  int T = 256;
  dim3 B(sqrt(T),sqrt(T));
  dim3 G((Nre + B.x - 1) / B.x, (Nim + B.y - 1) / B.y);

  // 5. record start event
  cudaEventRecord(start);

  //6. Launch the Kernel
  mandelbrotKernel<<< G,B >>>(Nre, Nim, cmin, dc, d_count);

  // 7. insert end record event in stream
  cudaEventRecord(end);

  // 8. copy from the GPU back to the host here
  cudaMemcpy(count, d_count, Nre * Nim * sizeof(float), cudaMemcpyDeviceToHost);

  // 9. print out elapsed time
  float elapsed;
  cudaEventSynchronize(end);
  cudaEventElapsedTime(&elapsed, start, end);
  elapsed /= 1000.; // convert to seconds

  printf("elapsed time: %g\n", elapsed);

  // 10. free arrays
  cudaFree(d_count);

  //clock_t end = clock(); //start time in CPU cycles

  // print elapsed time
  //printf("elapsed = %f\n", ((double)(end-start))/CLOCKS_PER_SEC);

  // output mandelbrot to ppm format image
  printf("Printing mandelbrot.ppm...");
  writeMandelbrot("mandelbrot.ppm", Nre, Nim, count, 0, 80);
  printf("done.\n");

  free(count);

  exit(0);
  return 0;
}


/* Output data as PPM file */
void saveppm(const char *filename, unsigned char *img, int width, int height){

  /* FILE pointer */
  FILE *f;

  /* Open file for writing */
  f = fopen(filename, "wb");

  /* PPM header info, including the size of the image */
  fprintf(f, "P6 %d %d %d\n", width, height, 255);

  /* Write the image data to the file - remember 3 byte per pixel */
  fwrite(img, 3, width*height, f);

  /* Make sure you close the file */
  fclose(f);
}



int writeMandelbrot(const char *fileName, int width, int height, float *img, int minI, int maxI){

  int n, m;
  unsigned char *rgb   = (unsigned char*) calloc(3*width*height, sizeof(unsigned char));

  for(n=0;n<height;++n){
    for(m=0;m<width;++m){
      int id = m+n*width;
      int I = (int) (768*sqrt((double)(img[id]-minI)/(maxI-minI)));

      // change this to change palette
      if(I<256)      rgb[3*id+2] = 255-I;
      else if(I<512) rgb[3*id+1] = 511-I;
      else if(I<768) rgb[3*id+0] = 767-I;
      else if(I<1024) rgb[3*id+0] = 1023-I;
      else if(I<1536) rgb[3*id+1] = 1535-I;
      else if(I<2048) rgb[3*id+2] = 2047-I;

    }
  }

  saveppm(fileName, rgb, width, height);

  free(rgb);
}
