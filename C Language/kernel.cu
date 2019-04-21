#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
#include <stdio.h>
#include <stdlib.h>
#include "cuComplex.h"
#include "fft_naive.h"

#define SIZE	(int)pow(2,10) //I can do 32 threads per block BUT I am not sure how exactly this number works

#ifndef M_PI
#define M_PI acos(-1.0)
#endif

int main()
{
	cuDoubleComplex * original_vector = (cuDoubleComplex *)calloc(SIZE, sizeof(cuDoubleComplex));
	for (int j = 0; j < SIZE; j++)
		original_vector[j] = make_cuDoubleComplex(j + 1, 0);

	cuDoubleComplex * d_init_vec;
	cudaMalloc(&d_init_vec, sizeof(cuDoubleComplex)); //initialise the CUDA subsystem

	cudaFree(d_init_vec);

	original_vector = fft_wrapper(original_vector, SIZE);

	printf("done");
	cudaFree(original_vector);
	return 0;
}
