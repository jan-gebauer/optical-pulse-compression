#include <stdio.h>
#include <time.h>


#ifndef M_PI
#define M_PI acos(-1.0)
#endif

__global__ void fftCUDA(cuDoubleComplex * sum_vector, cuDoubleComplex * wn_vector, cuDoubleComplex * fft_vector, int num_samples)
{
	int tid = threadIdx.x;
	int number_of_blocks = blockDim.x / 2;
	int step = 1;
	int target, source;
	sum_vector[blockIdx.x*num_samples + tid] = cuCmul(sum_vector[blockIdx.x*num_samples + tid], wn_vector[blockIdx.x*num_samples+tid]);
	while (number_of_blocks > 0)
	{
		if (tid < number_of_blocks)
		{
			target = blockIdx.x*num_samples + tid * step * 2;
			source = target + step;
			sum_vector[target] = cuCadd(sum_vector[target], sum_vector[source]);
			
		}

		step = step * 2;
		number_of_blocks = number_of_blocks / 2;

	}
	if (number_of_blocks == 0 && tid == 0)
	{
		fft_vector[blockIdx.x] = sum_vector[blockIdx.x*num_samples + tid];
	}
}

cuDoubleComplex powFunc(cuDoubleComplex base, int power)
{
	if (power == 0)
		return make_cuDoubleComplex(1, 0);
	else
	{
		cuDoubleComplex result = make_cuDoubleComplex(cuCreal(base), cuCimag(base));

		for (int i = 0; i < power - 1; i++)
		{
			result = cuCmul(result, base);
			printf("wn thing is %0.4f%+0.4fi\n", cuCreal(result), cuCimag(result));
		}

		return result;
	}
}

int * exponentFunc(int * exponent, int num_samples)
{
	int p = 0;
	while (p < num_samples*num_samples)
	{
		for (int k = 0; k < num_samples; k++)
		{
			for (int j = 0; j < num_samples; j++)
			{
				exponent[p] = j * k;
				p++;
			}
		}
	}

	return exponent;
}

cuDoubleComplex * wn_powered_func(cuDoubleComplex * cuda_wn_vector, int num_samples)
{
	double wn_exponent = -2 * M_PI / num_samples;
	int p = 0;
	for (int k = 0; k < num_samples; k++)
	{
		for (int j = 0; j < num_samples; j++)
		{
			cuda_wn_vector[p] = my_exp_pow(wn_exponent,j*k);
			p++;
		}
	}




	return cuda_wn_vector;
}

cuDoubleComplex * fft_conquer(cuDoubleComplex * original_vector, int num_samples, cuDoubleComplex * wn_vector)
{
	cuDoubleComplex * d_wn_vector;
	cuDoubleComplex * d_original_vector;
	cuDoubleComplex * d_fft_vector;
	cudaMalloc(&d_wn_vector, num_samples * num_samples * sizeof(cuDoubleComplex));
	cudaMalloc(&d_original_vector, num_samples * num_samples * sizeof(cuDoubleComplex));
	cudaMalloc(&d_fft_vector, num_samples * sizeof(cuDoubleComplex));

	cudaMemcpy(d_wn_vector, wn_vector, num_samples*num_samples * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
	for(int i = 0; i < num_samples; i++)
		cudaMemcpy(&d_original_vector[i*num_samples], original_vector, num_samples * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

	fftCUDA << <num_samples, num_samples >> > (d_original_vector, d_wn_vector, d_fft_vector, num_samples);

	cudaMemcpy(original_vector, d_fft_vector, num_samples*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

	cudaFree(d_wn_vector);
	cudaFree(d_original_vector);
	cudaFree(d_fft_vector);

	return original_vector;

}

cuDoubleComplex * fft_divide(cuDoubleComplex * original_vector, int num_samples, cuDoubleComplex * wn_vector)
{
	if (num_samples <= 32)
	{
		return fft_conquer(original_vector, 32, wn_vector);
	}
	else
	{
		cuDoubleComplex * factor = (cuDoubleComplex *)calloc(num_samples, sizeof(cuDoubleComplex));
		cuDoubleComplex * odd = (cuDoubleComplex *)calloc(num_samples / 2, sizeof(cuDoubleComplex));
		cuDoubleComplex * even = (cuDoubleComplex *)calloc(num_samples / 2, sizeof(cuDoubleComplex));
		cuDoubleComplex * comb = (cuDoubleComplex *)calloc(num_samples, sizeof(cuDoubleComplex));
		//Split the original vector
		int k = 0;
		for (int j = 0; j < num_samples; j++)
		{
			if (j % 2 == 0)
			{
				even[k] = original_vector[j];
			}
			else
			{
				odd[k] = original_vector[j];
				k++;
			}
		}

		even = fft_divide(even, num_samples / 2, wn_vector);
		odd = fft_divide(odd, num_samples / 2, wn_vector);

		//Calculate the factor
		for (int k = 0; k < num_samples; k++)
		{
			factor[k] = my_cexpf(make_cuDoubleComplex(0, -2 * M_PI*k / num_samples));
		}

		//Multiply out
		k = 0;
		while (k < num_samples)
		{
			for (int j = 0; j < num_samples / 2; j++)
			{
				original_vector[k] = cuCadd(even[j], cuCmul(odd[j], factor[k]));
				k++;
			}
		}

		free(factor);
		free(comb);
		free(odd);
		free(even);

		return original_vector;
	}
}

cuDoubleComplex * fft_wrapper(cuDoubleComplex * original_vector, int num_samples)
{
	cuDoubleComplex * wn_vector = (cuDoubleComplex *)calloc(num_samples*num_samples, sizeof(cuDoubleComplex));
	wn_vector = wn_powered_func(wn_vector, num_samples);
	original_vector = fft_conquer(original_vector, num_samples, wn_vector);

	return original_vector;
}