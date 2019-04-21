/*
	Recursive FFT and IFFT library
*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI acos(-1.0)
#endif

double complex * fft_conquer_func(double complex * original_vector, int num_samples)
{
	double complex * fft_vector = (double complex *)calloc(num_samples, sizeof(double complex));
	double complex wn = cexp(-2*M_PI*I/num_samples);
	int exponent;
	for(int k = 0; k < num_samples; k++)
	{
		for(int j = 0; j < num_samples; j++)
		{
			exponent = j*k;
			fft_vector[k] = fft_vector[k] + original_vector[j]*cpow(wn,exponent);
		}
		
	}
	//Copy the FFT vector to the passed pointer
	for(int j = 0; j < num_samples; j++)
	{
		*original_vector = *fft_vector;
		original_vector++;
		fft_vector++;
	}
	original_vector = original_vector - num_samples;
	fft_vector = fft_vector - num_samples;
	free(fft_vector);
	return original_vector;
}

double complex * fft_divide_func(double complex * original_vector, int num_samples)
{
	if(num_samples <= 2)
		return fft_conquer_func(original_vector, num_samples);
	else
	{	
		double complex factor[num_samples];
		double complex * odd = (double complex *)calloc(num_samples/2, sizeof(double complex));	
		double complex * even = (double complex *)calloc(num_samples/2, sizeof(double complex));
		double complex comb[num_samples];
		//Split the original vector
		int k = 0;
		for(int j = 0; j < num_samples; j++)
		{
			if(j % 2 == 0)
				even[k] = original_vector[j];
			else
			{
				odd[k] = original_vector[j]; 
				k++;
			}		
		}		
		even = fft_divide_func(even, num_samples/2);
		odd = fft_divide_func(odd, num_samples/2);
	
		//Calculate the factor
		for(int k = 0; k < num_samples; k++)
			factor[k] = cexp(-2*M_PI*I*k/num_samples);
		
		//Multiply out
		k = 0;
		while(k < num_samples)
		{
			for(int j = 0; j < num_samples/2; j++)
			{
				original_vector[k] = even[j] + odd[j]*factor[k];
				k++;
			}
		}		
		free(odd);
		free(even);	
		return original_vector;
	}	
}

double complex * fft_wrapper(double complex * original_vector, int num_samples)
{
	original_vector = fft_divide_func(original_vector, num_samples);
	return original_vector;
}

double complex * ifft_conquer_func(double complex * original_vector, int num_samples)
{
	double complex * ifft_vector = (double complex *)calloc(num_samples, sizeof(double complex));
	double complex wn = cexp(-2*M_PI*I/num_samples);
	int exponent;
	//int p = 0;
	for(int k = 0; k < num_samples; k++)
	{
		for(int j = 0; j < num_samples; j++)
		{
			exponent = -j*k;
			ifft_vector[k] = ifft_vector[k] + original_vector[j]*cpow(wn,exponent);
		}
	}
	
	for(int j = 0; j < num_samples; j++)
	{
		*original_vector = *ifft_vector;
		original_vector++;
		ifft_vector++;
	}
	original_vector = original_vector - num_samples;
	ifft_vector = ifft_vector - num_samples;
	free(ifft_vector);
	return original_vector;
}
	
double complex * ifft_divide_func(double complex * original_vector, int num_samples)
{
	if(num_samples <= 2)
	{
		return ifft_conquer_func(original_vector, num_samples);
	}
	else
	{		
		double complex factor[num_samples];
		double complex * odd = (double complex *)calloc(num_samples/2, sizeof(double complex));
		double complex * even = (double complex *)calloc(num_samples/2, sizeof(double complex));
		double complex comb[num_samples];
		
		//Split the original vector
		int k = 0;
		for(int j = 0; j < num_samples; j++)
		{
			if(j % 2 == 0)
			{
				even[k] = original_vector[j];
			}
			else
			{
				odd[k] = original_vector[j]; 
				k++;
			}		
		}
		
		even = ifft_divide_func(even, num_samples/2);
		odd = ifft_divide_func(odd, num_samples/2);
	
		//Calculate the factor
		for(int k = 0; k < num_samples; k++)
		{
			factor[k] = cexp(-2*M_PI*I*k/num_samples);
		}
		
		//Multiply out
		k = 0;
		while(k < num_samples)
		{
			for(int j = 0; j < num_samples/2; j++)
			{
				original_vector[k] = even[j] + odd[j]*factor[k];
				k++;
			}
		}	
		
		free(odd);
		free(even);	
		return original_vector;
	}
	
}

double complex * ifft_wrapper(double complex * original_vector, int num_samples)
{
	double complex * temp_vector = (double complex *)calloc(num_samples-1, sizeof(double complex));
	original_vector = ifft_divide_func(original_vector, num_samples);	
	for(int k = 0; k < num_samples; k++)
	{
		original_vector[k] = original_vector[k] * (1.0/num_samples);
	}
	
	int i = 1;
	for(int k = 1; k < num_samples/2; k++)
	{
		temp_vector[k-1] = original_vector[k];
		original_vector[k] = original_vector[num_samples - i];
		original_vector[num_samples - i] = temp_vector[k-1];
		i++;
	}
	
	free(temp_vector);	
	return original_vector;
}
