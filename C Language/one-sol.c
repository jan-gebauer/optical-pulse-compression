// Implementation of the whole SSFM code
//

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#include "fft-ifft-recursive.h"

#define NSAMPLES (int)pow(2,10)
#ifndef M_PI
#define M_PI acos(-1.0)
#endif

int main()
{
	
	double b2 = -10;
	double gamma = 0.1;
	double To = 10;
	double T = 20 * To;
	double C = 0;

	double N_order = 1;
	double Po = pow(N_order, 2);
	double Ld = pow(To, 2) / abs(b2);
	double Lnl = 1/(Po*gamma);

	double dz = Ld/100;;
	double L = 10*Ld;


	int z_size = L / dz;
	double z_vector[z_size];
	for(int i = 0; i < z_size; i++)
	{
		z_vector[i] = i*dz;
	}

										


	double Fs = (NSAMPLES - 1) / T;
	double df = 2*M_PI / T;
	double t[NSAMPLES];
	double f[NSAMPLES];
	for (int i = 0; i < NSAMPLES; i++)
	{
		t[i] = ((-NSAMPLES/2) + i)*(1/Fs);
		f[i] = ((-NSAMPLES/2) + i)*df;
	}	
	//Frequency vector shift
	for(int k = 0; k < NSAMPLES/2; k++)
	{
		double temp[NSAMPLES];
		
		temp[NSAMPLES-1] = f[0];
		for(int i = 0; i < NSAMPLES-1; i++)
			temp[i] = f[i+1];			
		
		for(int i = 0; i < NSAMPLES; i++)
			f[i] = temp[i];
	}


	printf("b2 is %0.2f, z_size is %d, ", b2,z_size);
	printf("dz is %0.2f, L is %0.2f, Ld is %0.2f\n", dz,L,Ld);

	double complex * simul_wave = (double complex *)calloc(NSAMPLES, sizeof(double complex));
	double complex * simul_wave_storage = (double complex *)calloc(z_size*NSAMPLES, sizeof(double complex));
	double complex * A = (double complex *)calloc(NSAMPLES, sizeof(double complex));
	
	for(int j = 0; j < NSAMPLES; j++)
	{
		A[j] = N_order*(1/ccosh(t[j]/To));
		simul_wave[j] = cpow(cabs(A[j]),2);	
	}
	


	
	
	
	//Find can be a simple loop instead, so can max
	/*fwhm1=find(abs(A)>abs(max(A)/2));
	fwhm1=length(fwhm1);
	*/
	//Loop here?
	
	printf("Gonna do an SSFM \n");

	double complex * D = (double complex *)calloc(NSAMPLES, sizeof(double complex));
	int total_track = 0;
	for(int i = 0; i < z_size; i++)
	{
		//Compute N
		double complex N[NSAMPLES];
		for(int j = 0; j < NSAMPLES; j++)
		{
			N[j] = cexp(I*gamma*dz*cpow(cabs(A[j]),2));
		}
		
		//Fourier transform of A
		A = fft_wrapper(A, NSAMPLES);
		
		//Compute D
		for(int j = 0; j < NSAMPLES; j++)
		{
			D[j] = cexp(I*(dz/2)*b2*cpow(f[j],2));
			D[j] = D[j]*A[j]; 
		}
		
		//Inverse Fourier Transform of D
		D = ifft_wrapper(D, NSAMPLES);		

		//A = D*N		
		for(int j = 0; j < NSAMPLES; j++)
		{
			A[j] = D[j]*N[j];
			simul_wave_storage[total_track] = cpow(cabs(A[j]),2);
			total_track++;
		}
	}

	printf("Gonna do a save\n");



	char * fs = fopen("D:\\Uni\\4thYear\\EG4013\\thesis-work\\MATLAB\\SSFM\\vals.csv", "w");
	fclose(fs);
	fs = fopen("D:\\Uni\\4thYear\\EG4013\\thesis-work\\MATLAB\\SSFM\\vals.csv", "a");
	
	printf("Original wave\n");
	for(int j = 0; j < NSAMPLES; j++)
	{
		fprintf(fs, "%f%+fi,", creal(simul_wave[j]), cimag(simul_wave[j]));		
	}
	fprintf(fs, "\n");
	
	for(int i = 0; i < NSAMPLES*z_size; i++)
	{
			if(i % NSAMPLES == 0)
				fprintf(fs, "\n");
			fprintf(fs, "%f%+fi,", creal(simul_wave_storage[i]), cimag(simul_wave_storage[i]));
	}
	
	free(D);
	free(simul_wave_storage);
	fclose(fs);
	return 0;
}


