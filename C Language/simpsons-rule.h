#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double simpsonRule(double * exp_vec, double * dis_vec, double tar, int arr_size)
{
	int a = 0;
	int b = arr_size-1;
	
	double amp_diffAB = exp_vec[b] - (exp_vec[b] - exp_vec[a])/2;
	double dis_diffAB = dis_vec[b] - dis_vec[a];
	double areaAB = amp_diffAB*dis_diffAB;
	
	double amp_diffAC = exp_vec[b/2] - (exp_vec[b/2] - exp_vec[a])/2;
	double dis_diffAC = dis_vec[b/2] - dis_vec[a];
	double areaAC = amp_diffAC*dis_diffAC;
	
	double amp_diffCB = exp_vec[b] - (exp_vec[b] - exp_vec[b/2])/2;
	double dis_diffCB = dis_vec[b] - dis_vec[b/2];
	double areaCB = amp_diffCB*dis_diffCB;
	
	double error = 0.1*absDouble(areaAB-areaAC-areaCB);
	
	if(error < tar)
	{
		printf("a is %0.4f, b is %0.4f\n",dis_vec[a],dis_vec[b]);
		return absDouble(areaAB);
	}

	
	else
	{
		double * left_exp = exp_vec;
		double * left_dis = dis_vec;
		double * right_exp = &exp_vec[b/2];
		double * right_dis = &dis_vec[b/2];
		
		double left = simpsonRule(left_exp, left_dis, tar, arr_size/2);
		double right = simpsonRule(right_exp, right_dis, tar, arr_size/2);
		
		left_exp = left_dis = right_exp = right_dis = NULL;
	
		return (left+right);
	}
}

double absDouble(double number)
{
	if (number < 0)
		return number*(-1);
	else
		return number;
}
