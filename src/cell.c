#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void vectorProduct(double *a2, double *a3, double *res)
{ 
	res[0] = (a2[1]*a3[2]-a3[1]*a2[2]);
	res[1] = -(a2[0]*a3[2]-a3[0]*a2[2]);
	res[2] = (a2[0]*a3[1]-a3[0]*a2[1]);
        return; 
}
	
double  cellVolume(double *a1, double *a2,  double *a3)
{
        double a_x, a_y, a_z; 
	double b_x, b_y, b_z; 
	double c_x, c_y, c_z; 
	double vol; 
	a_x = a1[0];
	a_y = a1[1];
	a_z = a1[2];
	b_x = a2[0];
	b_y = a2[1];
	b_z = a2[2];
	c_x = a3[0];
	c_y = a3[1];
	c_z = a3[2];
	vol = fabs(a_x*(b_y*c_z-c_y*b_z)-a_y*(b_x*c_z-c_x*b_z)+a_z*(b_x*c_y-c_x*b_y)); 
	return(vol); 
} 
