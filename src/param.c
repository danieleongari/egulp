#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "elements.h" 

#define MAXSTRL (128) 



void  InitParameters(double *qeqchi,  double *qeqmui, int *qeqids, int MaxElement, char *fname)
{
     FILE *fp; 
     int i, MaxParam, ElementNumber; 
     int el_id; 
     char tmpstr1[MAXSTRL]; 
     char tmpstr2[MAXSTRL]; 
     char tmpstr3[MAXSTRL];
   
     fp = fopen(fname, "r"); 
     printf("opening file with parameters %s\n", fname); 
     if(fp == NULL) {  
            printf("Error opening file %s\n", fname);  
	    exit(-1);  
     }   else { 	 
	 printf("parameter file %s is open\n", fname); 
	 fscanf(fp, "%s\n", tmpstr1);
	 MaxParam = atoi(tmpstr1); 
	 if((MaxParam < 1) || (MaxParam > MaxElement))  { 
	        printf("Number of parameters is %d...\n", MaxParam); 
		printf("should be in the (1..%d) range\n", MaxElement); 
		exit(-1); 
	 } 
/*	 printf("MaxElement = %d MaxParam=%d", MaxElement, MaxParam); */
	 for(i = 0; i   < MaxParam; i++) { 
		 fscanf(fp, "%s %s %s\n", tmpstr1, tmpstr2, tmpstr3);  
		 ElementNumber = atoi(tmpstr1);
		 el_id = element_index(ElementNumber, qeqids, MaxElement);
	 	 qeqchi[el_id] = atof(tmpstr2); 
		 qeqmui[el_id] = atof(tmpstr3); 
	 } 
	 printf("closing parameter file\n"); 
	 fclose(fp); 
     } 
     return; 
}

void PrintParameters(double *qeqchi,  double *qeqmui, double *qeqrad, int *qeqids, int MaxElement) 
{
       int i; 
       printf("Element  Electronegativity chi(eV)  Hardness 0.5J(eV) SIZE(A)\n");  
       for (i = 0; i < MaxElement; i++)  printf("%4d   %12.8f   %12.8f %12.8f\n", qeqids[i], qeqchi[i], qeqmui[i], qeqrad[i]);      
       return; 
} 

#undef MAXSTRL
