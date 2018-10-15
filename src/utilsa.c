#include <stdio.h>
#include <stdlib.h>


void allocate_double_array(double **ap, int n, char *errmsg)
{
     if (n < 1) {
         printf("utils.c: invalid array dimension\n");
	 exit(1); 
     }
     (*ap) = (double *) malloc(n*sizeof(double)); 
     if (!(*ap)) {
	printf(errmsg);
	exit(1); 
     } 
}

void allocate_float_array(float **ap, int n, char *errmsg)
{
     if (n < 1) {
         printf("utils.c: invalid array dimension\n");
         exit(1); 
     }
     (*ap) = (float *) malloc(n*sizeof(float)); 
     if (!(*ap)) {
        printf("utils.c: "); 
        printf(errmsg);
        exit(1); 
     } 
}

void allocate_int_array(int **ap, int n, char *errmsg)
{
     if (n < 1) {
         printf("utils.c: invalid array dimension\n");
         exit(1); 
     }
     (*ap) = (int *) malloc(n*sizeof(int)); 
     if (!(*ap)) {
        printf("utils.c: "); 
        printf(errmsg);
        exit(1); 
     } 
}

double **allocate_doubleptr_array(int n, char *errmsg)
{
    double **res;  
    if (n < 1) {
         printf("utils.c: invalid array dimension\n");
         exit(1); 
    }
    res = (double **) malloc(n*sizeof(double *)); 
    if (!(res)) {
        printf("utils.c: "); 
        printf(errmsg);
        exit(1); 
    } 
    return(res); 
}

void allocate_char_array(char **ap, int n, char *errmsg)
{
     if (n < 1) {
         printf("utils.c: invalid array dimension\n");
	 exit(1); 
     }
     (*ap) = (char *) malloc(n*sizeof(char)); 
     if (!(*ap)) {
	printf(errmsg);
	exit(1); 
     } 
}

char **allocate_charptr_array(int n, char *errmsg)
{
    char **res;  
    if (n < 1) {
         printf("utils.c: invalid array dimension\n");
         exit(1); 
    }
    res = (char **) malloc(n*sizeof(char *)); 
    if (!(res)) {
        printf("utils.c: "); 
        printf(errmsg);
        exit(1); 
    } 
    return(res); 
}

double sum_doubles(double *a, int n)
{
          int i; 
	  double res = 0.0; 
	  for(i = 0; i < n; i++) res = res + a[i]; 
          return(res); 
}
