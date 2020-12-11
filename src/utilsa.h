#ifndef _UTILSA_H_
#define _UTILSA_H_

void allocate_double_array(double **ap, int n, char *errmsg);
void allocate_float_array(float **ap, int n, char *errmsg);
void allocate_int_array(int **ap, int n, char *errmsg);
double **allocate_doubleptr_array(int n, char *errmsg);

char **allocate_charptr_array(int n, char *errmsg);
void allocate_char_array(char **ap, int n, char *errmsg);
double sum_doubles(double *a, int n);

#endif
