#ifndef _PSINV_H_
#define _PSINV_H_

void psinvmat(double *mat, int n_fil, int n_col, double *psinv_matrix,
              int ldebug);

/* ###example###
int main(void) {
        int n_fil, n_col;
        int i, j;
        double mat[6] =    {1.0, -1.0, 0.0, 1.0, 0.0, -1.0};
        double psinv[6] = {0.0,   0.0, 0.0, 0.0, 0.0,  0.0};
        printf("psinv test\n");
        printf("2 row x 3 column matrix\n");
        n_fil = 2;
        n_col = 3;
        psinvmat(mat, n_fil, n_col, psinv, 1);
        printf("pseudo inverse\n");
        for(i = 0; i < n_col; i++)
            for(j = 0; j < n_fil; j++)
                 printf("row %3d column %3d : %10.6f\n", i+1, j+1,
psinv[i*n_fil+j]);

        printf("3 row x 2 column matrix\n");
        n_fil = 3;
        n_col = 2;
        psinvmat(mat, n_fil, n_col, psinv, 1);
        printf("pseudo inverse\n");
        for(i = 0; i < n_col; i++)
            for(j = 0; j < n_fil; j++)
                 printf("row %3d column %3d : %10.6f\n", i+1, j+1,
psinv[i*n_fil+j]);


        return(0);
} */

#endif
