#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void psinv_small_row_large_column(double *mat, int n_fil, int n_col,
                                  double *psinv_matrix, int ldebug);
void psinv_large_row_small_column(double *mat, int n_fil, int n_col,
                                  double *psinv_matrix, int ldebug);

void psinvmat(double *mat, int n_fil, int n_col, double *psinv_matrix,
              int ldebug) {
  if (n_fil < n_col) {
    psinv_small_row_large_column(mat, n_fil, n_col, psinv_matrix, ldebug);
  } else {
    psinv_large_row_small_column(mat, n_fil, n_col, psinv_matrix, ldebug);
  }
  return;
}

void psinv_small_row_large_column(double *mat, int n_fil, int n_col,
                                  double *psinv_matrix, int ldebug) {
  unsigned i = 0;
  unsigned j = 0;
  gsl_matrix *gA = gsl_matrix_alloc(n_fil, n_col);
  for (i = 0; i < n_fil; i++)
    for (j = 0; j < n_col; j++) gsl_matrix_set(gA, i, j, mat[i * n_col + j]);

  /*       gsl_matrix_set(gA, 0, 0,  1.0);
         gsl_matrix_set(gA, 0, 1, -1.0);
         gsl_matrix_set(gA, 0, 2,  0.0);

         gsl_matrix_set(gA, 1, 0,  1.0);
         gsl_matrix_set(gA, 1, 1,  0.0);
         gsl_matrix_set(gA, 1, 2,  -1.0); */

  /* Computing the transpose of gA */
  gsl_matrix *gA_t = gsl_matrix_alloc(n_col, n_fil);
  gsl_matrix_transpose_memcpy(gA_t, gA);

  gsl_matrix *U = gsl_matrix_alloc(n_col, n_fil);
  gsl_matrix *V = gsl_matrix_alloc(n_fil, n_fil);
  gsl_vector *S = gsl_vector_alloc(n_fil);

  /* Computing the SVD of the transpose of A
  The matrix 'gA_t' will contain 'U' after the function is called */
  gsl_vector *work = gsl_vector_alloc(n_fil);
  gsl_linalg_SV_decomp(gA_t, V, S, work);
  gsl_vector_free(work);

  gsl_matrix_memcpy(U, gA_t);

  /* //Inverting S//
  //----------------------------------------------------------
  // Matrix 'S' is diagonal, so it is contained in a vector.
  // We operate to convert the vector 'S' into the matrix 'Sp'.
  //Then we invert 'Sp' to 'Spu'
  //----------------------------------------------------------*/
  gsl_matrix *Sp = gsl_matrix_alloc(n_fil, n_fil);
  gsl_matrix_set_zero(Sp);
  for (i = 0; i < n_fil; i++) gsl_matrix_set(Sp, i, i, gsl_vector_get(S, i));

  gsl_matrix *SI = gsl_matrix_calloc(n_fil, n_fil);

  for (i = 0; i < n_fil; i++) {
    if (ldebug == 1) printf("S [%d] = %12.8f\n", i, gsl_vector_get(S, i));
    if (gsl_vector_get(S, i) > 0.0000000001)
      gsl_matrix_set(SI, i, i, 1.0 / gsl_vector_get(S, i));
  }

  gsl_matrix *VT = gsl_matrix_alloc(n_fil, n_fil);
  gsl_matrix_transpose_memcpy(VT, V);
  /* // Tranpose of V */

  /*
          //THE PSEUDOINVERSE//
          //----------------------------------------------------------
          //Computation of the pseudoinverse of trans(A) as pinv(A) =
     U·inv(S).trans(V)   with trans(A) = U.S.trans(V)
          //----------------------------------------------------------
  */
  gsl_matrix *SIpVT = gsl_matrix_alloc(n_fil, n_fil);
  /* Calculating inv(S).trans(V) */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, SI, VT, 0.0, SIpVT);
  /* Calculating U inv(S).trans(V) */
  gsl_matrix *pinv = gsl_matrix_alloc(n_col, n_fil);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, SIpVT, 0.0, pinv);

  for (i = 0; i < n_col; i++)
    for (j = 0; j < n_fil; j++)
      psinv_matrix[i * n_fil + j] = gsl_matrix_get(pinv, i, j);

  if (ldebug == 1) {
    printf("pinv:\n");
    for (i = 0; i < n_col; i++)
      for (j = 0; j < n_fil; j++)
        printf("m(%d,%d) = %g\n", i, j, gsl_matrix_get(pinv, i, j));
    printf("\n");
  }

  gsl_matrix_free(gA);
  gsl_matrix_free(gA_t);

  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(S);

  gsl_matrix_free(Sp);
  gsl_matrix_free(SI);

  gsl_matrix_free(VT);
  gsl_matrix_free(SIpVT);
  gsl_matrix_free(pinv);
  return;
}

void psinv_large_row_small_column(double *mat, int n_fil, int n_col,
                                  double *psinv_matrix, int ldebug) {
  unsigned i = 0;
  unsigned j = 0;
  gsl_matrix *gA = gsl_matrix_alloc(n_fil, n_col);
  for (i = 0; i < n_fil; i++)
    for (j = 0; j < n_col; j++) gsl_matrix_set(gA, i, j, mat[i * n_col + j]);

  gsl_matrix *U = gsl_matrix_alloc(n_fil, n_col);
  gsl_matrix *V = gsl_matrix_alloc(n_col, n_col);
  gsl_vector *S = gsl_vector_alloc(n_col);

  /* Computing the SVD of the transpose of A
  The matrix 'gA_t' will contain 'U' after the function is called */
  gsl_vector *work = gsl_vector_alloc(n_col);
  gsl_linalg_SV_decomp(gA, V, S, work);
  gsl_vector_free(work);

  gsl_matrix_memcpy(U, gA);

  /* //Inverting S//
  //----------------------------------------------------------
  // Matrix 'S' is diagonal, so it is contained in a vector.
  // We operate to convert the vector 'S' into the matrix 'Sp'.
  //----------------------------------------------------------*/
  gsl_matrix *Sp = gsl_matrix_alloc(n_col, n_col);
  gsl_matrix_set_zero(Sp);
  for (i = 0; i < n_col; i++) gsl_matrix_set(Sp, i, i, gsl_vector_get(S, i));

  gsl_matrix *SI = gsl_matrix_calloc(n_col, n_col);

  for (i = 0; i < n_col; i++) {
    if (ldebug == 1) printf("S [%d] = %12.8f\n", i, gsl_vector_get(S, i));

    if (gsl_vector_get(S, i) > 0.0000000001)
      gsl_matrix_set(SI, i, i, 1.0 / gsl_vector_get(S, i));
  }

  gsl_matrix *UT = gsl_matrix_alloc(n_col, n_fil);
  gsl_matrix_transpose_memcpy(UT, U);
  /* // Tranpose of V */

  /*
          //THE PSEUDOINVERSE//
          //----------------------------------------------------------
          //Computation of the pseudoinverse of trans(A) as pinv(A) =
     U·inv(S).trans(V)   with trans(A) = U.S.trans(V)
          //----------------------------------------------------------
  */
  gsl_matrix *SIpUT = gsl_matrix_alloc(n_col, n_fil);
  /* Calculating inv(S).trans(U) */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, SI, UT, 0.0, SIpUT);
  /* Calculating V inv(S).trans(U) */
  gsl_matrix *pinv = gsl_matrix_alloc(n_col, n_fil);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, SIpUT, 0.0, pinv);

  for (i = 0; i < n_col; i++)
    for (j = 0; j < n_fil; j++)
      psinv_matrix[i * n_fil + j] = gsl_matrix_get(pinv, i, j);
  printf("\n");

  if (ldebug == 1) {
    printf("pinv:\n");
    for (i = 0; i < n_col; i++)
      for (j = 0; j < n_fil; j++)
        printf("m(%d,%d) = %g\n", i, j, gsl_matrix_get(pinv, i, j));
    printf("\n");
  }

  gsl_matrix_free(gA);
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_matrix_free(Sp);
  gsl_matrix_free(SI);

  gsl_matrix_free(UT);
  gsl_matrix_free(SIpUT);
  gsl_matrix_free(pinv);
  return;
}
