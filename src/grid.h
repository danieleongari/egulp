#ifndef _GRID_H_
#define _GRID_H_

typedef struct tGridData {
  /* cube stuff for  grid */
  double resolution[3];
  double vdw_factor_i;
  double vdw_factor_f;
  int use_vdw_factor;
  double offset;

  int gridn1;
  int gridn2;
  int gridn3;

  double gridh1[3];
  double gridh2[3];
  double gridh3[3];

  double dv3;

  int gridnl;

  int *gridL1;
  int *gridL2;
  int *gridL3;
  int *gridLxyz_inv;
} tGrid;

void InitGrid(tGrid *gridptr, double res1, double res2, double res3,
              double vdw_factor_i, double vdw_factor_f, int use_vdw_factor,
              double offset, double *av1, double *av2, double *av3, int natoms,
              int *atoms_num, double *atoms_xyz);

void WriteGridAsCube(char *fName, tGrid grid, int natoms, int *atoms_num,
                     double *atoms_xyz);
void WriteFunctionAsCube(char *fName, tGrid grid, int natoms, int *atoms_num,
                         double *atoms_xyz, double *func);

void DestroyGrid(tGrid *gridptr);

void InitGridFromCube(char *fName, tGrid *newgridptr, int natoms);
double IntegrateFromGrid(int nl, double *f, double dv3);

#endif
