#ifndef _GENPOT_H_
#define _GENPOT_H_

void    genpot(double *av1, double *av2, double *av3,  
  int natoms, int *atoms_num, double *atoms_charge, 
  double *atoms_xyz, double *J, 
  double *qeqrad, int *qeqids, int *qeqbasetype, int MaxEl, int *CalcMatrix); 
  
void   genpot_on_grid(double *av1, double *av2, double *av3,  
  int natoms, double *atoms_charge,  double *atoms_xyz, 
  int gridnl, double *pot, int *L1, int *L2, int *L3, double *h1, double *h2, double *h3); 
  

void   genpot_on_points(double *av1, double *av2, double *av3,  
  int nsource, double *source_charge,  double *source_xyz, 
  int npoints_to_gen, double *pot,  double *points_to_gen_xyz, int natoms); 
  

#endif
