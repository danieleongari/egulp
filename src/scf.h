#ifndef _SCF_H_
#define _SCF_H_

typedef struct tEGULPScf {
  int niter;
  int literate;
  int lconverged;
  int natoms;
  int current_iter;
  double *atoms_charge;
  double *atoms_charge_old;
  double qeqscfcrit;
  double MixParameter;
  double *energy;
  double *QEQ_JArray;
  double *QEQ_ChiArray;
} tSCF;

void InitEGULPScf(tSCF *SCF, double *atoms_charge, int *atom_num,
                  int NumberOfAtoms);
void PrintEGULPScf(tSCF *SCF);
void DestroyEGULPSCF(tSCF *SCF);
double getQEQEnergy(int natoms, double *atoms_charge, double *Jarray,
                    double *chi);
double getQEQPCEnergy(int natoms, double *atoms_charge, double *pcpot);

#endif
