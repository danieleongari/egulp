#include "eseqparam.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "elements.h"
#include "psinv.h"
#include "utilsa.h"

#define MAXSTRL (128)

void InitESEQParams(double *qeqchi, double *qeqmui, int *qeqids, int MaxElement,
                    char *fname, tESEQParams *eseqptr) {
  FILE *fp;
  int i, itype, jtype;
  int MaxParam;
  int ElementNumber;
  int ElNum1;
  int ElNum2;

  int MaxPairs;
  int el_id;
  char tmpstr1[MAXSTRL];
  char tmpstr2[MAXSTRL];
  char tmpstr3[MAXSTRL];
  char tmpstr4[MAXSTRL];

  fp = fopen(fname, "r");
  printf("opening file with parameters %s\n", fname);
  if (fp == NULL) {
    printf("Error opening file %s\n", fname);
    exit(-1);
  } else {
    printf("parameter file %s is open\n", fname);
    fscanf(fp, "%s\n", tmpstr1);
    MaxParam = atoi(tmpstr1);
    if ((MaxParam < 1) || (MaxParam > MaxElement)) {
      printf("Invalid number of parameters %d...\n", MaxParam);
      printf("should be in the (1..%d) range\n", MaxElement);
      exit(-1);
    }
    for (i = 0; i < MaxParam; i++) {
      fscanf(fp, "%s %s %s\n", tmpstr1, tmpstr2, tmpstr3);
      ElementNumber = atoi(tmpstr1);
      el_id = element_index(ElementNumber, qeqids, MaxElement);
      qeqchi[el_id] = atof(tmpstr2);
      qeqmui[el_id] = atof(tmpstr3);
    }
    fscanf(fp, "%s\n", tmpstr1);
    eseqptr->number_of_bond_types = atoi(tmpstr1);
    MaxPairs = (MaxElement) * (MaxElement + 1) / 2;
    if ((eseqptr->number_of_bond_types < 1) ||
        (eseqptr->number_of_bond_types > MaxPairs)) {
      printf("Invalid number of bond types %d...\n",
             eseqptr->number_of_bond_types);
      printf("should be in the (1..%d) range\n", MaxPairs);
      exit(-1);
    }

    allocate_int_array(&(eseqptr->bond_atom_types),
                       2 * eseqptr->number_of_bond_types,
                       "can not allocate bond_atom_types\n");
    allocate_double_array(&(eseqptr->elneg_bondcorr),
                          eseqptr->number_of_bond_types,
                          "can not allocate elneg_bondcorr\n");
    allocate_double_array(&(eseqptr->hard_bondcorr),
                          eseqptr->number_of_bond_types,
                          "can not allocate hard_bondcorr\n");

    /* the same parameters but in the geometry order */
    allocate_double_array(&(eseqptr->elneg_bondcorr_geom),
                          eseqptr->number_of_bond_types,
                          "can not allocate elneg_bondcorr_geom\n");
    allocate_double_array(&(eseqptr->hard_bondcorr_geom),
                          eseqptr->number_of_bond_types,
                          "can not allocate hard_bondcorr_geom\n");

    for (i = 0; i < eseqptr->number_of_bond_types; i++) {
      eseqptr->elneg_bondcorr_geom[i] = 0.0;
      eseqptr->hard_bondcorr_geom[i] = 0.0;
      fscanf(fp, "%s %s %s %s\n", tmpstr1, tmpstr2, tmpstr3, tmpstr4);
      ElNum1 = atoi(tmpstr1);
      if ((ElNum1 < 1) || (ElNum1 > MaxElement)) {
        printf("bond line %d:\n", i + 1);
        printf("invalid first element number %d\n", ElNum1);
        exit(-1);
      }

      ElNum2 = atoi(tmpstr2);
      if ((ElNum2 < 1) || (ElNum2 > MaxElement)) {
        printf("bond line %d:\n", i + 1);
        printf("invalid second element number %d\n", ElNum2);
        exit(-1);
      }
      if (ElNum1 < ElNum2) {
        eseqptr->bond_atom_types[2 * i + 0] = ElNum1;
        eseqptr->bond_atom_types[2 * i + 1] = ElNum2;
      } else {
        eseqptr->bond_atom_types[2 * i + 0] = ElNum2;
        eseqptr->bond_atom_types[2 * i + 1] = ElNum1;
      }
      eseqptr->elneg_bondcorr[i] = atof(tmpstr3);
      eseqptr->hard_bondcorr[i] = atof(tmpstr4);
    }

    for (itype = 0; itype < eseqptr->number_of_bond_types; itype++)
      for (jtype = itype + 1; jtype < eseqptr->number_of_bond_types; jtype++)
        if ((eseqptr->bond_atom_types[2 * itype] ==
             eseqptr->bond_atom_types[2 * jtype]) &&
            (eseqptr->bond_atom_types[2 * itype + 1] ==
             eseqptr->bond_atom_types[2 * jtype] + 1)) {
          printf("identical type number %d and %d\n", itype, jtype);
          exit(-1);
        }
    printf("closing parameter file\n");
    fclose(fp);
  }
  return;
}

void FreeESEQParams(tESEQParams *eseqptr) {
  free(eseqptr->bond_atom_types);
  free(eseqptr->elneg_bondcorr);
  free(eseqptr->hard_bondcorr);
  free(eseqptr->elneg_bondcorr_geom);
  free(eseqptr->hard_bondcorr_geom);

  return;
}

void InitTransferMatrices(tESEQParams *eseqptr, int ncharges, int nbonds,
                          int *nbonds_indexij) {
  int ibond, iatom;
  int atom_i, atom_j;
  allocate_double_array(&(eseqptr->TMatrix), ncharges * nbonds,
                        "can not allocate TMatrix\n");
  allocate_double_array(&(eseqptr->TMatrixPINV), ncharges * nbonds,
                        "can not allocate TMatrixPINV\n");
  for (iatom = 0; iatom < ncharges; iatom++) {
    for (ibond = 0; ibond < nbonds; ibond++) {
      atom_i = nbonds_indexij[2 * ibond + 0];
      atom_j = nbonds_indexij[2 * ibond + 1];
      if ((iatom + 1) == atom_i)
        eseqptr->TMatrix[iatom * nbonds + ibond] = 1.0;
      else if ((iatom + 1) == atom_j)
        eseqptr->TMatrix[iatom * nbonds + ibond] = -1.0;
      else
        eseqptr->TMatrix[iatom * nbonds + ibond] = 0.0;
    }
  }
  psinvmat(eseqptr->TMatrix, ncharges, nbonds, eseqptr->TMatrixPINV, 0);
  return;
}

void DestroyTransferMatrices(tESEQParams *eseqptr) {
  free(eseqptr->TMatrix);
  free(eseqptr->TMatrixPINV);
  return;
}

void BuildParamsGEOOrder(tESEQParams *eseqptr, int tbonds,
                         int *tbonds_indexij) {
  int geo_bondt;
  int inp_bondt;
  int l_bond_type_found;
  for (geo_bondt = 0; geo_bondt < tbonds; geo_bondt++) {
    l_bond_type_found = 0;
    for (inp_bondt = 0; inp_bondt < eseqptr->number_of_bond_types;
         inp_bondt++) {
      if ((eseqptr->bond_atom_types[2 * inp_bondt + 0] ==
           tbonds_indexij[2 * geo_bondt + 0]) &&
          (eseqptr->bond_atom_types[2 * inp_bondt + 1] ==
           tbonds_indexij[2 * geo_bondt + 1])) {
        l_bond_type_found = 1;
        eseqptr->elneg_bondcorr_geom[geo_bondt] =
            eseqptr->elneg_bondcorr[inp_bondt];
        eseqptr->hard_bondcorr_geom[geo_bondt] =
            eseqptr->hard_bondcorr[inp_bondt];
        break;
      }
    }
    if (l_bond_type_found != 1) {
      printf(
          "geometry bond of type (%d %d) is not found among input parameters\n",
          tbonds_indexij[2 * geo_bondt + 0], tbonds_indexij[2 * geo_bondt + 1]);
      exit(-1);
    }
  }
  printf("parameters in input order: type electronegativity hardness\n");
  for (inp_bondt = 0; inp_bondt < eseqptr->number_of_bond_types; inp_bondt++)
    printf("bond type %2d %2d %10.6f %10.6f\n",
           eseqptr->bond_atom_types[2 * inp_bondt + 0],
           eseqptr->bond_atom_types[2 * inp_bondt + 1],
           eseqptr->elneg_bondcorr[inp_bondt],
           eseqptr->hard_bondcorr[inp_bondt]);
  printf("parameters in geometry order: type electronegativity hardness\n");
  for (geo_bondt = 0; geo_bondt < tbonds; geo_bondt++)
    printf("bond type %2d %2d %10.6f %10.6f\n",
           tbonds_indexij[2 * geo_bondt + 0], tbonds_indexij[2 * geo_bondt + 1],
           eseqptr->elneg_bondcorr_geom[geo_bondt],
           eseqptr->hard_bondcorr_geom[geo_bondt]);
  return;
}

void BuildBondHCorrections(tESEQParams *eseqptr, int natoms, int nbonds,
                           int *nbonds_indexij, int *nbonds_type, int tbonds,
                           int *tbonds_indexij, double *BondH) {
  double *Jprime;
  double *JT;
  int ibond, jbond, icharge, jcharge;
  int i, k;
  int bond_type;
  double sm;
  double hardness;

  allocate_double_array(&Jprime, nbonds * nbonds, "can not allocate Jprime\n");
  printf("bond correction matrix...\n");
  for (ibond = 0; ibond < nbonds; ibond++) {
    bond_type = nbonds_type[ibond];
    for (jbond = 0; jbond < nbonds; jbond++) {
      if (ibond == jbond)
        Jprime[ibond * nbonds + jbond] = eseqptr->hard_bondcorr_geom[bond_type];
      else
        Jprime[ibond * nbonds + jbond] = 0.0;
    }
  }

  allocate_double_array(&JT, natoms * nbonds, "can not allocate JT\n");

  for (ibond = 0; ibond < nbonds; ibond++) {
    for (icharge = 0; icharge < natoms; icharge++) {
      sm = 0.0;
      for (k = 0; k < nbonds; k++)
        sm = sm + Jprime[ibond * nbonds + k] *
                      eseqptr->TMatrixPINV[k * natoms + icharge];
      JT[ibond * natoms + icharge] = sm;
    }
  }

  for (icharge = 0; icharge < natoms; icharge++) {
    for (jcharge = 0; jcharge < natoms; jcharge++) {
      sm = 0.0;
      for (k = 0; k < nbonds; k++)
        sm = sm + (eseqptr->TMatrixPINV[k * natoms + icharge]) *
                      JT[k * natoms + jcharge];
      BondH[icharge * natoms + jcharge] = sm;
    }
  }

  free(JT);
  free(Jprime);
  return;
}

/*     BuildBondECorrections(&eseq, geometry.natoms, geometry.nbonds,
   geometry.nbonds_indexij, geometry.nbonds_type, geometry.tbonds,
   geometry.tbonds_indexij, BondE);  */

void BuildBondECorrections(tESEQParams *eseqptr, int natoms, int nbonds,
                           int *nbonds_indexij, int *nbonds_type, int tbonds,
                           int *tbonds_indexij, double *BondE) {
  int icharge, ibonds;
  int bond_type;
  double electroneg;
  for (icharge = 0; icharge < natoms; icharge++) {
    BondE[icharge] = 0.0;
    for (ibonds = 0; ibonds < nbonds; ibonds++) {
      bond_type = nbonds_type[ibonds];
      electroneg = eseqptr->elneg_bondcorr_geom[bond_type];
      BondE[icharge] = BondE[icharge] +
                       electroneg * eseqptr->TMatrix[icharge * nbonds + ibonds];
    }
  }
  return;
}

#undef MAXSTRL
