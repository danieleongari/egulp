#include "geometry.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cell.h"
#include "utilsa.h"
#include "vdw.h"

/* included to map atoms into the central cell */

#define MAXSTRL (128)
#define LARGED (1000.0)
#define AUTOANG (0.5291772)
/* #####################
   defined at the end of
   the file
   ##################### */

double matrix_det(double *mat, int n);
void matrix_inverse(double *mat, int n, double *imat);

const char *atsym_nospace[] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",   "F",   "Ne",   "Na", "Mg",
    "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",  "Sc",  "Ti",   "V",  "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",  "As",  "Se",   "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",  "Rh",  "Pd",   "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",  "La",  "Ce",   "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",  "Tm",  "Yb",   "Lu", "Hf",
    "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",  "Tl",  "Pb",   "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",   "Np",  "Pu",   "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "O_N", "O_S", "O_SH", "H_O"};

void InitGeometry(tGeometry *geometry, char *fname) {
  FILE *fp;
  int i, a;
  char tmpstr1[MAXSTRL];
  char atoms_sym[3];
  char fline[100];
  char *tmpstr0;
  double tempd;

  fp = fopen(fname, "r");
  printf("opening geometry file %s\n", fname);
  if (fp == NULL) {
    printf("Error opening file %s\n", fname);
    exit(-1);
  } else {
    printf("geometry file %s is open\n", fname);

    while (fgets(fline, 100, fp) != NULL) {
      // Read the cell parameters from cif
      if (strstr(fline, "_cell_length_a") != NULL) {
        tmpstr0 = strtok(fline, " ");
        tmpstr0 = strtok(NULL, " ");
        geometry->cell_len[0] = atof(tmpstr0);
      }
      if (strstr(fline, "_cell_length_b") != NULL) {
        tmpstr0 = strtok(fline, " ");
        tmpstr0 = strtok(NULL, " ");
        geometry->cell_len[1] = atof(tmpstr0);
      }
      if (strstr(fline, "_cell_length_c") != NULL) {
        tmpstr0 = strtok(fline, " ");
        tmpstr0 = strtok(NULL, " ");
        geometry->cell_len[2] = atof(tmpstr0);
      }
      if (strstr(fline, "_cell_angle_alpha") != NULL) {
        tmpstr0 = strtok(fline, " ");
        tmpstr0 = strtok(NULL, " ");
        geometry->cell_ang[0] = atof(tmpstr0);
      }
      if (strstr(fline, "_cell_angle_beta") != NULL) {
        tmpstr0 = strtok(fline, " ");
        tmpstr0 = strtok(NULL, " ");
        geometry->cell_ang[1] = atof(tmpstr0);
      }
      if (strstr(fline, "_cell_angle_gamma") != NULL) {
        tmpstr0 = strtok(fline, " ");
        tmpstr0 = strtok(NULL, " ");
        geometry->cell_ang[2] = atof(tmpstr0);
        // All cell parameters are read: print them and convert to cell matrix
        printf("Lattice lenghts: A %10.7f B %10.7f C %10.7f\n",
               geometry->cell_len[0], geometry->cell_len[1],
               geometry->cell_len[2]);
        printf("Lattice angles: alpha %10.7f beta %10.7f gamma %10.7f\n",
               geometry->cell_ang[0], geometry->cell_ang[1],
               geometry->cell_ang[2]);

        tempd = (cos(geometry->cell_ang[0] * M_PI / 180.0) -
                 cos(geometry->cell_ang[2] * M_PI / 180.0) *
                     cos(geometry->cell_ang[1] * M_PI / 180.0)) /
                sin(geometry->cell_ang[2] * M_PI / 180.0);
        geometry->av1[0] = geometry->cell_len[0];
        geometry->av1[1] = 0.0;
        geometry->av1[2] = 0.0;
        geometry->av2[0] =
            geometry->cell_len[1] * cos(geometry->cell_ang[2] * M_PI / 180.0);
        geometry->av2[1] =
            geometry->cell_len[1] * sin(geometry->cell_ang[2] * M_PI / 180.0);
        geometry->av2[2] = 0.0;
        geometry->av3[0] =
            geometry->cell_len[2] * cos(geometry->cell_ang[1] * M_PI / 180.0);
        geometry->av3[1] = geometry->cell_len[2] * tempd;
        geometry->av3[2] =
            geometry->cell_len[2] *
            sqrt(1 - pow(cos(geometry->cell_ang[1] * M_PI / 180.0), 2.0) -
                 pow(tempd, 2.0));
        printf("Lattice vector #1 %10.7f %10.7f %10.7f\n", geometry->av1[0],
               geometry->av1[1], geometry->av1[2]);
        printf("Lattice vector #2 %10.7f %10.7f %10.7f\n", geometry->av2[0],
               geometry->av2[1], geometry->av2[2]);
        printf("Lattice vector #3 %10.7f %10.7f %10.7f\n", geometry->av3[0],
               geometry->av3[1], geometry->av3[2]);
        geometry->dv = cellVolume(geometry->av1, geometry->av2, geometry->av3);
        printf("unit cell volume is %10.7f\n", geometry->dv);
      }
      // Now look for the string _atom_site_fract_x, skip the _atom_site
      // section, read the number of atoms, close the file and initialize
      // arrays.
      if (strstr(fline, "_atom_site_fract_x") != NULL) break;
    }
    int countatoms = 0;
    while (fgets(fline, 100, fp) != NULL) {
      if (strstr(fline, "_atom_site") != NULL)  // stil a header
        strcpy(tmpstr1,
               fline);  // this will remember the last string before the atoms
      if ((strstr(fline, "_atom_site") == NULL) &&
          (fline != "\n"))  // an atom line
        countatoms++;
      if (fline == "loop_")  // may be the end of the atoms section
        break;
    }
    geometry->natoms = countatoms;
    printf("Number of atoms is %d\n", geometry->natoms);
    if (geometry->natoms < 1) {
      printf("Invalid number of atoms...exiting\n");
      exit(-1);
    }
    allocate_double_array(&(geometry->atoms_xyz), 3 * (geometry->natoms),
                          "can not allocate natoms xyz\n");
    allocate_double_array(&(geometry->atoms_fract), 3 * (geometry->natoms),
                          "can not allocate natoms fract\n");
    allocate_double_array(&(geometry->atoms_charge), geometry->natoms,
                          "can not allocate atom_charge\n");
    allocate_double_array(&(geometry->input_atoms_charge), geometry->natoms,
                          "can not allocate input_atom_charge\n");
    allocate_int_array(&(geometry->atoms_num), geometry->natoms,
                       "can not allocate atoms_num\n");
    allocate_int_array(&(geometry->atoms_base_num), geometry->natoms,
                       "can not allocate atoms_base_num\n");
    printf("Arrays for atoms allocated\n");
    for (i = 0; i < geometry->natoms; i++) {
      geometry->atoms_charge[i] = 0.0;
      geometry->input_atoms_charge[i] = 0.0;
    }

    printf("Rewind file %s to read the atoms after line: %s\n", fname, tmpstr1);
    rewind(fp);
    bool found_atoms = false;
    while (fgets(fline, 100, fp) != NULL) {
      // look for occurence of '_atomic_sites' in fline
      if (strstr(fline, tmpstr1) != NULL) {
        found_atoms = true;
        break;
      }
    }
    if (!found_atoms) {
      printf("Unable to find _atomic_sites tag\n");
      exit(-1);
    }

    printf("Reading atoms\n");
    for (i = 0; i < geometry->natoms; i++) {
      fgets(fline, 100, fp);
      tmpstr0 = strtok(fline, " ");
      tmpstr0 = strtok(NULL, " ");
      // strcpy(atoms_sym,tmpstr0);
      // gcc on MacOS doesn't like copying char * into char[2]
      strncpy(atoms_sym, tmpstr0, sizeof(atoms_sym));
      atoms_sym[2] = '\0';  // need to add null character manually
      // lookup the corrisponding atom number
      for (a = 0; a < 107; a++) {
        if (strcmp(atoms_sym, atsym_nospace[a]) == 0)
          geometry->atoms_num[i] = a + 1;
      }
      tmpstr0 = strtok(NULL, " ");
      geometry->atoms_fract[3 * i + 0] = atof(tmpstr0);
      tmpstr0 = strtok(NULL, " ");
      geometry->atoms_fract[3 * i + 1] = atof(tmpstr0);
      tmpstr0 = strtok(NULL, " ");
      geometry->atoms_fract[3 * i + 2] = atof(tmpstr0);
      geometry->input_atoms_charge[i] = 0.0;  // don't read it
      // convert in xyz coordinates
      geometry->atoms_xyz[3 * i + 0] =
          geometry->atoms_fract[3 * i + 0] * geometry->av1[0] +
          geometry->atoms_fract[3 * i + 1] * geometry->av1[1] +
          geometry->atoms_fract[3 * i + 2] * geometry->av1[2];
      geometry->atoms_xyz[3 * i + 1] =
          geometry->atoms_fract[3 * i + 0] * geometry->av2[0] +
          geometry->atoms_fract[3 * i + 1] * geometry->av2[1] +
          geometry->atoms_fract[3 * i + 2] * geometry->av2[2];
      geometry->atoms_xyz[3 * i + 2] =
          geometry->atoms_fract[3 * i + 0] * geometry->av3[0] +
          geometry->atoms_fract[3 * i + 1] * geometry->av3[1] +
          geometry->atoms_fract[3 * i + 2] * geometry->av3[2];
      printf("%s %5d %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n",
             atoms_sym, geometry->atoms_num[i], geometry->atoms_xyz[3 * i + 0],
             geometry->atoms_xyz[3 * i + 1], geometry->atoms_xyz[3 * i + 2],
             geometry->atoms_fract[3 * i + 0], geometry->atoms_fract[3 * i + 1],
             geometry->atoms_fract[3 * i + 2], geometry->input_atoms_charge[i]);
    }
    fclose(fp);
  }
  return;
}

void MapAtomsIntoCentralCell(tGeometry *geometry) {
  int iatom;
  int la, lb, lg;
  double xx, yy, zz;
  double alpha, beta, gamma;
  double malpha, mbeta, mgamma;
  double transfo[9];
  double itransfo[9];
  /* correct !!! transfo is transposed matrix of
 lattice vectors */
  transfo[0] = geometry->av1[0];
  transfo[1] = geometry->av2[0];
  transfo[2] = geometry->av3[0];

  transfo[3] = geometry->av1[1];
  transfo[4] = geometry->av2[1];
  transfo[5] = geometry->av3[1];
  transfo[6] = geometry->av1[2];
  transfo[7] = geometry->av2[2];
  transfo[8] = geometry->av3[2];
  /* get inverse transfo */
  matrix_inverse(transfo, 3, itransfo);
  for (iatom = 0; iatom < geometry->natoms; iatom++) {
    xx = geometry->atoms_xyz[3 * iatom + 0];
    yy = geometry->atoms_xyz[3 * iatom + 1];
    zz = geometry->atoms_xyz[3 * iatom + 2];
    printf("inp atom %3d: %15.12f %15.12f %15.12f\n", iatom + 1, xx, yy, zz);
    alpha = itransfo[0] * xx + itransfo[1] * yy + itransfo[2] * zz;
    beta = itransfo[3] * xx + itransfo[4] * yy + itransfo[5] * zz;
    gamma = itransfo[6] * xx + itransfo[7] * yy + itransfo[8] * zz;
    la = (alpha < 0.0) || (alpha >= 1.0);
    lb = (beta < 0.0) || (beta >= 1.0);
    lg = (gamma < 0.0) || (gamma >= 1.0);
    if (la || lb || lg) {
      printf("atom %3d is not mapped into the central cell\n", iatom + 1);
      printf("atom %3d inp natural coordinates %12.8f %12.8f %12.8f\n",
             iatom + 1, alpha, beta, gamma);
    }
    malpha = fmod(alpha, 1.0);
    if (malpha < 0.0) malpha = malpha + 1.0;
    mbeta = fmod(beta, 1.0);
    if (mbeta < 0.0) mbeta = mbeta + 1.0;
    mgamma = fmod(gamma, 1.0);
    if (mgamma < 0.0) mgamma = mgamma + 1.0;
    if (la || lb || lg)
      printf("atom %3d new natural coordinates %12.8f %12.8f %12.8f\n",
             iatom + 1, malpha, mbeta, mgamma);

    xx = geometry->atoms_xyz[3 * iatom + 0] =
        transfo[0] * malpha + transfo[1] * mbeta + transfo[2] * mgamma;
    yy = geometry->atoms_xyz[3 * iatom + 1] =
        transfo[3] * malpha + transfo[4] * mbeta + transfo[5] * mgamma;
    zz = geometry->atoms_xyz[3 * iatom + 2] =
        transfo[6] * malpha + transfo[7] * mbeta + transfo[8] * mgamma;
    printf("new atom %3d: %15.12f %15.12f %15.12f\n", iatom + 1, xx, yy, zz);
  }
  return;
}

void MapAtomsCheck(tGeometry *geometry) {
  int iatom;
  int la, lb, lg;
  int notmapped = 0;
  double xx, yy, zz;
  double alpha, beta, gamma;
  double malpha, mbeta, mgamma;
  double transfo[9];
  double itransfo[9];
  /* correct !!! transfo is transposed matrix of
 lattice vectors */

  transfo[0] = geometry->av1[0];
  transfo[1] = geometry->av2[0];
  transfo[2] = geometry->av3[0];
  transfo[3] = geometry->av1[1];
  transfo[4] = geometry->av2[1];
  transfo[5] = geometry->av3[1];
  transfo[6] = geometry->av1[2];
  transfo[7] = geometry->av2[2];
  transfo[8] = geometry->av3[2];
  /* get inverse transfo */
  matrix_inverse(transfo, 3, itransfo);
  for (iatom = 0; iatom < geometry->natoms; iatom++) {
    xx = geometry->atoms_xyz[3 * iatom + 0];
    yy = geometry->atoms_xyz[3 * iatom + 1];
    zz = geometry->atoms_xyz[3 * iatom + 2];
    alpha = itransfo[0] * xx + itransfo[1] * yy + itransfo[2] * zz;
    beta = itransfo[3] * xx + itransfo[4] * yy + itransfo[5] * zz;
    gamma = itransfo[6] * xx + itransfo[7] * yy + itransfo[8] * zz;
    la = (alpha < 0.0) || (alpha >= 1.0);
    lb = (beta < 0.0) || (beta >= 1.0);
    lg = (gamma < 0.0) || (gamma >= 1.0);
    if (la || lb || lg) {
      printf("atom %3d inp Cartesian coordinates %15.12f %15.12f %15.12f\n",
             iatom + 1, xx, yy, zz);
      printf("atom %3d inp natural coordinates %12.8f %12.8f %12.8f\n",
             iatom + 1, alpha, beta, gamma);
      printf("atom %3d is not mapped into the central cell\n", iatom + 1);
      if (la) {
        printf("natural coordinate 1: %15.12f\n", alpha);
      }
      if (lb) {
        printf("natural coordinate 2: %15.12f\n", beta);
      }
      if (lg) {
        printf("natural coordinate 3: %15.12f\n", gamma);
      }
      notmapped = notmapped + 1;
    }
  }
  if (notmapped > 0) {
    printf("%d atoms were not mapped into input cell...exiting...\n",
           notmapped);
    exit(-1);
  }
  return;
}

void DestroyGeometry(tGeometry *geometry) {
  free(geometry->atoms_xyz);
  free(geometry->atoms_charge);
  free(geometry->input_atoms_charge);
  free(geometry->atoms_num);
  free(geometry->atoms_base_num);
  return;
}

void DestroyBonds(tGeometry *geometry) {
  free(geometry->nbonds_indexij);
  free(geometry->tbonds_indexij);
  free(geometry->nbonds_type);
  return;
}

/* determines number of bonds and bond indices
   the bond is considered covalent if
   R_ij < 1.2*(R_i+R_j)
   R_i and R_j are obtained from setvdw/setcov */

/* At present, InitBonds atom-type atoms based on qeqbasetype */

void InitBonds(tGeometry *geometry) {
  int iatom, jatom, ibondt, ibond;
  int nbonds;
  int niatom, njatom;
  int icell, jcell, kcell;
  int elements[103];
  int num_distinct_el, btmax;
  int l_bond_type_found, lii, ljj;
  int imin, imax;
  double vdw[103];
  double rijmin, rij, rijvdw;
  double xatomi, yatomi, zatomi;
  double xatomj, yatomj, zatomj;
  int iatom_index, jatom_index;
  for (iatom = 0; iatom < 103; iatom++) elements[iatom] = 0;
  for (iatom = 0; iatom < geometry->natoms; iatom++) {
    niatom = geometry->atoms_base_num[iatom];
    if (elements[niatom - 1] == 0) elements[niatom - 1] = 1;
  }
  num_distinct_el = 0;
  for (iatom = 0; iatom < 103; iatom++)
    num_distinct_el = num_distinct_el + elements[iatom];
  geometry->natom_types = num_distinct_el;
  printf("InitBonds  atom type assignment is based on BASE TYPES\n");
  printf("number of atom types found=%d\n", geometry->natom_types);
  /* number of possible bonds: n*(n-1)/2+n */
  btmax = geometry->natom_types * (geometry->natom_types + 1) / 2;
  allocate_int_array(&(geometry->tbonds_indexij), 2 * btmax,
                     "can not allocate tbonds_indexij\n");
  for (ibondt = 0; ibondt < btmax; ibondt++)
    geometry->tbonds_indexij[2 * ibondt + 0] =
        geometry->tbonds_indexij[2 * ibondt + 1] = 0;
  geometry->tbonds = 0;

  /* SetVDWRADIUS(vdw);   */
  /* Covalent radius is in A? */
  SetCOVRADIUS(vdw);
  for (iatom = 0; iatom < 103; iatom++)
    vdw[iatom] = 1.2 * vdw[iatom] / 0.5291772;

  nbonds = 0;
  for (iatom = 0; iatom < geometry->natoms; iatom++) {
    xatomi = geometry->atoms_xyz[3 * iatom + 0];
    yatomi = geometry->atoms_xyz[3 * iatom + 1];
    zatomi = geometry->atoms_xyz[3 * iatom + 2];
    niatom = geometry->atoms_base_num[iatom];
    for (jatom = iatom + 1; jatom < geometry->natoms; jatom++) {
      njatom = geometry->atoms_base_num[jatom];
      rijmin = LARGED;
      for (icell = -1; icell <= 1; icell++) {
        for (jcell = -1; jcell <= 1; jcell++) {
          for (kcell = -1; kcell <= 1; kcell++) {
            xatomj = (geometry->atoms_xyz)[jatom * 3 + 0] +
                     icell * (geometry->av1[0]) + jcell * (geometry->av2[0]) +
                     kcell * (geometry->av3[0]);
            yatomj = (geometry->atoms_xyz)[jatom * 3 + 1] +
                     icell * (geometry->av1[1]) + jcell * (geometry->av2[1]) +
                     kcell * (geometry->av3[1]);
            zatomj = (geometry->atoms_xyz)[jatom * 3 + 2] +
                     icell * (geometry->av1[2]) + jcell * (geometry->av2[2]) +
                     kcell * (geometry->av3[2]);
            rij = sqrt((xatomi - xatomj) * (xatomi - xatomj) +
                       (yatomi - yatomj) * (yatomi - yatomj) +
                       (zatomi - zatomj) * (zatomi - zatomj));
            if (rij < 1.0e-8) {
              printf(
                  "atomj=%d maps into atomi=%d under primitive translation %d "
                  "%d %d\n",
                  iatom + 1, jatom + 1, icell, jcell, kcell);
              exit(-1);
            }
            if (rij < rijmin) rijmin = rij;
            rijvdw = (vdw[niatom - 1] + vdw[njatom - 1]) * (AUTOANG);
            if (rij < rijvdw) {
              nbonds = nbonds + 1;
              l_bond_type_found = 0;
              for (ibondt = 0; ibondt < geometry->tbonds; ibondt++) {
                lii = ((geometry->tbonds_indexij[2 * ibondt + 0] == niatom) &&
                       (geometry->tbonds_indexij[2 * ibondt + 1] == njatom));
                ljj = ((geometry->tbonds_indexij[2 * ibondt + 0] == njatom) &&
                       (geometry->tbonds_indexij[2 * ibondt + 1] == niatom));
                if (lii || ljj) {
                  l_bond_type_found = 1;
                  break;
                }
              }
              if (l_bond_type_found == 0) {
                imin = (niatom > njatom) ? njatom : niatom;
                imax = (niatom > njatom) ? niatom : njatom;
                geometry->tbonds_indexij[2 * geometry->tbonds + 0] = imin;
                geometry->tbonds_indexij[2 * geometry->tbonds + 1] = imax;
                geometry->tbonds = geometry->tbonds + 1;
              }
            }
          }
        }
      }
      if (fabs(rijmin - LARGED) < 1.0e-8) {
        printf("atomi=%d is too far from all other atoms\n", iatom + 1);
        exit(-1);
      }
    }
  }

  geometry->nbonds = nbonds;
  printf("%d bonds were determined from the input geometry\n",
         geometry->nbonds);
  printf("%d bond types were determined from the input geometry\n",
         geometry->tbonds);
  for (ibondt = 0; ibondt < geometry->tbonds; ibondt++)
    printf("bond type number %4d: atom types (%4d,%4d)\n", ibondt + 1,
           geometry->tbonds_indexij[2 * ibondt + 0],
           geometry->tbonds_indexij[2 * ibondt + 1]);

  if (geometry->nbonds < 1) {
    printf("invalid number of bonds %d\n", geometry->nbonds);
    exit(-1);
  }

  allocate_int_array(&(geometry->nbonds_indexij), 2 * geometry->nbonds,
                     "can not allocate nbonds_indexij\n");
  allocate_int_array(&(geometry->nbonds_type), geometry->nbonds,
                     "can not allocate nbonds_type\n");
  for (ibond = 0; ibond < geometry->nbonds; ibond++) {
    geometry->nbonds_indexij[2 * ibond + 0] = 0;
    geometry->nbonds_indexij[2 * ibond + 1] = 0;
    geometry->nbonds_type[ibond] = 0;
  }

  nbonds = 0;
  for (iatom = 0; iatom < geometry->natoms; iatom++) {
    xatomi = geometry->atoms_xyz[3 * iatom + 0];
    yatomi = geometry->atoms_xyz[3 * iatom + 1];
    zatomi = geometry->atoms_xyz[3 * iatom + 2];
    niatom = geometry->atoms_base_num[iatom];
    for (jatom = iatom + 1; jatom < geometry->natoms; jatom++) {
      njatom = geometry->atoms_base_num[jatom];
      rijmin = LARGED;
      for (icell = -1; icell <= 1; icell++) {
        for (jcell = -1; jcell <= 1; jcell++) {
          for (kcell = -1; kcell <= 1; kcell++) {
            xatomj = (geometry->atoms_xyz)[jatom * 3 + 0] +
                     icell * (geometry->av1[0]) + jcell * (geometry->av2[0]) +
                     kcell * (geometry->av3[0]);
            yatomj = (geometry->atoms_xyz)[jatom * 3 + 1] +
                     icell * (geometry->av1[1]) + jcell * (geometry->av2[1]) +
                     kcell * (geometry->av3[1]);
            zatomj = (geometry->atoms_xyz)[jatom * 3 + 2] +
                     icell * (geometry->av1[2]) + jcell * (geometry->av2[2]) +
                     kcell * (geometry->av3[2]);
            rij = sqrt((xatomi - xatomj) * (xatomi - xatomj) +
                       (yatomi - yatomj) * (yatomi - yatomj) +
                       (zatomi - zatomj) * (zatomi - zatomj));
            if (rij < rijmin) rijmin = rij;
            rijvdw = (vdw[niatom - 1] + vdw[njatom - 1]) * (AUTOANG);
            if (rij < rijvdw) {
              (geometry->nbonds_indexij)[2 * nbonds + 0] = iatom + 1;
              (geometry->nbonds_indexij)[2 * nbonds + 1] = jatom + 1;
              l_bond_type_found = 0;
              for (ibondt = 0; ibondt < geometry->tbonds; ibondt++) {
                lii = ((geometry->tbonds_indexij[2 * ibondt + 0] == niatom) &&
                       (geometry->tbonds_indexij[2 * ibondt + 1] == njatom));
                ljj = ((geometry->tbonds_indexij[2 * ibondt + 0] == njatom) &&
                       (geometry->tbonds_indexij[2 * ibondt + 1] == niatom));
                if (lii || ljj) {
                  l_bond_type_found = 1;
                  geometry->nbonds_type[nbonds] = ibondt;
                  break;
                }
              }
              if (l_bond_type_found == 0) {
                printf("bond type for atoms (%d,%d) is not found\n", iatom,
                       jatom);
                exit(-1);
              }
              nbonds = nbonds + 1;
            }
          }
        }
      }
    }
  }
  if (nbonds != geometry->nbonds) {
    printf("change in the number of bonds from old=%d to new=%d\n",
           geometry->nbonds, nbonds);
    exit(-1);
  }
  PrintBondInfo(geometry);
  return;
}

void PrintBondInfo(tGeometry *geometry) {
  int ibond;
  printf("printing bond info...\n");
  for (ibond = 0; ibond < geometry->nbonds; ibond++)
    printf("bond %10d: %4d %4d\n", ibond + 1,
           geometry->nbonds_indexij[2 * ibond + 0],
           geometry->nbonds_indexij[2 * ibond + 1]);
  return;
}

double matrix_det(double *mat, int n) {
  double res = 0;
  double a11, a12, a13;
  double a21, a22, a23;
  double a31, a32, a33;
  if (n == 2) {
    a11 = mat[2 * 0 + 0];
    a12 = mat[2 * 0 + 1];
    a21 = mat[2 * 1 + 0];
    a22 = mat[2 * 1 + 1];
    res = a11 * a22 - a12 * a21;
  } else if (n == 3) {
    a11 = mat[3 * 0 + 0];
    a12 = mat[3 * 0 + 1];
    a13 = mat[3 * 0 + 2];
    a21 = mat[3 * 1 + 0];
    a22 = mat[3 * 1 + 1];
    a23 = mat[3 * 1 + 2];
    a31 = mat[3 * 2 + 0];
    a32 = mat[3 * 2 + 1];
    a33 = mat[3 * 2 + 2];
    res = a11 * (a33 * a22 - a32 * a23) - a21 * (a33 * a12 - a32 * a13) +
          a31 * (a23 * a12 - a22 * a13);
  } else {
    printf("invalid matrix dimension in matrixdet\n");
    exit(-1);
  }
  return (res);
}

void matrix_inverse(double *mat, int n, double *imat) {
  double det;
  double a11, a12, a13;
  double a21, a22, a23;
  double a31, a32, a33;
  det = matrix_det(mat, n);
  if (fabs(det) < 1.0e-6) {
    printf("matrix determinant %10.6f is too small\n", det);
    exit(-1);
  }
  if (n == 2) {
    a11 = mat[2 * 0 + 0];
    a12 = mat[2 * 0 + 1];
    a21 = mat[2 * 1 + 0];
    a22 = mat[2 * 1 + 1];
    imat[2 * 0 + 0] = a22 / det;
    imat[2 * 0 + 1] = -a12 / det;
    imat[2 * 1 + 0] = -a21 / det;
    imat[2 * 1 + 1] = a11 / det;
  } else if (n == 3) {
    a11 = mat[3 * 0 + 0];
    a12 = mat[3 * 0 + 1];
    a13 = mat[3 * 0 + 2];
    a21 = mat[3 * 1 + 0];
    a22 = mat[3 * 1 + 1];
    a23 = mat[3 * 1 + 2];
    a31 = mat[3 * 2 + 0];
    a32 = mat[3 * 2 + 1];
    a33 = mat[3 * 2 + 2];

    imat[3 * 0 + 0] = (a33 * a22 - a32 * a23) / det;
    imat[3 * 0 + 1] = -(a33 * a12 - a32 * a13) / det;
    imat[3 * 0 + 2] = (a23 * a12 - a22 * a13) / det;
    imat[3 * 1 + 0] = -(a33 * a21 - a31 * a23) / det;
    imat[3 * 1 + 1] = (a33 * a11 - a31 * a13) / det;
    imat[3 * 1 + 2] = -(a23 * a11 - a21 * a13) / det;
    imat[3 * 2 + 0] = (a32 * a21 - a31 * a22) / det;
    imat[3 * 2 + 1] = -(a32 * a11 - a31 * a12) / det;
    imat[3 * 2 + 2] = (a22 * a11 - a21 * a12) / det;
  } else {
    printf("invalid matrix dimension in matrixdet\n");
    exit(-1);
  }

  return;
}

#undef MAXSTRL
#undef LARGED
#undef AUTOANG
