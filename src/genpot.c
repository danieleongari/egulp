#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cell.h"
#include "elements.h"
#include "geometry.h"
#include "interact.h"
#include "utilsa.h"

int prepareKSUM(double GMax, double eta, int NUMR, double *b1, double *b2,
                double *b3, double *bgridx, double *bgridy, double *bgridz,
                double *bgridcoef) {
  int GVEC = 0;
  double bvec[3] = {0.0, 0.0, 0.0};
  double bvecsz;
  double coef;
  int k1, k2, k3;
  /* for (k1=-NUMR; k1 <= NUMR; k1++) { */
  for (k1 = -NUMR; k1 <= 0; k1++) {
    for (k2 = -NUMR; k2 <= NUMR; k2++) {
      for (k3 = -NUMR; k3 <= NUMR; k3++) {
        if ((k1 == 0) && (k2 == 0) && (k3 == 0)) continue;
        bvec[0] = k1 * b1[0] + k2 * b2[0] + k3 * b3[0];
        bvec[1] = k1 * b1[1] + k2 * b2[1] + k3 * b3[1];
        bvec[2] = k1 * b1[2] + k2 * b2[2] + k3 * b3[2];
        bvecsz =
            sqrt(bvec[0] * bvec[0] + bvec[1] * bvec[1] + bvec[2] * bvec[2]);
        if (bvecsz < GMax) {
          GVEC = GVEC + 1;
          bgridx[GVEC - 1] = bvec[0];
          bgridy[GVEC - 1] = bvec[1];
          bgridz[GVEC - 1] = bvec[2];
          coef = exp((-0.25 / eta) * (bvecsz * bvecsz)) / (bvecsz * bvecsz);
          if (k1 < 0) coef = 2.0 * coef;
          bgridcoef[GVEC - 1] = coef;
        }
      }
    }
  }
  printf("preparing K sum: number of K-vectors in k sum is %d\n", GVEC);
  return (GVEC);
}

int principle_qnum(int anum) {
  int npqni = 0;
  if ((anum < 1) || (anum > 103)) {
    printf("invalid anum=%d in principle_qnum\n", anum);
    exit(-1);
  }

  if ((anum >= 1) && (anum <= 2)) {
    npqni = 1;
  } else if ((anum >= 3) && (anum <= 10)) {
    npqni = 2;
  } else if ((anum >= 11) && (anum <= 18)) {
    npqni = 3;
  } else if ((anum >= 19) && (anum <= 36)) {
    npqni = 4;
  } else if ((anum >= 37) && (anum <= 54)) {
    npqni = 5;
  } else if ((anum >= 55) && (anum <= 86)) {
    npqni = 6;
  } else
    npqni = 7;

  return (npqni);
}

double erf(double x) {
  int sign;
  const double a1 = 0.254829592;
  const double a2 = -0.284496736;
  const double a3 = 1.421413741;
  const double a4 = -1.453152027;
  const double a5 = 1.061405429;
  const double p = 0.3275911;
  double absx;
  double t;
  double y;
  double res;
  sign = 1;
  if (x < 0) {
    sign = -1;
  }
  absx = fabs(x);
  /* A&S formula 7.1.26 */
  t = 1.0 / (1.0 + p * absx);
  y = 1.0 -
      (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-absx * absx);
  res = ((double)(sign)) * y;
  return (res); /*erf(-x) = -erf(x) */
}

/* qmatrixelement(dx,dy,dz,av1,av2,av3, dvREAL, eta,  radmax, tweatpi,
    bgridx, bgridy, bgridz, bgridcoef, kvectors, &vijre, &vijself, &vijR); */

void qmatrixelement(double dx, double dy, double dz, double *av1, double *av2,
                    double *av3, double V, double eta, double radmax,
                    double tweatpi, double *bgridx, double *bgridy,
                    double *bgridz, double *bgridcoef, int kvectors,
                    double *vijre, double *vijself, double *vijR) {
  double sumre;
  int k;
  int lzdisp;
  int numr;
  int t1, t2, t3;
  double av1sz, av2sz, av3sz;
  double avecszmin;
  double tvecsz;
  double drx, dry, drz, drr;
  double tvec[3] = {0.0, 0.0, 0.0};

  sumre = 0.0;
  for (k = 0; k < kvectors; k++)
    sumre = sumre + cos(dx * bgridx[k] + dy * bgridy[k] + dz * bgridz[k]) *
                        bgridcoef[k];
  (*vijre) = (4.0 * M_PI / V) * sumre;

  lzdisp = 0;
  (*vijself) = 0.0;
  if (sqrt(dx * dx + dy * dy + dz * dz) < 1.0e-8) {
    (*vijself) = -tweatpi;
    lzdisp = 1;
  }

  av1sz = sqrt(av1[0] * av1[0] + av1[1] * av1[1] + av1[2] * av1[2]);
  av2sz = sqrt(av2[0] * av2[0] + av2[1] * av2[1] + av2[2] * av2[2]);
  av3sz = sqrt(av3[0] * av2[0] + av3[1] * av3[1] + av3[2] * av3[2]);
  avecszmin = av1sz;
  if (av2sz < avecszmin) avecszmin = av2sz;
  if (av3sz < avecszmin) avecszmin = av3sz;

  numr = ((int)(radmax / avecszmin)) + 2;

  sumre = 0.0;
  for (t1 = -numr; t1 <= numr; t1++) {
    for (t2 = -numr; t2 <= numr; t2++) {
      for (t3 = -numr; t3 <= numr; t3++) {
        if ((lzdisp == 1) && (t1 == 0) && (t2 == 0) && (t3 == 0)) continue;
        tvec[0] = t1 * av1[0] + t2 * av2[0] + t3 * av3[0];
        tvec[1] = t1 * av1[1] + t2 * av2[1] + t3 * av3[1];
        tvec[2] = t1 * av1[2] + t2 * av2[2] + t3 * av3[2];
        tvecsz =
            sqrt(tvec[0] * tvec[0] + tvec[1] * tvec[1] + tvec[2] * tvec[2]);
        drx = dx + tvec[0];
        dry = dy + tvec[1];
        drz = dz + tvec[2];
        drr = sqrt(drx * drx + dry * dry + drz * drz);
        if (drr < radmax)
          sumre = sumre + (1.0 / drr) * (1.0 - erf(sqrt(eta) * drr));
      }
    }
  }
  (*vijR) = sumre;
  return;
}

void genpot(double *av1, double *av2, double *av3, int natoms, int *atoms_num,
            double *atoms_charge, double *atoms_xyz, double *J, double *qeqrad,
            int *qeqids, int *qeqbasetype, int MaxEl, int *CalcMatrix) {
  double dvREAL, dvRECIP;
  double bv1[3];
  double bv2[3];
  double bv3[3];
  const double autoangs = 0.529177;
  const double angstoev = 14.3997584;
  const double accuracy = 8.0;
  const double rqeq = 15.0;
  const double cuts = 0.6;
  const double qeqlambda = 0.5;
  double rradmax;
  double radmax;
  double tweatpi;
  double bvecsz[3];
  double bvecszmin;
  int NUMR;
  int kvectors;
  double *bgridx;
  double *bgridy;
  double *bgridz;
  double *bgridcoef;
  int iatom, jatom, i;
  double dx, dy, dz;
  double vijre, vijself, vijR;
  int maxloop2[3] = {0, 0, 0};
  double avec1sz, avec2sz, avec3sz;
  double tvec[3] = {0.0, 0.0, 0.0};
  int iatom_ni, jatom_nj;
  int t1, t2, t3;
  double zeta_i, zeta_j;
  double dxx, dyy, dzz, dr;
  double gamma;
  double eta, accf;
  int count, iatom_index, jatom_index;

  dvREAL = cellVolume(av1, av2, av3);
  vectorProduct(av2, av3, bv1);
  vectorProduct(av3, av1, bv2);
  vectorProduct(av1, av2, bv3);
  for (i = 0; i <= 2; i++) {
    bv1[i] = bv1[i] * (2.0 * M_PI / dvREAL);
    bv2[i] = bv2[i] * (2.0 * M_PI / dvREAL);
    bv3[i] = bv3[i] * (2.0 * M_PI / dvREAL);
  }
  dvRECIP = cellVolume(bv1, bv2, bv3);
  printf("\n");
  printf("*** Entering genpot utility***\n");
  printf("Reciprocal lattice vector #1: %10.6f %10.6f %10.6f\n", bv1[0], bv1[1],
         bv1[2]);
  printf("Reciprocal lattice vector #2: %10.6f %10.6f %10.6f\n", bv2[0], bv2[1],
         bv2[2]);
  printf("Reciprocal lattice vector #3: %10.6f %10.6f %10.6f\n", bv3[0], bv3[1],
         bv3[2]);
  printf("Real unit cell volume: %10.6f\n", dvREAL);
  printf("Reciprocal unit cell volume: %10.6f\n", dvRECIP);

  accf = sqrt(log(pow(10.0, accuracy)));
  eta = pow(sqrt(natoms) * pow(M_PI, 3.0) / (dvREAL * dvREAL), 1.0 / 3.0);

  printf("QEq corrections: real cut-off rqeq=%10.6f\n", rqeq);
  printf("QEq corrections: self energy cut-off=%10.6f\n", cuts);
  printf("accf = %12.8f\n", accf);
  printf("eta = %12.8f\n", eta);

  rradmax = 2.0 * accf * sqrt(eta);
  printf("1 rradmax  = %12.8f\n", rradmax);
  if ((rradmax * rradmax) > 140.0 * eta) rradmax = sqrt(140.0 * eta);
  printf("2 rradmax  = %12.8f\n", rradmax);

  radmax = accf / sqrt(eta);
  printf("1 radmax = %12.8f\n", radmax);
  if ((radmax * radmax) > 85.0 / eta) radmax = sqrt(85.0 / eta);
  printf("2 radmax = %12.8f\n", radmax);
  tweatpi = 2.0 * sqrt(eta / M_PI);
  bvecsz[0] = sqrt(bv1[0] * bv1[0] + bv1[1] * bv1[1] + bv1[2] * bv1[2]);
  bvecsz[1] = sqrt(bv2[0] * bv2[0] + bv2[1] * bv2[1] + bv2[2] * bv2[2]);
  bvecsz[2] = sqrt(bv3[0] * bv3[0] + bv3[1] * bv3[1] + bv3[2] * bv3[2]);
  bvecszmin = bvecsz[0];
  if (bvecsz[1] < bvecszmin) bvecszmin = bvecsz[1];
  if (bvecsz[2] < bvecszmin) bvecszmin = bvecsz[2];
  NUMR = ((int)(rradmax / bvecszmin)) + 2;
  allocate_double_array(&(bgridx),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridx\n");
  allocate_double_array(&(bgridy),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridy\n");
  allocate_double_array(&(bgridz),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridz\n");
  allocate_double_array(&(bgridcoef),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridcoef\n");
  kvectors = prepareKSUM(rradmax, eta, NUMR, bv1, bv2, bv3, bgridx, bgridy,
                         bgridz, bgridcoef);

  for (iatom = 0; iatom < natoms; iatom++) {
    for (jatom = 0; jatom <= iatom; jatom++) {
      if (CalcMatrix[iatom * natoms + jatom] == 0) continue;
      dx = atoms_xyz[3 * jatom + 0] - atoms_xyz[3 * iatom + 0];
      dy = atoms_xyz[3 * jatom + 1] - atoms_xyz[3 * iatom + 1];
      dz = atoms_xyz[3 * jatom + 2] - atoms_xyz[3 * iatom + 2];
      qmatrixelement(dx, dy, dz, av1, av2, av3, dvREAL, eta, radmax, tweatpi,
                     bgridx, bgridy, bgridz, bgridcoef, kvectors, &vijre,
                     &vijself, &vijR);
      J[iatom * natoms + jatom] = (vijre + vijself + vijR);
      /* printf("iatom=%d jatom=%d  vijre=%10.6f vijself=%10.6f vijR=%10.6f\n",
        iatom, jatom,  vijre*angstoev, vijself*angstoev, vijR*angstoev); */
      if (!(iatom == jatom))
        J[jatom * natoms + iatom] = J[iatom * natoms + jatom];
    }
  }
  /* printf("interaction matrix without diagonal after Ewald\n");
   count = 0;
   for (iatom = 0; iatom < natoms; iatom++)
              for (jatom  = 0; jatom < natoms; jatom++) {
                      printf("%10.6f ",  J[iatom*natoms+jatom]*angstoev );
                      count = count + 1;
                      if(count %10 == 0) printf("\n");
               }  */
  /*#########################################################################
  # now J is computed using Ewald sums...
  # Need to modify for orbital interactions...
  #########################################################################*/
  avec1sz = sqrt(av1[0] * av1[0] + av1[1] * av1[1] + av1[2] * av1[2]);
  avec2sz = sqrt(av2[0] * av2[0] + av2[1] * av2[1] + av2[2] * av2[2]);
  avec3sz = sqrt(av3[0] * av3[0] + av3[1] * av3[1] + av3[2] * av3[2]);
  maxloop2[0] = ((int)(rqeq / avec1sz)) + 2;
  maxloop2[1] = ((int)(rqeq / avec2sz)) + 2;
  maxloop2[2] = ((int)(rqeq / avec3sz)) + 2;
  printf("orbital corrections loop %4d %4d %4d\n", maxloop2[0], maxloop2[1],
         maxloop2[2]);
  for (iatom = 0; iatom < natoms; iatom++) {
    iatom_index = element_index(atoms_num[iatom], qeqids, MaxEl);
    iatom_ni = principle_qnum(qeqbasetype[iatom_index]);
    /* iatom_ni = principle_qnum( atoms_num[iatom]); */
    /* zeta_i = 0.5*qeqlambda*(2.0*iatom_ni+1.0)/qeqrad[atoms_num[iatom]-1] ; */
    zeta_i = 0.5 * qeqlambda * (2.0 * iatom_ni + 1.0) / qeqrad[iatom_index];

    /* printf("iatom_ni=%d, zeta_i=%10.6f\n", iatom_ni, zeta_i); */
    if (atoms_num[iatom] == 1) zeta_i = zeta_i + atoms_charge[iatom] / autoangs;

    for (jatom = 0; jatom <= iatom; jatom++) {
      jatom_index = element_index(atoms_num[jatom], qeqids, MaxEl);
      jatom_nj = principle_qnum(qeqbasetype[jatom_index]);
      /* jatom_nj = principle_qnum( atoms_num[jatom])  ; */
      /* zeta_j = 0.5*qeqlambda*(2.0*jatom_nj+1.0)/qeqrad[atoms_num[jatom]-1];
       */
      zeta_j = 0.5 * qeqlambda * (2.0 * jatom_nj + 1.0) / qeqrad[jatom_index];
      /* printf("jatom_nj=%d, zeta_j=%10.6f\n", jatom_nj, zeta_j);  */
      if (atoms_num[jatom] == 1)
        zeta_j = zeta_j + atoms_charge[jatom] / autoangs;

      if (CalcMatrix[iatom * natoms + jatom] == 0) continue;

      dx = atoms_xyz[3 * jatom + 0] - atoms_xyz[3 * iatom + 0];
      dy = atoms_xyz[3 * jatom + 1] - atoms_xyz[3 * iatom + 1];
      dz = atoms_xyz[3 * jatom + 2] - atoms_xyz[3 * iatom + 2];
      /* printf("dx=%10.6f dy=%10.6f dz=%10.6f\n", dx, dy, dz);   */
      for (t1 = -maxloop2[0]; t1 <= maxloop2[0]; t1++) {
        for (t2 = -maxloop2[1]; t2 <= maxloop2[1]; t2++) {
          for (t3 = -maxloop2[2]; t3 <= maxloop2[2]; t3++) {
            tvec[0] = t1 * av1[0] + t2 * av2[0] + t3 * av3[0];
            tvec[1] = t1 * av1[1] + t2 * av2[1] + t3 * av3[1];
            tvec[2] = t1 * av1[2] + t2 * av2[2] + t3 * av3[2];
            dxx = dx + tvec[0];
            dyy = dy + tvec[1];
            dzz = dz + tvec[2];
            dr = sqrt(dxx * dxx + dyy * dyy + dzz * dzz);
            if ((dr > rqeq) || (dr < cuts)) continue;
            gamma = orb_interaction(iatom_ni, zeta_i, jatom_nj, zeta_j, dr);
            J[iatom * natoms + jatom] =
                J[iatom * natoms + jatom] + gamma - 1.0 / dr;
          }
        }
      }
      if (!(iatom == jatom))
        J[jatom * natoms + iatom] = J[iatom * natoms + jatom];
    }
  }

  for (iatom = 0; iatom < natoms; iatom++) {
    for (jatom = 0; jatom < natoms; jatom++) {
      if (CalcMatrix[iatom * natoms + jatom] == 0) continue;
      J[iatom * natoms + jatom] = J[iatom * natoms + jatom] * angstoev;
    }
  }
  /*      printf("interaction matrix without diagonal is\n");
         count = 0;
         for (iatom = 0; iatom < natoms; iatom++)
                  for (jatom  = 0; jatom < natoms; jatom++) {
                          printf("%10.6f ",  J[iatom*natoms+jatom] );
                          count = count + 1;
                          if(count %10 == 0) printf("\n");
                   }  */

  /*      printf("genpot is done...\n");  */
  free(bgridx);
  free(bgridy);
  free(bgridz);
  free(bgridcoef);
  return;
}

void qmatrixelement_grid(double dx, double dy, double dz, double *av1,
                         double *av2, double *av3, double avecszmin, double V,
                         double eta, double radmax, double *bgridx,
                         double *bgridy, double *bgridz, double *bgridcoef,
                         int kvectors, double *vijre, double *vijR) {
  double sumre;
  int k;
  int numr;
  int t1, t2, t3;
  double tvecsz;
  double drx, dry, drz, drr;
  double tvec[3] = {0.0, 0.0, 0.0};

  sumre = 0.0;
  for (k = 0; k < kvectors; k++)
    sumre = sumre + cos(dx * bgridx[k] + dy * bgridy[k] + dz * bgridz[k]) *
                        bgridcoef[k];
  (*vijre) = (4.0 * M_PI / V) * sumre;

  numr = ((int)(radmax / avecszmin)) + 2;

  sumre = 0.0;
  for (t1 = -numr; t1 <= numr; t1++) {
    for (t2 = -numr; t2 <= numr; t2++) {
      for (t3 = -numr; t3 <= numr; t3++) {
        tvec[0] = t1 * av1[0] + t2 * av2[0] + t3 * av3[0];
        tvec[1] = t1 * av1[1] + t2 * av2[1] + t3 * av3[1];
        tvec[2] = t1 * av1[2] + t2 * av2[2] + t3 * av3[2];
        tvecsz =
            sqrt(tvec[0] * tvec[0] + tvec[1] * tvec[1] + tvec[2] * tvec[2]);
        drx = dx + tvec[0];
        dry = dy + tvec[1];
        drz = dz + tvec[2];
        drr = sqrt(drx * drx + dry * dry + drz * drz);
        if (drr < radmax)
          sumre = sumre + (1.0 / drr) * (1.0 - erf(sqrt(eta) * drr));
      }
    }
  }
  (*vijR) = sumre;
  return;
}

void genpot_on_grid(double *av1, double *av2, double *av3, int natoms,
                    double *atoms_charge, double *atoms_xyz, int gridnl,
                    double *pot, int *L1, int *L2, int *L3, double *h1,
                    double *h2, double *h3) {
  double dvREAL, dvRECIP;
  double bv1[3];
  double bv2[3];
  double bv3[3];
  const double autoangs = 0.529177;
  const double angstoev = 14.3997584;
  const double accuracy = 8.0;
  const double rqeq = 15.0;
  const double cuts = 0.6;
  const double qeqlambda = 0.5;
  double rradmax;
  double radmax;
  double tweatpi;
  double bvecsz[3];
  double bvecszmin;
  int NUMR;
  int kvectors;
  double *bgridx;
  double *bgridy;
  double *bgridz;
  double *bgridcoef;
  int iatom, jatom, i;
  double dx, dy, dz;
  double vijre, vijR;
  int maxloop2[3] = {0, 0, 0};
  double avec1sz, avec2sz, avec3sz;
  double tvec[3] = {0.0, 0.0, 0.0};
  int iatom_ni, jatom_nj;
  int t1, t2, t3;
  double zeta_i, zeta_j;
  double dxx, dyy, dzz, dr;
  double gamma;
  double eta, accf;
  double xx, yy, zz;
  int count, ip;
  double av1sz, av2sz, av3sz;
  double avecszmin;

  dvREAL = cellVolume(av1, av2, av3);
  vectorProduct(av2, av3, bv1);
  vectorProduct(av3, av1, bv2);
  vectorProduct(av1, av2, bv3);
  for (i = 0; i <= 2; i++) {
    bv1[i] = bv1[i] * (2.0 * M_PI / dvREAL);
    bv2[i] = bv2[i] * (2.0 * M_PI / dvREAL);
    bv3[i] = bv3[i] * (2.0 * M_PI / dvREAL);
  }
  dvRECIP = cellVolume(bv1, bv2, bv3);
  printf("\n");
  printf("*** Entering genpot_on_grid utility***\n");
  printf("Reciprocal lattice vector #1: %10.6f %10.6f %10.6f\n", bv1[0], bv1[1],
         bv1[2]);
  printf("Reciprocal lattice vector #2: %10.6f %10.6f %10.6f\n", bv2[0], bv2[1],
         bv2[2]);
  printf("Reciprocal lattice vector #3: %10.6f %10.6f %10.6f\n", bv3[0], bv3[1],
         bv3[2]);
  printf("Real unit cell volume: %10.6f\n", dvREAL);
  printf("Reciprocal unit cell volume: %10.6f\n", dvRECIP);

  accf = sqrt(log(pow(10.0, accuracy)));
  eta = pow(sqrt(natoms) * pow(M_PI, 3.0) / (dvREAL * dvREAL), 1.0 / 3.0);

  printf("real cut-off rqeq=%10.6f\n", rqeq);
  printf("self energy cut-off=%10.6f\n", cuts);
  printf("accf = %12.8f\n", accf);
  printf("eta = %12.8f\n", eta);

  rradmax = 2.0 * accf * sqrt(eta);
  printf("1 rradmax  = %12.8f\n", rradmax);
  if ((rradmax * rradmax) > 140.0 * eta) rradmax = sqrt(140.0 * eta);
  printf("2 rradmax  = %12.8f\n", rradmax);

  radmax = accf / sqrt(eta);
  printf("1 radmax = %12.8f\n", radmax);
  if ((radmax * radmax) > 85.0 / eta) radmax = sqrt(85.0 / eta);
  printf("2 radmax = %12.8f\n", radmax);
  tweatpi = 2.0 * sqrt(eta / M_PI);
  bvecsz[0] = sqrt(bv1[0] * bv1[0] + bv1[1] * bv1[1] + bv1[2] * bv1[2]);
  bvecsz[1] = sqrt(bv2[0] * bv2[0] + bv2[1] * bv2[1] + bv2[2] * bv2[2]);
  bvecsz[2] = sqrt(bv3[0] * bv3[0] + bv3[1] * bv3[1] + bv3[2] * bv3[2]);
  bvecszmin = bvecsz[0];
  if (bvecsz[1] < bvecszmin) bvecszmin = bvecsz[1];
  if (bvecsz[2] < bvecszmin) bvecszmin = bvecsz[2];
  NUMR = ((int)(rradmax / bvecszmin)) + 2;
  allocate_double_array(&(bgridx),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridx\n");
  allocate_double_array(&(bgridy),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridy\n");
  allocate_double_array(&(bgridz),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridz\n");
  allocate_double_array(&(bgridcoef),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridcoef\n");
  kvectors = prepareKSUM(rradmax, eta, NUMR, bv1, bv2, bv3, bgridx, bgridy,
                         bgridz, bgridcoef);
  av1sz = sqrt(av1[0] * av1[0] + av1[1] * av1[1] + av1[2] * av1[2]);
  av2sz = sqrt(av2[0] * av2[0] + av2[1] * av2[1] + av2[2] * av2[2]);
  av3sz = sqrt(av3[0] * av2[0] + av3[1] * av3[1] + av3[2] * av3[2]);
  avecszmin = av1sz;
  if (av2sz < avecszmin) avecszmin = av2sz;
  if (av3sz < avecszmin) avecszmin = av3sz;

  for (ip = 0; ip < gridnl; ip++) {
    xx = L1[ip] * h1[0] + L2[ip] * h2[0] + L3[ip] * h3[0];
    yy = L1[ip] * h1[1] + L2[ip] * h2[1] + L3[ip] * h3[1];
    zz = L1[ip] * h1[2] + L2[ip] * h2[2] + L3[ip] * h3[2];
    pot[ip] = 0.0;
    for (iatom = 0; iatom < natoms; iatom++) {
      dx = xx - atoms_xyz[3 * iatom + 0];
      dy = yy - atoms_xyz[3 * iatom + 1];
      dz = zz - atoms_xyz[3 * iatom + 2];
      qmatrixelement_grid(dx, dy, dz, av1, av2, av3, avecszmin, dvREAL, eta,
                          radmax, bgridx, bgridy, bgridz, bgridcoef, kvectors,
                          &vijre, &vijR);
      pot[ip] = pot[ip] + (vijre + vijR) * atoms_charge[iatom];
    }
  }

  /* potential is returned in eV */
  for (ip = 0; ip < gridnl; ip++) pot[ip] = pot[ip] * angstoev;

  free(bgridx);
  free(bgridy);
  free(bgridz);
  free(bgridcoef);
  return;
}

/*
   av1 - lattice vector #1
   av2 - lattice vector #2
   av3 - lattice vector #3
   nsource - number of source charges
   source_charge - "nsource" source charges used to construct potential
   source_xyz - 3x"nsource" coordinates of source charges
   npoints_to_gen - number of points to generate potential in...
   points_to_gen_xyz - coordinates of points to generate potential in...
   natoms - number of atoms in the system which is used to determine
            parameters of the Ewald summation...
*/

void genpot_on_points(double *av1, double *av2, double *av3, int nsource,
                      double *source_charge, double *source_xyz,
                      int npoints_to_gen, double *pot,
                      double *points_to_gen_xyz, int natoms) {
  double dvREAL, dvRECIP;
  double bv1[3];
  double bv2[3];
  double bv3[3];
  const double autoangs = 0.529177;
  const double angstoev = 14.3997584;
  const double accuracy = 8.0;
  const double rqeq = 15.0;
  const double cuts = 0.6;
  const double qeqlambda = 0.5;
  double rradmax;
  double radmax;
  double tweatpi;
  double bvecsz[3];
  double bvecszmin;
  int NUMR;
  int kvectors;
  double *bgridx;
  double *bgridy;
  double *bgridz;
  double *bgridcoef;
  int iatom, jatom, i, isource;
  double dx, dy, dz;
  double vijre, vijR;
  int maxloop2[3] = {0, 0, 0};
  double avec1sz, avec2sz, avec3sz;
  double tvec[3] = {0.0, 0.0, 0.0};
  int iatom_ni, jatom_nj;
  int t1, t2, t3;
  double zeta_i, zeta_j;
  double dxx, dyy, dzz, dr;
  double gamma;
  double eta, accf;
  double xx, yy, zz;
  int count, ip;
  double av1sz, av2sz, av3sz;
  double avecszmin;

  dvREAL = cellVolume(av1, av2, av3);
  vectorProduct(av2, av3, bv1);
  vectorProduct(av3, av1, bv2);
  vectorProduct(av1, av2, bv3);
  for (i = 0; i <= 2; i++) {
    bv1[i] = bv1[i] * (2.0 * M_PI / dvREAL);
    bv2[i] = bv2[i] * (2.0 * M_PI / dvREAL);
    bv3[i] = bv3[i] * (2.0 * M_PI / dvREAL);
  }
  dvRECIP = cellVolume(bv1, bv2, bv3);
  printf("\n");
  printf("*** Entering genpot_on_grid utility***\n");
  printf("Reciprocal lattice vector #1: %10.6f %10.6f %10.6f\n", bv1[0], bv1[1],
         bv1[2]);
  printf("Reciprocal lattice vector #2: %10.6f %10.6f %10.6f\n", bv2[0], bv2[1],
         bv2[2]);
  printf("Reciprocal lattice vector #3: %10.6f %10.6f %10.6f\n", bv3[0], bv3[1],
         bv3[2]);
  printf("Real unit cell volume: %10.6f\n", dvREAL);
  printf("Reciprocal unit cell volume: %10.6f\n", dvRECIP);

  accf = sqrt(log(pow(10.0, accuracy)));
  eta = pow(sqrt(natoms) * pow(M_PI, 3.0) / (dvREAL * dvREAL), 1.0 / 3.0);

  printf("real cut-off rqeq=%10.6f\n", rqeq);
  printf("self energy cut-off=%10.6f\n", cuts);
  printf("accf = %12.8f\n", accf);
  printf("eta = %12.8f\n", eta);

  rradmax = 2.0 * accf * sqrt(eta);
  printf("1 rradmax  = %12.8f\n", rradmax);
  if ((rradmax * rradmax) > 140.0 * eta) rradmax = sqrt(140.0 * eta);
  printf("2 rradmax  = %12.8f\n", rradmax);

  radmax = accf / sqrt(eta);
  printf("1 radmax = %12.8f\n", radmax);
  if ((radmax * radmax) > 85.0 / eta) radmax = sqrt(85.0 / eta);
  printf("2 radmax = %12.8f\n", radmax);
  tweatpi = 2.0 * sqrt(eta / M_PI);
  bvecsz[0] = sqrt(bv1[0] * bv1[0] + bv1[1] * bv1[1] + bv1[2] * bv1[2]);
  bvecsz[1] = sqrt(bv2[0] * bv2[0] + bv2[1] * bv2[1] + bv2[2] * bv2[2]);
  bvecsz[2] = sqrt(bv3[0] * bv3[0] + bv3[1] * bv3[1] + bv3[2] * bv3[2]);
  bvecszmin = bvecsz[0];
  if (bvecsz[1] < bvecszmin) bvecszmin = bvecsz[1];
  if (bvecsz[2] < bvecszmin) bvecszmin = bvecsz[2];
  NUMR = ((int)(rradmax / bvecszmin)) + 2;
  allocate_double_array(&(bgridx),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridx\n");
  allocate_double_array(&(bgridy),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridy\n");
  allocate_double_array(&(bgridz),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridz\n");
  allocate_double_array(&(bgridcoef),
                        (2 * NUMR + 1) * (2 * NUMR + 1) * (2 * NUMR + 1),
                        "can not allocate bgridcoef\n");
  kvectors = prepareKSUM(rradmax, eta, NUMR, bv1, bv2, bv3, bgridx, bgridy,
                         bgridz, bgridcoef);
  av1sz = sqrt(av1[0] * av1[0] + av1[1] * av1[1] + av1[2] * av1[2]);
  av2sz = sqrt(av2[0] * av2[0] + av2[1] * av2[1] + av2[2] * av2[2]);
  av3sz = sqrt(av3[0] * av2[0] + av3[1] * av3[1] + av3[2] * av3[2]);
  avecszmin = av1sz;
  if (av2sz < avecszmin) avecszmin = av2sz;
  if (av3sz < avecszmin) avecszmin = av3sz;

  for (ip = 0; ip < npoints_to_gen; ip++) {
    xx = points_to_gen_xyz[3 * ip + 0];
    yy = points_to_gen_xyz[3 * ip + 1];
    zz = points_to_gen_xyz[3 * ip + 2];
    pot[ip] = 0.0;
    for (isource = 0; isource < nsource; isource++) {
      dx = xx - source_xyz[3 * isource + 0];
      dy = yy - source_xyz[3 * isource + 1];
      dz = zz - source_xyz[3 * isource + 2];
      qmatrixelement_grid(dx, dy, dz, av1, av2, av3, avecszmin, dvREAL, eta,
                          radmax, bgridx, bgridy, bgridz, bgridcoef, kvectors,
                          &vijre, &vijR);
      pot[ip] = pot[ip] + (vijre + vijR) * source_charge[isource];
    }
  }

  /* potential is returned in eV */
  for (ip = 0; ip < npoints_to_gen; ip++) pot[ip] = pot[ip] * angstoev;

  free(bgridx);
  free(bgridy);
  free(bgridz);
  free(bgridcoef);
  return;
}
