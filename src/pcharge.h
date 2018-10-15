#ifndef _PCHARGE_H_
#define _PCHARGE_H_

typedef struct tPointChargeData {
    int npoints; 
    double *points_xyz; 
    double *points_charge; 
    double *atatom_pot; 
    double *atpoint_pot; 
} tPointCharge;

void InitPointCharge(tPointCharge *pcptr, double *av1, double *av2, double *av3, int natoms, double *atoms_xyz, char *fName ); 
void DestroyPointCharge( tPointCharge *pcptr); 

#endif
