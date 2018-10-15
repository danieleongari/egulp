#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

typedef struct tGeometryData {
    double av1[3];            
    double av2[3];            
    double av3[3];
    int natoms; 
    double *atoms_xyz; 
    int *atoms_num; 
    double *atoms_charge; 
    double dv; 
    int nbonds; 
    int *nbonds_indexij; 
} tGeometry;

void InitGeometry(tGeometry *geometry, char *fname);
void InitBonds(tGeometry *geometry); 

void DestroyGeometry(tGeometry *geometry, int do_bonds); 
double  cellVolume(double *a1, double *a2,  double *a3);

void vectorProduct(double *a2, double *a3, double *res); 

#endif
