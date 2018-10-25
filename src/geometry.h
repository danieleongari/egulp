#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

typedef struct tGeometryData {
    double av1[3];            
    double av2[3];            
    double av3[3];
    double cell_len[3];
    double cell_ang[3];
    int natoms; 
    double *atoms_xyz; 
    double *atoms_fract;
    int *atoms_num; 
    int *atoms_base_num; 
    double *atoms_charge; 
    double *input_atoms_charge; 
    double dv; 
   /* number of distinct atomic numbers - set in bonds */ 
    int natom_types; 
    
   /* number of bonds */ 
    int nbonds; 
   /* bond indices for atoms */  
    int *nbonds_indexij;
   /* bond indices for atom types */   
    int *nbonds_type; 
    /*  k: 1..nbonds 
          nbonds_indexij[2*k+0] 
          nbonds_indexij[2*k+1] 
	  will return 1..natoms 
	  nbonds_type[k] will return 0..tbonds-1 */ 
   /* number of bond types */  
    int tbonds; 
   /* bond type indices  */   
    int *tbonds_indexij; 

} tGeometry;

void InitGeometry(tGeometry *geometry, char *fname);


void MapAtomsIntoCentralCell(tGeometry *geometry); 
void MapAtomsCheck(tGeometry *geometry) ; 

void InitBonds(tGeometry *geometry); 
void PrintBondInfo(tGeometry *geometry); 

void DestroyGeometry(tGeometry *geometry);
void DestroyBonds(tGeometry *geometry);  


#endif
