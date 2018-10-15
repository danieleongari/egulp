#ifndef _ESEQPARAM_H_
#define _ESEQPARAM_H_

typedef struct tESEQParamsData {
      /* these are input parameters */ 
       int number_of_bond_types; 
       int *bond_atom_types;             /* for each bond C-C, C-H, for example */ 
       double *elneg_bondcorr;          /* electronegativity correction */ 
       double *hard_bondcorr;           /* hardness correction */
       /* parameters to be used in calcs */ 
        
       double *TMatrix; 
       double *TMatrixPINV; 	
       double *elneg_bondcorr_geom; 
       double *hard_bondcorr_geom; 
}  tESEQParams; 

void InitESEQParams(double *qeqchi,  double *qeqmui, int *qeqids, int MaxElement, char *fname, tESEQParams *eseqptr);

void InitTransferMatrices(tESEQParams *eseqptr, int ncharges, int nbonds,  int *nbonds_indexij);
void DestroyTransferMatrices(tESEQParams *eseqptr);  

void BuildParamsGEOOrder(tESEQParams *eseqptr,  int tbonds, int *tbonds_indexij);
void BuildBondHCorrections(tESEQParams *eseqptr,  int natoms, int nbonds,  int *nbonds_indexij, 
            int *nbonds_type,  int tbonds, int *tbonds_indexij, double *BondH); 
      /* ELECTRONEGATIVITY BOND CORRECTIONS E = T*c */

void BuildBondECorrections(tESEQParams *eseqptr, int natoms, int nbonds,  int *nbonds_indexij, 
            int *nbonds_type, int tbonds,  int *tbonds_indexij, double *BondE); 


void FreeESEQParams(tESEQParams *eseqptr); 


#endif
