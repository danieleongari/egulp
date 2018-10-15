#include  "scf.h"
#include "utilsa.h"


void InitEGULPScf(tSCF *SCF, double *atoms_charge, int *atom_num, int NumberOfAtoms)
{
    int i, j; 
    SCF->niter = 1; 
    SCF->literate = 0; 
    SCF->lconverged = 0; 
    SCF->current_iter = 0; 
    SCF->natoms = NumberOfAtoms;
    allocate_double_array(&(SCF->atoms_charge), SCF->natoms, "can not allocate atom_charge in SCF\n");      
    allocate_double_array(&(SCF->atoms_charge_old), SCF->natoms, "can not allocate atom_charge_old in SCF\n");  
    for (i = 0; i < NumberOfAtoms; i++) { 
           SCF->atoms_charge[i] = atoms_charge[i];
	   SCF->atoms_charge_old[i] =  atoms_charge[i];
    }
      	       
    for (i = 0; i <  NumberOfAtoms; i++) {
          if(atom_num[i] == 1) { 
	        SCF->literate = 1; 
		SCF->niter = 50; 
                 break;  
           } 
    } 

    SCF->qeqscfcrit =  1.0e-5; 
/*    original mixing parameter ... */     
/*    SCF->MixParameter = 0.75;
      SCF->qeqscfcrit =  1.0e-6;   */ 
    SCF->MixParameter = 0.3; 
    allocate_double_array(&(SCF->energy), SCF->niter, "can not allocate energy in SCF\n");      
    for(i = 0; i < SCF->niter; i++) SCF->energy[i]  = 0.0 ; 
    allocate_double_array(&(SCF->QEQ_JArray), (SCF->natoms)*(SCF->natoms), "can not allocate energy in SCF\n");      
    allocate_double_array(&(SCF->QEQ_ChiArray), SCF->natoms, "can not allocate energy in SCF\n");      
    for(i = 0; i < SCF->natoms; i++) { 
         SCF->QEQ_ChiArray[i] = 0.0; 
	 for(j = 0; j < SCF->natoms; j++) SCF->QEQ_JArray[i*SCF->natoms+j] = 0.0 ; 
     }	 
    return;  
} 

void PrintEGULPScf(tSCF *SCF) 
{
        if(SCF->literate == 1) { 
	      printf("hydrogen is present\n"); 
	      printf("iterative solution will be attempted\n"); 
	      printf("convergence criteria = %10.7f\n", SCF->qeqscfcrit ); 
	      printf("mixing parameter = %10.7f\n", SCF->MixParameter); 
	}
	printf("literate = %d\n", SCF->literate); 
	printf("max number of iterations is %d\n", SCF->niter); 
	printf("number of atoms in SCF is %d\n", SCF->natoms); 
        return; 
} 

void DestroyEGULPSCF(tSCF *SCF)
{
       free(SCF->atoms_charge); 
       free(SCF->atoms_charge_old); 
       free(SCF->energy); 
       free(SCF->QEQ_JArray); 
       free(SCF->QEQ_ChiArray); 
       return; 
}

double getQEQEnergy(int natoms, double *atoms_charge, double *JArray, double *ChiArray) 
{
            int i, j; 
	    double res = 0; 
	    for(i = 0; i < natoms; i++) { 
	          res = res + ChiArray[i]*atoms_charge[i]; 
		  for (j = 0; j < natoms; j++) res = res + 0.5*JArray[i*natoms+j]*atoms_charge[i]*atoms_charge[j]; 
            }
	    return (res) ;  
}

double getQEQPCEnergy(int natoms, double *atoms_charge, double *pcpot)
{
            int i;
	    double res = 0; 
	    for(i = 0; i < natoms; i++) res = res + pcpot[i]*atoms_charge[i]; 
	    return (res) ;  
}
