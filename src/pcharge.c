#include <math.h>
#include <stdlib.h> 
#include <stdio.h>

#include "utilsa.h"
#include "pcharge.h" 
#define MAXSTRL (128)

void InitPointCharge(tPointCharge *pcptr, double *av1, double *av2, double *av3, int natoms, double *atoms_xyz, char *fName  )
{
	int i; 
	FILE *fp; 
        char tmpstr1[MAXSTRL]; 
        char tmpstr2[MAXSTRL]; 
        char tmpstr3[MAXSTRL];	
        char tmpstr4[MAXSTRL];	
	
        fp = fopen(fName, "r"); 
        printf("opening point charge file %s for reading\n", fName); 
        if(fp != NULL) {  
               printf("point charge file %s is open\n", fName);
               fscanf(fp, "%s\n", tmpstr1);
               pcptr->npoints = atoi(tmpstr1);
               printf("number of point charges is %d\n", pcptr->npoints); 
               if( pcptr->npoints < 1)  { 
                      printf("Invalid number of point charges %d...exiting\n", pcptr->npoints); 
                      exit(-1); 
               }
               allocate_double_array(&(pcptr->points_xyz), 3*(pcptr->npoints), "can not allocate points_xyz\n"); 
               allocate_double_array(&(pcptr->points_charge), pcptr->npoints, "can not allocate points_charge\n");  
               allocate_double_array(&(pcptr->atatom_pot), natoms, "can not allocate atatom_pot\n");     
	       allocate_double_array(&(pcptr->atpoint_pot), pcptr->npoints, "can not allocate atpoint_pot\n");     
               for(i = 0; i < natoms; i++) pcptr->atatom_pot[i] = 0.0; 
	       for(i = 0; i < pcptr->npoints; i++) pcptr->atpoint_pot[i] = 0.0; 
	       for(i = 0; i < pcptr->npoints; i++) pcptr->points_charge[i] = 0.0; 
               for(i = 0; i < pcptr->npoints; i++) {
                      fscanf(fp, "%s %s %s %s\n", tmpstr1, tmpstr2, tmpstr3, tmpstr4);
                      pcptr->points_xyz[3*i+0] = atof(tmpstr1);
		      pcptr->points_xyz[3*i+1] = atof(tmpstr2); 
		      pcptr->points_xyz[3*i+2] = atof(tmpstr3);
		      pcptr->points_charge[i] =  atof(tmpstr4);
                      printf("point charge %5d: %10.7f %10.7f %10.7f %10.7f\n", i+1, 
                                pcptr->points_xyz[3*i+0],                          
                                pcptr->points_xyz[3*i+1], 
                                pcptr->points_xyz[3*i+2],
				pcptr->points_charge[i] );  
	       } 	
	 } else { 
	              printf("can not open point charge file %s\n", fName); 
		      exit(-1);  
	 }             
         return;   
} 


void DestroyPointCharge( tPointCharge *pcptr) 
{ 
	 free(pcptr->points_xyz); 
	 free(pcptr->points_charge); 
	 free(pcptr->atatom_pot); 
	 free(pcptr->atpoint_pot); 
         return; 
} 


#undef MAXSTRL 
