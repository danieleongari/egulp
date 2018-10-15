#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "econfig.h"

void check_resolution(double dhsz);
void check_vdw_factor(double vdw_factor); 


void InitEconfigDefaults(tConfigure *confptr) 
{
       confptr->dh1sz  = 0.25; 
       confptr->dh2sz  = 0.25;
       confptr->dh3sz  = 0.25; 
       confptr->vdw_factor_i = 1.0; 
       confptr->vdw_factor_f = 2.0;
       confptr->use_vdw_factor = 1; 
       confptr->offset = 3.0;  
       confptr->build_grid = 0; 
       confptr->build_grid_from_scratch = 1; 
       confptr->save_grid = 0; 
       strcpy(confptr->input_grid_file, "none");
       strcpy(confptr->output_grid_file, "none"); 
       strcpy(confptr->output_pot_file, "none"); 
       confptr->calculate_pot_diff = 0; 
       confptr->calculate_pot = 0;        
       confptr->skip_everything = 0; 
       confptr->point_charges_present = 0; 
       confptr->include_pceq = 0;
       /* imethod=0 standard QEQ 
       imethod=1 split charge */  
       confptr->imethod = 0 ; 
}

void InitEconfig(tConfigure *confptr, char *fname)
{
     FILE *fp; 
     int i; 
     char tmpstr1[128]; 
     char tmpstr2[128]; 
     char tmpstr3[128];
     char tmpstr4[128];
     char tmpstr5[128];
     char tmpstr6[128];
     char tmpstr7[128];
     char tmpstr8[128];
     char tmpstr9[128];
     char tmpstr10[128];
   
     fp = fopen(fname, "r"); 
     printf("opening config file %s\n", fname); 
     if(fp == NULL) {  
            printf("Error opening file %s\n", fname);  
	    exit(-1);  
     }   else { 	 
	 printf("config file %s is open\n", fname); 
	 fscanf(fp, "%s %s", tmpstr1, tmpstr2);  
	 confptr->build_grid = atoi(tmpstr2);
	 printf("build_grid is read as %d\n", confptr->build_grid);
	 if(! ((confptr->build_grid == 0) || (confptr->build_grid == 1))   ) {
            printf("invalid build_grid=%d. Should be 1(True) or 0(False)\n", confptr->build_grid);  
	    exit(-1);  
	 }
	 fscanf(fp, "%s %s %s %s %s %s %s %s %s %s", tmpstr1, tmpstr2, confptr->input_grid_file, tmpstr4, tmpstr5, tmpstr6, tmpstr7,  tmpstr8, tmpstr9, tmpstr10);  
	 confptr->build_grid_from_scratch = atoi(tmpstr2);
         printf("build_grid_from_scratch is read as %d\n", confptr->build_grid_from_scratch);
	 if(! ((confptr->build_grid_from_scratch == 0) || (confptr->build_grid_from_scratch == 1))   ) {
             printf("invalid build_grid_from_scratch=%d. Should be 1(True) or 0(False)\n", confptr->build_grid_from_scratch);  
	     exit(-1);  
	 }
	 
	 printf("checking out resolution #1\n"); 
	 confptr->dh1sz = atof(tmpstr4); 
	 check_resolution(confptr->dh1sz); 
	 
	 printf("checking out resolution #2\n"); 
         confptr->dh2sz = atof(tmpstr5);
	 check_resolution(confptr->dh2sz); 
	 
	 printf("checking out resolution #3\n"); 
         confptr->dh3sz = atof(tmpstr6);
	 check_resolution(confptr->dh3sz); 

	 printf("checking out vdw_factor_i\n"); 
         confptr->vdw_factor_i = atof(tmpstr7); 
	 check_vdw_factor(confptr->vdw_factor_i); 

	 printf("checking out vdw_factor_f\n"); 
         confptr->vdw_factor_f = atof(tmpstr8); 
	 check_vdw_factor(confptr->vdw_factor_f); 

         confptr->use_vdw_factor = atoi(tmpstr9); 
         printf("use_vdw_factor is read as %d\n", confptr->use_vdw_factor);	 
	 if(! ((confptr->use_vdw_factor == 0) || (confptr->use_vdw_factor == 1))   ) {
             printf("invalid use_vdw_factor=%d. Should be 1(True) or 0(False)\n", confptr->use_vdw_factor);  
	     exit(-1);  
	 }

         confptr->offset = atoi(tmpstr10); 
         printf("offset is read as %10.6f\n", confptr->offset);	 
	 if(confptr->offset == -1) printf("All region minus vdw will be used if use_vdw_factor==0\n"); 

	 if ( (confptr->build_grid == 1) && (confptr->build_grid_from_scratch == 0)  )  { 
	 	printf("grid will be loaded from file %s\n", confptr->input_grid_file); 
		printf("resolution vdw_factors and offsets will be read from %s, not from the configure file\n", confptr->input_grid_file); 
	 }
	 
	 fscanf(fp, "%s %s %s", tmpstr1, tmpstr2, confptr->output_grid_file);   	
         confptr->save_grid = atoi(tmpstr2); 
	 if(! ((confptr->save_grid == 0) || (confptr->save_grid == 1))   ) {
             printf("invalid save_grid=%d. Should be 1(True) or 0(False)\n", confptr->save_grid);  
	     exit(-1);  
	 }
	 
	 if ( (confptr->build_grid == 1) && (confptr->save_grid == 1)  )   printf("grid will be saved into file %s\n", confptr->output_grid_file); 
	 
	 fscanf(fp, "%s %s", tmpstr1, tmpstr2);
	 confptr->calculate_pot_diff = atoi(tmpstr2); 
	 if(! ((confptr->calculate_pot_diff  == 0) || (confptr->calculate_pot_diff  == 1))   ) {
             printf("invalid calculate_pot_diff =%d. Should be 1(True) or 0(False)\n", confptr->calculate_pot_diff );  
	     exit(-1);  
	 }
	 if (confptr->calculate_pot_diff  == 1) printf("electrostatic potential based on the output-input charges will be constructed\n"); 
	 else printf("electrostatic potential based on the output-input charges will not be constructed\n"); 

	 fscanf(fp, "%s %s %s", tmpstr1, tmpstr2, confptr->output_pot_file);
	 confptr->calculate_pot = atoi(tmpstr2); 
	 if(!  ( (confptr->calculate_pot  == 1) || (confptr->calculate_pot  == 2) || (confptr->calculate_pot  == 0) )    ) {
             printf("invalid calculate_pot =%d. Should be (1,2 - calculate input,output or 0 - do not calculate)\n", confptr->calculate_pot );  
	     exit(-1);  
	 }
	 if (confptr->calculate_pot == 1) printf("electrostatic potential based on the input  charges will be constructed\n"); 
	 else if (confptr->calculate_pot == 2) printf("electrostatic potential based on the output  charges will be constructed\n");  
	 else printf("electrostatic potential based on the input  or output charges will not be constructed\n"); 
         if ((confptr->calculate_pot == 1)  || (confptr->calculate_pot == 2)) printf("electrostatic potential will be saved in file %s\n", confptr->output_pot_file); 
	  	 
	 fscanf(fp, "%s %s", tmpstr1, tmpstr2);
	 confptr->skip_everything = atoi(tmpstr2); 
	 printf("skip_everything is read as %d\n", confptr->skip_everything);
	 if(! ((confptr->skip_everything   == 0) || (confptr->skip_everything  == 1))   ) {
             printf("invalid skip_everything  =%d. Should be 1(True) or 0(False)\n", confptr->skip_everything );  
	     exit(-1);  
	 }
	 if (confptr->skip_everything  == 1) printf("calculation of charges and potentials will be skipped\n"); 
	 else printf("calculation of charges and potentials will not be skipped\n"); 

	 fscanf(fp, "%s %s", tmpstr1, tmpstr2);
	 confptr->point_charges_present = atoi(tmpstr2); 
	 printf("point_charges_present is read as %d\n", confptr->point_charges_present);
	 if(! ((confptr->point_charges_present   == 0) || (confptr->point_charges_present  == 1))   ) {
             printf("invalid point_charges_present  =%d. Should be 1(True) or 0(False)\n", confptr->point_charges_present );  
	     exit(-1);  
	 }
	 if (confptr->point_charges_present  == 1) printf("point charges are present...\n"); 
	 else printf("point charges are not present\n"); 

	 fscanf(fp, "%s %s", tmpstr1, tmpstr2);
	 confptr->include_pceq = atoi(tmpstr2); 
	 printf("include_pceq is read as %d\n", confptr->include_pceq);
	 if(! ((confptr->include_pceq   == 0) || (confptr->include_pceq  == 1))   ) {
             printf("invalid include_pceq  =%d. Should be 1(True) or 0(False)\n", confptr->include_pceq );  
	     exit(-1);  
	 }
	 if (confptr->include_pceq  == 1) { 
	        printf("response to point charges will be computed\n"); 
		printf("equations will be modified to account for point charges\n"); 
		printf("Vq term will be reported separately\n"); 
	 } else { 
	        printf("response to point charges will not be computed\n"); 	 
		printf("equations will not be modified to account for point charges\n"); 
		printf("Vq term will be reported separately\n"); 
	 }	 
	 
	 fscanf(fp, "%s %s", tmpstr1, tmpstr2);
	 confptr->imethod = atoi(tmpstr2); 
	 printf("imethod is read as %d\n", confptr->imethod);
	 if(! ((confptr->imethod   == 0) || (confptr->imethod  == 1))   ) {
             printf("invalid imethod =%d. Should be 1(SPLIT) or 0(QEQ)\n", confptr->imethod );  
	     exit(-1);  
	 }
	 if (confptr->imethod  == 1) printf("SPLIT CHARGE CALCULATION\n"); 
	 else if (confptr->imethod  == 0) printf("QEQ CALCULATION\n"); 
	 else  printf("INVALID METHOD OPTION...\n"); 
 	 fclose(fp); 
     } 
     return; 
}

void DestroyEconfig(tConfigure conf)
{
        return; 
}

void check_resolution(double dhsz)
{
	if( (dhsz < 0.05) || (dhsz > 0.9)) { 
	        printf("resolution %10.6f is in invalid range\n", dhsz); 
		printf("valid range is from 0.05 to 0.9 A\n"); 
		exit(-1); 
	} else printf("valid resolution %10.6f\n", dhsz); 
	return; 
}

void check_vdw_factor(double vdw_factor)
{
	if  ((vdw_factor < 0.5) ||   (vdw_factor > 3.0)) { 
	        printf("vdw_factor %10.6f is in invalid range\n", vdw_factor); 
		printf("the valid range is 0.5..3.0\n"); 
		exit(-1); 
	} else printf("valid vdw_factor %10.6f\n", vdw_factor); 
	return; 
} 
