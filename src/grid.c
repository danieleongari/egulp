#include <math.h>
#include <stdlib.h> 
#include <stdio.h>

#include "utilsa.h"
#include "grid.h" 
#include "vdw.h" 
#include "cell.h"


#define AUTOANG_GRID (0.5291772)

#define MAXSTRL (128)
#define PATTERN (6)
#define PATTERNSTR1 ("%s %s %s %s\n")
#define PATTERNSTR2 ("%6d%12.6f%12.6f%12.6f\n")
#define PATTERNSTR3 ("%s %s %s %s %s %s\n")
#define VALIDPOINT (1.0)
#define INVALIDPOINT (100.0)

void ReadPATTERNSTR1(FILE **fp, char *str1, char *str2, char *str3, char *str4, 
  int *i1, double *d1, double *d2, double *d3) 
{  
	 fscanf(*fp, PATTERNSTR1, str1, str2, str3, str4);  
         (*i1) = atoi(str1); 
         (*d1) = atof(str2); 
         (*d2) = atof(str3); 
         (*d3) = atof(str4); 
         return;    
} 	 

void InitGrid(tGrid *gridptr, double res1, double res2, double res3, 
   double vdw_factor_i, double vdw_factor_f, int use_vdw_factor, double offset, 
   double *av1, double *av2, double *av3, int natoms,  
   int *atoms_num, double *atoms_xyz)
{
       /* ESP resolution is wired in... */         

	double av1sz, av2sz, av3sz;   
	double h1sz, h2sz, h3sz; 
	int  i, j, k, ntot, iatom, count;
	int icell, jcell, kcell; 
	double xx, yy, zz; 
	double xatom, yatom, zatom;
	double dist;  
	double vdw[103]; 
        double r1vdw, r2vdw; 
	int point_is_outside_flag; 

	SetVDWRADIUS(vdw); 
	
	av1sz = sqrt( (av1)[0]*(av1)[0]+(av1)[1]*(av1)[1]+(av1)[2]*(av1)[2]  );
	av2sz = sqrt( (av2)[0]*(av2)[0]+(av2)[1]*(av2)[1]+(av2)[2]*(av2)[2]  );
	av3sz = sqrt( (av3)[0]*(av3)[0]+(av3)[1]*(av3)[1]+(av3)[2]*(av3)[2]  );
	
	gridptr->resolution[0] = res1; 
	gridptr->resolution[1] = res2; 
	gridptr->resolution[2] = res3; 
        
	gridptr->vdw_factor_i = vdw_factor_i;
	gridptr->vdw_factor_f = vdw_factor_f;
	gridptr->use_vdw_factor = use_vdw_factor;
	if(! ((gridptr->use_vdw_factor == 0) || (gridptr->use_vdw_factor == 1)) ) { 
	     printf("grid.c: invalid use_vdw_factor %d\n", gridptr->use_vdw_factor); 
	     exit(-1); 
	}
	gridptr->offset = offset;
	if(gridptr->offset == -1) gridptr->offset = 1000.0; 
	 
        
	gridptr->gridn1 = (int) floor(av1sz/gridptr->resolution[0]); 
	gridptr->gridn2 = (int) floor(av2sz/gridptr->resolution[1]); 
	gridptr->gridn3 = (int) floor(av3sz/gridptr->resolution[2]); 
	
	for(i = 0; i < 3; i++) { 
		 gridptr->gridh1[i] = av1[i]/gridptr->gridn1; 
		 gridptr->gridh2[i] = av2[i]/gridptr->gridn2; 
		 gridptr->gridh3[i] = av3[i]/gridptr->gridn3; 
	}
	
	h1sz = sqrt( (gridptr->gridh1)[0]*(gridptr->gridh1)[0]+(gridptr->gridh1)[1]*(gridptr->gridh1)[1]+(gridptr->gridh1)[2]*(gridptr->gridh1)[2]  );
	h2sz = sqrt( (gridptr->gridh2)[0]*(gridptr->gridh2)[0]+(gridptr->gridh2)[1]*(gridptr->gridh2)[1]+(gridptr->gridh2)[2]*(gridptr->gridh2)[2]  );
	h3sz = sqrt( (gridptr->gridh3)[0]*(gridptr->gridh3)[0]+(gridptr->gridh3)[1]*(gridptr->gridh3)[1]+(gridptr->gridh3)[2]*(gridptr->gridh3)[2]  );
     
	 															printf("real space grid will be constructed....\n"); 
	 printf("%3d steps in A along lattice vector #1 %10.7f %10.7f %10.7f size %10.7f\n", 
	     gridptr->gridn1, gridptr->gridh1[0],  gridptr->gridh1[1], gridptr->gridh1[2], h1sz );  
	 printf("%3d steps in A along lattice vector #2 %10.7f %10.7f %10.7f size %10.7f\n", 
	     gridptr->gridn2, gridptr->gridh2[0],  gridptr->gridh2[1], gridptr->gridh2[2], h2sz );  
	 printf("%3d steps in A along lattice vector #3 %10.7f %10.7f %10.7f size %10.7f\n", 
	     gridptr->gridn3, gridptr->gridh3[0],  gridptr->gridh3[1], gridptr->gridh3[2], h3sz ); 
	 ntot = (gridptr->gridn1)*(gridptr->gridn2)*(gridptr->gridn3); 
	 printf("total number of points is %d\n", ntot); 
	 gridptr->dv3 =  cellVolume(gridptr->gridh1, gridptr->gridh2,  gridptr->gridh3);
	 printf("volume element is %10.7f\n", gridptr->dv3); 

	 allocate_int_array(&(gridptr->gridLxyz_inv), ntot, "can not allocate Lxyz_inv\n");
	 for(i = 0; i < ntot; i++)  (gridptr->gridLxyz_inv)[i] = 1;      	 
	 
	 for (i = 0; i <  gridptr->gridn1; i++) { 
	 for (j = 0; j <  gridptr->gridn2; j++) { 
	 for (k = 0; k < gridptr->gridn3; k++) { 
	        point_is_outside_flag = 1;
		xx = i*(gridptr->gridh1)[0]+j*(gridptr->gridh2)[0]+k*(gridptr->gridh3)[0]; 
		yy = i*(gridptr->gridh1)[1]+j*(gridptr->gridh2)[1]+k*(gridptr->gridh3)[1];  
		zz = i*(gridptr->gridh1)[2]+j*(gridptr->gridh2)[2]+k*(gridptr->gridh3)[2];  
		for (iatom = 0; iatom < natoms; iatom++) { 
		       r1vdw = (vdw[atoms_num[iatom]-1])*(AUTOANG_GRID)*(gridptr->vdw_factor_i); 
		       switch(gridptr->use_vdw_factor) {  
                            case  0: { 
			        r2vdw = r1vdw  + gridptr->offset; 
				break; 
			    }
			    case  1: { 
			         r2vdw = (vdw[atoms_num[iatom]-1])*(AUTOANG_GRID)*(gridptr->vdw_factor_f); 
				 break; 
			    }
			    default: { 
			         printf("invalid use_vdw_factor\n"); 
				 exit(-1);  
			    } 	 
		       }   

		       for(icell = -1; icell <= 1; icell++)  { 
		       for(jcell = -1; jcell <= 1; jcell++)  { 
		       for(kcell = -1; kcell <= 1; kcell++) {
		             xatom = (atoms_xyz)[iatom*3+0]+icell*(av1[0])+jcell*(av2[0])+kcell*(av3[0]);
                             yatom = (atoms_xyz)[iatom*3+1]+icell*(av1[1])+jcell*(av2[1])+kcell*(av3[1]); 
			     zatom = (atoms_xyz)[iatom*3+2]+icell*(av1[2])+jcell*(av2[2])+kcell*(av3[2]); 
			     dist = sqrt(  (xx-xatom)*(xx-xatom)+(yy-yatom)*(yy-yatom)+(zz-zatom)*(zz-zatom) ); 
			     /* add another goto statment? */
			     if(dist < r2vdw)  point_is_outside_flag = 0; 
			     if ( dist < r1vdw)   { 
			     	       gridptr->gridLxyz_inv[i*(gridptr->gridn2)*(gridptr->gridn3)+j*(gridptr->gridn3)+k] = 0; 
			               goto out_label_now;
		             }
		       }
		       }
		       }
		}
		 
		if(point_is_outside_flag == 1) gridptr->gridLxyz_inv[i*(gridptr->gridn2)*(gridptr->gridn3)+j*(gridptr->gridn3)+k] = 0; 
		out_label_now:
	           ;         
	  }
	  }
	  }
	  
	 gridptr->gridnl = 0; 
	 for (i = 0; i <  gridptr->gridn1; i++) { 
	 for (j = 0; j <  gridptr->gridn2; j++) { 
	 for (k = 0; k < gridptr->gridn3; k++) { 
	      if((gridptr->gridLxyz_inv)[i*(gridptr->gridn2)*(gridptr->gridn3)+j*(gridptr->gridn3)+k] == 1) gridptr->gridnl = gridptr->gridnl + 1; 
	  }
	  }
	  }

	  printf("resolution was set to %10.7f %10.7f %10.7f A\n",  
	       gridptr->resolution[0], gridptr->resolution[1], gridptr->resolution[2]); 
	  printf("the actual resolution might be different\n");      
	  printf("final grid vdw_factor was set to %10.7f\n", gridptr->vdw_factor_f); 
	  printf("initial grid vdw factor was set to %10.7f\n", gridptr->vdw_factor_i);
	  if(gridptr->use_vdw_factor == 0) printf("grid was constructed using fixed offset %10.6f\n", gridptr->offset); 
	  if(gridptr->use_vdw_factor == 1) printf("grid was constructed using final vdw factor %10.6f\n", gridptr->vdw_factor_f);  
	  printf("total number of points %d\n", ntot); 
	  printf("total number of simulation points %d\n", gridptr->gridnl); 
	  printf("allocating L1 L2 L3...\n"); 
	  allocate_int_array(&(gridptr->gridL1), gridptr->gridnl, "can not allocate gridL1\n");
	  allocate_int_array(&(gridptr->gridL2), gridptr->gridnl, "can not allocate gridL2\n");
	  allocate_int_array(&(gridptr->gridL3), gridptr->gridnl, "can not allocate gridL3\n");
	  printf("L1  L2 L3 are allocated\n"); 
	  count = 0; 
	  for (i = 0; i <  gridptr->gridn1; i++) { 
	  for (j = 0; j <  gridptr->gridn2; j++) { 
	  for (k = 0; k < gridptr->gridn3; k++) { 
               if( gridptr->gridLxyz_inv[i*(gridptr->gridn2)*(gridptr->gridn3)+j*(gridptr->gridn3)+k] == 1) { 
	                gridptr->gridL1[count] = i; 
			gridptr->gridL2[count] = j; 
			gridptr->gridL3[count] = k; 
			count = count + 1; 
               } 
	 }
	 }
	 } 
	 if(count != gridptr->gridnl) { 
	 	printf("error in number of simulation points: count=%d,gridnl=%d\n", count,  gridptr->gridnl); 
		exit(-1); 
	 } else printf("total number of simulation points %d\n", count); 
         return;   
} 

void DestroyGrid(tGrid *gridptr) 
{ 
	 free(gridptr->gridLxyz_inv); 
	 free(gridptr->gridL1); 
	 free(gridptr->gridL2); 
	 free(gridptr->gridL3); 
         return; 
} 

void WriteGridAsCube(char *fName, tGrid grid, int natoms,  int *atoms_num, double *atoms_xyz)
{
    int ixi, iyi, izi;
    int iatom;  
    int ifun_val; 
    char str[50]; 
    FILE *fp; 
    fp = fopen(fName, "w"); 
    printf("writing Gaussian cube %s\n", fName); 
    if(fp == NULL) { 
         printf("can not open file %s for writing\n", fName); 
	 exit(-1);   
    }

    fprintf(fp, "resolution(A) %10.6f %10.6f %10.6f\n",  grid.resolution[0], grid.resolution[1], grid.resolution[2]); 
    fprintf(fp, "vdw_factor_i %10.6f  vdw_factor_f %10.6f  filename %s simpoints %d use_vdw_factor %d offset %10.6f\n", grid.vdw_factor_i, grid.vdw_factor_f,  fName, grid.gridnl, grid.use_vdw_factor, grid.offset); 
    fprintf(fp, "%6d %11.6f %11.6f %11.6f\n", natoms, 0.0, 0.0, 0.0); 
    fprintf(fp, "%6d %11.6f %11.6f %11.6f\n", grid.gridn1, 
        grid.gridh1[0]/(AUTOANG_GRID), 
	grid.gridh1[1]/(AUTOANG_GRID), 
	grid.gridh1[2]/(AUTOANG_GRID)
     ); 
    fprintf(fp, "%6d %11.6f %11.6f %11.6f\n", grid.gridn2, 
        grid.gridh2[0]/(AUTOANG_GRID), 
	grid.gridh2[1]/(AUTOANG_GRID), 
	grid.gridh2[2]/(AUTOANG_GRID)  
     ); 
    fprintf(fp, "%6d %11.6f %11.6f %11.6f\n", grid.gridn3, 
        grid.gridh3[0]/(AUTOANG_GRID), 
	grid.gridh3[1]/(AUTOANG_GRID), 
	grid.gridh3[2]/(AUTOANG_GRID) 
     ); 
    for(iatom = 0; iatom < natoms; iatom++)
    	fprintf(fp, "%6d %11.6f %11.6f %11.6f %11.6f\n", 
	  atoms_num[iatom], 0.0,  
	  atoms_xyz[3*iatom+0]/(AUTOANG_GRID), 
	  atoms_xyz[3*iatom+1]/(AUTOANG_GRID), 
	  atoms_xyz[3*iatom+2]/(AUTOANG_GRID)  ); 
	 
    
    for (ixi = 0; ixi < grid.gridn1; ixi++) {
    for (iyi = 0; iyi < grid.gridn2; iyi++) {
    for (izi = 0; izi < grid.gridn3; izi++) {
    	  ifun_val = grid.gridLxyz_inv[ixi*(grid.gridn2)*(grid.gridn3)+iyi*(grid.gridn3)+izi]; 
	  if(ifun_val == 1) sprintf(str, "%13.6E ",  VALIDPOINT);
	  else sprintf(str, "%13.6E ",  INVALIDPOINT);  
	  fprintf(fp, str);  
	  if((izi % 6) == 5) fprintf(fp, "\n"); 
    }
          if( (grid.gridn3 % 6) != 0) fprintf(fp, "\n");  
    }
    }
    fclose(fp); 
    return;  
}

void WriteFunctionAsCube(char *fName, tGrid grid, int natoms,  int *atoms_num, double *atoms_xyz, double *func)
{
    int ixi, iyi, izi;
    int iatom, lpoint;  
    int ifun_val; 
    char str[50]; 
    FILE *fp; 
    fp = fopen(fName, "w"); 
    printf("writing Gaussian cube %s with a function\n", fName); 
    if(fp == NULL) { 
         printf("can not open file %s for writing\n", fName); 
	 exit(-1);   
    }

    fprintf(fp, "resolution(A) %10.6f %10.6f %10.6f\n",  grid.resolution[0], grid.resolution[1], grid.resolution[2]); 
    fprintf(fp, "vdw_factor_i %10.6f  vdw_factor_f %10.6f  filename %s simpoints %d use_vdw_factor %d offset %10.6f\n", grid.vdw_factor_i, grid.vdw_factor_f,  fName, grid.gridnl, grid.use_vdw_factor, grid.offset); 
    fprintf(fp, "%6d %11.6f %11.6f %11.6f\n", natoms, 0.0, 0.0, 0.0); 
    fprintf(fp, "%6d %11.6f %11.6f %11.6f\n", grid.gridn1, 
        grid.gridh1[0]/(AUTOANG_GRID), 
	grid.gridh1[1]/(AUTOANG_GRID), 
	grid.gridh1[2]/(AUTOANG_GRID)
     ); 
    fprintf(fp, "%6d %11.6f %11.6f %11.6f\n", grid.gridn2, 
        grid.gridh2[0]/(AUTOANG_GRID), 
	grid.gridh2[1]/(AUTOANG_GRID), 
	grid.gridh2[2]/(AUTOANG_GRID)  
     ); 
    fprintf(fp, "%6d %11.6f %11.6f %11.6f\n", grid.gridn3, 
        grid.gridh3[0]/(AUTOANG_GRID), 
	grid.gridh3[1]/(AUTOANG_GRID), 
	grid.gridh3[2]/(AUTOANG_GRID) 
     ); 
    for(iatom = 0; iatom < natoms; iatom++)
    	fprintf(fp, "%6d %11.6f %11.6f %11.6f %11.6f\n", 
	  atoms_num[iatom], 0.0,  
	  atoms_xyz[3*iatom+0]/(AUTOANG_GRID), 
	  atoms_xyz[3*iatom+1]/(AUTOANG_GRID), 
	  atoms_xyz[3*iatom+2]/(AUTOANG_GRID)  ); 
	 
    lpoint = 0; 
    for (ixi = 0; ixi < grid.gridn1; ixi++) {
    for (iyi = 0; iyi < grid.gridn2; iyi++) {
    for (izi = 0; izi < grid.gridn3; izi++) {
    	  ifun_val = grid.gridLxyz_inv[ixi*(grid.gridn2)*(grid.gridn3)+iyi*(grid.gridn3)+izi]; 
	  if(ifun_val == 1) {
	         sprintf(str, "%13.6E ",  func[lpoint]); 
		  lpoint = lpoint + 1;   
	 } 	  else sprintf(str, "%13.6E ",  1000.0);  
	  fprintf(fp, str);
	  if((izi % 6) == 5) fprintf(fp, "\n"); 
    }
          if( (grid.gridn3 % 6) != 0) fprintf(fp, "\n");  
    }
    }
    fclose(fp); 
    return;  
}

void InitGridFromCube(char *fName, tGrid *newgridptr, int natoms)
{
     double h1sz, h2sz, h3sz; 
     char tmpstr1[MAXSTRL]; 
     char tmpstr2[MAXSTRL]; 
     char tmpstr3[MAXSTRL];
     char tmpstr4[MAXSTRL]; 
     char tmpstr5[MAXSTRL]; 
     char tmpstr6[MAXSTRL]; 
     char tmpstr7[MAXSTRL]; 
     char tmpstr8[MAXSTRL]; 
     char tmpstr9[MAXSTRL]; 
     char tmpstr10[MAXSTRL];      
     char tmpstr11[MAXSTRL]; 
     char tmpstr12[MAXSTRL];      
     
     char tmpstr_array[PATTERN][MAXSTRL];
     int NumOfAtoms, iatom, icoord; 
     double cxx, cyy, czz; 
     int ixi, iyi, izi, i, j, k, itmp; 
     int ntot, full_lines, rem, count; 
     int *atoms_num; 
     double *atoms_xyz; 
     FILE *fp; 

     allocate_int_array(&(atoms_num), natoms, "can not allocate atoms_num in InitGridFromCube\n");
     allocate_double_array(&(atoms_xyz), 3*natoms, "can not allocate atoms_xyz in InitGridFromCube\n"); 
      
     fp = fopen(fName, "r"); 
     printf("opening cube file %s for reading\n", fName); 
     if(fp == 0) {  
             printf("Error opening file %s\n", fName);  
	     exit(-1); 
      }    
      printf("cube file %s is open\n", fName);

      fscanf(fp, "%s %s %s %s\n", tmpstr1, tmpstr2, tmpstr3, tmpstr4);
      newgridptr->resolution[0] = atof(tmpstr2); 
      newgridptr->resolution[1] = atof(tmpstr3); 
      newgridptr->resolution[2] = atof(tmpstr4); 
      printf("resolution (A)   %10.6f %10.6f %10.6f\n", newgridptr->resolution[0], newgridptr->resolution[1], newgridptr->resolution[2]); 
/*   fprintf(fp, "vdw_factor_i %10.6f  vdw_factor_f %10.6f  filename %s simpoints %d use_vdw_factor %d offset %10.6f\n", grid.vdw_factor_i, grid.vdw_factor_f,  fName, grid.gridnl, grid.use_vdw_factor, grid.offset); */
      fscanf(fp, "%s %s %s %s %s %s %s %s %s %s %s %s\n",  tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12); 
      newgridptr->vdw_factor_i = atof(tmpstr2); 
      printf("vdw_factor_i %10.6f\n", newgridptr->vdw_factor_i); 

      newgridptr->vdw_factor_f = atof(tmpstr4); 
      printf("vdw_factor_f %10.6f\n", newgridptr->vdw_factor_f); 
      
      itmp = atof(tmpstr8); 
      printf("number of simulation points is read as %d\n", itmp); 
      
      newgridptr->use_vdw_factor = atoi(tmpstr10);
      newgridptr->offset = atof(tmpstr12);
      printf("use_vdw_factor = %d offset = %10.6f\n",   newgridptr->use_vdw_factor, newgridptr->offset);  	   
      /* read number of atoms and origin */ 
      ReadPATTERNSTR1(&fp, tmpstr1, tmpstr2, tmpstr3, tmpstr4,  &NumOfAtoms, &cxx, &cyy, &czz);
      /* read number of points in e_1, e_2, e_3 and unit vectors */ 
      ReadPATTERNSTR1(&fp, tmpstr1, tmpstr2, tmpstr3, tmpstr4, &(newgridptr->gridn1), &(newgridptr->gridh1[0]), &(newgridptr->gridh1[1]), &(newgridptr->gridh1[2]) ); 	 
      ReadPATTERNSTR1(&fp, tmpstr1, tmpstr2, tmpstr3, tmpstr4, &(newgridptr->gridn2), &(newgridptr->gridh2[0]), &(newgridptr->gridh2[1]), &(newgridptr->gridh2[2]) );
      ReadPATTERNSTR1(&fp, tmpstr1, tmpstr2, tmpstr3, tmpstr4, &(newgridptr->gridn3), &(newgridptr->gridh3[0]), &(newgridptr->gridh3[1]), &(newgridptr->gridh3[2]) ); 	 
      for(icoord = 0; icoord < 3; icoord++) { 
	      newgridptr->gridh1[icoord] = (newgridptr->gridh1[icoord])*AUTOANG_GRID ; 
	      newgridptr->gridh2[icoord] = (newgridptr->gridh2[icoord])*AUTOANG_GRID ; 
	      newgridptr->gridh3[icoord] = (newgridptr->gridh3[icoord])*AUTOANG_GRID ; 
      }
       newgridptr->dv3 =  cellVolume(newgridptr->gridh1, newgridptr->gridh2,  newgridptr->gridh3);
       printf("%d atoms will be read in\n", NumOfAtoms); 
       if(NumOfAtoms != natoms)  { 
	      printf("mismatched number of atoms:\n"); 
	      printf("%d atoms from cube vs. %d atoms from geometry\n", NumOfAtoms, natoms); 
	} 
        
	
	for(iatom = 0; iatom < NumOfAtoms; iatom++) {
        	 fscanf(fp, "%s %s %s %s %s\n", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5);  
                 atoms_num[iatom] = atoi(tmpstr1); 
                 atoms_xyz[iatom*3+0] = atof(tmpstr3)*(AUTOANG_GRID); 
                 atoms_xyz[iatom*3+1] = atof(tmpstr4)*(AUTOANG_GRID);
		 atoms_xyz[iatom*3+2] = atof(tmpstr5)*(AUTOANG_GRID);  
         }

         full_lines = (newgridptr->gridn3)/(PATTERN); 
         rem = (newgridptr->gridn3)%(PATTERN); 
	 h1sz = sqrt( (newgridptr->gridh1)[0]*(newgridptr->gridh1)[0]+(newgridptr->gridh1)[1]*(newgridptr->gridh1)[1]+(newgridptr->gridh1)[2]*(newgridptr->gridh1)[2]  );
	 h2sz = sqrt( (newgridptr->gridh2)[0]*(newgridptr->gridh2)[0]+(newgridptr->gridh2)[1]*(newgridptr->gridh2)[1]+(newgridptr->gridh2)[2]*(newgridptr->gridh2)[2]  );
	 h3sz = sqrt( (newgridptr->gridh3)[0]*(newgridptr->gridh3)[0]+(newgridptr->gridh3)[1]*(newgridptr->gridh3)[1]+(newgridptr->gridh3)[2]*(newgridptr->gridh3)[2]  );
	 printf("full_lines=%d rem=%d\n", full_lines, rem); 
	 printf("%3d steps along lattice vector #1 %10.7f %10.7f %10.7f size %10.7f\n", 
	     newgridptr->gridn1, newgridptr->gridh1[0],  newgridptr->gridh1[1], newgridptr->gridh1[2], h1sz );  
	 printf("%3d steps along lattice vector #2 %10.7f %10.7f %10.7f size %10.7f\n", 
	     newgridptr->gridn2, newgridptr->gridh2[0],  newgridptr->gridh2[1], newgridptr->gridh2[2], h2sz );  
	 printf("%3d steps along lattice vector #3 %10.7f %10.7f %10.7f size %10.7f\n", 
	     newgridptr->gridn3, newgridptr->gridh3[0],  newgridptr->gridh3[1], newgridptr->gridh3[2], h3sz ); 
	 ntot = (newgridptr->gridn1)*(newgridptr->gridn2)*(newgridptr->gridn3); 
	 printf("total number of points is %d\n", ntot); 
         printf("volume element is %10.7f\n", newgridptr->dv3); 

	 allocate_int_array(&(newgridptr->gridLxyz_inv), ntot, "can not allocate Lxyz_inv\n");    
	/* printf("entering the loop...\n"); */
	 newgridptr->gridnl = 0; 
         for(ixi = 0; ixi < newgridptr->gridn1; ixi++) {  
         for(iyi = 0; iyi < newgridptr->gridn2; iyi++) { 
	          for(izi = 0; izi < newgridptr->gridn3; izi++)  newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+izi] = 0;   		
                  for(i = 0; i < full_lines; i++) {  
		         fscanf(fp, PATTERNSTR3, &(tmpstr_array[0][0]), &(tmpstr_array[1][0]), &(tmpstr_array[2][0]), 
			                                      &(tmpstr_array[3][0]), &(tmpstr_array[4][0]), &(tmpstr_array[5][0]) ); 
			for(j = 0; j < 	PATTERN; j++) 			      
		              if( atof(&(tmpstr_array[j][0])) == VALIDPOINT) { 
			          newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+i*PATTERN+j]  = 1; 
				  (newgridptr->gridnl)++; 
		              }  
		   }	         		  
	       	   
		   switch(rem) { 
                       case 0:  
		           break; 
                       case 1:  
		            fscanf(fp, "%s\n", tmpstr1); 
		            if( atof(tmpstr1) == VALIDPOINT) { 
			                   newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+0]  = 1; 
                                           (newgridptr->gridnl)++; 
                            }        
                            break; 
	               case 2: 
		            fscanf(fp, "%s %s\n", tmpstr1, tmpstr2); 
			     if( atof(tmpstr1) == VALIDPOINT) { 
			                   newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+0]  = 1; 
                                           (newgridptr->gridnl)++; 
                             }  
			     if( atof(tmpstr2) == VALIDPOINT) { 
			             newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+1]  = 1; 
                                    (newgridptr->gridnl)++;  
                            }  
                            break; 
	               case 3:     
	                    fscanf(fp, "%s %s %s\n", tmpstr1, tmpstr2, tmpstr3); 
			    if( atof(tmpstr1) == VALIDPOINT) 
			    { 
			            newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+0]  = 1; 
                                    (newgridptr->gridnl)++; 
                             }  
			    if( atof(tmpstr2) == VALIDPOINT) 
			    { 
			            newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+1]  = 1; 
                                    (newgridptr->gridnl)++; 
                             } 
			     if( atof(tmpstr3) == VALIDPOINT) 
			     { 
			             newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+2]  = 1; 
                                     (newgridptr->gridnl)++; 
                              }    
                              break; 
	               case 4:     
	                    fscanf(fp, "%s %s %s %s\n", tmpstr1, tmpstr2, tmpstr3, tmpstr4); 
			    if( atof(tmpstr1) == VALIDPOINT) {
			            newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+0]  = 1; 
                                    (newgridptr->gridnl)++; 
                            } 
			    if( atof(tmpstr2) == VALIDPOINT) {
			            newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+1]  = 1; 
                                    (newgridptr->gridnl)++; 
                            }
			    if( atof(tmpstr3) == VALIDPOINT) { 
			             newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+2]  = 1; 
                                     (newgridptr->gridnl)++; 
                            }
			    if( atof(tmpstr4) == VALIDPOINT) { 
			             newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+3]  = 1; 
                                    (newgridptr->gridnl)++; 
                            } 
	                    break;
	               case 5:     
	                    fscanf(fp, "%s %s %s %s %s\n", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5); 
			    if( atof(tmpstr1) == VALIDPOINT) { 
			              newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+0]  = 1; 
				       (newgridptr->gridnl)++; 
			    }  	       
			    if( atof(tmpstr2) == VALIDPOINT) { 
			              newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+1]  = 1; 
				      (newgridptr->gridnl)++;
			    } 	       
			    if( atof(tmpstr3) == VALIDPOINT) { 
			              newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+2]  = 1; 
				       (newgridptr->gridnl)++; 
			    } 	       
			    if( atof(tmpstr4) == VALIDPOINT) { 
			              newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+3]  = 1; 
				       (newgridptr->gridnl)++; 
		            } 		       
			    if( atof(tmpstr5) == VALIDPOINT) { 
			              newgridptr->gridLxyz_inv[ixi*(newgridptr->gridn2)*(newgridptr->gridn3)+iyi*(newgridptr->gridn3)+full_lines*PATTERN+4]  = 1; 
				     (newgridptr->gridnl)++; 
		            } 		       
                            break;
                       default:
                           printf("rem error in InitGridFromCube\n"); 
                           exit(-1); 
                    } 		   
     } 
     }
      printf("total number of simulation points %d\n", newgridptr->gridnl); 
      allocate_int_array(&(newgridptr->gridL1), newgridptr->gridnl, "can not allocate gridL1\n");
      allocate_int_array(&(newgridptr->gridL2), newgridptr->gridnl, "can not allocate gridL2\n");
      allocate_int_array(&(newgridptr->gridL3), newgridptr->gridnl, "can not allocate gridL3\n");
      
      count = 0; 
      for (i = 0; i <  newgridptr->gridn1; i++) { 
      for (j = 0; j <  newgridptr->gridn2; j++) { 
      for (k = 0; k < newgridptr->gridn3; k++) { 
               if( newgridptr->gridLxyz_inv[i*(newgridptr->gridn2)*(newgridptr->gridn3)+j*(newgridptr->gridn3)+k] > 0) { 
	                newgridptr->gridL1[count] = i; 
			newgridptr->gridL2[count] = j; 
			newgridptr->gridL3[count] = k; 
			count = count + 1; 
               } 
	}
	}
	} 

     free(atoms_num); 
     free(atoms_xyz); 
     
     fclose(fp);  
     return; 
}	 

double IntegrateFromGrid(int nl, double *f, double dv3)
{
        double res = 0.0; 
	int i; 
	for(i = 0; i < nl; i++) res = res + f[i]; 
	res = res*dv3; 
        return(res); 
} 

#undef AUTOANG_GRID
#undef MAXSTRL 
#undef PATTERN
#undef PATTERNSTR1
#undef PATTERNSTR2
#undef PATTERNSTR3
#undef VALIDPOINT 
#undef INVALIDPOINT

