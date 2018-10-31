#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
 #include <gsl/gsl_linalg.h>

#include "scf.h" 
#include "geometry.h"
#include  "utilsa.h"
#include "genpot.h"
#include "param.h"
#include "grid.h"
#include "econfig.h"
#include "pcharge.h"
#include "eseqparam.h"
#include "elements.h"

const double eVtoau = 27.2214;
const double angstoev = 14.3997584;
const double autoangs =  0.529177;

tSCF scf; 
tGeometry geometry;
tGrid grid;
tConfigure config; 
tPointCharge pcharge; 
tESEQParams eseq; 
 
double qtot = 0.0; 
double *CMatrix; 
double *DMatrix; 
double *JJ; 
int *RecalcMatrix; 
int *AllcalcMatrix; 
double *BondH; 
double *BondE; 

const char *atsym[] = { "H ","He","Li","Be","B ","C ","N ","O ","F ", "Ne","Na","Mg","Al","Si","P ","S ","Cl","Ar", 
               "K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co", "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", 
               "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh", "Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", 
               "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu", "Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf", 
               "Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl", "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th", 
               "Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es", "Fm","Md","No","Lr", "O_N", "O_S", "O_SH", "H_O"}; 

/* Paramter chi in eV from RGIII/TableI (Rappe-Goddard the 3rd)	 
 modules.f90 of gulp except for Cu (4.200/4.220 SHOLL)*/        
double qeqchi[MAXELEMENTS] = {4.528,  9.660, 3.006, 4.877,5.110,5.343,6.899,8.741,10.874,11.04, 
                                 2.843,  3.951, 4.060, 4.168,5.463,6.928,8.564,9.465,2.421,3.231,
				 3.395,  3.470, 3.650, 3.415,3.325,3.760,4.105,4.465, 4.200,5.106,
				 3.641,  4.051, 5.188, 6.428,7.790,8.505,2.331,3.024,3.830,3.400, 
                                 3.550,  3.465, 3.290, 3.575,3.975,4.320,4.436,5.034,3.506,3.987,
				 4.899,  5.816, 6.822, 7.595,2.183,2.814,2.8355,2.744,2.858,2.8685,
				 2.881,2.9115,2.8785,3.1665,3.018,3.0555,3.127,3.1865,3.2514,3.2889, 
				 2.9629,3.700,5.100,  4.630,3.960,5.140, 5.000,4.790,4.894,6.270,
				 3.200,3.900, 4.690,4.210,4.750,5.370,2.00,2.843,2.835,3.175, 
                                 2.985,3.341,3.549,3.243,2.9895,2.8315, 3.1935,3.197,3.333,3.400,
				 3.470,3.475, 3.500,
				  8.741, 8.741, 8.741, 4.528};       

/*
These look like 0.5J in eV from RGIII/TableI 
and, therefore, should be multiplied by 2 */
	  
double qeqmui[MAXELEMENTS] = {6.9452,14.92,2.386,4.443, 4.750,5.063,5.880,6.682,7.474,10.55, 
                                   2.296,3.693,3.590,3.487,4.000,4.486, 4.946,6.355,1.920,2.880,
				   3.080,3.380,3.410,3.865,4.105,4.140,4.175,4.205, 4.220,4.285,
				   3.160,3.438,3.809,4.131,4.425,5.715,1.846,2.440,2.810,3.550, 
                                   3.380,3.755,3.990,4.015,4.005,4.000,3.134,3.957,2.896,3.124,
				   3.342,3.526,3.762,4.975,1.711,2.396,2.7415,2.692,2.564,2.6205,
				   2.673,2.7195,2.7875,2.9745,2.834,2.8715,2.891,2.9145,2.9329,2.965, 
				   2.4629,3.400,2.850,3.310,3.920,3.630,4.000,4.430,2.586,4.160,
				   2.900,3.530, 3.740,4.210,4.750,5.370,2.000,2.434,2.835,2.905, 
                                   2.905,2.853,2.717,2.819,3.0035,3.1895, 3.0355,3.101,3.089,3.10,
				   3.110,3.175, 3.200,
				   6.682, 6.682, 6.682, 6.9452};

/* These are R "sizes" from RGIII/Table I in A */ 
double qeqrad[MAXELEMENTS]  = {0.371,1.300,1.557,1.240, 0.822,0.759,0.715,0.669,0.706,1.768, 
                                  2.085,1.500,1.201,1.176,1.102,1.047,0.994,2.108,2.586,2.000,
				  1.750,1.607,1.470,1.402,1.533,1.393,1.406,1.398, 1.434,1.400,
				  1.211,1.189,1.204,1.224, 1.141,2.270,2.770,2.415,1.998,1.758, 
                                  1.603,1.530,1.500,1.500,1.509,1.544, 1.622,1.600,1.404,1.354,
				  1.404,1.380,1.333,2.459,2.984,2.442,2.071,1.925, 2.007,2.007,
				  2.000,1.978,2.227,1.968, 1.954,1.934,1.925,1.915,2.000,2.158, 
                                  1.896,1.759,1.605,1.538,1.600,1.700, 1.866,1.557,1.618,1.600,
				  1.530,1.444,1.514,1.480,1.470,2.200, 2.30,2.200,2.108,2.018, 
                                  1.800,1.713,1.800,1.840,1.942,1.900, 1.900,1.900,1.900,1.900,
				  1.900,1.900, 1.900, 
				  0.669, 0.669, 0.669, 0.371};

				   	     
int qeqids[MAXELEMENTS] = {1,  2,   3,   4,  5,   6,   7,   8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                            21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 
					    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 
					    61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 
					    81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
					    101, 102, 103,
					    800, 801, 802, 1001 }; 

/* used to identify principle quantum number */ 
int qeqbasetype[MAXELEMENTS] = {1,  2,   3,   4,  5,   6,   7,   8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                            21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 
					    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 
					    61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 
					    81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
					    101, 102, 103,
					    8, 8, 8, 1}; 

int main(int argc, char *argv[])
{
      FILE *outdat;
      FILE *outcif;
      FILE *outxyz;
      FILE *outen;
      int i, iatom, jatom, ibond, status; 
      int iatom_index, jatom_index; 
      double zetah0;
      double rjfac; 
      double qdiff; 
      double *pot, *func;
      double *atom_charge_diff, *charges;  
      double chargesum, pot_sum, pot_fitness, func_sum; 
      double fenergy, pcenergy; 
      double func_min, func_max; 
 
      gsl_matrix_view CMatrixView; 
      gsl_vector_view DVectorView; 
      gsl_vector *QVector; 
      gsl_permutation * p; 
      printf("egulp version 2.0.2\n");
      printf("number of command line parameters is %d\n", argc);
      if(argc != 4) { 
            printf("one of the input files is missing...Run as\n"); 
	    printf("egulp <geometry> <parameters> <configure>\n"); 
	    exit(-1); 
      } else { 
      	    printf("input parameters are assumed in order\n"); 
            printf("egulp <geometry> <parameters> <configure>\n"); 
      } 

      InitGeometry(&geometry, &(argv[1][0]) );
      SetBaseAtomType(&(geometry.atoms_num[0]),  &(geometry.atoms_base_num[0]),  
              geometry.natoms,  qeqids, qeqbasetype, MAXELEMENTS); 
      MapAtomsIntoCentralCell(&geometry);
      /* configure ... */  
      InitEconfigDefaults(&config); 
      InitEconfig(&config, &(argv[3][0])); 
      if(config.imethod == 0) { 
            InitParameters(qeqchi,  qeqmui, qeqids, MAXELEMENTS, &(argv[2][0]));
            PrintParameters(qeqchi,  qeqmui, qeqrad, qeqids, MAXELEMENTS); 
      } else if (config.imethod == 1) { 
            InitBonds(&geometry);  
	    InitESEQParams(qeqchi,  qeqmui, qeqids, MAXELEMENTS,  &(argv[2][0]), &eseq);
	    if(eseq.number_of_bond_types < geometry.tbonds) { 
	            printf("number of input bond types %d is smaller than %d\n", eseq.number_of_bond_types, geometry.tbonds); 
		    exit(-1); 
	    } 
	    
            printf("allocating bond hardness matrix and electronegativity corrections...\n"); 
            allocate_double_array(&(BondH), (geometry.natoms)*(geometry.natoms), "can not allocate BondH\n"); 
            allocate_double_array(&(BondE), (geometry.natoms), "can not allocate BondE\n"); 
            printf("allocating and constructing transfer matrices...\n"); 
            InitTransferMatrices(&eseq, geometry.natoms, geometry.nbonds,  geometry.nbonds_indexij);
            printf("transfer matrix...\n"); 
            for(iatom = 0; iatom < geometry.natoms; iatom++) 
                  for(ibond = 0; ibond < geometry.nbonds; ibond++) 
	                 printf("%3d %3d %10.6f\n", iatom, ibond, eseq.TMatrix[iatom*geometry.nbonds+ibond]); 
            printf("inverse transfer matrix...\n"); 
            for(ibond = 0; ibond < geometry.nbonds; ibond++) 
                  for(iatom = 0; iatom < geometry.natoms; iatom++) 
	                 printf("%3d %3d %10.6f\n", ibond, iatom, eseq.TMatrixPINV[ibond*geometry.natoms+iatom]); 
	     BuildParamsGEOOrder(&eseq,  geometry.tbonds, geometry.tbonds_indexij); 
	    /* HARDNESS BOND CORRECTIONS = (T^{-1})^T J' T^{-1}) */
             BuildBondHCorrections(&eseq, geometry.natoms, geometry.nbonds,  geometry.nbonds_indexij, 
                   geometry.nbonds_type, geometry.tbonds, geometry.tbonds_indexij, BondH); 
            /* ELECTRONEGATIVITY BOND CORRECTIONS E = T*c */
	     BuildBondECorrections(&eseq, geometry.natoms, geometry.nbonds,  geometry.nbonds_indexij, 
                   geometry.nbonds_type, geometry.tbonds, geometry.tbonds_indexij, BondE); 

      } else { 
      	    printf("invalid charge method...\n");
	    printf("use 0 for QEQ and 1 for SPLITQ\n"); 
	    exit(-1);  
      }
      
      if( (config.build_grid == 1) && (config.build_grid_from_scratch == 1)   ) 
            InitGrid(&grid, config.dh1sz, config.dh2sz, config.dh3sz, 
	    config.vdw_factor_i, config.vdw_factor_f, config.use_vdw_factor, config.offset, geometry.av1, geometry.av2, geometry.av3,  
	          geometry.natoms, geometry.atoms_base_num, geometry.atoms_xyz); 
      else if ( (config.build_grid == 1) && (config.build_grid_from_scratch == 0) ) 
            InitGridFromCube(config.input_grid_file, &grid, geometry.natoms);
      else printf("grid will not be constructed\n");  	    
       
       if (  (config.build_grid == 1) &&  ((config.build_grid_from_scratch == 1) || (config.build_grid_from_scratch == 0)) && (config.save_grid == 1)   ) 
            WriteGridAsCube(config.output_grid_file,  grid, geometry.natoms,  geometry.atoms_base_num,  geometry.atoms_xyz);  

       if(config.skip_everything == 1) { 
             if(config.build_grid == 1) DestroyGrid(&grid);  
             DestroyGeometry(&geometry);
	     if(config.imethod == 1) { 
	          DestroyBonds(&geometry); 
		  free(BondH); 
		  free(BondE); 
                  DestroyTransferMatrices(&eseq);
                  FreeESEQParams(&eseq); 
	     } 	  
             DestroyEconfig(config); 
	     return(0); 
       }	    
      
      InitEGULPScf(&scf, geometry.atoms_charge, geometry.atoms_num, geometry.natoms); 
      PrintEGULPScf(&scf);  

      if(config.point_charges_present == 1) { 
             InitPointCharge(&pcharge, geometry.av1, geometry.av2, geometry.av3,  geometry.natoms,  geometry.atoms_xyz, "pointcharge.dat"  ); 
	     genpot_on_points(geometry.av1, geometry.av2, geometry.av3,  pcharge.npoints, pcharge.points_charge,  pcharge.points_xyz, 
                geometry.natoms, pcharge.atatom_pot,  geometry.atoms_xyz,  geometry.natoms); 

      }
      printf("allocating ...\n"); 
      allocate_double_array(&(CMatrix), (geometry.natoms+1)*(geometry.natoms+1), "can not allocate CMatrix\n"); 
      allocate_double_array(&(JJ), (geometry.natoms)*(geometry.natoms), "can not allocate JJ\n");    
      allocate_double_array(&(DMatrix), (geometry.natoms+1), "can not allocate DMatrix\n");      
      allocate_int_array(&(RecalcMatrix), (geometry.natoms)*(geometry.natoms), "can not allocate RecalcMatrix\n");          
      allocate_int_array(&(AllcalcMatrix), (geometry.natoms)*(geometry.natoms), "can not allocate AllcalcMatrix\n");          
      printf("allocating...done\n"); 
      
      printf("allocating GSL objects...\n"); 
      QVector = gsl_vector_alloc (geometry.natoms+1 ); 
      p = gsl_permutation_alloc (geometry.natoms+1); 
      printf("allocating...done\n"); 

      for(iatom = 0; iatom < (geometry.natoms+1); iatom++) { 
            DMatrix[iatom] = 0.0; 
	    for(jatom = 0; jatom < (geometry.natoms+1); jatom++) 
	         CMatrix[iatom*(geometry.natoms+1)+jatom] = 0.0; 
      } 	    
      
      for(iatom = 0; iatom < geometry.natoms; iatom++)  { 
          for(jatom = 0; jatom < geometry.natoms; jatom++)  { 
	       JJ[iatom*geometry.natoms+jatom] = 0.0; 
	       if((geometry.atoms_num[iatom] == 1) || (geometry.atoms_num[jatom]==1)) RecalcMatrix[iatom*geometry.natoms+jatom] = 1; 
	       else RecalcMatrix[iatom*geometry.natoms+jatom] = 0; 
	       AllcalcMatrix[iatom*geometry.natoms+jatom] = 1;
	  } 
      }	        
      if(scf.current_iter != 0) { 
           printf("current iteration outside the loop is not 0...\n"); 
	   exit(-1);  
      }
      
     
      while((scf.current_iter < scf.niter) && (scf.lconverged == 0)) { 
	       
	             
	     if(scf.current_iter == 0) 	   {   
	              genpot(geometry.av1, geometry.av2, geometry.av3,  geometry.natoms, 
	              geometry.atoms_num, scf.atoms_charge, geometry.atoms_xyz, JJ, qeqrad, qeqids, qeqbasetype, MAXELEMENTS,  AllcalcMatrix); 
		      for(iatom = 0; iatom < geometry.natoms; iatom++) { 
		      	      iatom_index = element_index(geometry.atoms_num[iatom], qeqids, MAXELEMENTS); 
		              scf.QEQ_ChiArray[iatom] = qeqchi[iatom_index];
		       }      
	     } else   	      
	              genpot(geometry.av1, geometry.av2, geometry.av3,  geometry.natoms, 
	              geometry.atoms_num, scf.atoms_charge, geometry.atoms_xyz, JJ, qeqrad, qeqids, qeqbasetype, MAXELEMENTS,  RecalcMatrix); 
	       
	       for(iatom = 0; iatom < geometry.natoms; iatom++) { 
		      iatom_index = element_index(geometry.atoms_num[iatom], qeqids, MAXELEMENTS); 
	              DMatrix[iatom] = (-1.0)*qeqchi[iatom_index]; 
		      if(config.imethod == 1) DMatrix[iatom] = DMatrix[iatom]+BondE[iatom]; 
		      if((config.point_charges_present == 1) && (config.include_pceq == 1)) 
		              DMatrix[iatom] = DMatrix[iatom] - pcharge.atatom_pot[iatom] ; 
                }
               /*electronegativity corrections go in here, add to DMatrix */                  
	       DMatrix[geometry.natoms] = qtot; 

	       for(iatom = 0; iatom < geometry.natoms; iatom++)  {
		      for (jatom  = 0; jatom < geometry.natoms; jatom++)  {
			    CMatrix[iatom*(geometry.natoms+1)+jatom] = JJ[iatom*geometry.natoms+jatom]; 
			    scf.QEQ_JArray[iatom*geometry.natoms+jatom] = JJ[iatom*geometry.natoms+jatom]; 
		      }
	       }
	       	
	       for(iatom = 0; iatom < geometry.natoms; iatom++)  { 
		      if(geometry.atoms_num[iatom] == 1) { 
			    zetah0 = 0.75/(qeqrad[0]/0.529177); 
			    rjfac = 1.0+scf.atoms_charge[iatom]/zetah0; 
		    	    CMatrix[iatom*(geometry.natoms+1)+iatom] = 
			    CMatrix[iatom*(geometry.natoms+1)+iatom] + 2.0*qeqmui[0]*rjfac; 
		       } else { 
		             iatom_index = element_index(geometry.atoms_num[iatom], qeqids, MAXELEMENTS); 
		             CMatrix[iatom*(geometry.natoms+1)+iatom] = 
			     CMatrix[iatom*(geometry.natoms+1)+iatom] + 2.0*qeqmui[iatom_index]; 
		       } 	
		       scf.QEQ_JArray[iatom*geometry.natoms+iatom] = CMatrix[iatom*(geometry.natoms+1)+iatom] ;
     	
		       CMatrix[iatom*(geometry.natoms+1)+geometry.natoms] = 1.0; 
		       CMatrix[geometry.natoms*(geometry.natoms+1)+iatom]= 1.0; 	
	      } 	
              CMatrix[geometry.natoms*(geometry.natoms+1)+geometry.natoms] = 0.0; 
	      if(config.imethod == 1) { 
	              /* ADDING BOND HARDNESS CORRECTIONS */ 
	              for(iatom = 0; iatom < geometry.natoms; iatom++)  { 
		      for (jatom  = 0; jatom < geometry.natoms; jatom++)  {  
		              CMatrix[iatom*(geometry.natoms+1)+jatom] = CMatrix[iatom*(geometry.natoms+1)+jatom] + BondH[iatom*geometry.natoms+jatom]; 
		      }  	      
	              }  
              }	      

	      CMatrixView = gsl_matrix_view_array (CMatrix, geometry.natoms+1, geometry.natoms+1); 
	      DVectorView = gsl_vector_view_array (DMatrix, geometry.natoms+1);
              gsl_linalg_LU_decomp (&CMatrixView.matrix, p, &status);
              gsl_linalg_LU_solve (&CMatrixView.matrix, p, &DVectorView.vector, QVector);	
	      
   	      /* QVector = CMatrix.I * DMatrix 
	   invert matrix ... */    
	      printf("iteration #%d\n", scf.current_iter+1);  
	      for(iatom  = 0; iatom < geometry.natoms; iatom++) 
	             scf.atoms_charge[iatom] = gsl_vector_get(QVector, iatom); 
              if(scf.literate == 1) { 
		  qdiff = 0.0;
		  for (iatom = 0; iatom < geometry.natoms; iatom++)
			qdiff = qdiff + fabs(scf.atoms_charge[iatom]-scf.atoms_charge_old[iatom]); 
		  qdiff = qdiff/geometry.natoms; 
		  if (qdiff < scf.qeqscfcrit) scf.lconverged = 1; 
		  else scf.lconverged = 0; 	  
		  printf("Cycle: %4d Qdiff: %12.8f\n", scf.current_iter + 1, qdiff); 
		  if (scf.lconverged == 0) {  
			for(iatom  = 0; iatom < geometry.natoms; iatom++) { 
				scf.atoms_charge[iatom] = scf.MixParameter *scf.atoms_charge[iatom] + 
				    (1.0-scf.MixParameter )*scf.atoms_charge_old[iatom]; 
				scf.atoms_charge_old[iatom] =  scf.atoms_charge[iatom];
		 	} 	
		  } 		
	      }
	      scf.current_iter = scf.current_iter + 1; 	
	      fenergy = scf.energy[scf.current_iter-1] = getQEQEnergy(geometry.natoms, scf.atoms_charge, scf.QEQ_JArray, scf.QEQ_ChiArray) ; 
	      if (config.point_charges_present == 1) pcenergy = getQEQPCEnergy(geometry.natoms, scf.atoms_charge, pcharge.atatom_pot); 
	      else pcenergy = 0.0; 
      }   
      if( (scf.lconverged == 1) || ((scf.lconverged == 0) && (scf.literate == 0))) { 
            printf("SCF CONVERGED\n"); 
      } else { 
            printf("SCF NOT CONVERGED\n"); 
      }
      printf("Final charges from Qeq:\n"); 
      for (iatom = 0; iatom < geometry.natoms; iatom++) { 
	        geometry.atoms_charge[iatom] = scf.atoms_charge[iatom];  
	        printf("%4d %4d %12.7f\n", iatom+1, geometry.atoms_num[iatom], geometry.atoms_charge[iatom]); 
      }	
      
      /* write charges to file */ 
      outdat = fopen("charges.dat", "w"); 
      if(outdat != NULL) {
              for (iatom = 0; iatom < geometry.natoms; iatom++) 
	      fprintf(outdat, "%4d %4d %12.7f\n", iatom+1, geometry.atoms_num[iatom], geometry.atoms_charge[iatom]); 
	      fclose(outdat); 
      }	
      /* write .cif file with charges*/
      outcif = fopen("charges.cif", "w"); 
      fprintf(outcif, "data_crystal\n");
      fprintf(outcif, "\n");
      fprintf(outcif, "_cell_length_a %12.7f\n", geometry.cell_len[0]);
      fprintf(outcif, "_cell_length_b %12.7f\n", geometry.cell_len[1]);
      fprintf(outcif, "_cell_length_c %12.7f\n", geometry.cell_len[2]);
      fprintf(outcif, "_cell_angle_alpha %12.7f\n", geometry.cell_ang[0]);
      fprintf(outcif, "_cell_angle_beta %12.7f\n", geometry.cell_ang[1]);
      fprintf(outcif, "_cell_angle_gamma %12.7f\n", geometry.cell_ang[2]); 
      fprintf(outcif, "\n"); 
      fprintf(outcif, "_symmetry_space_group_name_Hall 'P 1'\n"); 
      fprintf(outcif, "_symmetry_space_group_name_H-M  'P 1'\n");
      fprintf(outcif, "\n");
      fprintf(outcif, "loop_\n");
      fprintf(outcif, "_symmetry_equiv_pos_as_xyz\n");
      fprintf(outcif, " 'x,y,z' \n");
      fprintf(outcif, "\n");
      fprintf(outcif, "loop_\n");
      fprintf(outcif, "_atom_site_label\n");
      fprintf(outcif, "_atom_site_type_symbol\n");
      fprintf(outcif, "_atom_site_fract_x\n");
      fprintf(outcif, "_atom_site_fract_y\n");
      fprintf(outcif, "_atom_site_fract_z\n");
      fprintf(outcif, "_atom_site_charge\n");
      if(outcif != NULL) {
              for (iatom = 0; iatom < geometry.natoms; iatom++) 
              fprintf(outcif, "%s  %s  %12.7f %12.7f %12.7f %12.7f\n",  
                atsym[geometry.atoms_base_num[iatom]-1], 
                atsym[geometry.atoms_base_num[iatom]-1], 
                geometry.atoms_fract[3*iatom+0], 
                geometry.atoms_fract[3*iatom+1], 
                geometry.atoms_fract[3*iatom+2],   
                geometry.atoms_charge[iatom]); 
                fclose(outcif); 
      }
      /* write .xyz file with charges*/
      outxyz = fopen("charges.xyz", "w"); 
      fprintf(outxyz, "%d\n\n",  geometry.natoms);  
      if(outxyz != NULL) {
              for (iatom = 0; iatom < geometry.natoms; iatom++) 
	      fprintf(outxyz, "%s %12.7f %12.7f %12.7f %12.7f\n",  
                atsym[geometry.atoms_base_num[iatom]-1], 
	        geometry.atoms_xyz[3*iatom+0], 
		geometry.atoms_xyz[3*iatom+1], 
		geometry.atoms_xyz[3*iatom+2],   
	        geometry.atoms_charge[iatom]); 
	        fclose(outxyz); 
      }	
      outen = fopen("energy.dat", "w");
      if(outen != NULL) { 
      	    fprintf(outen, "%d\n", scf.current_iter); 
           for(i = 0; i < scf.current_iter; i++) fprintf(outen, "%3d %12.7f  (EV)\n",  i+1, scf.energy[i]);  
	   fprintf(outen, "FINAL CONFIGURATIONAL %12.7f  (EV)\n",  fenergy);  
	   fprintf(outen, "DUE TO POINT CHARGES %12.7f  (EV)\n", pcenergy); 
	   fprintf(outen, "TOTAL %12.7f (EV)\n", fenergy+pcenergy); 
	   fclose(outen); 
      }  
      if(config.point_charges_present == 1) { 
	    genpot_on_points(geometry.av1, geometry.av2, geometry.av3,  geometry.natoms, geometry.atoms_charge,  
	       geometry.atoms_xyz, pcharge.npoints, pcharge.atpoint_pot,  pcharge.points_xyz,  geometry.natoms); 
	    outen = fopen("pchout.dat", "w");
            if(outen != NULL) { 
      	         fprintf(outen, "%d\n", pcharge.npoints); 
                 for(i = 0; i < pcharge.npoints; i++) fprintf(outen, "%3d %12.7f  %12.7f\n",  i+1, pcharge.points_charge[i], pcharge.atpoint_pot[i]);  
	         fclose(outen); 
            } 		 
      }  
	       
      if( (config.build_grid ==  1) && (config.calculate_pot_diff == 1) ) { 
      	   printf("calculating potential difference due to charges: ouput-input\n");  
           allocate_double_array(&pot, grid.gridnl, "egulp: can not allocate pot\n");      
	   allocate_double_array(&atom_charge_diff, geometry.natoms, "egulp: can not allocate atom_charge_diff\n");
	   for(i = 0; i < grid.gridnl; i++) pot[i] = 0.0; 
	   for(iatom = 0; iatom < geometry.natoms; iatom++)  atom_charge_diff[iatom] = geometry.atoms_charge[iatom]-geometry.input_atoms_charge[iatom];  
	   chargesum = sum_doubles(atom_charge_diff, geometry.natoms); 
	   printf("total charge %10.7f\n",  chargesum); 
	   for(iatom = 0; iatom < geometry.natoms; iatom++) atom_charge_diff[iatom]=atom_charge_diff[iatom]-chargesum/geometry.natoms; 
	   printf("total charge %10.7f\n",   sum_doubles(atom_charge_diff, geometry.natoms) ); 

           genpot_on_grid(geometry.av1, geometry.av2, geometry.av3,  geometry.natoms,  atom_charge_diff,  
	                  geometry.atoms_xyz, grid.gridnl, pot,  grid.gridL1,  grid.gridL2,  grid.gridL3, grid.gridh1,  grid.gridh2,  grid.gridh3); 

	   printf("number of grid points %d\n", grid.gridnl); 
	   printf("volume element is %10.7f\n", grid.dv3); 
	   pot_sum = sum_doubles(pot, grid.gridnl); 
           for(i = 0; i < grid.gridnl; i++) pot[i] = pot[i] - pot_sum/grid.gridnl;
	   pot_fitness = 0.0; 
	   for(i = 0; i < grid.gridnl; i++) pot_fitness = pot_fitness + fabs(pot[i]); 
	   
	   printf("integral potential error per MOF %10.7f  (eV)\n",  pot_fitness*grid.dv3); 
	   printf("potential error per grid point      %10.7f  (eV)\n",  pot_fitness/grid.gridnl);
           outdat = fopen("potential.dat", "w"); 
           if(outdat != NULL) {
	       fprintf(outdat, "PER MOF   %12.7f EV\n",  pot_fitness*grid.dv3); 
	       fprintf(outdat, "PER POINT %12.7f EV\n",  pot_fitness/grid.gridnl); 
	       fclose(outdat); 
           }	
	   free(pot); 
	   free(atom_charge_diff);   
      }	   
      
      if( (config.build_grid ==  1) && ( (config.calculate_pot == 1) ||  (config.calculate_pot == 2) )  ) { 
	   if (config.calculate_pot == 1) printf("calculating potential due to input charges\n"); 
	   if (config.calculate_pot == 2) printf("calculating potential due to output charges\n"); 
           allocate_double_array(&func, grid.gridnl, "egulp: can not allocate func\n");      
	   allocate_double_array(&charges, geometry.natoms, "egulp: can not allocate charges\n");
	   for(i = 0; i < grid.gridnl; i++) func[i] = 0.0; 
	   if (config.calculate_pot == 1) for(iatom = 0; iatom < geometry.natoms; iatom++)  charges[iatom] = geometry.input_atoms_charge[iatom];  
	   if (config.calculate_pot == 2) for(iatom = 0; iatom < geometry.natoms; iatom++)  charges[iatom] = geometry.atoms_charge[iatom]; 

	   chargesum = sum_doubles(charges, geometry.natoms); 
	   printf("total charge %10.7f\n",  chargesum); 
	   for(iatom = 0; iatom < geometry.natoms; iatom++) charges[iatom]= charges[iatom]-chargesum/geometry.natoms; 
	   printf("total charge %10.7f\n",   sum_doubles(charges, geometry.natoms) ); 

           genpot_on_grid(geometry.av1, geometry.av2, geometry.av3,  geometry.natoms,  charges,  
	                  geometry.atoms_xyz, grid.gridnl, func,  grid.gridL1,  grid.gridL2,  grid.gridL3, grid.gridh1,  grid.gridh2,  grid.gridh3); 

	   printf("number of grid points %d\n", grid.gridnl); 
	   printf("volume element is %10.7f\n", grid.dv3); 
	   func_sum = sum_doubles(func, grid.gridnl); 
           for(i = 0; i < grid.gridnl; i++) func[i] = func[i] - func_sum/grid.gridnl;
           func_min = 1.0e10; 
	   func_max = -1.0e10; 
	   for(i = 0; i < grid.gridnl; i++) {
	   	if (func[i] > func_max) func_max = func[i]; 
		if (func[i] < func_min)  func_min = func[i]; 
	   } 

           outdat = fopen("range.dat", "w"); 
           if(outdat != NULL) {
	       fprintf(outdat, "calculate_pot %d min %10.6f max %10.6f\n",  config.calculate_pot, func_min, func_max ); 
	       fclose(outdat); 
           }	
	   /* printf("calling WriteFunctionAsCube...\n");  */ 
	   WriteFunctionAsCube(config.output_pot_file, grid, geometry.natoms, geometry.atoms_base_num, geometry.atoms_xyz, func);
	   printf("out of WriteFunctionAsCube...\n"); 

	   free(func); 
	   free(charges);   
      }	   
      
      if(config.build_grid == 1) { 
          printf("destroying grid....\n"); 
          DestroyGrid(&grid);  
      }  
      printf("destroying geometry....\n"); 
      DestroyGeometry(&geometry);
      if(config.imethod == 1) { 
	     DestroyBonds(&geometry); 
             free(BondH); 
             free(BondE); 
             DestroyTransferMatrices(&eseq);
             FreeESEQParams(&eseq); 
      } 	        
      DestroyEconfig(config); 
      DestroyEGULPSCF(&scf); 
      if(config.point_charges_present == 1) DestroyPointCharge(&pcharge); 
      
      free(CMatrix); 
      free(DMatrix); 
      free(JJ); 
      free(AllcalcMatrix); 
      free(RecalcMatrix); 
      
      gsl_permutation_free (p);
      gsl_vector_free (QVector);      
      
      return(0); 
}





