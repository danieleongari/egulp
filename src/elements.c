
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int  element_index(int id, int *qeqids, int maxel)
{
     int id_res;
     int i;  
     id_res = -1; 
     for(i = 0; i <  maxel; i++) { 
           if(id == qeqids[i]) id_res = i; 
	   else continue;  
     }    
     
     if(id_res == -1) { 
            printf("element %d is not found by element_index\n", id); 
	    exit(-1); 
     } 
     
     return(id_res); 
}

void SetBaseAtomType(int *atom_nums, int *atom_base_nums, int natoms, int *qeqids, int *qeqbasetype, int MaxEl)
{ 
      int iatom;  
      int iatom_index;  
      printf("Atomic numbers and base numbers\n"); 
      for (iatom = 0; iatom < natoms; iatom++) { 
            iatom_index = element_index(atom_nums[iatom], qeqids,  MaxEl);
	    atom_base_nums[iatom] = qeqbasetype[iatom_index]; 
	    printf("%4d %4d %4d\n", iatom+1, atom_nums[iatom], atom_base_nums[iatom]); 
      }
      return; 
}
	 
