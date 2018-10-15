#ifndef _ELEMENTS_H_
#define _ELEMENTS_H_

#define MAXELEMENTS (107)

int  element_index(int id, int *qeqids, int maxel); 
void SetBaseAtomType(int *atom_nums, int *atom_base_nums, int natoms, int *qeqids, int *qeqbasetype, int MaxEl); 
#endif
