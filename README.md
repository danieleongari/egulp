# eGULP

 * Efficient Qeq algorithm for periodic systems. Original method by [Rappé et al.](http://pubs.acs.org/doi/abs/10.1021/j100161a070)

 * Consider using one of the following parameters: 
 - GMP.param (from UFF) 
 - MEPO.param from [Kadantsev et al.](http://pubs.acs.org/doi/10.1021/jz401479k)
 - GULP.param same as GMP except for Cu and Ce

 * The main modifications of 25 Oct 18, allows for reading cif files and printig charge.cif,
   but please consider that the parsing of the cif file occurr in the following way (and therefore not all cif formats are compatible):
   1) reads cell info from _cell_length_a, _cell_length_b, _cell_length_c, _cell_angle_alpha, _cell_angle_beta, _cell_angle_gamma
   2) only P1 symmetry is compatible 
   3) only fractional coordinates of atoms are compatible
   4) starts to read the atoms from the first line after _atom_site_fract_x that does not contains "_atom_site"
   5) for each atom line in the cif file, the program reads:
       column 1: ignored
       column 2: atom name
       column 3: x fractional coordinate
       column 4: y fractional coordinate
       column 5: z fractional coordinate
       column 6+: if present, ignored
