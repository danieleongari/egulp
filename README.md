# eGULP

Efficient Qeq algorithm for periodic systems. Original method by [Rappé et al.](http://pubs.acs.org/doi/abs/10.1021/j100161a070)

Credit to [Kadantsev et al.](http://pubs.acs.org/doi/10.1021/jz401479k)

## Installation

Prerequisites:
 * C compiler (either icc or gcc)
 * GNU scientific library (GSL)
   for ubuntu: `sudo apt-get install libgsl0-dev`
   for macos: `sudo port install gsl`


```bash
cd src
# edit Makefile to select compiler (gcc.arch or intel.arch)
make
```

## Usage
```
egulp HKUST.cif GMP.param configure.input
```
prints files `charge.cif`,`charge.dat`,`charge.xyz` and `energy.dat`.

 * Consider using one of the following parameters:
   1) [GMP.param](data/GMP.param) (from UFF as implemented in [Openbabel](https://github.com/openbabel/openbabel/blob/master/data/qeq.txt))
   2) [MEPO.param](data/MEPO.param) from [Kadantsev et al.](http://pubs.acs.org/doi/10.1021/jz401479k)
   3) [GULP.param](data/GULP.param) same as GMP except for Cu and Ce

 * This version of the code allows reading cif files and printig charge.cif, but please consider that the parsing of the cif file is done in the following way (and therefore not all cif files are compatible):
   1) it reads cell info from _cell_length_a, _b, _c, _cell_angle_alpha, _beta, _gamma
   2) only P1 symmetry is compatible
   3) only fractional coordinates of atoms are compatible
   4) it starts to read the atoms from the first line after _atom_site_fract_x that does not contains "_atom_site"
   5) for each atom line in the cif file, it reads:

     | column |    read                  |
     |:------:|:------------------------ |
     |    1   | ignored                  |
     |    2   | atom name                |
     |    3   | x fractional coordinate  |
     |    4   | y fractional coordinate  |
     |    5   | z fractional coordinate  |
     |    6+  | if present, ignored      |

 * File `configure.input` allows to specify extra settings:
 ```
build_grid 0
build_grid_from_scratch 1 none 0.25 0.25 0.25 1.0 2.0 0 0.3
save_grid 0 grid.cube
calculate_pot_diff 0
calculate_pot 0 repeat.cube
skip_everything 0
point_charges_present 0
include_pceq 0
imethod 0
 ```
which stand for:
 ```
build_grid ............ {bool}
build_grid_from_scratch {bool}
                         {grid_input name}
                         {dx, dy and dz grid spacing}
                         {vdw_min vdw_max}
                         {bool, consider only gridpoint between vdw_min and vdw_max}
save_grid ............. {bool}
                         {cube filename}
calculate_pot_diff .... {bool}
calculate_pot ......... {0:False, 1: Input pot, 2: output pot}
                         {cube filename}
skip_everything ....... {bool, to skip charge and pot calculation}
point_charges_present . {bool, to start the method from initial charges}
include_pceq .......... {bool}
imethod ............... {0: Qeq, 1: SPLIT-Qeq}
 ```

## Development

Before making changes to the code, please install the pre-commit hooks (for automatic code formatting and linting):
```
pip install pre-commit
pre-commit install
```
