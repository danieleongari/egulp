#!/bin/bash
rm -rf output/
mkdir output
../../egulp HKUST1.cif GMP.param configure.input > output/egulp.log
mv charges.cif charges.dat charges.xyz energy.dat output/
