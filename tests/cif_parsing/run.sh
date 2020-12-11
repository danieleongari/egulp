#!/bin/bash
# Test parsing a couple of CIF files
set -e

for cif in *.cif; do
    echo "## Running $cif"
    ../../src/egulp $cif GMP.param configure.input > $cif.log
    echo "## Success $cif"
done
