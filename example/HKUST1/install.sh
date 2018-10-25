#!/bin/bash

echo '*** You need icc: this install.sh is working for deneb.epfl.ch ***'

module load intel
cd ../../src
make
mv egulp ..
