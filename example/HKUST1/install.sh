#!/bin/bash

echo '*** You need icc: this install.sh is working for deneb.epfl.ch'

module load intel
cd ../../1_original_download
make
rm *o
mv egulppot ..
