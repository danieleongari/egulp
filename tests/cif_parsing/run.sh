#!/bin/bash
# test that we can parse cif file in both formats
set -e

egulp manage_crystal.cif GMP.param configure.input
egulp pycifrw.cif GMP.param configure.input
