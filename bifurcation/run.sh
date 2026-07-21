#!/usr/bin/env bash

source ../src_hemocell/loadHemoCell.sh

mpirun -np 16 ./bifurcation config.xml 2>&1 | tee out.txt
