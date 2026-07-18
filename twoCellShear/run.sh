#!/usr/bin/env bash

source ../src_hemocell/loadHemoCell.sh

mpirun -np 2 ./twoCellShear config.xml 2>&1 | tee out.txt