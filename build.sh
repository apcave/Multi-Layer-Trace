#!/bin/bash
source ./venv/bin/activate

cd build
make clean
make all
cd ..


python scripts/steel_water.py
#python scripts/test.py