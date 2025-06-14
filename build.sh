#!/bin/bash
source ./venv/bin/activate

cd build
make all
cd ..

python scripts/test.py