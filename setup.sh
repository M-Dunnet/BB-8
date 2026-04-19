#!/bin/bash

cwd=$(pwd)

# Resolve dependencies
echo 'Pip installables (scipy, numpy, mappy, edlib, networkx, pyabpoa, Cython)'
python3 -m pip install --upgrade pip
python3 -m pip install scipy numpy mappy edlib networkx pyabpoa Cython

echo 'conk'
python3 -m pip install wheel setuptools
git clone https://github.com/M-Dunnet/conk
cd conk && make
cd $cwd

echo 'blat'
mkdir blat
cd blat
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat
chmod +x blat
cd $cwd

echo 'Done'
