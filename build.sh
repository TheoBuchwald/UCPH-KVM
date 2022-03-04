#!/bin/bash

pip install numpy --user
pip install requests --user
pip install matplotlib --user
python3 setup.py install --user

# Checking for ase in .bashrc:
if grep -Fq "#module load ase" ~/.bashrc
then
 echo 'Ase is not loaded: This will need to be loaded'
elif grep -Fq "module load ase" ~/.bashrc
 then
  :
else:
 echo 'module load ase/3.19.0' >> ~/.bashrc
 echo 'Ase was not loaded or in you bashrc. It has now been added, please restart you terminal'
fi
