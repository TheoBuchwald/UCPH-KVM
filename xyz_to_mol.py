import numpy as np
import sys
from collections import Counter #For number of unique elements
import dependencies.chemical_information as ci

#Read Here!
#Basis
basis = sys.argv[2]
RIbasis = sys.argv[3]
charge = sys.argv[4]

#.xyzfile
molfile = sys.argv[1]

namesmol = []
molx, moly, molz = [], [], []

#Convert label into atomic number

#Read in xyz coordinates
with open(molfile, "r") as f:
    lines = f.readlines()
    for i in range(2, len(lines)):
       x = lines[i].split()
       namesmol.append(x[0])
       molx.append(float(x[1]))
       moly.append(float(x[2]))
       molz.append(float(x[3]))

#Lines to be inserted in .mol file
lines_mol =[]

#Starting lines
lines_mol.append('ATOMBASIS\n')
lines_mol.append('./'+molfile+'\n')
lines_mol.append('Hej Theo\n')
lines_mol.append('Atomtypes='+str(len(set(namesmol)))+f' Charge={charge} NoSymmetry Angstrom\n')

unique_indices = list(Counter(namesmol).keys())

#Write coordinates and basis set
for i in unique_indices:
    ind = unique_indices.index(i)
    count = list(Counter(namesmol).values())[ind]
    lines_mol.append(f'  {ci.AtomicInformation(i).atomnr():.4f}     {count} Bas={basis} Aux={RIbasis}\n')
    for j, atom in enumerate(namesmol):
        if atom == i:
            lines_mol.append(''.join([atom.ljust(2),' ',f"{molx[j]:.9f}".rjust(14),' ', f"{moly[j]:.9f}".rjust(19), ' ',f"{molz[j]:.9f}".rjust(19) ,'\n']))

#Output is .mol
with open(molfile[:-4] + '_' + basis + '.mol','w') as f:
    f.writelines(lines_mol)
