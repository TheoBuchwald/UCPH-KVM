import numpy as np
import sys
from collections import Counter #For number of unique elements

def get_atm_number(i):
    dic = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92}
    return dic[i]

#Read Here!
#Basis
basis = sys.argv[2]
RIbasis = sys.argv[3]
charge = sys.argv[4]

#.xyzfile
molxyz = sys.argv[1]

namesmol = []
molx, moly, molz = [], [], []

#Convert label into atomic number

#Read in xyz coordinates
with open(molxyz, "r") as f:
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
lines_mol.append('./'+molxyz+'\n')
lines_mol.append('Hej Theo\n')
lines_mol.append('Atomtypes='+str(len(set(namesmol)))+f' Charge={charge} NoSymmetry Angstrom\n')

unique_indices = list(Counter(namesmol).keys())

#Write coordinates and basis set
for i in unique_indices:
    ind = unique_indices.index(i)
    count = list(Counter(namesmol).values())[ind]
    lines_mol.append(f'  {get_atm_number(i):.4f}     {count} Bas={basis} Aux={RIbasis}\n')
    for j in range(len(namesmol)):
        if namesmol[j] == i:
            lines_mol.append(''.join([namesmol[j].ljust(2),' ',f"{molx[j]:.9f}".rjust(14),' ', f"{moly[j]:.9f}".rjust(19), ' ',f"{molz[j]:.9f}".rjust(19) ,'\n']))

#Output is .mol
with open(molxyz[:-4] + '_' + basis + '.mol','w') as f:
    f.writelines(lines_mol)
