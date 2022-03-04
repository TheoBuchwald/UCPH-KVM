#Make Gaussian style input file from xyz file
#Use the BSE API to obtain basis sets not already in Gaussian

import numpy as np
import argparse
from Kurt import structures as struct
from Kurt import chemical_information as ci
import fnmatch as fn

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Use to make Gaussian input files from a given .xyz file.
Basis sets not implemented in Gaussian are imported from BSE: https://www.basissetexchange.org
NB! When specifying Pople-style basis sets with polarization functions, parentheses must be escaped by \\
Stars can still be used, however.''', epilog='''For help contact
        Theo Juncker von Buchwald
        fnc970@alumni.ku.dk
        Magnus Bukhave Johansen
        qhw298@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs=1, help='The file(s) to convert into .com files', metavar='.xyz file')
    parser.add_argument('calc', type=str,nargs=1,help='Keywords for the calculation. If there are spaces, have quotes around the entire thing.')

    CalculationGroup = parser.add_argument_group('Calculation options')
    CalculationGroup.add_argument('--charge', default=[0], nargs=1, type=int, help='Include to specify charge - 0 if not included')
    CalculationGroup.add_argument('--mult', default=[1], nargs=1, type=int, help='Include to specify multiplicity - 1 if not included')
    CalculationGroup.add_argument('--basis', default=['pc-1'], nargs=1, type=str, help='Include to specify basis set - pc-1 if not included')
    CalculationGroup.add_argument('--method', default=['cam-b3lyp'], nargs=1, type=str, help='Include to specify method for calculation - CAM-B3LYP if not included')
    CalculationGroup.add_argument('--cpu', default=[8], nargs=1, type=int, help='Include to specify the amount of cpu cores - 8 if not included')
    CalculationGroup.add_argument('--mem', default=[8], nargs=1, type=int, help='Include to specify the amount of memory in GB - 8 if not included')
    CalculationGroup.add_argument('--chk',action='store_false', help='Include to insert a checkpoint file with the same name as the .xyz file')

    args = parser.parse_args()

    infile = args.infile[0]
    calc = args.calc[0]

    charge = args.charge[0]
    multiplicity = args.mult[0]
    basis = args.basis[0]
    method = args.method[0]
    ncpus = args.cpu[0]
    mem = args.mem[0]
    chk = args.chk


    #Get molecular information
    namesmol = np.array([])
    molxyz = struct.Molecule(np.empty((0, 3)))
    with open(infile, 'r') as f:
        lines = f.readlines()
        for i in range(2, len(lines)):
            x = lines[i].split()
            namesmol = np.append(namesmol, x[0])
            molxyz.molecule = np.vstack([molxyz.molecule, np.array([float(x[1]), float(x[2]), float(x[3])])])

    #Check if basis set is in Gaussian already (Only checks the most common)
    basis_sets_gauss = ["6-31*[Gg]*","*[Cc][Cc]-[Pp][Vv]*","[Ss][Tt][Oo]-3[Gg]","3-21*[Gg]*", "6-21*[Gg]*", "4-31*[Gg]*",  "[Ll][Aa][Nn][Ll]2*"]
    BSE = True

    for i in basis_sets_gauss:
        if fn.filter([basis],i):
            BSE = False
            break

    if BSE:
        #Use the BSE API to get the basis set
        print("Need to get basis set from BSE!")
        BasisSet = ci.BasisSet()
        try:
            basis_mol = BasisSet.GenerateBasisSet('gaussian94', basis, namesmol, SupressHeader=True)
        except RuntimeError:
            print("Failed to get basis set from BSE. Please check the spelling, upper-/lowercase is not important")
            exit()
        basis_name = basis
        basis = "GEN"

    #Build .com file
    filename = infile[:-4]
    lines_to_add = []
    if not chk:
        lines_to_add.append(f"%chk={filename}.chk\n")
    lines_to_add.append(f"%mem={mem}GB\n")
    lines_to_add.append(f"%nprocshared={ncpus}\n")
    lines_to_add.append(f"# {calc} {method}/{basis} \n")
    lines_to_add.append('\n')
    if BSE:
        lines_to_add.append(f"Hej Magnus - {filename}.xyz - {basis_name}\n")
    else:
        lines_to_add.append(f"Hej Magnus - {filename}.xyz\n")
    lines_to_add.append('\n')
    lines_to_add.append(f"{charge} {multiplicity}\n")
    for i in range(len(namesmol)):
        lines_to_add.append(''.join([namesmol[i],' ',f"{molxyz.molecule[i,0]:.9f}".rjust(14),' ', f"{molxyz.molecule[i,1]:.9f}".rjust(19), ' ',f"{molxyz.molecule[i,2]:.9f}".rjust(19) ,'\n']))
    lines_to_add.append('\n')

    if BSE:
        lines_to_add.append(basis_mol+'\n')

    with open(filename+'.com','w') as f:
        f.writelines(lines_to_add)

    print(filename+'.com')

