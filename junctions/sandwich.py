#Makes a molecule sandwiched between two hemispheres

import numpy as np
from collections import Counter #For number of unique elements
import argparse
from Kurt import structures as struct
from Kurt import chemical_information as ci

# ------------------------------------ INPUTS ------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''A script to make junction consisting of a molecule and nanoparticles

    To use the following must be given:
        xyz-file   atom   atom   diameter''', epilog='''For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='.xyz file')
    parser.add_argument('atom1', type=int, nargs=1, help='Atom 1 that should be aligned between the nanoparticles')
    parser.add_argument('atom2', type=int, nargs=1, help='Atom 2 that should be aligned between the nanoparticles')
    parser.add_argument('diameter', type=float, nargs=1, help='Diameter of the nanoparticle')

    CrystalGroup = parser.add_argument_group('Nanoparticles')
    CrystalGroup.add_argument('-au', action='store_true', help='Include to make gold nanoparticles')
    CrystalGroup.add_argument('-ag', action='store_true', help='Include to make silver nanoparticles')
    CrystalGroup.add_argument('-cu', action='store_true', help='Include to make copper nanoparticles')
    CrystalGroup.add_argument('-tio2', action='store_true', help='Include to make titanium dioxide nanoparticles')
    CrystalGroup.add_argument('-nacl', action='store_true', help='Include to make salt nanoparticles')
    CrystalGroup.add_argument('-pd', action='store_true', help='Include to make palladium nanoparticles')
    CrystalGroup.add_argument('-pt', action='store_true', help='Include to make platinum nanoparticles')
    CrystalGroup.add_argument('-cosb3', action='store_true', help='Include to make CoSb3 nanoparticles')

    CalculationGroup = parser.add_argument_group('Calculation options')
    CalculationGroup.add_argument('--charge', default=[0], nargs=1, type=int, help='Include to specify charge - 0 if not included')
    CalculationGroup.add_argument('--basis', default=['pc-1'], nargs=1, type=str, help='Include to specify basis set - pc-1 if not included')

    AdditionalCommandsGroup = parser.add_argument_group('Additional commands')
    AdditionalCommandsGroup.add_argument('--outwards', action='store_false', help='Include to turn the nanoparticles outwards')
    AdditionalCommandsGroup.add_argument('--noxyz', action='store_false', help='Include to TURN OFF creation of .xyz files of the nanoparticle system')
    AdditionalCommandsGroup.add_argument('-l', '--linenumber', action='store_true', help='Include to use the linenumber of the atoms in the xyz file instead')

    args = parser.parse_args()

    Arguments = {   #Only crystal structures
        'Au' : args.au,
        'Ag' : args.ag,
        'Cu' : args.cu,
        'Pt' : args.pt,
        'Pd' : args.pd,
        'TiO2' : args.tio2,
        'NaCl' : args.nacl,
        'CoSb3' : args.cosb3
    }

    if all(value == False for value in Arguments.values()):
        Arguments['Au'] = True

    input_files = args.infile
    atom1 = args.atom1[0]
    atom2 = args.atom2[0]
    diameter = args.diameter[0]

    charge = args.charge[0]
    basis = args.basis[0]

    inwards = args.outwards
    linenumber = args.linenumber
    returnxyz = args.noxyz

    if linenumber:
        atom1 -=3
        atom2 -=3
    else:
        atom1 -=1
        atom2 -=1

    for molfile in input_files:
        namesmol = np.array([])
        molxyz = struct.Molecule(np.empty((0, 3)))
        with open(molfile, 'r') as f:
            lines = f.readlines()
            for i in range(2, len(lines)):
                x = lines[i].split()
                namesmol = np.append(namesmol, x[0])
                molxyz.molecule = np.vstack([molxyz.molecule, np.array([float(x[1]), float(x[2]), float(x[3])])])

        #We get the axis between the two points to be parallel to the x-axis and then translate it such that it lies in the x-axis
        #Get normalized vector defining the axis
        v1 = molxyz.molecule[atom1, :] - molxyz.molecule[atom2, :]
        v1 *= 1 / np.sqrt(np.dot(v1, v1))

        #Calculate angle in order to be parallel to the x-axis
        theta = np.arccos(np.dot(v1, np.array([1, 0, 0])))

        #The direction vector for the rotation in orthogonal to both v1 and the x-axis. Note that the normalization is crucial.
        dir_vec = np.cross(v1, np.array([1, 0, 0]))
        dir_vec *= 1/np.sqrt(np.dot(dir_vec, dir_vec))

        molxyz.get_rotation_matrix(molxyz.molecule[atom1], dir_vec, theta)
        molxyz.rotateMolecule(atom1)

        molxyz.__len__()

        mol_min = molxyz.min()
        mol_max = molxyz.max()

        index_min = molxyz.index_min()
        index_max = molxyz.index_max()

        crystal_structures = [item[0] for item in Arguments.items() if item[1] == True]
        for crystal_structure in crystal_structures:
            NP = struct.NanoParticle(crystal_structure)

            NP.setInwards(inwards)
            NP.setDiameter(diameter)

            atoms_symbol, atoms_pos = NP.makeNanoparticle(diameter)

            left, right, left_symbols, right_symbols = NP.makeSandwich(molxyz, namesmol)

            lines_to_add = []
            if returnxyz:
                # Build .xyz file
                lines_to_add.append(str(left[:, 0].size+right[:, 0].size+molxyz.molecule[:, 0].size)+'\n')
                lines_to_add.append('\n')
                for i in range(left[:, 0].size):
                    lines_to_add.append(''.join([f"{atoms_symbol[i]}", ' ', f"{left[i, 0]:.6f}", ' ', f"{left[i, 1]:.6f}", ' ', f"{left[i, 2]:.6f}" , '\n']))
                for i in range(right[:, 0].size):
                    lines_to_add.append(''.join([f"{atoms_symbol[i]}", ' ', f"{right[i, 0]:.6f}", ' ', f"{right[i, 1]:.6f}", ' ', f"{right[i, 2]:.6f}" , '\n']))

                for i in range(molxyz.molecule[:, 0].size):
                    lines_to_add.append(''.join([namesmol[i], ' ', f"{molxyz.molecule[i, 0]:.6f}", ' ', f"{molxyz.molecule[i, 1]:.6f}", ' ', f"{molxyz.molecule[i, 2]:.6f}" , '\n']))
                with open(molfile[:-4] + f'_charge_{charge}_{crystal_structure}.xyz', 'w') as f:
                    f.writelines(lines_to_add)

            lines_mol =[]

            lines_mol.append('ATOMBASIS\n')
            lines_mol.append('./'+molfile+'\n')
            lines_mol.append(f'Hej - The distance between NPs is {NP.distance:.4f} AA\n')
            lines_mol.append('Atomtypes='+str(len(set(namesmol)))+f' Charge={charge} NoSymmetry Angstrom\n')

            unique_indices = list(Counter(namesmol).keys())

            for i in unique_indices:
                ind = unique_indices.index(i)
                count = list(Counter(namesmol).values())[ind]
                lines_mol.append(f'  {ci.AtomicInformation(i).getAtomnr():.4f}     {count} Bas={basis}\n')
                for j in range(len(namesmol)):
                    if namesmol[j] == i:
                        lines_mol.append(''.join([namesmol[j].ljust(2), ' ', f"{molxyz.molecule[j, 0]:.9f}".rjust(14), ' ', f"{molxyz.molecule[j, 1]:.9f}".rjust(19), ' ', f"{molxyz.molecule[j, 2]:.9f}".rjust(19) , '\n']))

            with open(molfile[:-4] + f'_charge_{charge}_{crystal_structure}.mol', 'w') as f:
                f.writelines(lines_mol)

            #Build .pol file
            lines_pol = []

            lines_pol.append('AA\n')
            lines_pol.append(str(left[:, 0].size*2)+'   0 1 1\n')

            n = 1

            for i in range(left[:, 0].size):
                if left_symbols[i] == 'O':
                    lines_pol.append(''.join([str(n), f" {left[i, 0]:.6f} ", f"{left[i, 1]:.6f} ", f"{left[i, 2]:.6f} ", "-2.000000 ", f"{ci.AtomicInformation(left_symbols[i]).polarizability():.6f}", '\n']))
                elif left_symbols[i] == 'Ti':
                    lines_pol.append(''.join([str(n), f" {left[i, 0]:.6f} ", f"{left[i, 1]:.6f} ", f"{left[i, 2]:.6f} ", "4.000000 ", f"{ci.AtomicInformation(left_symbols[i]).polarizability():.6f}", '\n']))
                else:
                    lines_pol.append(''.join([str(n), f" {left[i, 0]:.6f} ", f"{left[i, 1]:.6f} ", f"{left[i, 2]:.6f} ", "0.000000 ", f"{ci.AtomicInformation(left_symbols[i]).polarizability():.6f}", '\n']))
            n += 1

            for i in range(right[:, 0].size):
                if right_symbols[i] == 'O':
                    lines_pol.append(''.join([str(n), f" {right[i, 0]:.6f} ", f"{right[i, 1]:.6f} ", f"{right[i, 2]:.6f} ", "-2.000000 ", f"{ci.AtomicInformation(right_symbols[i]).polarizability():.6f}", '\n']))
                elif right_symbols[i] == 'Ti':
                    lines_pol.append(''.join([str(n), f" {right[i, 0]:.6f} ", f"{right[i, 1]:.6f} ", f"{right[i, 2]:.6f} ", "4.000000 ", f"{ci.AtomicInformation(right_symbols[i]).polarizability():.6f}", '\n']))
                else:
                    lines_pol.append(''.join([str(n), f" {right[i, 0]:.6f} ", f"{right[i, 1]:.6f} ", f"{right[i, 2]:.6f} ", "0.000000 ", f"{ci.AtomicInformation(right_symbols[i]).polarizability():.6f}", '\n']))
            n += 1

            with open(molfile[:-4]+f'_charge_{charge}_{crystal_structure}.pot', 'w') as f:
                f.writelines(lines_pol)

            print("Outputs:")
            if returnxyz:
                print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}.xyz')
            print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}.pot')
            print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}.mol')
