#Makes a molecule sandwiched between two hemispheres
#Then makes Gaussian16 input files for both

import numpy as np
import argparse
import dependencies.structures as struct
import dependencies.chemical_information as ci
import fnmatch as fn

#INPUTS HERE
#------------------------------

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
    CalculationGroup.add_argument('--basis', default=['pc-1'], nargs=1, type=str, help='Include to specify basis set of the molecular atoms - pc-1 if not included')
    CalculationGroup.add_argument('--NPbasis', default=['LANL2DZ'], nargs=1, type=str, help='Include to specify basis set of the nanoparticle atoms - LANL2DZ if not included')
    CalculationGroup.add_argument('--ECPbasis', default=['LANL2'], nargs=1, type=str, help='Include to specify electronic core potential basis set for nanoparticle atoms - LANL2 if not included')
    CalculationGroup.add_argument('--method', default=['cam-b3lyp'], nargs=1, type=str, help='Include to specify method for calculation - CAM-B3LYP if not included')
    CalculationGroup.add_argument('--cpu', default=[16], nargs=1, type=int, help='Include to specify the amount of cpu cores - 16 if not included')
    CalculationGroup.add_argument('--mem', default=[16], nargs=1, type=int, help='Include to specify the amount of memory in GB - 16 if not included')

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

    input_files = args.infile
    atom1 = args.atom1[0]
    atom2 = args.atom2[0]
    diameter = args.diameter[0]

    charge = args.charge[0]
    basis_mol = args.basis[0]
    basis_NP = args.NPbasis[0]
    basis_ECP = args.ECPbasis[0]
    method = args.method[0]
    ncpus=args.cpu[0]
    mem=args.mem[0]

    inwards = args.outwards
    returnxyz = args.noxyz
    linenumber = args.linenumber

    if linenumber:
        atom1 -=3
        atom2 -=3
    else:
        atom1 -=1
        atom2 -=1

    #Check if basis set is in Gaussian already (Only checks the most common)
    basis_sets_gauss = ["6-31*[Gg]*","*[Cc][Cc]-[Pp][Vv]*","[Ss][Tt][Oo]-3[Gg]","3-21*[Gg]*", "6-21*[Gg]*", "4-31*[Gg]*",  "[Ll][Aa][Nn][Ll]2*"]
    BSE_mol = True
    BSE_NP = True
    BSE_ECP = True

    for i in basis_sets_gauss:
        if fn.filter([basis_mol],i):
            BSE_mol = False
            break
    for i in basis_sets_gauss:
        if fn.filter([basis_NP],i):
            BSE_NP = False
            break
    for i in basis_sets_gauss:
        if fn.filter([basis_ECP],i):
            BSE_ECP = False
            break

    # Getting basis set from Basis set exchange
    BasisSet = ci.BasisSet()

    for molfile in input_files:
        namesmol = np.array([])
        molxyz = struct.Molecule(np.empty((0, 3)))
        with open(molfile, 'r') as f:
            lines = f.readlines()
            for i in range(2, len(lines)):
                x = lines[i].split()
                namesmol = np.append(namesmol, x[0])
                molxyz.molecule = np.vstack([molxyz.molecule, np.array([float(x[1]), float(x[2]), float(x[3])])])

        if BSE_mol:
            basis_mol = BasisSet.GenerateBasisSet('gaussian94', basis_mol, namesmol)


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

        #Rescaling to get index in array

        #Crystal Structure (Au, Ag, Cu, TiO2, NaCl, CoSb3, Pt, Pd)
        crystal_structures = [item[0] for item in Arguments.items() if item[1] == True]
        for crystal_structure in crystal_structures:
            NP = struct.NanoParticle(crystal_structure)

            NP.setInwards(inwards)
            NP.setDiameter(diameter)

            atoms_symbol, atoms_pos = NP.makeNanoparticle(diameter)

            left, right, left_symbols, right_symbols = NP.makeSandwich(molxyz, namesmol)

            atmtype = NP.atomtypes

            if BSE_NP:
                basis_NP = BasisSet.GenerateBasisSet('gaussian94', basis_NP, atmtype)
            if BSE_ECP:
                basis_ECP = BasisSet.GenerateBasisSet('gaussian94', basis_ECP, atmtype)

            if returnxyz:
                #Build .xyz files
                for j in ['left','right']:
                    lines_to_add = []
                    lines_to_add.append(str(left[:,0].size+molxyz.molecule[:,0].size)+'\n')
                    lines_to_add.append('\n')

                    for i in range(len(left[:,0])):
                        if j == 'left':
                            lines_to_add.append(''.join([atoms_symbol[i],' ',f"{left[i,0]:.6f}",' ', f"{left[i,1]:.6f}", ' ',f"{left[i,2]:.6f}" ,'\n']))
                        else:
                            lines_to_add.append(''.join([atoms_symbol[i],' ',f"{right[i,0]:.6f}",' ', f"{right[i,1]:.6f}", ' ',f"{right[i,2]:.6f}" ,'\n']))

                    for i in range(len(molxyz.molecule[:,0])):
                        lines_to_add.append(''.join([namesmol[i],' ',f"{molxyz.molecule[i,0]:.6f}",' ', f"{molxyz.molecule[i,1]:.6f}", ' ',f"{molxyz.molecule[i,2]:.6f}" ,'\n']))

                    with open(molfile[:-4] + f'_charge_{charge}_{crystal_structure}_'+j+'.xyz','w') as f:
                        f.writelines(lines_to_add)


            for j in ['left','right']:
                filename = molfile[:-4] + f'_charge_{charge}_{crystal_structure}_'+j
                lines_com = []
                lines_com.append(f"%mem={mem}GB\n")
                lines_com.append(f"%nprocshared={ncpus}\n")
                lines_com.append(f"%chk={filename}.chk\n")
                lines_com.append(f"# {method}/GEN PSEUDO=READ scf=qc pop=full iop(3/33=1)\n")
                lines_com.append('\n')
                lines_com.append(f"Hej Theo - {molfile}-{j}\n")
                lines_com.append('\n')
                if charge == 0 or abs(charge) == 2:
                    lines_com.append(f"{charge} 1 {charge} 1 0 1\n")
                else:
                    lines_com.append(f"{charge} 2 {charge} 2 0 1\n")
                for i in range(len(namesmol)):
                    lines_com.append(''.join([namesmol[i],'(Fragment=1)',' ',f"{molxyz.molecule[i,0]:.9f}",' ', f"{molxyz.molecule[i,1]:.9f}", ' ',f"{molxyz.molecule[i,2]:.9f}" ,'\n']))
                if j == 'left':
                    for i in range(len(left[:,0])):
                        lines_com.append(''.join([atoms_symbol[i],'(Fragment=2)',' ',f"{left[i,0]:.9f}",' ', f"{left[i,1]:.9f}", ' ',f"{left[i,2]:.9f}" ,'\n']))
                else:
                    for i in range(len(right[:,0])):
                        lines_com.append(''.join([atoms_symbol[i],'(Fragment=2)',' ',f"{right[i,0]:.9f}",' ', f"{right[i,1]:.9f}", ' ',f"{right[i,2]:.9f}" ,'\n']))
                lines_com.append('\n')
                for i in atmtype:
                    if i in set(namesmol):
                        pass
                    else:
                        lines_com.append(i+'\n')
                        lines_com.append(basis_NP+'\n')
                        lines_com.append('****'+'\n')
                lines_com.append(basis_mol+'\n')
                lines_com.append('\n')
                for i in atmtype:
                    lines_com.append(i+' 0\n')
                    lines_com.append(basis_ECP+'\n')
                    lines_com.append('\n')
                lines_com.append('\n')
                with open(filename+'.com','w') as f:
                    f.writelines(lines_com)

            if returnxyz:
                print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}_left.xyz')
                print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}_right.xyz')
            print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}_left.com')
            print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}_right.com')
