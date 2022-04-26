
import argparse
from ...KurtGroup.Kurt import xyz
from ...KurtGroup.Kurt import chemical_information as ci

def generateDaltonInputFileText(XYZ: xyz.xyz_to, charge: int) -> str:
    """Makes the text for a Dalton input file
    Args:
        charge (int): The charge of the molecule
    Returns:
        str: Returns the file text
    """
    if not hasattr(XYZ, 'atoms'):
        print(f'''You need to process the xyz file using the class function processXYZ before generating the file:
Exiting program''')
        exit()

    XYZ.input_filename = XYZ.filename.replace('.xyz', '.mol')

    unique_atoms = set(XYZ.atoms[:,0])
    XYZ.filetext = f'''ATOMBASIS
./{XYZ.filename}
Atomtypes={len(unique_atoms)} Charge={charge} NoSymmetry Angstrom
'''

    for unique_atom in unique_atoms:
        count = list(XYZ.atoms[:,0]).count(unique_atom)
        if not XYZ.BSE:
            XYZ.filetext += f'  {ci.AtomicInformation(unique_atom).atomnr():.4f}     {count} Bas={XYZ.basis}\n'
        else:
            BasisSet = ci.BasisSet()
            try:
                basis_mol = BasisSet.AtomBasisSet(XYZ.program, XYZ.basis, unique_atom, SupressHeader=True)
            except RuntimeError:
                print(f'''Failed to get basis set from BSE. Please check the spelling, upper-/lowercase IS important
The problem may also be that the basis set does not exist for {unique_atom}''')
                exit()

            BlockTypes = ['s functions', 'p functions', 'd functions', 'f functions', 'g functions', 'h functions', 'i functions', 'j functions', 'k functions']
            blocks = 0
            for j in BlockTypes:
                if j in basis_mol:
                    blocks += 1
            Block = f'{blocks}'
            for j in BlockTypes[:blocks]:
                Block += f' {basis_mol.count(j)}'
            XYZ.filetext += f'  Charge={ci.AtomicInformation(unique_atom).atomnr():.4f}     Atoms={count}     Blocks={Block}\n'
            basis_mol = basis_mol.replace('H','').split('\n')[5:-2]
            basis_mol = [i for i in basis_mol if 'functions' not in i]
        for j, atom in enumerate(XYZ.atoms[:,0]):
            if atom == unique_atom:
                XYZ.filetext += f'{atom: <2} {XYZ.atoms[j,1]: >14.9} {XYZ.atoms[j,2]: >19.9} {XYZ.atoms[j,3]: >19.9}\n'
        if XYZ.BSE:
            for j in basis_mol:
                XYZ.filetext += f'{j}\n'
    return XYZ.filetext

def writeInputfile(XYZ: xyz.xyz_to):
    """Writes the file text to an input file with similar name as the xyz file
    """
    with open(XYZ.input_filename, 'w') as input_file:
        input_file.writelines(XYZ.filetext)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''A script to convert xyz files to mol files for DALTON''', epilog='''For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='.xyz file')
    parser.add_argument('--charge', default=[0], nargs=1, type=int, help='Include to specify charge - 0 if not included')
    parser.add_argument('--basis', default=['pc-1'], nargs=1, type=str, help='Include to specify basis set of the molecular atoms - pc-1 if not included')
    # parser.add_argument('--RIbasis', nargs=1, type=str, help='Include to specify basis set of the molecular atoms - RI-BASIS if not included')

    args = parser.parse_args()

    input_files = args.infile
    basis = args.basis[0]
    charge = args.charge[0]

    for input_file in input_files:
        A = xyz.xyz_to('Dalton', input_file)
        A.processXYZ()
        A.setBasis(basis)
        filetext = generateDaltonInputFileText(A, charge)
        writeInputfile(A)

if __name__ == '__main__':
    main()