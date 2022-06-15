#Make Gaussian style input file from xyz file
#Use the BSE API to obtain basis sets not already in Gaussian

import argparse
import sys

sys.path.append('../')

import KurtGroup.Kurt.chemical_information as ci
import KurtGroup.Kurt.xyz as xyz

def generateGaussianInputFileText(XYZ: xyz.xyz_to, charge: int) -> str:
    """Makes the text for a Gaussian input file

    Args:
        charge (int): The charge of the molecule

    Returns:
        str: Returns the file text
    """
    if not hasattr(XYZ, 'atoms'):
        print(f'''You need to process the xyz file using the class function processXYZ before generating the file:
Exiting program''')
        exit()

    basis_name = XYZ.basis
    if XYZ.BSE:
        BasisSet = ci.BasisSet()
        try:
            basis_mol = BasisSet.GenerateBasisSet('gaussian94', XYZ.basis, XYZ.atoms[:,0], SupressHeader=True)
        except RuntimeError:
            print("Failed to get basis set from BSE. Please check the spelling, upper-/lowercase is not important")
            exit()
        basis_name = XYZ.basis
        XYZ.basis = "GEN"

    filename_no_ext = XYZ.filename.replace('.xyz', '')
    XYZ.input_filename = f'{filename_no_ext}.com'

    if charge % 2 == 0:
        multiplicity = 1
    else:
        multiplicity = 2

    if XYZ.method.upper() == 'DFT':
        method = XYZ.functional
    else:
        method = XYZ.method

    XYZ.filetext = f'''%chk={filename_no_ext}.chk
%mem={XYZ.mem}GB
%nprocshared={XYZ.ncpus}
# {XYZ.calc} {method}/{XYZ.basis}

'''
    if XYZ.BSE:
        XYZ.filetext += f'{filename_no_ext}.xyz - {basis_name}\n'
    else:
        XYZ.filetext += f'{filename_no_ext}.xyz\n'

    XYZ.filetext += f'\n{charge} {multiplicity}\n'

    for atom in XYZ.atoms:
        XYZ.filetext += f'{atom[0]} {atom[1]: >14.9} {atom[2]: >19.9} {atom[3]: >19.9}\n'

    XYZ.filetext += '\n'

    if XYZ.BSE:
        XYZ.filetext += '****\n'
        XYZ.filetext += f'{basis_mol}\n'
        XYZ.filetext += '\n'

    XYZ.basis = basis_name

    return XYZ.filetext

def writeInputfile(XYZ: xyz.xyz_to):
    """Writes the file text to an input file with similar name as the xyz file
    """
    with open(XYZ.input_filename, 'w') as input_file:
        input_file.writelines(XYZ.filetext)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Use to make Gaussian input files from a given .xyz file.
Basis sets not implemented in Gaussian are imported from BSE: https://www.basissetexchange.org
NB! When specifying Pople-style basis sets with polarization functions, parentheses must be escaped by \\
Stars can still be used, however.''', epilog='''For help contact
        Theo Juncker von Buchwald
        fnc970@alumni.ku.dk
        Magnus Bukhave Johansen
        qhw298@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs=1, help='The file to convert into a .com file', metavar='.xyz file')
    parser.add_argument('calc', type=str,nargs=1,help='Keywords for the calculation. If there are spaces, have quotes around the entire thing.')

    CalculationGroup = parser.add_argument_group('Calculation options')
    CalculationGroup.add_argument('--charge', default=[0], nargs=1, type=int, help='Include to specify charge - 0 if not included')
    CalculationGroup.add_argument('--basis', default=['pc-1'], nargs=1, type=str, help='Include to specify basis set - pc-1 if not included')
    CalculationGroup.add_argument('--method', default=['cam-b3lyp'], nargs=1, type=str, help='Include to specify method for calculation - CAM-B3LYP if not included')
    CalculationGroup.add_argument('--cpu', default=[8], nargs=1, type=int, help='Include to specify the amount of cpu cores - 8 if not included')
    CalculationGroup.add_argument('--mem', default=[8], nargs=1, type=int, help='Include to specify the amount of memory in GB - 8 if not included')

    args = parser.parse_args()

    infile = args.infile[0]
    calc = args.calc[0]

    charge = args.charge[0]
    basis = args.basis[0]
    method = args.method[0]
    ncpus = args.cpu[0]
    mem = args.mem[0]


    A = xyz.xyz_to('Gaussian94', infile, calculation_type=calc, ncpus=ncpus, memory=mem)
    A.processXYZ()
    if A.checkFunctional(method):
        functional = method
        method = 'DFT'
        A.setMethod(method, functional)
    else:
        A.setMethod(method)
    A.setBasis(basis)
    filetext = generateGaussianInputFileText(A, charge)
    writeInputfile(A)

if __name__ == '__main__':
    main()
