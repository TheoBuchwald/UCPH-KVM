#Make Gaussian style input file from xyz file
#Use the BSE API to obtain basis sets not already in Gaussian

import argparse
from Kurt import xyz

def main():
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
    filetext = A.generateGaussianInputFileText(charge)
    A.writeInputfile()

if __name__ == '__main__':
    main()