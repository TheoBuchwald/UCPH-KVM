
import argparse
from Kurt import xyz

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''A script to convert xyz files to mol files for DALTON''', epilog='''For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='.xyz file')
    parser.add_argument('--charge', default=[0], nargs=1, type=int, help='Include to specify charge - 0 if not included')
    parser.add_argument('--basis', default=['pc-1'], nargs=1, type=str, help='Include to specify basis set of the molecular atoms - pc-1 if not included')
    parser.add_argument('--RIbasis', nargs=1, type=str, help='Include to specify basis set of the molecular atoms - RI-BASIS if not included')

    args = parser.parse_args()

    input_files = args.infile
    basis = args.basis[0]
    charge = args.charge[0]

    for input_file in input_files:
        A = xyz.xyz_to('Dalton', input_file)
        A.processXYZ()
        A.setBasis(basis)
        filetext = A.generateDaltonInputFileText(charge)
        A.writeInputfile()

if __name__ == '__main__':
    main()
