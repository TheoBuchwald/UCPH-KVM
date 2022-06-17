# Imports
# -------
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''A script to convert xyz files to inp files for ORCA, specialized for calculations using diffuse basis sets''', epilog='''For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk
    Magnus Bukhave Johansen
    qhw298@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='.xyz file')
    parser.add_argument('calctype',type=str, nargs=1,choices=['opt','exc'],help='The type of calculation to do. opt is geometry optimzation and exc is excitation energies.',metavar="calculation type")
    parser.add_argument('--charge', default=[0], nargs=1, type=int, help='Include to specify charge - 0 if not included')
    parser.add_argument('--mem', default=[4800], nargs=1, type=int, help='Include to specify the amount of memory in MB pr. core - 4800 if not included')
    parser.add_argument('--noaug', default=['pc-2'], nargs=1, type=str, help='Include to specify the basis set without diffuse functions - default is pc-2')
    parser.add_argument('--aug', default=['aug-pc-2'], nargs=1, type=str, help='Include to specify the basis set with diffuse functions - default is aug-pc-2')
    parser.add_argument('--method', default=['M062X'], nargs=1, type=str, help='Include to specify the method used for the geometry optimization - default is M06-2X')

    args = parser.parse_args()

    input_files = args.infile
    noaug = args.noaug[0]
    aug = args.aug[0]
    method = args.method[0]
    keyword = args.calctype

    if method == 'HF':
        pre = 'RHF'
    else:
        pre = 'RKS'

    if noaug != 'pc-2' and aug == 'aug-pc-2':
        aug = 'aug-' + noaug

    memory = f'{args.mem[0]}' #per mpi process in MB

    charge = args.charge[0]
    if charge % 2 == 0:
        multiplicity = 1
    else:
        multiplicity = 2

    # Driver part of the script
    # -------------------------

    #xyzfile = 'test.xyz'
    for xyzfile in input_files:
        with open(xyzfile) as thefile:
            content=thefile.readlines()

        name, ext = os.path.splitext(xyzfile)

  
        # Read xyz coordinates
        # -------------------

        for xyzfile in input_files:
            token = []
            x = []
            y = []
            z = []

            c = 1
            for line in content:
                if c > 2:
                    line = line.strip().split()
                    token.append(line[0])
                    x.append(float(line[1]))
                    y.append(float(line[2]))
                    z.append(float(line[3]))
                c += 1

            with open(name + '.inp', 'w') as slutfil:
            # Writing orca input file
            # -------------------------

            #Add more jobtypes if needed

                slutfil.write(f"! {pre} {noaug} LSD VeryTightSCF\n")
                slutfil.write("%method\nend\n")
                slutfil.write('%maxcore '+memory+'\n')
                slutfil.write(r'%base "temp"' +'\n')
                slutfil.write('* xyz '+str(charge)+' '+str(multiplicity)+'\n')

                for i in range(len(x)):
                    slutfil.write(token[i])
                    slutfil.write('  ')
                    slutfil.write('%f' % x[i])
                    slutfil.write('  ')
                    slutfil.write('%f' % y[i])
                    slutfil.write('  ')
                    slutfil.write('%f' % z[i] + '\n')

                slutfil.write('*'+'\n'+'\n')
                slutfil.write('$new_job\n')
                slutfil.write(f"! {pre} {aug} BP86 VeryTightSCF TRAH MOREAD\n")
                slutfil.write("%method\nend\n")
                slutfil.write('%maxcore '+memory+'\n')
                slutfil.write(r'%moinp "temp.gbw"' + '\n')
                slutfil.write(r'%base "temp2"' + '\n')
                slutfil.write(r"%scf GuessMode CMatrix" + '\n')    
                slutfil.write(r"  sthresh  1e-6" + '\n')
                slutfil.write('end'+'\n')
                slutfil.write('* xyz '+str(charge)+' '+str(multiplicity)+'\n')

                for i in range(len(x)):
                    slutfil.write(token[i])
                    slutfil.write('  ')
                    slutfil.write('%f' % x[i])
                    slutfil.write('  ')
                    slutfil.write('%f' % y[i])
                    slutfil.write('  ')
                    slutfil.write('%f' % z[i] + '\n')

                slutfil.write('*'+'\n'+'\n')
                slutfil.write('$new_job\n')
                if keyword == 'opt':
                    slutfil.write(f"! {pre} {aug} {method} VeryTightSCF TRAH MOREAD TightOpt Freq\n")
                else:
                    slutfil.write(f"! {pre} {aug} {method} VeryTightSCF TRAH MOREAD\n")
                    slutfil.write(f"%TDDFT  NROOTS 15\n    TRIPLETS TRUE\n    TDA FALSE\nend\n")
                slutfil.write("%method\nend\n")
                slutfil.write('%maxcore '+memory+'\n')
                slutfil.write(r'%moinp "temp2.gbw"'+'\n')
                slutfil.write(r'%scf sthresh 1e-6' +'\n')
                slutfil.write('end'+'\n')
                slutfil.write('* xyz '+str(charge)+' '+str(multiplicity)+'\n')

                for i in range(len(x)):
                    slutfil.write(token[i])
                    slutfil.write('  ')
                    slutfil.write('%f' % x[i])
                    slutfil.write('  ')
                    slutfil.write('%f' % y[i])
                    slutfil.write('  ')
                    slutfil.write('%f' % z[i] + '\n')

                slutfil.write('*'+'\n'+'\n')

