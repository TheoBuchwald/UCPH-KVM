# Imports
# -------
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''A script to convert xyz files to inp files for ORCA''', epilog='''For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='.xyz file')
    parser.add_argument('keyword', nargs=1, type=int, help='Include to specify keyword string', choices=range(1,11))
    parser.add_argument('--charge', default=0, nargs=1, type=int, help='Include to specify charge - 0 if not included')
    parser.add_argument('--mem', default=4800, nargs=1, type=int, help='Include to specify the amount of memory in MB pr. core - 4800 if not included')
    parser.add_argument('--extra1', action='store_true')
    parser.add_argument('--extra2', action='store_true')

    args = parser.parse_args()

    input_files = args.infile[0]
    jobtype = args.keyword[0]

    memory = f'{args.mem}' #per mpi process in MB

    charge = args.charge
    if charge % 2:
        multiplicity = 1
    else:
        multiplicity = 2

    extra = args.extra1
    extra2 = args.extra2

    keyword_string = { 1: """! RKS 6-311+G* M062X
! Opt Freq
%method
Grid 7
end""",

    2: """! RKS 6-31G** B3LYP VeryTightSCF TightOPT
! Opt Freq
%method
end""",

    3: """! RKS 6-311G* M062X TightSCF TightOPT SlowConv OptTS
! Freq
%method
Grid 7
end
%geom
Calc_Hess true
end""",

    4: """! aug-cc-pVDZ MP2 VeryTightSCF NoFrozenCore
%pal nprocs 32
end
%mp2 density unrelaxed
end
%elprop dipole true
end""",

    5: """! aug-cc-pVDZ CCSD VeryTightSCF NoFrozenCore
%pal nprocs 32
end
%mdci density unrelaxed
end
%elprop dipole true
end""",

    6: """! RI-MP2 VeryTightSCF NoFrozenCore
%pal nprocs 32
end
%basis basis "aug-cc-pVDZ/C"
newGTO H "cc-pVDZ/C" end
end
%mp2 density unrelaxed
end
%elprop dipole true
end""",

    7: """! CCSD VeryTightSCF NoFrozenCore
%pal nprocs 32
end
%basis basis "aug-cc-pVDZ"
newGTO H "cc-pVDZ" end
end
%mdci density unrelaxed
end
%elprop dipole true
end""",

    8: """! RHF DLPNO-CCSD cc-pVDZ cc-pVDZ/C TightSCF TIGHTOPT TightPNO NOFROZENCORE
! Opt NumGrad
%pal nprocs 16
end""",

    9: """! 6-311++G** TightSCF MP2
%mp2
density unrelaxed
natorbs true
end""",

    10: """! PAL4
! RKS pc-1 CAM-B3LYP VeryTightSCF TightOPT
! Opt Freq
%method
end
%elprop
Polar 1
end"""}

    # Driver part of the script
    # -------------------------

    #xyzfile = 'test.xyz'
    for xyzfile in input_files:
        with open(xyzfile) as thefile:
            content=thefile.readlines()

        name, ext = os.path.splitext(xyzfile)

        extra_calc="""$new_job
! RKS 6-311+G* M062X TightSCF TightOpt
! Opt Freq
%method
Grid 7
end
%maxcore 4000
* xyzfile 0 1"""


        extra_calc2 = '''$new_job
! 6-311++G** TightSCF
! moread
%moinp "'''+name+'''.mp2nat"
%maxcore 8000
%casscf
nel 8
norb 8
PTMethod FIC_CASPT2K
end'''


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

                slutfil.write(keyword_string[jobtype] + '\n')

                slutfil.write('%maxcore '+memory+'\n')
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

                if extra:
                    slutfil.write(extra_calc+'\n')
                if extra2:
                    slutfil.write(extra_calc2+'\n')
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
