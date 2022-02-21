# Imports
# -------
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''A script to convert xyz files to molpro files''', epilog='''For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='.xyz file')
    parser.add_argument('keyword', nargs=1, type=int, help='Include to specify keyword string', choices=range(1,4))

    args = parser.parse_args()

    input_files = args.infile[0]
    keyword = args.keyword[0]

    # Keyword specification
    # ---------------------

    memory = '4800' #per mpi process in Mwords

    keyword_string = {1: """basis=6-311++G(D,P)
hf
{mp2
core,0
cphf,maxit=100
natorb,record=2200}
{casscf
start,record=2200,type=natural
maxiter,40
closed,29
occ,37
wf,66,1,0}
{rs2c,shift=0.1;
closed,29
occ,37
wf,66,1,0}""",
    2: """basis=cc-pVTZ
hf
ccsd
optg
""",
    3: """basis={
default,VTZ-F12
set,df
default,VTZ-F12/mp2fit   !density fitting basis
set,jk
default,VTZ-F12/jkfit    !density fitting basis for Fock and exchange matrices
set,ri
default,VTZ-F12/optri    !ri cabs basis
}
hf
ccsd(t)-f12,df_basis=df,df_basis_exch=jk,ri_basis=ri"""}

    # Driver part of the script
    # -------------------------

    for xyzfile in input_files:
        with open(xyzfile) as thefile:
            content=thefile.readlines()

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

        name, ext = os.path.splitext(xyzfile)
        with open(name + '.inp', 'w') as slutfil:

        # Writing molpro input file
        # -------------------------

            slutfil.write('***,'+name+'\n')
            slutfil.write('memory,'+memory+',m'+'\n')
            slutfil.write('file,1,'+name+'.int'+'\n')
            slutfil.write('file,2,'+name+'.wfu'+'\n')
            slutfil.write('angstrom'+'\n'+'geometry={'+'\n')

            for i in range(len(x)):
                slutfil.write(token[i])
                slutfil.write('  ')
                slutfil.write('%f' % x[i])
                slutfil.write('  ')
                slutfil.write('%f' % y[i])
                slutfil.write('  ')
                slutfil.write('%f' % z[i] + '\n')

            slutfil.write('}'+'\n'+'\n')
            slutfil.write(keyword_string[keyword])
