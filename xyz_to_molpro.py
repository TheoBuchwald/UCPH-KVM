# Imports
# -------
import numpy as np
import sys
import os

# Keyword specification
# ---------------------

memory = '4800' #per mpi process in Mwords

keyword_string = """basis=6-311++G(D,P)
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
wf,66,1,0}"""

keyword_string2 = """basis=cc-pVTZ
hf
ccsd
optg
"""

keyword_string3 = """basis={
default,VTZ-F12
set,df
default,VTZ-F12/mp2fit   !density fitting basis
set,jk
default,VTZ-F12/jkfit    !density fitting basis for Fock and exchange matrices
set,ri
default,VTZ-F12/optri    !ri cabs basis
}
hf
ccsd(t)-f12,df_basis=df,df_basis_exch=jk,ri_basis=ri"""

# Driver part of the script
# -------------------------

xyzfile = sys.argv[1]

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
    slutfil.write(keyword_string2)
