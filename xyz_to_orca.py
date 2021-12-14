# Imports
# -------
import numpy as np
import sys
import os


# choices
# -------
extra = False
extra2 = False

# Driver part of the script
# -------------------------

xyzfile = sys.argv[1]
jobtype = sys.argv[2]
#xyzfile = 'test.xyz'

with open(xyzfile) as thefile:
    content=thefile.readlines()

name, ext = os.path.splitext(xyzfile)

# Keyword specification
# ---------------------

memory = '4800' #per mpi process in MB

multiplicity = 1
charge = 0

keyword_string = """! RKS 6-311+G* M062X
! Opt Freq
%method
Grid 7
end"""

keyword_string3 = """! RKS 6-31G** B3LYP VeryTightSCF TightOPT
! Opt Freq
%method
end"""

keyword_string2 = """! RKS 6-311G* M062X TightSCF TightOPT SlowConv OptTS
! Freq
%method
Grid 7
end
%geom
Calc_Hess true
end"""

keyword_string4 = """! aug-cc-pVDZ MP2 VeryTightSCF NoFrozenCore
%pal nprocs 32
end
%mp2 density unrelaxed
end
%elprop dipole true
end"""

keyword_string6 = """! aug-cc-pVDZ CCSD VeryTightSCF NoFrozenCore
%pal nprocs 32
end
%mdci density unrelaxed
end
%elprop dipole true
end"""

keyword_string7 = """! RI-MP2 VeryTightSCF NoFrozenCore
%pal nprocs 32
end
%basis basis "aug-cc-pVDZ/C"
newGTO H "cc-pVDZ/C" end
end
%mp2 density unrelaxed
end
%elprop dipole true
end"""

keyword_string8 = """! CCSD VeryTightSCF NoFrozenCore
%pal nprocs 32
end
%basis basis "aug-cc-pVDZ"
newGTO H "cc-pVDZ" end
end
%mdci density unrelaxed
end
%elprop dipole true
end"""

keyword_string9 = """! RHF DLPNO-CCSD cc-pVDZ cc-pVDZ/C TightSCF TIGHTOPT TightPNO NOFROZENCORE
! Opt NumGrad
%pal nprocs 16
end"""

keyword_string5 = """! 6-311++G** TightSCF MP2
%mp2
density unrelaxed
natorbs true
end"""

keyword_string10 = """! PAL4
! RKS pc-1 CAM-B3LYP VeryTightSCF TightOPT
! Opt Freq
%method
end
%elprop
Polar 1
end"""

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

slutfil = open(name+'.inp', 'w')

# Writing orca input file
# -------------------------

#Add more jobtypes if needed

if jobtype=='1':
    slutfil.write(keyword_string+'\n')
if jobtype=='2':
    slutfil.write(keyword_string2+'\n')
if jobtype=='3':
    slutfil.write(keyword_string3+'\n')
if jobtype=='4':
    slutfil.write(keyword_string4+'\n')
if jobtype=='5':
    slutfil.write(keyword_string5+'\n')
if jobtype=='6':
    slutfil.write(keyword_string6+'\n')
if jobtype=='7':
    slutfil.write(keyword_string7+'\n')
if jobtype=='8':
    slutfil.write(keyword_string8+'\n')
if jobtype=='9':
    slutfil.write(keyword_string9+'\n')
    memory = '7000'
if jobtype=='10':
    slutfil.write(keyword_string10+'\n')

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
