#Makes a molecule sandwiched between two hemispheres
#Then makes Gaussian16 input files for both

import sys #For argv
import numpy as np
from collections import Counter #For number of unique elements
from ase.visualize.plot import plot_atoms
from ase.lattice.cubic import FaceCenteredCubic
from ase.visualize import view
from ase.build import fcc111, bcc111, hcp0001, diamond100
from ase.io import read, write
from ase.spacegroup import crystal

#INPUTS HERE
#------------------------------
#Inputs from command line
#.xyzfile of molecule
molfile = sys.argv[1]

#choose between using line number in xyz file or atom position in list (e.g. molden)
linenumber = False

#Line number of atoms in molxyz to lie on the connecting line of the NPs
atom1 = int(sys.argv[2])
atom2 = int(sys.argv[3])

#Charge
charge = int(sys.argv[4])

#Crystal Structure (Au, Ag, Cu, TiO2, NaCl, CoSb3, Pt, Pd)
crystal_structure = sys.argv[5]

#Computational
ncpus=16
mem=16

#Basis
basis_mol = '''H     0
S    3   1.00
      0.122520D+02           0.228220D-01
      0.186870D+01           0.155640D+00
      0.418210D+00           0.488980D+00
S    1   1.00
      0.106100D+00           0.100000D+01
P    1   1.00
      0.100000D+01           0.100000D+01
****
C     0
S    5   1.00
      0.125260D+04           0.557360D-02
      0.188570D+03           0.414960D-01
      0.428390D+02           0.182630D+00
      0.118180D+02           0.461290D+00
      0.355670D+01           0.449310D+00
S    5   1.00
      0.188570D+03          -0.274400D-03
      0.428390D+02          -0.255830D-02
      0.118180D+02          -0.333750D-01
      0.355670D+01          -0.873050D-01
      0.542580D+00           0.534150D+00
S    1   1.00
      0.160580D+00           0.100000D+01
P    3   1.00
      0.914260D+01           0.444990D-01
      0.192980D+01           0.231080D+00
      0.525220D+00           0.512270D+00
P    1   1.00
      0.136080D+00           0.100000D+01
D    1   1.00
      0.800000D+00           0.100000D+01
****
N     0
S    5   1.00
      0.174220D+04           0.546030D-02
      0.262210D+03           0.407060D-01
      0.595850D+02           0.180380D+00
      0.164890D+02           0.460090D+00
      0.499450D+01           0.451730D+00
S    5   1.00
      0.262210D+03          -0.301040D-03
      0.595850D+02          -0.255110D-02
      0.164890D+02          -0.349000D-01
      0.499450D+01          -0.874290D-01
      0.786730D+00           0.526600D+00
S    1   1.00
      0.228370D+00           0.100000D+01
P    3   1.00
      0.131320D+02           0.455150D-01
      0.281080D+01           0.237400D+00
      0.763510D+00           0.513920D+00
P    1   1.00
      0.195600D+00           0.100000D+01
D    1   1.00
      0.900000D+00           0.100000D+01
****
S     0
S    8   1.00
      0.555100D+05           0.588990D-03
      0.832670D+04           0.454600D-02
      0.189560D+04           0.233350D-01
      0.536770D+03           0.899650D-01
      0.174670D+03           0.254860D+00
      0.625780D+02           0.439980D+00
      0.235730D+02           0.315490D+00
      0.589080D+01           0.328260D-01
S    8   1.00
      0.832670D+04          -0.204800D-04
      0.189560D+04          -0.454290D-04
      0.536770D+03          -0.179370D-02
      0.174670D+03          -0.972520D-02
      0.625780D+02          -0.663990D-01
      0.235730D+02          -0.855020D-01
      0.589080D+01           0.495880D+00
      0.218520D+01           0.649880D+00
S    8   1.00
      0.189560D+04           0.834610D-05
      0.536770D+03           0.983060D-04
      0.174670D+03           0.100400D-02
      0.625780D+02           0.544090D-02
      0.235730D+02           0.114000D-01
      0.589080D+01          -0.951340D-01
      0.218520D+01          -0.247350D+00
      0.452820D+00           0.639730D+00
S    1   1.00
      0.162020D+00           0.100000D+01
P    6   1.00
      0.395450D+03           0.462660D-02
      0.928340D+02           0.358980D-01
      0.292660D+02           0.148540D+00
      0.105950D+02           0.358780D+00
      0.402720D+01           0.454170D+00
      0.154840D+01           0.205450D+00
P    6   1.00
      0.928340D+02           0.109340D-03
      0.292660D+02          -0.163580D-02
      0.105950D+02          -0.506320D-02
      0.402720D+01          -0.239600D-01
      0.154840D+01           0.103880D+00
      0.461340D+00           0.580630D+00
P    1   1.00
      0.131310D+00           0.100000D+01
D    1   1.00
      0.550000D+00           0.100000D+01'''
basis_NP = 'LANL2DZ'
ecp_NP = 'LANL2'

#Method
method = 'cam-b3lyp'

#Sphere diameter
diameter = 13

#Point both inwards or both outwards
inwards = True

#Return .xyzfile of both cases
returnxyz = True

# ------------------------------------ FUNCTIONS ------------------------------------

#Add new elements to first five functions if introducing new atoms
def atm_type(i):
    dic = {'Au' : ['Au'], 'Cu' : ['Cu'], 'Ag' : ['Ag'], 'TiO2' : ['Ti','O'], 'Pd' : ['Pd'], 'Pt' : ['Pt'], 'NaCl' : ['Na','Cl'], 'CoSb3' : ['Co','Sb']}
    if i in dic:
        return dic[i]
    else:
        print(i, "Is not in the Atm Type dictionary")
        exit()

def lattice_length(i):
				#In Ångstrom
    dic = {'Cu' : 3.597, 'Pd' : 3.859, 'Ag' : 4.079, 'Pt' : 3.912, 'Au' : 4.065, 'TiO2' : [4.6,2.95], 'NaCl' : 5.64, 'CoSb3' : 9.04}
    if i in dic:
        return dic[i]
    else:
        print(i, "Is not in the Lattice Length dictionary")
        exit()

#Rocksalt = Halite
def lattice_type(i):
    dic = {'Cu' : 'FCC', 'Pd' : 'FCC', 'Ag' : 'FCC', 'Pt' : 'FCC', 'Au' : 'FCC', 'TiO2' : 'Rutile', 'NaCl' : 'Rocksalt', 'CoSb3' : 'Skutterudite'}
    if i in dic:
        return dic[i]
    else:
        print(i, "Is not in the Lattice Type dictionary")
        exit()

def vdW_radius(i):
				#In pm
    dic = {'H' : 1.09, 'C' : 1.7, 'O' : 1.52, 'Na' : 2.27, 'S' : 1.8, 'Cl' : 1.75, 'Cu' : 1.4, 'Pd' : 1.63, 'Ag' : 1.72, 'Sb' : 2.06, 'Pt' : 1.75, 'Au' : 1.66, 'TiO2' : 1.52,
           'NaCl' : 2.27, 'CoSb3' : 2.06, 'N' : 1.55}
				#https://periodictable.com/Properties/A/VanDerWaalsRadius.v.html
    if i in dic:
        return dic[i]
    else:
        print(i, "Is not in the Van der Waals Radius dictionary")
        exit()

def polarizability(i):
				#In au^3
    dic = {'H' : 4.50114, 'C' : 8.465, 'O' : 16.15551282, 'Cu' : 33.742, 'Ag' : 49.9843, 'Au' : 31.04, 'Ti' : 1.254517, 'Pt' : 42.515}
				#Au - J. Phys. Chem. C, 114 (48) (2010), pp. 20870-20876
				#Ag - M. Pereiro, D. Baldomir, Structure and static response of small silver clusters to an external electric field, arXiv preprint physics, 2007, 0702238.
				#Cu - Phys. Rev. A, 99 (1) (2019), p. 012503, J. Chem. Phys., 117 (7) (2002), pp. 3208-3218, J. Chem. Phys., 120 (22) (2004), pp. 10450-10454
				#Ti - Phys. Rev. B 71, 085418
				#O  - Phys. Rev. B 71, 085418
    if i in dic:
        return dic[i]
    else:
        print(i, "Is not in the Polarizability dictionary")
        exit()

def get_atm_number(i): #Contains elements up to U
    dic = {'H' : 1, 'He' : 2, 'Li' : 3, 'Be' : 4, 'B' : 5, 'C' : 6, 'N' : 7, 'O' : 8, 'F' : 9, 'Ne' : 10, 'Na' : 11, 'Mg' : 12, 'Al' : 13, 'Si' : 14, 'P' : 15, 'S' : 16, 'Cl' : 17, 'Ar' : 18, 'K' : 19, 'Ca' : 20, 'Sc' : 21, 'Ti' : 22, 'V' : 23, 'Cr' : 24, 'Mn':  25, 'Fe' : 26, 'Co' : 27, 'Ni' : 28, 'Cu' : 29, 'Zn' : 30, 'Ga' : 31, 'Ge' : 32, 'As' : 33, 'Se' : 34, 'Br' : 35, 'Kr' : 36, 'Rb' : 37, 'Sr' : 38, 'Y' : 39, 'Zr' : 40, 'Nb' : 41, 'Mo' : 42, 'Tc' : 43, 'Ru' : 44, 'Rh' : 45, 'Pd' : 46, 'Ag' : 47, 'Cd' : 48, 'In' : 49, 'Sn' : 50, 'Sb' : 51, 'Te' : 52, 'I' : 53, 'Xe' : 54, 'Cs' : 55, 'Ba' : 56, 'La' : 57, 'Ce' : 58, 'Pr' : 59, 'Nd' : 60, 'Pm' : 61, 'Sm' : 62, 'Eu' : 63, 'Gd' : 64, 'Tb' : 65, 'Dy' : 66, 'Ho' : 67, 'Er' : 68, 'Tm' : 69, 'Yb' : 70, 'Lu' : 71, 'Hf' : 72, 'Ta' : 73, 'W' : 74, 'Re' : 75, 'Os' : 76, 'Ir' : 77, 'Pt' : 78, 'Au' : 79, 'Hg' : 80, 'Tl':  81, 'Pb' : 82, 'Bi' : 83, 'Po' : 84, 'At' : 85, 'Rn':  86, 'Fr' : 87, 'Ra' : 88, 'Ac' : 89, 'Th' : 90, 'Pa' : 91, 'U' : 92, 'Np' : 93, 'Pu' : 94, 'Am' : 95, 'Cm' : 96, 'Bk' : 97, 'Cf' : 98, 'Es' : 99, 'Fm' : 100, 'Md' : 101, 'No' : 102, 'Lr' : 103, 'Rf' : 104, 'Db' : 105, 'Sg' : 106, 'Bh' : 107, 'Hs' : 108, 'Mt' : 109, 'Ds' : 110, 'Rg' : 111, 'Cn' : 112, 'Nh' : 113, 'Fl' : 114, 'Mc' : 115, 'Lv' : 116, 'Ts' : 117, 'Og' : 118}
    if i in dic:
        return dic[i]
    else:
        print(i, "Is not in the Get Atom Number dictionary")
        exit()

def make_sphere(atom_array,sphere_radii):
    shiftx = max(atom_array[:,0])
    shifty = max(atom_array[:,1])
    shiftz = max(atom_array[:,2])
    atom_array[:,0] = atom_array[:,0] - 0.5*shiftx
    atom_array[:,1] = atom_array[:,1] - 0.5*shifty
    atom_array[:,2] = atom_array[:,2] - 0.5*shiftz
    sphere = np.empty((0,3))
    for i in range(len(atom_array)):
        if sum(atom_array[i,:]**2)**0.5 < sphere_radii:
            sphere = np.vstack([sphere,atom_array[i,:]])
    return sphere

def spherical(radius, mol_max, mol_min, inwards, vacuum_dist_left, vacuum_dist_right, atoms, symbols):
    sphere_coor = make_sphere(atoms,radius)

    left_hemisphere = np.empty((0,3))
    left_symbols = np.array([])
    for i in range(len(sphere_coor)):
        if sphere_coor[i,0] < 0:
            left_hemisphere = np.vstack([left_hemisphere,sphere_coor[i,:]])
            left_symbols = np.append(left_symbols,symbols[i])
    right_hemisphere = np.copy(left_hemisphere) * np.array([-1,1,1])
    right_symbols = left_symbols
    if inwards == False:
        left_hemisphere[:,0] = left_hemisphere[:,0] + mol_min - vacuum_dist_left
        right_hemisphere[:,0] = right_hemisphere[:,0] + mol_max + vacuum_dist_right
    else:
        left_hemisphere[:,0] = left_hemisphere[:,0] + mol_max + radius + vacuum_dist_left
        right_hemisphere[:,0] = right_hemisphere[:,0] + mol_min - radius - vacuum_dist_right
    return left_hemisphere, right_hemisphere, left_symbols, right_symbols

def slab(mol_max, mol_min, lat_length, vacuum_dist_left, vacuum_dist_right, atoms, symbols):
    shifty = max(atoms[:,1])
    shiftz = max(atoms[:,2])
    atoms[:,1] -= 0.5*shifty
    atoms[:,2] -= 0.5*shiftz
    left_slab = np.copy(atoms)
    right_slab = np.copy(atoms) * np.array([-1,1,1])
    left_symbols = np.copy(symbols)
    right_symbols = np.copy(symbols)
    left_slab[:,0] = left_slab[:,0] + mol_min - vacuum_dist_left - lat_length[0]
    right_slab[:,0] = right_slab[:,0] + mol_max + vacuum_dist_right + lat_length[0]
    return left_slab,right_slab,left_symbols,right_symbols

def get_rot_matrix(point_for_rot,direc_vec,thet):
    #Takes a point on the rotational axis, the normalized direction vector and the angle
    #The rotation matrix is 4x4 with the fourth dimension translating the origin to (a,b,c) and back
    a,b,c = point_for_rot
    u,v,w = direc_vec
    u2,v2,w2 = u*u, v*v,w*w
    cosT = np.cos(thet)
    sinT = np.sin(thet)
    m11 = u2 + (v2 + w2) * cosT
    m12 = u*v * (1-cosT) - w*sinT
    m13 = u*w * (1-cosT) + v*sinT
    m14 = (a*(v2 + w2) - u*(b*v + c*w)) * (1-cosT) + (b*w - c*v)*sinT
    m21 = u*v * (1-cosT) + w*sinT
    m22 = v2 + (u2 + w2) * cosT
    m23 = v*w * (1-cosT) - u*sinT
    m24 = (b*(u2 + w2) - v*(a*u + c*w))*(1-cosT) + (c*u - a*w)*sinT
    m31 = u*w * (1-cosT) - v*sinT
    m32 = v*w * (1-cosT) + u*sinT
    m33 = w2 + (u2 + v2) * cosT
    m34 = (c*(u2 + v2) - w*(a*u + b*v))*(1-cosT) + (a*v - b*u)*sinT

    return np.array([m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34])

def rotate_point_new(point, rot_mat):
    #Essentially matrix multiplication
    p = np.zeros(3)
    x,y,z = point
    m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34 = rot_mat
    p[0] = m11*x + m12*y + m13*z + m14
    p[1] = m21*x + m22*y + m23*z + m24
    p[2] = m31*x + m32*y + m33*z + m34

    return p

# ------------------------------------ SETUP ------------------------------------

atmtype = atm_type(crystal_structure)

#Rescaling to get index in array

if linenumber:
    atom1 -=3
    atom2 -=3
else:
    atom1 -=1
    atom2 -=1

namesmol = np.array([])
molxyz = np.empty((0,3))
rotxyz = np.empty((0,3))

# ------------------------------------ ROTATING MOLECULE ------------------------------------

lines_to_add = []

with open(molfile, "r") as f:
    lines = f.readlines()
    for i in range(2, len(lines)):
       x = lines[i].split()
       namesmol = np.append(namesmol,x[0])
       molxyz = np.vstack([molxyz,np.array([float(x[1]),float(x[2]),float(x[3])])])
#We get the axis between the two points to be parallel to the x-axis and then translate it such that it lies in the x-axis
#Get normalized vector defining the axis
v1 = molxyz[atom1,:] - molxyz[atom2,:]
v1 *= 1 / np.sqrt(np.dot(v1,v1))

#Calculate angle in order to be parallel to the x-axis
theta = np.arccos(np.dot(v1,np.array([1,0,0])))

#The direction vector for the rotation in orthogonal to both v1 and the x-axis. Note that the normalization is crucial.
dir_vec = np.cross(v1,np.array([1,0,0]))
dir_vec *= 1/np.sqrt(np.dot(dir_vec,dir_vec))

#Get elements of rotation matrix and keep using them to convert every point
rot_mat_elems = get_rot_matrix(molxyz[atom1,:],dir_vec,theta)

#Rotate all atoms
for i in range(len(molxyz)):
    p = rotate_point_new(molxyz[i,:],rot_mat_elems)
    rotxyz = np.vstack([rotxyz,np.array([p[0],p[1],p[2]])])

#Place atom1 in (0,0,0) and move the rest accordingly
for i in range(3):
    rotxyz[:,i] -= rotxyz[atom1,i]

#Place furthest atom in (0,y,z) and move the rest accordingly
rotxyz[:,0] -= np.amax(rotxyz[:,0])

#Get length of molecule in the x-direction
lengthmol = rotxyz[:,0].max() - rotxyz[:,0].min()

#Place molecule at correct x-position
molxyz = rotxyz

# ------------------------------------ CREATING NANOPARTICLES ------------------------------------

mol_min = min(molxyz[:,0])
mol_max = max(molxyz[:,0])

index_min = np.argmin(molxyz[:,0])
index_max = np.argmax(molxyz[:,0])

radius = diameter*0.5


lat_length = lattice_length(crystal_structure)
lat_type = lattice_type(crystal_structure)

vdW_electrode = vdW_radius(crystal_structure)
vdW_left = vdW_radius(namesmol[index_min])
vdW_right = vdW_radius(namesmol[index_max])

vacuum_dist_left = vdW_left + vdW_electrode
vacuum_dist_right = vdW_right + vdW_electrode

if type(lat_length) == float:
    dim = int(np.ceil(diameter/lat_length)) * 2
else:
    dim1 = int(np.ceil(diameter/lat_length[0]))
    dim2 = int(np.ceil(diameter/lat_length[1]))

if lat_type == 'FCC':
    atoms = fcc111(symbol=atmtype[0],size=(dim,dim,dim),a=lat_length,orthogonal=True)
elif lat_type == 'BCC':
    atoms = bcc111(symbol=atmtype[0],size=(dim,dim,4.5*dim),a=lat_length,orthogonal=True)
elif lat_type == 'HCP':
    atoms = hcp0001(symbol=atmtype[0],size=(dim,dim,3*dim),a=lat_length[0],c=lat_length[1],orthogonal=True)
elif lat_type == 'Diamond':
    atoms = diamond100(symbol=atmtype[0],size=(dim,dim,3*dim),a=lat_length,orthogonal=True)
elif lat_type == 'Rutile':
    a = lat_length[0]
    c = lat_length[1]
    atoms = crystal([atmtype[0], atmtype[1]], basis=[(0, 0, 0), (0.3, 0.3, 0.0)], spacegroup=136, cellpar=[a, a, c, 90, 90, 90],size=(1,dim1,dim2))
elif lat_type == 'Rocksalt':
    a = lat_length
    atoms = crystal([atmtype[0], atmtype[1]], basis=[(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225, cellpar=[a, a, a, 90, 90, 90],size=(2,dim1,dim2))
elif lat_type == 'Skutterudite':
    a = 9.04
    atoms = crystal([atmtype[0], atmtype[1]], basis=[(0.25, 0.25, 0.25), (0.0, 0.335, 0.158)], spacegroup=204, cellpar=[a, a, a, 90, 90, 90], size=(2,dim1,dim2))

atoms_symbol = atoms.get_chemical_symbols()
atoms = atoms.get_positions()

if type(lat_length) == float:
    left, right, left_symbols, right_symbols = spherical(radius, mol_max, mol_min, inwards, vacuum_dist_left, vacuum_dist_right, atoms, atoms_symbol)
else:
    left, right, left_symbols, right_symbols = slab(mol_max, mol_min, lat_length, vacuum_dist_left, vacuum_dist_right, atoms, atoms_symbol)

# ------------------------------------ CREATING OUTPUT FILES ------------------------------------

if returnxyz:
    #Build .xyz files
    for j in ['left','right']:
      lines_to_add = []
      lines_to_add.append(str(left[:,0].size+molxyz[:,0].size)+'\n')
      lines_to_add.append('\n')

      for i in range(len(left[:,0])):
         if j == 'left':
           lines_to_add.append(''.join([atoms_symbol[i],' ',f"{left[i,0]:.6f}",' ', f"{left[i,1]:.6f}", ' ',f"{left[i,2]:.6f}" ,'\n']))
         else:
           lines_to_add.append(''.join([atoms_symbol[i],' ',f"{right[i,0]:.6f}",' ', f"{right[i,1]:.6f}", ' ',f"{right[i,2]:.6f}" ,'\n']))

      for i in range(len(molxyz[:,0])):
         lines_to_add.append(''.join([namesmol[i],' ',f"{molxyz[i,0]:.6f}",' ', f"{molxyz[i,1]:.6f}", ' ',f"{molxyz[i,2]:.6f}" ,'\n']))

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
     lines_com.append(''.join([namesmol[i],'(Fragment=1)',' ',f"{molxyz[i,0]:.9f}",' ', f"{molxyz[i,1]:.9f}", ' ',f"{molxyz[i,2]:.9f}" ,'\n']))
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
           lines_com.append('*********'+'\n')
#   lines_com.append(''.join([f"{x} " for x in list(set(namesmol))])+'\n')
   lines_com.append(basis_mol+'\n')
   lines_com.append('*********'+'\n')
   lines_com.append('\n')
   for i in atmtype:
       lines_com.append(i+' 0\n')
       lines_com.append(ecp_NP+'\n')
       lines_com.append('\n')
   lines_com.append('\n')
   with open(filename+'.com','w') as f:
    f.writelines(lines_com)


if returnxyz:
  print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}_left.xyz')
  print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}_right.xyz')
print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}_left.com')
print(molfile[:-4]+f'_charge_{charge}_{crystal_structure}_right.com')
