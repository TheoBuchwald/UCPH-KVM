
import numpy as np
from ase.build import fcc110, bcc110, hcp0001, diamond100
from ase.spacegroup import crystal
import Kurt.chemical_information as ci
# from . import chemical_information as ci

class NanoParticle():
    def __init__(self, nanoparticle: str) -> None:
        self.structure = nanoparticle

        Atoms = {
            'Au' : ['Au'], 'Cu' : ['Cu'], 'Ag' : ['Ag'], 'TiO2' : ['Ti', 'O'],
            'Pd' : ['Pd'], 'Pt' : ['Pt'], 'NaCl' : ['Na', 'Cl'], 'CoSb3' : ['Co', 'Sb']
        }

        Length = { # In Ångstrom
            'Cu' : 3.597, 'Pd' : 3.859, 'Ag' : 4.079, 'Pt' : 3.912, 'Au' : 4.065,
            'TiO2' : [4.6, 2.95], 'NaCl' : 5.64, 'CoSb3' : 9.04
        }

        Type = {
            'Cu' : 'FCC', 'Pd' : 'FCC', 'Ag' : 'FCC', 'Pt' : 'FCC', 'Au' : 'FCC',
            'TiO2' : 'Rutile', 'NaCl' : 'Rocksalt', 'CoSb3' : 'Skutterudite'
        }

        VdW = {
            'H' : 1.09, 'C' : 1.7, 'O' : 1.52, 'Na' : 2.27, 'S' : 1.8, 'Cl' : 1.75,
            'Cu' : 1.4, 'Pd' : 1.63, 'Ag' : 1.72, 'Sb' : 2.06, 'Pt' : 1.75,
            'Au' : 1.66, 'TiO2' : 1.52, 'NaCl' : 2.27, 'CoSb3' : 2.06, 'N' : 1.55
        }

        Errors = 0

        try:
            self.atomtypes = Atoms[self.structure]
        except KeyError:
            print(f'{self.structure} has not been implemented in Atoms dictionary')
            Errors += 1

        try:
            self.lattice_length = Length[self.structure]
        except KeyError:
            print(f'{self.structure} has not been implemented in Length dictionary')
            Errors += 1

        try:
            self.lattice_type = Type[self.structure]
        except KeyError:
            print(f'{self.structure} has not been implemented in Type dictionary')
            Errors += 1

        try:
            self.vdw = VdW[self.structure]
        except KeyError:
            print(f'{self.structure} has not been implemented in Van der Waals dictionary')
            Errors += 1

        if Errors > 0:
            print(f'There were {Errors} errors for nanoparticle {nanoparticle}')
            raise ValueError

    def setInwards(self, set_to) -> None:
        self.inwards = set_to

    def setDiameter(self, dia) -> None:
        self.diameter = dia
        self.radius = self.diameter / 2

    def setRadius(self, rad) -> None:
        self.radius = rad
        self.diameter = self.radius * 2

    def makeNanoparticle(self, diameter):
        if type(self.lattice_length) == float:
            dim = int(np.ceil(diameter/self.lattice_length)) * 2
        else:
            dim1 = int(np.ceil(diameter/self.lattice_length[0]))
            dim2 = int(np.ceil(diameter/self.lattice_length[1]))

        if self.lattice_type == 'FCC':
            self.atoms = fcc110(symbol=self.atomtypes[0], size=(dim, dim, dim), a=self.lattice_length, orthogonal=True)
        elif self.lattice_type == 'BCC':
            self.atoms = bcc110(symbol=self.atomtypes[0], size=(dim, dim, 4.5*dim), a=self.lattice_length, orthogonal=True)
        elif self.lattice_type == 'HCP':
            self.atoms = hcp0001(symbol=self.atomtypes[0], size=(dim, dim, 3*dim), a=self.lattice_length[0], c=self.lattice_length[1], orthogonal=True)
        elif self.lattice_type == 'Diamond':
            self.atoms = diamond100(symbol=self.atomtypes[0], size=(dim, dim, 3*dim), a=self.lattice_length, orthogonal=True)
        elif self.lattice_type == 'Rutile':
            a = self.lattice_length[0]
            c = self.lattice_length[1]
            self.atoms = crystal([self.atomtypes[0], self.atomtypes[1]], basis=[(0, 0, 0), (0.3, 0.3, 0.0)], spacegroup=136, cellpar=[a, a, c, 90, 90, 90], size=(2, dim1, dim2))
        elif self.lattice_type == 'Rocksalt':
            a = self.lattice_length
            self.atoms = crystal([self.atomtypes[0], self.atomtypes[1]], basis=[(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225, cellpar=[a, a, a, 90, 90, 90], size=(2, dim1, dim2))
        elif self.lattice_type == 'Skutterudite':
            a = 9.04
            self.atoms = crystal([self.atomtypes[0], self.atomtypes[1]], basis=[(0.25, 0.25, 0.25), (0.0, 0.335, 0.158)], spacegroup=204, cellpar=[a, a, a, 90, 90, 90], size=(2, dim1, dim2))
        else:
            print(f'The crystal structure type {self.lattice_type} is not implemented')
            raise ValueError

        self.atoms_symbols = self.atoms.get_chemical_symbols()
        self.atoms_pos = self.atoms.get_positions()

        return self.atoms_symbols, self.atoms_pos

    def makeSandwich(self, molecule: list, molecule_symbols: list) -> None:

        vdW_left = ci.AtomicInformation(molecule_symbols[molecule.index_min()]).VdW()
        vdW_right = ci.AtomicInformation(molecule_symbols[molecule.index_max()]).VdW()

        vacuum_dist_left = vdW_left + self.vdw
        vacuum_dist_right = vdW_right + self.vdw

        self.distance = molecule.__len__() + vacuum_dist_left + vacuum_dist_right

        if type(self.lattice_length) == float:
            left, right, left_symbols, right_symbols = self.spherical(molecule.max(), molecule.min(), vacuum_dist_left, vacuum_dist_right)
        else:
            left, right, left_symbols, right_symbols = self.slab(molecule.max(), molecule.min(), vacuum_dist_left, vacuum_dist_right)

        self.left = left
        self.right = right
        self.left_symbols = left_symbols
        self.right_symbols = right_symbols

        return self.left, self.right, self.left_symbols, self.right_symbols

    def make_sphere(self, atom_array, sphere_radii) -> list:
        shiftx = max(atom_array[:, 0])
        shifty = max(atom_array[:, 1])
        shiftz = max(atom_array[:, 2])
        atom_array[:, 0] = atom_array[:, 0] - 0.5*shiftx
        atom_array[:, 1] = atom_array[:, 1] - 0.5*shifty
        atom_array[:, 2] = atom_array[:, 2] - 0.5*shiftz
        sphere = np.empty((0, 3))
        for i in range(len(atom_array)):
            if sum(atom_array[i, :]**2)**0.5 < sphere_radii:
                sphere = np.vstack([sphere, atom_array[i, :]])
        return sphere

    def spherical(self, mol_max, mol_min, vacuum_dist_left, vacuum_dist_right):
        sphere_coor = self.make_sphere(self.atoms_pos, self.radius)

        left_hemisphere = np.empty((0, 3))
        left_symbols = np.array([])
        for i in range(len(sphere_coor)):
            if sphere_coor[i, 0] < 0:
                left_hemisphere = np.vstack([left_hemisphere, sphere_coor[i, :]])
                left_symbols = np.append(left_symbols, self.atoms_symbols[i])
        right_hemisphere = np.copy(left_hemisphere) * np.array([-1, 1, 1])
        right_symbols = left_symbols
        if self.inwards:
            left_hemisphere[:, 0] = left_hemisphere[:, 0] + mol_max + self.radius + vacuum_dist_left
            right_hemisphere[:, 0] = right_hemisphere[:, 0] + mol_min - self.radius - vacuum_dist_right
        else:
            left_hemisphere[:, 0] = left_hemisphere[:, 0] + mol_min - vacuum_dist_left
            right_hemisphere[:, 0] = right_hemisphere[:, 0] + mol_max + vacuum_dist_right
        return left_hemisphere, right_hemisphere, left_symbols, right_symbols

    def slab(self, mol_max, mol_min, vacuum_dist_left, vacuum_dist_right):
        shifty = max(self.atoms_pos[:, 1])
        shiftz = max(self.atoms_pos[:, 2])
        self.atoms_pos[:, 1] -= 0.5*shifty
        self.atoms_pos[:, 2] -= 0.5*shiftz
        left_slab = np.copy(self.atoms_pos)
        right_slab = np.copy(self.atoms_pos) * np.array([-1, 1, 1])
        left_symbols = np.copy(self.atoms_symbols)
        right_symbols = np.copy(self.atoms_symbols)
        left_slab[:, 0] = left_slab[:, 0] + mol_min - vacuum_dist_left - 2 * self.lattice_length[0]
        right_slab[:, 0] = right_slab[:, 0] + mol_max + vacuum_dist_right + 2 * self.lattice_length[0]
        return left_slab, right_slab, left_symbols, right_symbols


class Molecule():
    def __init__(self, molecule: list) -> None:
        self.molecule = molecule

    def __len__(self) -> float:
        return self.molecule[:, 0].max() - self.molecule[:, 0].min()

    def min(self) -> float:
        return min(self.molecule[:, 0])

    def max(self) -> float:
        return max(self.molecule[:, 0])

    def index_min(self) -> float:
        return np.argmin(self.molecule[:, 0])

    def index_max(self) -> float:
        return np.argmax(self.molecule[:, 0])

    def get_rotation_matrix(self, point_for_rot, direc_vec, theta) -> None:
        # Takes a point on the rotational axis, the normalized direction vector and the angle
        # The rotation matrix is 4x4 with the fourth dimension translating the origin to (a, b, c) and back
        a, b, c = point_for_rot
        u, v, w = direc_vec
        u2, v2, w2 = u*u, v*v, w*w
        cosT = np.cos(theta)
        sinT = np.sin(theta)
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

        self.rotationMatrix = np.array([m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34])

    def rotated_point(self, point) -> None:
        #Essentially matrix multiplication
        p = np.zeros(3)
        x, y, z = point
        m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34 = self.rotationMatrix
        p[0] = m11*x + m12*y + m13*z + m14
        p[1] = m21*x + m22*y + m23*z + m24
        p[2] = m31*x + m32*y + m33*z + m34

        self.rotationElements = p

    def rotateMolecule(self, alignmentAtom) -> None:
        rotxyz = np.empty((0, 3))
        for atom in self.molecule:
            self.rotated_point(atom)
            rotxyz = np.vstack([rotxyz, np.array([self.rotationElements[0], self.rotationElements[1], self.rotationElements[2]])])
        for i in range(3):
            rotxyz[:, i] -= rotxyz[alignmentAtom, i]
        rotxyz[:, 0] -= np.amax(rotxyz[:, 0])

        self.molecule = rotxyz

    def moveMolecule(self, translation) -> None:
        self.molecule += translation

