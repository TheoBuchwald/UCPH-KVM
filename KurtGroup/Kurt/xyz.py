
import Kurt.chemical_information as ci
# from . import chemical_information as ci
import numpy as np
import json
import os


def checkBasis(program: str, basis: str) -> bool:
    """Checks whether a basis set exits for a given program

    Args:
        program (str): Program to check
        functional (str): Basis set to check

    Returns:
        bool: Returns true or false
    """
    BASE_DIR = os.path.dirname(__file__)

    with open(f'{BASE_DIR}/Info.json') as info:
        info = json.load(info)
    if basis in info[program]['basis sets']:
        return True
    elif ci.BasisSet.CheckBasisSet(basis):
        return True
    return False

def checkFunctional(program: str, functional: str) -> bool:
    """Checks whether a functional exits for a given program

    Args:
        program (str): Program to check
        functional (str): Functional to check

    Returns:
        bool: Returns true or false
    """
    BASE_DIR = os.path.dirname(__file__)

    with open(f'{BASE_DIR}/Info.json') as info:
        info = json.load(info)

    return functional in info[program]['functionals']


class xyz_to:
    def __init__(self, program: str, xyz_file: str = None, memory: int = 8, ncpus: int = 8, calculation_type: str = 'Opt') -> None:
        """xyz_to ismade to make it easy to generate input files for quantum chemistry programs

        Args:
            program (str): The program to generate the input file for
            xyz_file (str, optional): The xyz file used for the input file. Defaults to None.
        """
        self.BSE = False
        self.__base_dir__ = os.path.dirname(__file__)
        self.program = program
        self.basis = 'pc-1'
        self.method = 'DFT'
        self.functional = 'B3LYP'
        self.mem = memory
        self.ncpus = ncpus
        self.calc = calculation_type

        with open(f'{self.__base_dir__}/Info.json') as info:
            self.info = json.load(info)

        if xyz_file:
            self.readXYZ(xyz_file)

    def setBasis(self, basis: str) -> None:
        """Set basis set

        Args:
            basis (str): Basis set you wish to use
        """
        if self.checkBasis(basis):
            self.basis = basis
        else:
            print(f'''{basis} was not found in the standard library of {self.program}:
Checking Basis Set Exchange''')
            if ci.BasisSet.CheckBasisSet(self.program, basis):
                self.BSE = True
                self.basis = basis
                print('Basis set found on the Basis Set Exchange\n')
            else:
                print(f'''{basis} does not seem to be available for {self.program}:
Upper and lower case characters ARE important, so you may wish to check that they are written correctly
Exiting program''')
                exit()

    def checkBasis(self, basis: str) -> bool:
        """Checks if the basis set exists for the given program

        Args:
            basis (str): Basis set to check

        Returns:
            bool: Returns true or false
        """
        return basis.upper() in self.info[self.program.upper()]['BASIS SETS']

    def setMethod(self, method: str, functional: str = None) -> None:
        """Set method (HF, MP2, DFT and so on...)

        Args:
            method (str): Method you wish to use
            functional (str, optional): Functional incase of DFT
        """
        self.method = method.upper()
        if self.method.upper() == 'DFT':
            self.setFunctional(functional)

    def setFunctional(self, functional: str) -> None:
        """Functional to use. Only relevant if using DFT method

        Args:
            functional (str): Functional uo wish to use
        """
        if self.checkFunctional(functional):
            self.functional = functional
        else:
            print(f'''{functional} does not seem to be available for {self.program}:
Exiting program''')
            exit()

    def checkFunctional(self, functional: str) -> bool:
        """Checks if the functional exists for the program

        Args:
            functional (_type_): The functional to check

        Returns:
            bool: Returns true or false
        """
        return functional.upper() in self.info[self.program.upper()]['FUNCTIONALS']

    def readXYZ(self, xyz_file: str) -> None:
        """Reads the xyz file

        Args:
            xyz_file (str): xyz file to read
        """
        self.filename = xyz_file
        with open(xyz_file, "r") as xyz:
            self.xyz_file = xyz.readlines()

    def processXYZ(self) -> None:
        """Read through the xyz file and store the contents in self.atoms
        """
        if not self.xyz_file:
            print('''xyz file has not been loaded (read):
Exiting program''')
            exit()
        self.atoms = np.empty((len(self.xyz_file[2:]), 4), dtype=object)

        for linenumber, line in enumerate(self.xyz_file[2:]):
            self.atoms[linenumber] = np.array(line.split())

    def generateDaltonInputFileText(self, charge: int) -> str:
        """Makes the text for a Dalton input file

        Args:
            charge (int): The charge of the molecule

        Returns:
            str: Returns the file text
        """
        if not hasattr(self, 'atoms'):
            print(f'''You need to process the xyz file using the class function processXYZ before generating the file:
Exiting program''')
            exit()

        self.input_filename = self.filename.replace('.xyz', '.mol')

        unique_atoms = set(self.atoms[:,0])
        self.filetext = f'''ATOMBASIS
./{self.filename}

Atomtypes={len(unique_atoms)} Charge={charge} NoSymmetry Angstrom
'''

        for unique_atom in unique_atoms:
            count = list(self.atoms[:,0]).count(unique_atom)
            if not self.BSE:
                self.filetext += f'  {ci.AtomicInformation(unique_atom).atomnr():.4f}     {count} Bas={self.basis}\n'
            else:
                BasisSet = ci.BasisSet()
                try:
                    basis_mol = BasisSet.AtomBasisSet(self.program, self.basis, unique_atom, SupressHeader=True)
                except RuntimeError:
                    print(f'''Failed to get basis set from BSE. Please check the spelling, upper-/lowercase IS important
The problem may also be that the basis set does not exist for {unique_atom}''')
                    exit()

                BlockTypes = ['s functions', 'p functions', 'd functions', 'f functions', 'g functions', 'h functions', 'i functions', 'j functions', 'k functions']
                blocks = 0
                for j in BlockTypes:
                    if j in basis_mol:
                        blocks += 1
                Block = f'{blocks}'
                for j in BlockTypes[:blocks]:
                    Block += f' {basis_mol.count(j)}'
                self.filetext += f'  Charge={ci.AtomicInformation(unique_atom).atomnr():.4f}     Atoms={count}     Blocks={Block}\n'
                basis_mol = basis_mol.replace('H','').split('\n')[5:-2]
                basis_mol = [i for i in basis_mol if 'functions' not in i]
            for j, atom in enumerate(self.atoms[:,0]):
                if atom == unique_atom:
                    self.filetext += f'{atom: <2} {self.atoms[j,1]: >14.9} {self.atoms[j,2]: >19.9} {self.atoms[j,3]: >19.9}\n'
            if self.BSE:
                for j in basis_mol:
                    self.filetext += f'{j}\n'
        return self.filetext

    def generateGaussianInputFileText(self, charge: int) -> str:
        """Makes the text for a Gaussian input file

        Args:
            charge (int): The charge of the molecule

        Returns:
            str: Returns the file text
        """
        if not hasattr(self, 'atoms'):
            print(f'''You need to process the xyz file using the class function processXYZ before generating the file:
Exiting program''')
            exit()

        basis_name = self.basis
        if self.BSE:
            BasisSet = ci.BasisSet()
            try:
                basis_mol = BasisSet.GenerateBasisSet('gaussian94', self.basis, self.atoms[:,0], SupressHeader=True)
            except RuntimeError:
                print("Failed to get basis set from BSE. Please check the spelling, upper-/lowercase is not important")
                exit()
            basis_name = self.basis
            self.basis = "GEN"

        filename_no_ext = self.filename.replace('.xyz', '')
        self.input_filename = f'{filename_no_ext}.com'

        if charge % 2 == 0:
            multiplicity = 1
        else:
            multiplicity = 2

        if self.method.upper() == 'DFT':
            method = self.functional
        else:
            method = self.method

        self.filetext = f'''%chk={filename_no_ext}.chk
%mem={self.mem}GB
%nprocshared={self.ncpus}
# {self.calc} {method}/{self.basis}

'''
        if self.BSE:
            self.filetext += f'{filename_no_ext}.xyz - {basis_name}\n'
        else:
            self.filetext += f'{filename_no_ext}.xyz\n'

        self.filetext += f'\n{charge} {multiplicity}\n'

        for atom in self.atoms:
            self.filetext += f'{atom[0]} {atom[1]: >14.9} {atom[2]: >19.9} {atom[3]: >19.9}\n'

        self.filetext += '\n'

        if self.BSE:
            self.filetext += '****\n'
            self.filetext += f'{basis_mol}\n'
            self.filetext += '\n'

        self.basis = basis_name

        return self.filetext

    def writeInputfile(self):
        """Writes the file text to an input file with similar name as the xyz file
        """
        with open(self.input_filename, 'w') as input_file:
            input_file.writelines(self.filetext)
