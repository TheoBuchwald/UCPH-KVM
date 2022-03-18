
from . import chemical_information as ci
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

    def setRIBasis(self, RIbasis: str) -> None:
        """Set basis set

        Args:
            basis (str): Basis set you wish to use
        """
        if self.checkBasis(RIbasis):
            self.RIbasis = RIbasis
        else:
            print(f'''{RIbasis} was not found in the standard library of {self.program}:
Checking Basis Set Exchange''')
            if ci.BasisSet.CheckBasisSet(self.program, RIbasis):
                self.BSE = True
                self.RIbasis = RIbasis
                print('Basis set found on the Basis Set Exchange\n')
            else:
                print(f'''{RIbasis} does not seem to be available for {self.program}:
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
