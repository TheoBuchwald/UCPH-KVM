
import subprocess
import numpy as np
from typing import List
from chemical_information import AtomicInformation

def Forward_search_last(file: str, text: str, error: str, quiet: bool = False) -> int:
    """Searches from the beggining of the file given to the end where it returns the linenumber of the last occurence

    Args:
        file (str): The file to search in
        text (str): The text string to search for
        error (str): If no occurences were found it will print 'No [error] could be found in [file]
        err (bool, optional): Whether or not to print error message if no occurences are found. Defaults to True.

    Returns:
        (int): Linenumber of last occurence
    """
    ps1 = subprocess.run(['grep', '-nT', text, file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = subprocess.run(['tail', '-n1'], input=ps1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    res = out.stdout
    if len(res) == 0:
        if not quiet:
            with open("collect_data.log", "a") as logfile:
                logfile.write(f'No {error} could be found in {file}\n')
        return 'NaN'
    res = res.split()[0]
    res = str(res).split(':')
    return int(res[0].replace('b\'','').replace('\'','')) - 1

def Forward_search_after_last(file: str, text1: str, text2: str, lines: int, error: str, quiet: bool = False) -> int:
    """Searches from beggining of file for last occurence of [text1] and in the following [lines] after for [text2]

    Args:
        file (str): File to search in
        text1 (str): From the last occurence of this this function will search
        text2 (str): This is what will be found in the lines following [text1]
        lines (int): How many lines after [text1] should the function search for [text2]
        error (str): If no occurences were found it will print 'No [error] could be found in [file]
        err (bool, optional): Whether or not to print error message if no occurences are found. Defaults to True.. Defaults to True.

    Returns:
        (int): Linenumber of [text2] occurence
    """
    ps1 = subprocess.run(['grep', '-nTA', f'{lines}', text1, file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps2 = subprocess.run(['tail', '-n', f'{lines + 1}'], input=ps1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = subprocess.run(['grep', text2], input=ps2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    res = out.stdout
    if len(res) == 0:
        if not quiet:
            with open("collect_data.log", "a") as logfile:
                logfile.write(f'No {error} could be found in {file}\n')
        return 'NaN'
    res = res.split()[0]
    res = str(res).split('-')
    return int(res[0].replace('b\'','').replace('\'','')) - 1

def Backward_search_last(file: str, text: str, filelength: int, error: str, quiet: bool = False) -> int:
    """Finds the last occurence of a text string in a file by searching from the end

    Args:
        file (str): File to search in
        text (str): Text string to look for
        filelength (int): The length of the file
        error (str): If no occurences were fount it will print 'No [error] could be found in [file]
        quiet (bool, optional): Whether or not to print error message. Defaults to False.

    Returns:
        (int): Linenumber of last occurence
    """
    ps1 = subprocess.run(['tac', file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = subprocess.run(['grep', '-nTm1', text], input=ps1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    res = out.stdout
    if len(res) == 0:
        if not quiet:
            with open("collect_data.log", "a") as logfile:
                logfile.write(f'No {error} could be found in {file}\n')
        return 'NaN'
    res = res.split()[0]
    res = str(res).split(':')
    return filelength - int(res[0].replace('b\'','').replace('\'',''))

def Forward_search_first(file: str, text: str, error: str, quiet: bool = False) -> int:
    """Searches from beginning of file and finds the first occurence of [text]

    Args:
        file (str): File to search in
        text (str): Text to look for
        error (str): If no occurences were found it will print 'No [error] could be found in [file]
        err (bool, optional): Whether or not to print error message if no occurences are found. Defaults to True.. Defaults to True.

    Returns:
        (int): Linenumber of first occurence
    """
    out = subprocess.run(['grep', '-nTm1', text, file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    res = out.stdout
    if len(res) == 0:
        if not quiet:
            with open("collect_data.log", "a") as logfile:
                logfile.write(f'No {error} could be found in {file}\n')
        return 'NaN'
    res = res.split()[0]
    res = str(res).split(':')
    return int(res[0].replace('b\'','').replace('\'','')) - 1

def Forward_search_all(file: str, text: str, error: str, quiet: bool = False) -> list:
    """Searches from beggining of file to end of file finding all occurences of [text]

    Args:
        file (str): File to search in
        text (str): Text to look for
        error (str): If no occurences were found it will print 'No [error] could be found in [file]
        err (bool, optional): Whether or not to print error message if no occurences are found. Defaults to True.. Defaults to True.

    Returns:
        (list): List of the linenumbers of all occurences
    """
    ps1 = subprocess.run(['grep', '-nT', text, file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = subprocess.run(['awk', '{print $1}'], input=ps1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    res = out.stdout
    if len(res) == 0:
        if not quiet:
            with open("collect_data.log", "a") as logfile:
                logfile.write(f'No {error} could be found in {file}\n')
        return 'NaN'
    res = str(res).split('\\n')
    return [int(val.replace('b\'','').replace('\'','').replace(':','')) - 1 for val in res[:-1]]

def CheckForOnlyNans(array: list) -> bool:
    """Function for checking if an array is fille only with the value 'NaN'

    Args:
        array (list): Array that should be looked through

    Returns:
        (bool): Returns a False/True based on whether or not the array the array only consists of 'NaN'
    """
    for i in array:
        if i != 'NaN':
           return False
    return True

def WriteToFile(filename : str, lines : list) -> None:
    """ Function for writing out to a file

    Args:
        lines (list) : List of lines that should be written in the file
    """
    with open(filename,'w') as wrt:
        wrt.writelines(lines)

def GenerateXYZ(lines : list, filename : str , start : int, end : int, lab_loc : int, transform : bool = False) -> None:
    """ Function for generating and writing out XYZ file from imput

    Args:
        lines (list): Lines in an input file
        filename (str): Filename for geometry file
        start, end (int): Starting and ending linenumber of the final geometry in the file
        lab_loc (int): Location of label in line
        transform (bool): Transforms atomic number into label, if needed
    """
    lines_to_add = []
    lines_to_add.append(str(end-(start))+ '\n')
    lines_to_add.append('\n')
    for line in lines[start:end]:
        words = line.split()
        if transform:
            atm = AtomicInformation(int(words[lab_loc]))
            lines_to_add.append(''.join([atm.atom.ljust(2),' ',words[-3].rjust(10),' ', words[-2].rjust(15), ' ',words[-1].rjust(15) ,'\n']))
        else:
            lines_to_add.append(''.join([words[lab_loc].ljust(2),' ',words[-3].rjust(10),' ', words[-2].rjust(15), ' ',words[-1].rjust(15) ,'\n']))
    WriteToFile(filename,lines_to_add)


class OutputType:
    def __init__(self, filename: str, *, Quiet: bool = False, Temperature: float = 298.15):
        self.filename = filename

        with open(self.filename,'r') as read:
            lines = read.readlines()[:100]

        AMS = False
        GAUSS = False
        for line in lines:
            if "Amsterdam Modeling Suite (AMS)" in line:
                AMS = True
            if "Gaussian, Inc.  All Rights Reserved." in line:
                GAUSS = True

        # The output file is determined to be of one of the following types

        # File type = ORCA
        if '* O   R   C   A *' in lines[4]:
            self.extract = OrcaExtract(self.filename, Quiet=Quiet, Temperature=Temperature)
            self.input = 'ORCA'

        # File type = DALTON
        elif '*************** Dalton - An Electronic Structure Program ***************' in lines[3]:
            self.extract = DaltonExtract(self.filename, Quiet=Quiet, Temperature=Temperature)
            self.input = 'DALTON'

        # File type = DIRAC
        elif 'DIRAC' in lines[0]:
            self.extract = DiracExtract(self.filename, Quiet=Quiet, Temperature=Temperature)
            self.input = 'DIRAC'

        # File type = GAUSSIAN
        elif GAUSS:
            self.extract = GaussianExtract(self.filename, Quiet=Quiet, Temperature=Temperature)
            self.input = 'GAUSSIAN'

        # File type = LSDALTON
        elif '**********  LSDalton - An electronic structure program  **********' in lines[2]:
            self.extract = LSDaltonExtract(self.filename, Quiet=Quiet, Temperature=Temperature)
            self.input = 'LSDALTON'

        # File type = VELOXCHEM
        elif 'VELOXCHEM' in lines[2]:
            self.extract = VeloxExtract(self.filename, Quiet=Quiet, Temperature=Temperature)
            self.input = 'VELOXCHEM'

        # File type = AMS
        elif AMS:
            self.extract = AMSExtract(self.filename, Quiet=Quiet, Temperature=Temperature)
            self.input = 'Amsterdam Modeling Suite'

        # File type = Q-Chem
        elif 'Welcome to Q-Chem' in lines[0]:
            self.extract = QChemExtract(self.filename, Quiet=Quiet, Temperature=Temperature)
            self.input = 'Q-Chem'

        # File type not implemented
        else:
            self.extract = UnknownExtract()
            self.input = 'Unknown'
            if not Quiet:
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"The output file {self.filename} is not of a known format\n")

    def getEnergy(self) -> float:
        try:
            return self.extract.tot_energy
        except AttributeError:
            try:
                self.extract._Energy()
                return self.extract.tot_energy
            except AttributeError: ...

    def getZeroPointVibrationalEnergy(self) -> float:
        try:
            return self.extract.zpv
        except AttributeError:
            try:
                self.extract._ZPV()
                return self.extract.zpv
            except AttributeError: ...

    def getEnthalpy(self) -> float:
        try:
            return self.extract.enthalpy
        except AttributeError:
            try:
                self.extract.tot_energy
            except AttributeError:
                try:
                    self.extract._Energy()
                except AttributeError: return
            try:
                self.extract.freq
            except AttributeError:
                try:
                    self.extract._Frequencies()
                except AttributeError: return
            try:
                self.extract._Enthalpy()
                return self.extract.enthalpy
            except AttributeError: ...

    def getEntropy(self) -> float:
        try:
            return self.entropy
        except AttributeError:
            try:
                self.extract.freq
            except AttributeError:
                try:
                    self.extract._Frequencies()
                except AttributeError: return
            try:
                self.extract._Entropy()
                return self.extract.entropy
            except AttributeError: ...

    def getGibbsFreeEnergy(self) -> float:
        try:
            return self.extract.gibbs
        except AttributeError:
            try:
                self.extract.tot_energy
            except AttributeError:
                try:
                    self.extract._Energy()
                except AttributeError: return
            try:
                self.extract.freq
            except AttributeError:
                try:
                    self.extract._Frequencies()
                except AttributeError: return
            try:
                self.extract.enthalpy
            except AttributeError:
                try:
                    self.extract._Enthalpy()
                except AttributeError: return
            try:
                self.extract.entropy
            except AttributeError:
                try:
                    self.extract._Entropy()
                except AttributeError: return
            try:
                self.extract._Gibbs()
                return self.extract.gibbs
            except AttributeError: ...

    def getDipoleMoment(self) -> List[float]:
        try:
            return [self.extract.dipolex, self.extract.dipoley, self.extract.dipolez, self.extract.total_dipole]
        except AttributeError:
            try:
                self.extract._Dipole_moments()
                return [self.extract.dipolex, self.extract.dipoley, self.extract.dipolez, self.extract.total_dipole]
            except AttributeError: ...

    def getPolarizability(self) -> List[float]:
        try:
            return [self.extract.polx, self.extract.poly, self.extract.polz, self.extract.iso_polar]
        except AttributeError:
            try:
                self.extract._Polarizabilities()
                return [self.extract.polx, self.extract.poly, self.extract.polz, self.extract.iso_polar]
            except AttributeError: ...

    def getExcitationEnergies(self) -> List[float]:
        try:
            return self.extract.exc_energies
        except AttributeError:
            try:
                self.extract._Excitation_energies()
                return self.extract.exc_energies
            except AttributeError: ...

    def getOscillatorStrengths(self) -> List[float]:
        try:
            return self.extract.osc_strengths
        except AttributeError:
            try:
                self.extract.exc_energies
            except AttributeError:
                try:
                    self.extract._Excitation_energies()
                except AttributeError: return
            try:
                self.extract._Oscillator_strengths()
                return self.extract.osc_strengths
            except AttributeError: ...

    def getRotationalStrengths(self) -> List[float]:
        try:
            return self.extract.rot_strengths
        except AttributeError:
            try:
                self.extract.exc_energies
            except AttributeError:
                try:
                    self.extract._Excitation_energies()
                except AttributeError: return
            try:
                self.extract._Rotational_strengths()
                return self.extract.rot_strengths
            except AttributeError: ...

    def getFrequencies(self) -> List[float]:
        try:
            return self.extract.freq
        except AttributeError:
            try:
                self.extract._Frequencies()
                return self.extract.freq
            except AttributeError: ...

    def getPartitionFunction(self) -> float:
        try:
            return self.extract.qTotal
        except AttributeError:
            try:
                self.extract.freq
            except AttributeError:
                try:
                    self.extract._Frequencies()
                except AttributeError: return
            try:
                self.extract._PartitionFunctions()
                return self.extract.qTotal
            except AttributeError: ...

    def getCPUTime(self) -> List[float]:
        try:
            return [self.extract.total_cpu_time, self.extract.wall_cpu_time]
        except AttributeError:
            try:
                self.extract._CPUS()
                return [self.extract.total_cpu_time, self.extract.wall_cpu_time]
            except AttributeError: ...

    def getOptimizedGeometry(self) -> None:
        self.extract._Optimized_Geometry()


class Constants:
    def __init__(self) -> None:
        self.ev_to_au = 0.036749405469679
        self.inv_cm_to_au = 1/219474.63068
        self.trans_const_fac = 1.5625517342018425307E+22 #Molar value assuming 1 bar standard pressure
        self.rot_lin_const = 20.83661793 #Assuming rigid, linear rotor and T>>Rotational temperature and rotational constant in GHz
        self.rot_poly_const = 168.5837766 #Assuming rigid, polyatomic rotor and T>>Rotational temperature and rotational constant in GHz
        self.vib_const = 3.157750419E+05 #Assuming harmonic oscillator and frequency in au
        self.gas_constant = 8.31446261815324E-03 # In kJ/(mol*K)
        self.s_trans_const = 0.3160965065 #Assuming 1 bar standard pressure and molar
        self.au_to_kJmol = 2625.4996394799
        self.bohr_to_ao = 0.529177249
        self.debye_to_au = 0.393456


class UnknownExtract:
    def __init__(self) -> None: ...


class VeloxExtract:
    def __init__(self, filename: str, *, Quiet: bool = False, Temperature: float = 298.15) -> None:
        self.filename = filename
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self) -> None:
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _Energy(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Total Energy', 'final energy', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-2])
            return
        self.tot_energy = 'NaN'

    def _Dipole_moments(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Ground-State Dipole Moment', 'dipole moment', quiet=self.quiet)
        if isinstance(linenumber, int):
            linenumber += 3
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[linenumber].split()[-4]), float(self.lines[linenumber+1].split()[-4]), float(self.lines[linenumber+2].split()[-4]), float(self.lines[linenumber+3].split()[-4])
            return
        self.dipolex, self.dipoley, self.dipolez, self.total_dipole = 'NaN'

    #NB! Only static polarizability considered
    def _Polarizabilities(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Polarizability (w=0.0000)', 'polarizability', quiet=self.quiet)
        if isinstance(linenumber, int):
            #Need to perform diagonalization
            PolarizabilityTensor = np.array([[float(self.lines[linenumber+3].split()[1]), float(self.lines[linenumber+3].split()[2]), float(self.lines[linenumber+3].split()[3])],
                                             [float(self.lines[linenumber+4].split()[1]), float(self.lines[linenumber+4].split()[2]), float(self.lines[linenumber+4].split()[3])],
                                             [float(self.lines[linenumber+5].split()[1]), float(self.lines[linenumber+5].split()[2]), float(self.lines[linenumber+5].split()[3])]])
            PolarizabilityEigenvalues = np.linalg.eigh(PolarizabilityTensor)[0]
            self.polx, self.poly, self.polz = PolarizabilityEigenvalues
            self.iso_polar = PolarizabilityEigenvalues.mean()
            return
        self.polx = self.poly = self.polz = self.iso_polar = 'NaN'

    def _Optimized_Geometry(self) -> None:
        start = Forward_search_last(self.filename, 'Molecular Geometry', 'geometry', quiet=self.quiet)
        if start != "NaN":# and end != "NaN":
            start += 5
            for i, line in enumerate(self.lines[start:]):
                if len(line.strip()) == 0:
                    end = start + i
                    break
            #Offset for going into actual coordinate list
            #Which position in the line is the atom label / number at
            label_location = 0
            OptGeomFilename = self.filename[:-4] + "_opt.xyz"
            GenerateXYZ(self.lines, OptGeomFilename, start, end, label_location)
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write("Final geometry has been saved to " + OptGeomFilename + "\n")


class AMSExtract:
    def __init__(self, filename: str, *, Quiet: bool = False, Temperature: float = 298.15) -> None:
        self.filename = filename
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()
        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self) -> None:
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _Energy(self) -> None:
        linenumber = Forward_search_last(self.filename, "Energy (hartree)", "final energy", quiet=self.quiet)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        self.tot_energy = 'NaN'

    def _Dipole_moments(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Dipole Moment', 'dipole moment', quiet=self.quiet)
        if isinstance(linenumber, int):
            linenumber += 3
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[linenumber].split()[-3])*self.constants.debye_to_au, float(self.lines[linenumber].split()[-2])*self.constants.debye_to_au, float(self.lines[linenumber].split()[-1])*self.constants.debye_to_au, float(self.lines[linenumber+1].split()[-1])*self.constants.debye_to_au
            return
        self.dipolex, self.dipoley, self.dipolez, self.total_dipole = 'NaN'

    def _Optimized_Geometry(self) -> None:
        start = Forward_search_last(self.filename, 'Formula:', 'geometry', quiet=self.quiet)
        if start != "NaN":
            start += 3
            #Offset for going into actual coordinate list
            for i, line in enumerate(self.lines[start:]):
                if len(line.strip()) == 0:
                    end = start + i
                    break
            #Which position in the line is the atom label / number at
            label_location = 1
            OptGeomFilename = self.filename[:-4] + "_opt.xyz"
            GenerateXYZ(self.lines, OptGeomFilename, start, end, label_location)
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write("Final geometry has been saved to " + OptGeomFilename + "\n")


class GaussianExtract:
    def __init__(self, filename: str, *, Quiet: bool = False, Temperature: float = 298.15) -> None:
        self.filename = filename
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self) -> None:
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _CPUS(self) -> None:
        linenumber = Backward_search_last(self.filename, 'Job cpu time:', self.end, 'CPU time', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.total_cpu_time = float(self.lines[linenumber].split()[3])*24*60 + float(self.lines[linenumber].split()[5])*60 + float(self.lines[linenumber].split()[7]) + float(self.lines[linenumber].split()[9])/60
            self.wall_cpu_time = float(self.lines[linenumber+1].split()[2])*24*60 + float(self.lines[linenumber+1].split()[4])*60 + float(self.lines[linenumber+1].split()[6]) + float(self.lines[linenumber+1].split()[8])/60
            return
        self.total_cpu_time = 'NaN'
        self.wall_cpu_time = 'NaN'

    def _Energy(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Sum of electronic and zero-point Energies=', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1]) - float(self.lines[linenumber-4].split()[-2])
            return
        linenumber = Forward_search_last(self.filename, 'SCF Done:', 'final energy', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[4])
            return
        self.tot_energy = 'NaN'

    def _ZPV(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Sum of electronic and zero-point Energies=', 'ZPV energy', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.zpv = float(self.lines[linenumber].split()[-1])
            return
        self.zpv = 'NaN'

    def _Dipole_moments(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Electric dipole moment (input orientation):', 'dipole moments', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[linenumber+4].split()[1].replace('D','E')), float(self.lines[linenumber+5].split()[1].replace('D','E')), float(self.lines[linenumber+6].split()[1].replace('D','E')), float(self.lines[linenumber+3].split()[1].replace('D','E'))
            return
        self.dipolex = self.dipoley = self.dipolez = self.total_dipole = 'NaN'

    def _Polarizabilities(self) -> None:
        linenumber = ['NaN', 'NaN', 'NaN', 'NaN']
        searchwords = [' xx ', ' yy ', ' zz ', ' iso ']
        for i in range(len(searchwords)):
            linenumber[i] = Forward_search_after_last(self.filename, 'Dipole polarizability, Alpha (input orientation).', searchwords[i], 15, 'polarizabilities', quiet=self.quiet)
        if linenumber != ['NaN', 'NaN', 'NaN', 'NaN']:
            self.polx, self.poly, self.polz, self.iso_polar = float(self.lines[linenumber[0]].split()[1].replace('D','E')), float(self.lines[linenumber[1]].split()[1].replace('D','E')), float(self.lines[linenumber[2]].split()[1].replace('D','E')), float(self.lines[linenumber[3]].split()[1].replace('D','E'))
            return
        self.polx = self.poly = self.polz = self.iso_polar = 'NaN'

    def _Frequencies(self) -> None:
        self.freq = []
        linenumbers = Forward_search_all(self.filename, 'Frequencies --', 'frequencies', quiet=self.quiet)
        if isinstance(linenumbers, list):
            for i in linenumbers:
                for j in self.lines[i].split()[2:]:
                    self.freq.append(float(j)* self.constants.inv_cm_to_au)
        if len(self.freq) == 0:
            self.freq = ['NaN']

    def _Excitation_energies(self) -> None:
        self.exc_energies = []
        linenumber = Forward_search_last(self.filename, 'Excitation energies and oscillator strengths:', 'excitation energies', quiet=True)
        if isinstance(linenumber, int):
            linenumbers = Forward_search_all(self.filename, 'Excited State', 'excitation energies', quiet=self.quiet)
            linenumbers = [i for i in linenumbers if i > linenumber]
            for i in linenumbers:
                self.exc_energies.append(float(self.lines[i].split()[4])* self.constants.ev_to_au)
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN']

    def _Oscillator_strengths(self) -> None:
        self.osc_strengths = []
        linenumber = Forward_search_last(self.filename, 'Excitation energies and oscillator strengths:', 'oscillator strengths', quiet=True)
        if isinstance(linenumber, int):
            linenumbers = Forward_search_all(self.filename, 'Excited State', 'oscillator strengths', quiet=self.quiet)
            linenumbers = [i for i in linenumbers if i > linenumber]
            for i in linenumbers:
                for j in self.lines[i].split():
                    if 'f=' in j:
                        self.osc_strengths.append(float(j.replace('f=','')))
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN']

    def _RotationalConsts(self) -> None:
        self.rots = []
        linenumbers = Forward_search_all(self.filename, 'Rotational constants (GHZ):', 'rotational constants', quiet=self.quiet)
        for i in self.lines[linenumbers[-1]].split()[3:]:
            self.rots.append(float(i))
        self.rots = np.array(self.rots)
        self.rots = self.rots[self.rots != 0.0]

    def _Mass(self) -> None:
        self.mass = 0.0
        linenumber = Forward_search_last(self.filename, 'Molecular mass', 'molecular mass', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.mass = float(self.lines[linenumber].split()[2])

    def _SymmetryNumber(self):
        self.symnum = 0
        linenumber = Forward_search_last(self.filename, 'Rotational symmetry number', 'rotational symmetry number', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.symnum = int(self.lines[linenumber].split()[-1].replace('.',''))

    def _Multiplicity(self) -> None:
        self.multi = 0
        linenumber = Forward_search_first(self.filename, 'Multiplicity', 'multiplicity', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.multi = int(self.lines[linenumber].split()[-1])

    def _PartitionFunctions(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping partition function calculation\n")
            self.qTotal = 'NaN'
            return
        self._RotationalConsts()
        self._Mass()
        self._SymmetryNumber()
        self._Multiplicity()
        self.qT = self.constants.trans_const_fac * self.mass ** (1.5) * self.T ** (2.5)
        if len(self.rots) == 1:
            self.qR = self.constants.rot_lin_const * self.T / (self.symnum * self.rots[0])
        else:
            self.qR = self.constants.rot_poly_const * self.T ** (1.5) / ( self.symnum * np.prod(np.array(self.rots)) ** (0.5))
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.qV = np.prod(1 / (1 - np.exp( - self.constants.vib_const * realfreq / self.T)))
        self.qE = self.multi #Good approximation for most closed-shell molecules
        self.qTotal = self.qT*self.qR*self.qV*self.qE

    def _Enthalpy(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping partition function calculation\n")
            self.enthalpy = 'NaN'
        self._RotationalConsts()
        self.E_T = 3/2 * self.T * self.constants.gas_constant
        if len(self.rots) == 1:
            self.E_R = self.T * self.constants.gas_constant
        else:
            self.E_R = 3/2 * self.T * self.constants.gas_constant
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.E_V = self.constants.gas_constant * np.sum(self.constants.vib_const * realfreq * (1/2 + 1 / (np.exp(self.constants.vib_const * realfreq /  self.T ) - 1)))
        self.E_e = 0 #Good approximation for most closed-shell molecules
        self.enthalpy = (self.E_T+self.E_R+self.E_V+self.constants.gas_constant *  self.T ) / self.constants.au_to_kJmol + self.tot_energy

    def _Entropy(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping partition function calculation\n")
            self.entropy = 'NaN'
            return
        self._RotationalConsts()
        self._Mass()
        self._SymmetryNumber()
        self._Multiplicity()
        self.S_T = self.constants.gas_constant * np.log(self.constants.s_trans_const * self.mass ** 1.5 * self.T ** 2.5)
        if len(self.rots) == 1:
            self.S_R = self.constants.gas_constant * np.log(self.constants.rot_lin_const * self.T / (self.symnum * self.rots[0]))
        else:
            self.S_R = self.constants.gas_constant * (3/2 + np.log(self.constants.rot_poly_const * self.T ** (1.5) / ( self.symnum * np.prod(np.array(self.rots)) ** (0.5))))
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.S_V = self.constants.gas_constant * np.sum(self.constants.vib_const * realfreq / self.T / (np.exp(self.constants.vib_const * realfreq /  self.T ) - 1) - np.log(1-np.exp(-self.constants.vib_const * realfreq /  self.T )))
        self.S_E = self.constants.gas_constant * np.log(self.multi) #Good approximation for most closed-shell molecules
        self.entropy = self.S_T+self.S_R+self.S_V+self.S_E

    def _Gibbs(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping free energy energy calculation\n")
            self.gibbs = 'NaN'
            return
        self.gibbs = self.enthalpy - self.T*self.entropy / self.constants.au_to_kJmol

    def _Optimized_Geometry(self) -> None:
        start = Forward_search_last(self.filename, 'Standard orientation', 'geometry', quiet=self.quiet)
        # end = Forward_search_after_last(self.filename, 'Standard orientation', 'Rotational constants', 200, "end of geometry", quiet=self.quiet)
        if start != "NaN":# and end != "NaN":
            #Offset for going into actual coordinate list
            start += 5
            for i, line in enumerate(self.lines[start:]):
                if '---------------------------------------------------------------------' in line:
                    end = start + i
                    break
            #Which position in the line is the atom label / number at
            label_location = 1
            OptGeomFilename = self.filename[:-4] + "_opt.xyz"
            GenerateXYZ(self.lines, OptGeomFilename, start, end, label_location, transform = True)
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write("Final geometry has been saved to " + OptGeomFilename + "\n")


class OrcaExtract:
    def __init__(self, filename: str, *, Quiet: bool = False, Temperature: float = 298.15) -> None:
        self.filename = filename
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self) -> None:
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _CPUS(self) -> None:
        linenumber = Backward_search_last(self.filename, 'Sum of individual times         ...', self.end, 'CPU time', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.wall_cpu_time = float(self.lines[linenumber].split()[-2])
            linenumber2 = Forward_search_last(self.filename, '%pal nprocs', 'CPU count', quiet=True)
            if isinstance(linenumber2, int):
                self.total_cpu_time = self.wall_cpu_time * int(self.lines[linenumber2].split()[-1])
            linenumber3 = Forward_search_last(self.filename, 'PAL', 'CPU count', quiet=self.quiet)
            if isinstance(linenumber3, int):
                self.total_cpu_time = self.wall_cpu_time * int(self.lines[linenumber3].split()[-1][3:])
            return
        self.total_cpu_time = 'NaN'
        self.wall_cpu_time = 'NaN'

    def _Energy(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Electronic energy', 'Final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-2])
            return
        linenumber = Forward_search_last(self.filename, 'FINAL SINGLE POINT ENERGY', 'Final energy', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        self.tot_energy = 'NaN'

    def _ZPV(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Electronic energy', 'ZPV energy', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.zpv = float(self.lines[linenumber].split()[-2]) + float(self.lines[linenumber+1].split()[-4])
            return
        self.zpv = 'NaN'

    def _Enthalpy(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping partition function calculation\n")
            self.enthalpy = 'NaN'
            return
        self._RotationalConsts()
        self.E_T = 3/2 * self.T * self.constants.gas_constant
        if len(self.rots) == 1:
            self.E_R = self.T * self.constants.gas_constant
        else:
            self.E_R = 3/2 * self.T * self.constants.gas_constant
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.E_V = self.constants.gas_constant * np.sum(self.constants.vib_const * realfreq * (1/2 + 1 / (np.exp(self.constants.vib_const * realfreq /  self.T ) - 1)))
        self.E_e = 0 #Good approximation for most closed-shell molecules
        self.enthalpy = (self.E_T+self.E_R+self.E_V+self.constants.gas_constant *  self.T ) / self.constants.au_to_kJmol + self.tot_energy

    def _Gibbs(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping free energy energy calculation\n")
            self.gibbs = 'NaN'
            return
        self.gibbs = self.enthalpy - self.T*self.entropy / self.constants.au_to_kJmol

    def _Dipole_moments(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Total Dipole Moment', 'dipole moment', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[linenumber].split()[-3]), float(self.lines[linenumber].split()[-2]), float(self.lines[linenumber].split()[-1]), float(self.lines[linenumber+2].split()[-1])
            return
        self.dipolex, self.dipoley, self.dipolez, self.total_dipole = 'NaN'

    def _Polarizabilities(self) -> None:
        linenumber = Forward_search_after_last(self.filename, 'THE POLARIZABILITY TENSOR', "'diagonalized tensor:'", 10, 'polarizability', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.polx, self.poly, self.polz, self.iso_polar = float(self.lines[linenumber+1].split()[0]), float(self.lines[linenumber+1].split()[1]), float(self.lines[linenumber+1].split()[2]), float(self.lines[linenumber+7].split()[-1])
            return
        self.polx = self.poly = self.polz = self.iso_polar = 'NaN'

    def _Excitation_energies(self) -> None:
        self.exc_energies = []
        linenumbers = Forward_search_all(self.filename, 'STATE ', 'excitation energies', quiet=self.quiet)
        if isinstance(linenumbers, list):
            for i in linenumbers:
                self.exc_energies.append(float(self.lines[i].split()[3]))
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN']

    def _Oscillator_strengths(self) -> None:
        self.osc_strengths = []
        linenumber = Forward_search_last(self.filename, 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS', 'oscillator strengths', quiet=self.quiet)
        if isinstance(linenumber, int):
            for i in range(len(self.exc_energies)):
                if len(self.lines[linenumber+5+i].split()) > 6:
                    self.osc_strengths.append(float(self.lines[linenumber+5+i].split()[3]))
                else:
                    self.osc_strengths.append('NaN')
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN']

    def _Frequencies(self) -> None:
        self.freq = []
        linenumber = Forward_search_last(self.filename, "VIBRATIONAL FREQUENCIES", 'frequencies', quiet=self.quiet)
        if isinstance(linenumber, int):
            for j in self.lines[linenumber+7: self.end]:
                if ": " and " 0.00 " in j:
                    pass
                elif ": " in j and not " 0.00 " in j:
                    self.freq.append(float(j.split()[1])* self.constants.inv_cm_to_au)
                else:
                    break
        if len(self.freq) == 0:
            self.freq = ['NaN']

    def _RotationalConsts(self) -> None:
        self.rots = []
        linenumbers = Forward_search_first(self.filename, 'Rotational constants in MHz', 'rotational constants', quiet=self.quiet)
        for i in self.lines[linenumbers].split()[-3:]:
            self.rots.append(float(i))
        self.rots = np.array(self.rots) * 1E-3
        self.rots = self.rots[self.rots != 0.0]

    def _Mass(self) -> None:
        self.mass = 0.0
        linenumber = Forward_search_last(self.filename, 'Total Mass', 'molecular mass', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.mass = float(self.lines[linenumber].split()[-2])

    def _SymmetryNumber(self) -> None:
        self.symnum = 0
        linenumber = Forward_search_last(self.filename, 'Symmetry Number', 'rotational symmetry number', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.symnum = int(self.lines[linenumber].split()[-1])

    def _Multiplicity(self) -> None:
        self.multi = 0
        linenumber = Forward_search_first(self.filename, 'Multiplicity', 'multiplicity', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.multi = int(self.lines[linenumber].split()[-1])

    def _PartitionFunctions(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping partition function calculation\n")
            self.qTotal = 'NaN'
            return
        self._RotationalConsts()
        self._Mass()
        self._Multiplicity()
        self._SymmetryNumber()
        self.qT = self.constants.trans_const_fac * self.mass ** (1.5) * self.T ** (2.5)
        if len(self.rots) == 1:
            self.qR = self.constants.rot_lin_const * self.T / (self.symnum * self.rots[0])
        else:
            self.qR = self.constants.rot_poly_const * self.T ** (1.5) / ( self.symnum * np.prod(np.array(self.rots)) ** (0.5))
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.qV = np.prod(1 / (1 - np.exp( - self.constants.rot_poly_const * realfreq /  self.T )))
        self.qE = self.multi #Good approximation for most closed-shell molecules
        self.qTotal = self.qT*self.qR*self.qV*self.qE

    def _Entropy(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping partition function calculation\n")
            self.entropy = 'NaN'
            return
        self._RotationalConsts()
        self._Mass()
        self._Multiplicity()
        self._SymmetryNumber()
        self.S_T = self.constants.gas_constant * np.log(self.constants.s_trans_const * self.mass ** 1.5 * self.T ** 2.5)
        if len(self.rots) == 1:
            self.S_R = self.constants.gas_constant * np.log(self.constants.rot_lin_const * self.T / (self.symnum * self.rots[0]))
        else:
            self.S_R = self.constants.gas_constant * (3/2 + np.log(self.constants.rot_poly_const * self.T ** (1.5) / ( self.symnum * np.prod(np.array(self.rots)) ** (0.5))))
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.S_V = self.constants.gas_constant * np.sum(self.constants.vib_const * realfreq / self.T / (np.exp(self.constants.vib_const * realfreq /  self.T ) - 1) - np.log(1-np.exp(-self.constants.vib_const * realfreq /  self.T )))
        self.S_E = self.constants.gas_constant * np.log(self.multi) #Good approximation for most closed-shell molecules
        self.entropy = self.S_T+self.S_R+self.S_V+self.S_E

    def _Optimized_Geometry(self) -> None:
        start = Forward_search_last(self.filename, 'CARTESIAN COORDINATES (ANGSTROEM)', 'geometry', quiet=self.quiet)
        if start != "NaN":
            #Offset for going into actual coordinate list
            start += 2
            for i, line in enumerate(self.lines[start:]):
                if'----------------------------' in line:
                    end = start + i - 1
                    break
            #Which position in the line is the atom label / number at
            label_location = 0
            OptGeomFilename = self.filename[:-4] + "_opt.xyz"
            GenerateXYZ(self.lines, OptGeomFilename, start, end, label_location)
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write("Final geometry has been saved to " + OptGeomFilename + "\n")


class DaltonExtract:
    def __init__(self, filename: str, NeededArguments: dict = None, Quiet: bool = False, Temperature: float = 298.15) -> None:
        self.filename = filename
        self.NeededArguments = NeededArguments
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self) -> None:
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _Complex_propagator(self) -> None:
        linenumbers = Forward_search_all(self.filename, 'Averaged value', 'polarizability with damping', quiet=self.quiet)
        if isinstance(linenumbers, list):
            self.complex_propagator = []
            for i in linenumbers:
                self.complex_propagator.append([float(self.lines[i].split()[-3]), float(self.lines[i].split()[-2]), float(self.lines[i].split()[-1])])
            return
        self.complex_propagator = 'NaN'

    def _CPUS(self) -> None:
        linenumber = Backward_search_last(self.filename, 'Total CPU  time used in DALTON:', self.end, 'CPU time', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.total_cpu_time = 0.
            self.wall_cpu_time = 0.
            total_time = self.lines[linenumber].split()[6:]
            pr_time = self.lines[linenumber+1].split()[6:]
            for i, time_value in enumerate(total_time[-2::-2]):
                if i*2 == 0:
                    self.total_cpu_time += float(time_value) / 60
                elif i*2 == 2:
                    self.total_cpu_time += float(time_value)
                elif i*2 == 4:
                    self.total_cpu_time += float(time_value) * 60
                elif i*2 == 6:
                    self.total_cpu_time += float(time_value) * 60 * 24
                else:
                    with open("collect_data.log", "a") as logfile:
                        logfile.write('''It was not expected that DALTON would print anything larger than days in the total CPU time
This will not be accounted for when printing the CPU time. The result will therefore not be correct
Please contact a maintainer of the script ot have this updated\n''')
            for i, time_value in enumerate(pr_time[-2::-2]):
                if i*2 == 0:
                    self.wall_cpu_time += float(time_value) / 60
                elif i*2 == 2:
                    self.wall_cpu_time += float(time_value)
                elif i*2 == 4:
                    self.wall_cpu_time += float(time_value) * 60
                elif i*2 == 6:
                    self.wall_cpu_time += float(time_value) * 60 * 24
                else:
                    with open("collect_data.log", "a") as logfile:
                        logfile.write('''It was not expected that DALTON would print anything larger than days in the total CPU time
This will not be accounted for when printing the CPU time. The result will therefore not be correct
Please contact a maintainer of the script ot have this updated\n''')
            return
        self.wall_cpu_time = 'NaN'
        self.total_cpu_time = 'NaN'

    def _Energy(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Total .*  energy:', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        linenumber = Forward_search_last(self.filename, '@    Final .* energy:', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        linenumber = Forward_search_last(self.filename, '@ Energy at final geometry is', 'final energy', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-2])
            return
        self.tot_energy = 'NaN'

    def _ZPV(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Total Molecular Energy', 'zero-point energy', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.zpv = float(self.lines[linenumber+5].split()[1])
            return
        self.zpv = 'NaN'

    def _Dipole_moments(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Dipole moment components', 'dipole moment', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[linenumber+5].split()[1]), float(self.lines[linenumber+6].split()[1]), float(self.lines[linenumber+7].split()[1]), float(self.lines[linenumber-3].split()[0])
            return
        self.dipolex = self.dipoley = self.dipolez = self.total_dipole = 'NaN'

    def _Polarizabilities(self) -> None:
        linenumber = Forward_search_last(self.filename, 'SECOND ORDER PROPERTIES', 'polarizabilities', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.polx, self.poly, self.polz = float(self.lines[linenumber+2].split()[-1]), float(self.lines[linenumber+5].split()[-1]), float(self.lines[linenumber+7].split()[-1])
            self.iso_polar = (self.polx + self.poly + self.polz)/3.
            return
        self.polx = self.poly = self.polz = self.iso_polar = 'NaN'

    def _Excitation_energies(self) -> None:
        self.exc_energies = []
        self.exc_type = None
        linenumber = Forward_search_last(self.filename, '@  Oscillator strengths are dimensionless.', 'excitation energies', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.exc_type = '.EXCITA'
            for i in self.lines[linenumber+5: self.end]:
                if "@ "in i:
                    self.exc_energies.append(float(i.split()[3])* self.constants.ev_to_au)
                else:
                    break
        linenumber = Forward_search_last(self.filename, '|  sym. | Exci.  |        CCSD       Excitation energies            | ||T1||  |', 'excitation energies', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.exc_type = '.EXCITA_sym'
            for i in self.lines[linenumber+4: self.end]:
                if '------' in i:
                    continue
                elif len(i.split()) > 1:
                    self.exc_energies.append(float(i.split()[5]))
                else:
                    break
        linenumber = Forward_search_last(self.filename, '@                  Oscillator and Scalar Rotational Strengths', 'excitation energies', quiet=self.quiet)
        if isinstance(linenumber, int) and self.exc_type != '.EXCITA':
            self.exc_type = '.ECD'
            for i in self.lines[linenumber+9: self.end]:
                if "@  " in i:
                    self.exc_energies.append(float(i.split()[3])* self.constants.ev_to_au)
                else:
                    break
        linenumbers = Forward_search_all(self.filename, '@ Excitation energy', 'excitation energies', quiet=self.quiet)
        if isinstance(linenumbers, list):
            self.exc_type = 'MCTDHF'
            for i in linenumbers:
                self.exc_energies.append(float(self.lines[i].split()[-2]))
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN']

    def _Oscillator_strengths(self) -> None:
        self.osc_strengths = []
        if self.exc_type == '.EXCITA':
            linenumber = Forward_search_last(self.filename, '@  Oscillator strengths are dimensionless.', 'oscillator strengths', quiet=self.quiet)
            if isinstance(linenumber, int):
                for i in self.lines[linenumber+5: self.end]:
                    if "@ " in i:
                        self.osc_strengths.append(float(i.split()[-1]))
                    else:
                        break
        elif self.exc_type == ".EXCITA_sym":
            linenumber = Forward_search_last(self.filename, '|  sym. | Exci.  |        CCSD       Length   Gauge Oscillator Strength       |', 'excitation energies', quiet=self.quiet)
            if isinstance(linenumber, int):
                for i in self.lines[linenumber+4: self.end]:
                    if '------' in i:
                        continue
                    elif len(i.split()) > 1:
                        self.osc_strengths.append(float(i.split()[-4]))
                    else:
                        break
        elif self.exc_type == 'MCTDHF':
            linenumbers = Forward_search_all(self.filename, '@ Excitation energy', 'oscillator strengths', quiet=self.quiet)
            if isinstance(linenumbers, list):
                for i in linenumbers:
                    osc = 0
                    for j in self.lines[i:i+15]:
                        if '@ Oscillator strength' in j:
                            osc += float(j.split()[5])**2
                    self.osc_strengths.append(osc**0.5)
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN']

    def _Rotational_strengths(self) -> None:
        self.rot_strengths = []
        if self.exc_type in ('.EXCITA', '.ECD'):
            linenumber = Forward_search_last(self.filename, '@                  Oscillator and Scalar Rotational Strengths', 'rotational strengths', quiet=self.quiet)
            if isinstance(linenumber, int):
                for i in self.lines[linenumber+9: self.end]:
                    if "@  " in i:
                        self.rot_strengths.append(float(i.split()[-2]))
                    else:
                        break
        if len(self.rot_strengths) == 0:
            self.rot_strengths = ['NaN']

    def _Frequencies(self) -> None:
        self.freq = []
        linenumber = Forward_search_last(self.filename, 'Vibrational Frequencies and IR Intensities', 'frequencies', quiet=self.quiet)
        if isinstance(linenumber, int):
            for i in self.lines[linenumber+7: self.end]:
                if len(i.split()) < 1:
                    break
                self.freq.append(float(i.split()[3]))
        if len(self.freq) == 0:
            self.freq = ['NaN']

    def _RotationalConsts(self) -> None:
        self.rots = []
        linenumbers = Forward_search_last(self.filename, 'Rotational constants', 'rotational constants', quiet=self.quiet)
        for i in self.lines[linenumbers+7].split()[:-1]:
            self.rots.append(float(i))
        self.rots = np.array(self.rots) * 1E-3
        self.rots = self.rots[self.rots != 0.0]

    def _Mass(self) -> None:
        self.mass = 0.0
        linenumber = Forward_search_last(self.filename, 'Total mass:', 'molecular mass')
        if isinstance(linenumber, int):
            self.mass = float(self.lines[linenumber].split()[-2])

    #Symmetry checking not implemented by default in Dalton
    #def _SymmetryNumber(self):
    #
    #    self.symnum = 0
    #    linenumber = Forward_search_last(self.file, 'Symmetry Number', 'rotational symmetry number')
    #    if isinstance(linenumber, int):
    #        self.symnum = int(self.lines[linenumber].split()[-1])

    def _Multiplicity(self) -> None:
        self.multi = 0
        linenumber = Forward_search_last(self.filename, 'Spatial symmetry', 'multiplicity', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.multi = int(self.lines[linenumber].split()[2])

    def _PartitionFunctions(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping partition function calculation\n")
            self.qTotal = 'NaN'
            return
        self._RotationalConsts()
        self._Mass()
        self._Multiplicity()
        self.qT = self.constants.trans_const_fac * self.mass ** (1.5) * self.T ** (2.5)
        #Rotational does not give the same as Dalton, due to a correction from the assymmetric top being applied: 10.1063/1.1748490
        if len(self.rots) == 1:
            self.qR = self.constants.rot_lin_const * self.T / (self.rots[0])
        else:
            self.qR = self.constants.rot_poly_const * self.T ** (1.5) / (np.prod(self.rots) ** (0.5))
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.qV = np.prod(1 / (1 - np.exp( - self.constants.vib_const * realfreq /  self.T )))
        self.qE = self.multi #Good approximation for most closed-shell molecules
        self.qTotal = self.qT*self.qR*self.qV*self.qE

    def _Entropy(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping partition function calculation\n")
            self.entropy = 'NaN'
            return
        self._RotationalConsts()
        self._Mass()
        self._Multiplicity()
        self.S_T = self.constants.gas_constant * np.log(self.constants.s_trans_const * self.mass ** 1.5 * self.T ** 2.5)
        if len(self.rots) == 1:
            self.S_R = self.constants.gas_constant * np.log(self.constants.rot_lin_const * self.T / (self.rots[0]))
        else:
            self.S_R = self.constants.gas_constant * (3/2 + np.log(self.constants.rot_poly_const * self.T ** (1.5) / ( np.prod(self.rots) ** (0.5))))
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.S_V = self.constants.gas_constant * np.sum(self.constants.vib_const * realfreq / self.T / (np.exp(self.constants.vib_const * realfreq /  self.T ) - 1) - np.log(1-np.exp(-self.constants.vib_const * realfreq /  self.T )))
        self.S_E = self.constants.gas_constant * np.log(self.multi) #Good approximation for most closed-shell molecules
        self.entropy = self.S_T+self.S_R+self.S_V+self.S_E

    def _Enthalpy(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping partition function calculation\n")
            self.enthalpy = 'NaN'
            return
        self._RotationalConsts()
        self.E_T = 3/2 * self.T * self.constants.gas_constant
        if len(self.rots) == 1:
            self.E_R = self.T * self.constants.gas_constant
        else:
            self.E_R = 3/2 * self.T * self.constants.gas_constant
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.E_V = self.constants.gas_constant * np.sum(self.constants.vib_const * realfreq * (1/2 + 1 / (np.exp(self.constants.vib_const * realfreq /  self.T ) - 1)))
        self.E_e = 0 #Good approximation for most closed-shell molecules
        self.enthalpy = (self.E_T+self.E_R+self.E_V+self.constants.gas_constant *  self.T ) / self.constants.au_to_kJmol + self.tot_energy

    def _Gibbs(self) -> None:
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f"No frequencies found in {self.filename}, skipping free energy energy calculation\n")
            self.gibbs = 'NaN'
            return
        self.gibbs = self.enthalpy - self.T*self.entropy / self.constants.au_to_kJmol

    def _Optimized_Geometry(self) -> None:
        start = Forward_search_last(self.filename, 'Final geometry (xyz format; angstrom)', 'final geometry', quiet=self.quiet)
        if start != "NaN":
            #Offset for going into actual coordinate list
            start += 5
            end = start + int(self.lines[start-2])
            #Which position in the line is the atom label / number at
            label_location = 0
            OptGeomFilename = self.filename[:-4] + "_opt.xyz"
            GenerateXYZ(self.lines, OptGeomFilename, start, end, label_location)
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write("Final geometry has been saved to " + OptGeomFilename + "\n")
        else:
            start = Forward_search_last(self.filename, 'Cartesian Coordinates', 'initial geometry', quiet=self.quiet)
            if start != "NaN":
                start += 4
                end = start + int(int(self.lines[start-1].split(' ')[-1])/3)
                lines_to_add = []
                lines_to_add.append(str(end-(start))+ '\n')
                lines_to_add.append('\n')
                for line in self.lines[start:end]:
                    words = line.split()
                    lines_to_add.append(''.join([words[0].ljust(2),' ',f"{float(words[-7]) * self.constants.bohr_to_ao:.7f}".rjust(20),' ', f"{float(words[-4]) * self.constants.bohr_to_ao:.7f}".rjust(25), ' ',f"{float(words[-1]) * self.constants.bohr_to_ao:.7f}".rjust(25) ,'\n']))
                OptGeomFilename = self.filename[:-4] + "_opt.xyz"
                WriteToFile(OptGeomFilename,lines_to_add)
                if not(self.quiet):
                    with open("collect_data.log", "a") as logfile:
                        logfile.write("Initial geometry has been saved to " + OptGeomFilename + "\n")


class DiracExtract:
    def __init__(self, filename: str, NeededArguments: dict = None, Quiet: bool = False, Temperature: float = 298.15) -> None:
        self.filename = filename
        self.NeededArguments = NeededArguments
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self) -> None:
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _Energy(self) -> None:
        linenumber = Forward_search_last(self.filename, '@ Total .*  energy:', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        linenumber = Forward_search_last(self.filename, 'SCF energy:', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        self.tot_energy = 'NaN'

    def _Excitation_energies(self) -> None:
        self.exc_energies = []
        self.exc_type = None
        linenumbers = Forward_search_all(self.filename, 'Full light-matter interaction: Isotropic case', 'excitation energies', quiet=True)
        if isinstance(linenumbers, list):
            self.exc_type = '.EXCITA'
            for linenumber in linenumbers:
                for i in self.lines[linenumber+13: self.end]:
                    if len(i.split()) > 0:
                        self.exc_energies.append(float(i.split()[1]))
                    else:
                        break
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN']

    def _Oscillator_strengths(self) -> None:
        self.osc_strengths = []
        if self.exc_type == '.EXCITA':
            linenumbers = Forward_search_all(self.filename, 'Full light-matter interaction: Isotropic case', 'excitation energies', quiet=True)
            if isinstance(linenumbers, list):
                for linenumber in linenumbers:
                    for i in self.lines[linenumber+13: self.end]:
                        if len(i.split()) > 0:
                            self.osc_strengths.append(float(i.split()[-1]))
                        else:
                            break
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN']


class LSDaltonExtract:
    def __init__(self, filename: str, NeededArguments: dict = None, Quiet: bool = False, Temperature: float = 298.15) -> None:
        self.filename = filename
        self.NeededArguments = NeededArguments
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self) -> None:
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _CPUS(self) -> None:
        linenumber = Backward_search_last(self.filename, '>>>  CPU Time used in LSDALTON is', self.end, 'CPU time', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.total_cpu_time = 0.
            self.wall_cpu_time = 0.
            total_time = self.lines[linenumber].split()[7:]
            pr_time = self.lines[linenumber+1].split()[7:]
            for i, time_value in enumerate(total_time[-2::-2]):
                if i*2 == 0:
                    self.total_cpu_time += float(time_value) / 60
                elif i*2 == 2:
                    self.total_cpu_time += float(time_value)
                elif i*2 == 4:
                    self.total_cpu_time += float(time_value) * 60
                elif i*2 == 6:
                    self.total_cpu_time += float(time_value) * 60 * 24
                else:
                    with open("collect_data.log", "a") as logfile:
                        logfile.write('''It was not expected that LSDALTON would print anything larger than days in the total CPU time
This will not be accounted for when printing the CPU time. The result will therefore not be correct
Please contact a maintainer of the script ot have this updated\n''')
            for i, time_value in enumerate(pr_time[-2::-2]):
                if i*2 == 0:
                    self.wall_cpu_time += float(time_value) / 60
                elif i*2 == 2:
                    self.wall_cpu_time += float(time_value)
                elif i*2 == 4:
                    self.wall_cpu_time += float(time_value) * 60
                elif i*2 == 6:
                    self.wall_cpu_time += float(time_value) * 60 * 24
                else:
                    with open("collect_data.log", "a") as logfile:
                        logfile.write('''It was not expected that LSDALTON would print anything larger than days in the total CPU time
This will not be accounted for when printing the CPU time. The result will therefore not be correct
Please contact a maintainer of the script ot have this updated\n''')
            return
        self.total_cpu_time = 'NaN'
        self.wall_cpu_time = 'NaN'

    def _Energy(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Total .*  energy:', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        linenumber = Forward_search_last(self.filename, '@    Final .* energy:', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        linenumber = Forward_search_last(self.filename, '@ Energy at final geometry is', 'final energy', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-2])
            return
        self.tot_energy = 'NaN'

    def _Energy(self) -> None:
        linenumber = Forward_search_last(self.filename, 'ENERGY SUMMARY', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            for i in self.lines[linenumber+3:self.end]:
                if 'E: ' in i:
                    self.tot_energy = float(i.split()[-1])
                else:
                    return
        linenumber = Forward_search_last(self.filename, 'Final .* energy:', 'final energy', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        self.tot_energy = 'NaN'

    def _Dipole_moments(self) -> None:
        linenumber = Forward_search_last(self.filename, 'Permanent dipole moment', 'dipole moment', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[linenumber+9].split()[1]), float(self.lines[linenumber+10].split()[1]), float(self.lines[linenumber+11].split()[1]), float(self.lines[linenumber+3].split()[0])
            return
        self.dipolex = self.dipoley = self.dipolez = self.total_dipole = 'NaN'

    def _Polarizabilities(self) -> None:
        linenumber = Forward_search_last(self.filename, '*          POLARIZABILITY TENSOR RESULTS (in a.u.)          *', 'polarizability', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.polx, self.poly, self.polz, self.iso_polar = float(self.lines[linenumber+10].split()[-3]), float(self.lines[linenumber+11].split()[-2]), float(self.lines[linenumber+12].split()[-1]), float(self.lines[linenumber+14].split()[-1])
            return
        self.polx = self.poly = self.polz = self.iso_polar = 'NaN'

    def _Excitation_energies(self) -> None:
        self.exc_energies = []
        linenumber = Forward_search_last(self.filename, '*                   ONE-PHOTON ABSORPTION RESULTS (in a.u.)                  *', 'excitation energies', quiet=self.quiet)
        if isinstance(linenumber, int):
            for i in range(linenumber+8,self.end):
                if len(self.lines[i].split()) < 1:
                    break
                self.exc_energies.append(float(self.lines[i].split()[0]))
        else:
            linenumber = Forward_search_last(self.filename, 'excitation energies', 'excitation energies', quiet=self.quiet)
            if isinstance(linenumber, int):
                for i in range(linenumber+4,self.end):
                    if len(self.lines[i].split()) < 1:
                        break
                    try:
                        self.exc_energies.append(float(self.lines[i].split()[1]))
                    except ValueError:
                        break
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN']

    def _Oscillator_strengths(self) -> None:
        self.osc_strengths = []
        linenumber = Forward_search_last(self.filename, '*                   ONE-PHOTON ABSORPTION RESULTS (in a.u.)                  *', 'oscillator strengths', quiet=self.quiet)
        if isinstance(linenumber, int):
            for i in range(len(self.exc_energies)):
                self.osc_strengths.append(float(self.lines[linenumber+8+i].split()[-1]))
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN']

    def _Optimized_Geometry(self) -> None:
        start = Forward_search_last(self.filename, 'Final geometry', 'geometry', quiet=self.quiet)
        if start != "NaN":
            #Offset for going into actual coordinate list
            start += 6
            atoms_in_molecule = int(int(self.lines[start-3].split()[-1])/3)
            end = start + atoms_in_molecule*4
            lines_to_add = []
            lines_to_add.append(str(atoms_in_molecule)+ '\n')
            lines_to_add.append('\n')
            for i, _ in enumerate(self.lines[start:end:4]):
                linenr = start + i*4
                # print(i, start-end)
                lines_to_add.append(''.join([self.lines[linenr].split()[1].ljust(2),' ', f"{float(self.lines[linenr].split()[-1]) * self.constants.bohr_to_ao:.7f}".rjust(20),' ', f"{float(self.lines[linenr+1].split()[-1]) * self.constants.bohr_to_ao:.7f}".rjust(25),' ', f"{float(self.lines[linenr+2].split()[-1]) * self.constants.bohr_to_ao:.7f}".rjust(25), '\n']))
            OptGeomFilename = self.filename[:-4] + "_opt.xyz"
            WriteToFile(OptGeomFilename,lines_to_add)
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write("Final geometry has been saved to " + OptGeomFilename + "\n")
        else:
            start = Forward_search_last(self.filename, 'PRINTING THE MOLECULE.INP FILE', 'initial geometry', quiet=self.quiet)
            if start != "NaN":
                start += 7
                current_line = start
                atoms_in_molecule = 0
                lines_to_add = []
                lines_to_add.append(str(atoms_in_molecule) + '\n')
                lines_to_add.append('\n')
                for _ in range(int(self.lines[start-1].split()[0].split('=')[-1])):
                    current_atom = AtomicInformation(int(self.lines[current_line].split()[0].split('.')[0]))
                    for i in range(int(self.lines[current_line].split()[1])):
                        lines_to_add.append(''.join([current_atom.atom.ljust(2),' ', f"{float(self.lines[current_line+i+1].split()[1]):.7f}".rjust(20),' ', f"{float(self.lines[current_line+i+1].split()[2]):.7f}".rjust(25),' ', f"{float(self.lines[current_line+i+1].split()[3]):.7f}".rjust(25), '\n']))
                    else:
                        current_line += i+2
                        atoms_in_molecule += i+1
                lines_to_add[0] = f"{atoms_in_molecule}\n"
            OptGeomFilename = self.filename[:-4] + "_opt.xyz"
            WriteToFile(OptGeomFilename,lines_to_add)
            if not(self.quiet):
                with open("collect_data.log", "a") as logfile:
                    logfile.write("Final geometry has been saved to " + OptGeomFilename + "\n")


class QChemExtract:
    def __init__(self, filename: str, NeededArguments: dict = None, Quiet: bool = False, Temperature: float = 298.15) -> None:
        self.filename = filename
        self.NeededArguments = NeededArguments
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self) -> None:
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _CPUS(self) -> None:
        linenumber = Backward_search_last(self.filename, 'Total job time:', self.end, 'CPU time', quiet=self.quiet)
        if isinstance(linenumber, int):
            self.total_cpu_time = float(self.lines[linenumber].split()[4].split('s')[0])/60
            self.wall_cpu_time = float(self.lines[linenumber].split()[3].split('s')[0])/60
            return
        self.total_cpu_time = 'NaN'
        self.wall_cpu_time = 'NaN'

    def _Energy(self) -> None:
        linenumber = Forward_search_last(self.filename, 'total energy:', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        linenumber = Forward_search_last(self.filename, 'MP2 energy:', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        linenumber = Forward_search_last(self.filename, 'SCF energy:', 'final energy', quiet=True)
        if isinstance(linenumber, int):
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        self.tot_energy = 'NaN'

    def _Excitation_energies(self) -> None:
        self.exc_energies = []
        linenumbers = Forward_search_all(self.filename, 'Energy GAP', 'excitation energies', quiet=self.quiet)
        if isinstance(linenumbers, list):
            for i in linenumbers:
                self.exc_energies.append(float(self.lines[i].split()[3]))
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN']

    def _Oscillator_strengths(self) -> None:
        self.osc_strengths = []
        linenumbers = Forward_search_all(self.filename, 'Oscillator strength', 'oscillator strengths', quiet=self.quiet)
        if isinstance(linenumbers, list):
            for i in linenumbers:
                self.osc_strengths.append(float(self.lines[i].split()[-1]))
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN']


