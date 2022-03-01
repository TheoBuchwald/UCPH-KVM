
import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from types import FunctionType
from colorama import Style, Fore, init

init(autoreset=True)

def Forward_search_last(file: str, text: str, error: str, err: bool =True, quiet: bool = False):
    """Searches from the beggining of the file given to the end where it returns the linenumber of the last occurence

    Args:
        file (str): The file to search in
        text (str): The text string to search for
        error (str): If no occurences were found it will print 'No [error] could be found in [file]
        err (bool, optional): Whether or not to print error message if no occurences are found. Defaults to True.

    Returns:
        (int): Linenumber of last occurence
    """
    ps1 = subprocess.run(['grep', '-nT', text, file], capture_output=True)
    out = subprocess.run(['tail', '-n1'], input=ps1.stdout, capture_output=True)
    res = out.stdout
    if len(res) == 0:
        if not(quiet):
            if err:
                print(f'{Fore.RED}No {error}{Style.RESET_ALL} could be found in {Style.BRIGHT}{Fore.YELLOW}{file}')
        return 'NaN'
    res = res.split()[0]
    res = str(res).split(':')
    return int(res[0].replace('b\'','').replace('\'','')) - 1

def Forward_search_after_last(file: str, text1: str, text2: str, lines: int, error: str, err: bool=True, quiet: bool = False):
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
    ps1 = subprocess.run(['grep', '-nTA', f'{lines}', text1, file], capture_output=True)
    ps2 = subprocess.run(['tail', '-n', f'{lines + 1}'], input=ps1.stdout, capture_output=True)
    out = subprocess.run(['grep', text2], input=ps2.stdout, capture_output=True)
    res = out.stdout
    if len(res) == 0:
        if not(quiet):
            if err:
                print(f'{Fore.RED}No {error}{Style.RESET_ALL} could be found in {Style.BRIGHT}{Fore.YELLOW}{file}')
        return 'NaN'
    res = res.split()[0]
    res = str(res).split('-')
    return int(res[0].replace('b\'','')) - 1

def Forward_search_first(file: str, text: str, error: str, err: bool=True, quiet: bool = False):
    """Searches from beginning of file and finds the first occurence of [text]

    Args:
        file (str): File to search in
        text (str): Text to look for
        error (str): If no occurences were found it will print 'No [error] could be found in [file]
        err (bool, optional): Whether or not to print error message if no occurences are found. Defaults to True.. Defaults to True.

    Returns:
        (int): Linenumber of first occurence
    """
    out = subprocess.run(['grep', '-nTm1', text, file], capture_output=True)
    res = out.stdout
    if len(res) == 0:
        if not(quiet):
            if err:
                print(f'{Fore.RED}No {error}{Style.RESET_ALL} could be found in {Style.BRIGHT}{Fore.YELLOW}{file}')
        return 'NaN'
    res = res.split()[0]
    res = str(res).split(':')
    return int(res[0].replace('b\'','')) - 1

def Forward_search_all(file: str, text: str, error: str, err: bool=True, quiet: bool = False):
    """Searches from beggining of file to end of file finding all occurences of [text]

    Args:
        file (str): File to search in
        text (str): Text to look for
        error (str): If no occurences were found it will print 'No [error] could be found in [file]
        err (bool, optional): Whether or not to print error message if no occurences are found. Defaults to True.. Defaults to True.

    Returns:
        (list): List of the linenumbers of all occurences
    """
    ps1 = subprocess.run(['grep', '-nT', text, file], capture_output=True)
    out = subprocess.run(['awk', '{print $1}'], input=ps1.stdout, capture_output=True)
    res = out.stdout
    if len(res) == 0:
        if not(quiet):
            if err:
                print(f'{Fore.RED}No {error}{Style.RESET_ALL} could be found in {Style.BRIGHT}{Fore.YELLOW}{file}')
        return 'NaN'
    res = str(res).split('\\n')
    return [int(val.replace('b\'','').replace(':','')) - 1 for val in res[:-1]]

def Resize(array: list):
    """Takes an array of arrays with varying sizes and resizes them to the same size.
    This is done by appending 'NaN' to the smaller lists.
    If the first value says 'Not implemented' the remainder of the array will be filled with 'Not implemented

    Args:
        array (list): Array of arrays with varying sizes
    """
    max_size = 0
    for arr in array:
        if len(arr) > max_size:
            max_size = len(arr)

    for i, arr in enumerate(array):
        if arr == ['Not implemented']:
            array[i] *= max_size
        else:
            array[i] += ['NaN'] * (max_size - len(arr))

def CheckForOnlyNans(array: list):
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

def Data_Extraction(infile, Needed_Values: dict, NeededArguments, quiet: bool = False, Temperature: float = 298.15):
    Extracted_values = dict()

    infile = Output_type(str(infile), NeededArguments, quiet, Temperature)

    Extract_data(quiet, Needed_Values, infile.filename, infile.filetype, infile.input)  #Extracting data
    infile.filetype.__delattr__('lines') #Removing the file text from memory

    dict_keys = [*infile.filetype.__dict__.keys()]
    collection_dict = dict()

    for i in dict_keys[1:]: #Collecting the data in dictionaries
        collection_dict[i] =  infile.filetype.__dict__[i]

    Extracted_values[infile.filename] = collection_dict

    return Extracted_values

def Extract_data(suppressed: bool, Wanted_Values: dict, infile: str, file_text: dict, input_type: str):
    for i in Wanted_Values:
        try:
            method = getattr(type(file_text),i)
            method(file_text)
        except AttributeError:
            if not(suppressed):
                print(f'{Style.BRIGHT}{Fore.YELLOW}{infile}{Style.RESET_ALL}: {Fore.RED}{i}{Style.RESET_ALL} has not been implemented for {Style.BRIGHT}{Fore.CYAN}{input_type}')

def Check_if_Implemented(input_file: str, Set_of_values: dict, Extracted_values: dict):
    for infile in input_file:
        for key in Set_of_values:
            for val in Set_of_values[key]:
                try:
                    Extracted_values[infile][val]
                except KeyError:
                    Extracted_values[infile][val] = ['Not implemented']

def Collect_and_sort_data(input_file: str, Set_of_values: dict, Extracted_values: dict):
    Final_arrays = dict()

    for key in Set_of_values:
        for val in Set_of_values[key]:
            Final_arrays[val] = []
            for infile in input_file:
                Final_arrays[val].append(Extracted_values[infile][val])
    return Final_arrays

def Downsizing_variable_arrays(Outputs: dict, Variable_arrays: dict, count: int, Final_arrays: dict):
    for item in Variable_arrays.items():
        if item[1] > 0:
            for val in Outputs[item[0]]:
                for file in range(0,count):
                    Final_arrays[val][file] = Final_arrays[val][file][0:item[1]]

def Create_Header(Header_text: dict, Set_of_values: dict, Final_arrays: dict):
    header = ['File']
    for key in Set_of_values.keys():
        for val in Set_of_values[key]:
            if len(Final_arrays[val][0]) > 1:
                for i in range(len(Final_arrays[val][0])):
                    header.append(f'{Header_text[val]} {i+1}')
            else:
                header.append(Header_text[val])
    return header

def Fill_output_array(Set_of_values: dict, array_input: dict, count: int, Final_arrays: dict, output_array: list):
    if count == 1:
        output_array[1,0] = array_input[0][0]
    else:
        output_array[1:,0] = np.array(np.concatenate(array_input))

    col = 1
    for key in Set_of_values.keys():
        for val in Set_of_values[key]:
            output_array[1:,col:col+len(np.array(Final_arrays[val][0]))] = np.array(Final_arrays[val])
            col += len(np.array(Final_arrays[val][0]))

def UVVIS_Spectrum(t: list, l: list, f: list, k: float, sigmacm: float):
    lambda_tot = np.zeros(len(t))
    for x in range(1,len(t)):
        lambda_tot[x] = sum((k/sigmacm)*f*np.exp(-4*np.log(2)*((1/t[x]-1/l)/(1E-7*sigmacm))**2))
    return lambda_tot

def Make_uvvis_spectrum(input_file: list, suppressed: bool, UVVIS_Spectrum: FunctionType, count: int, Final_arrays: dict, UVVIS: str, Extracted_Values: dict, SAVE: bool = True):
    # A LOT OF PLOT SETUP
    rc('text', usetex=True)
    xlabel_font = ylabel_font = title_font = 16
    plt.rc('font', size=12) # x/y axis font size
    N=1000 # number of calculated points in curve
    NA=6.02214199*10**23 #avogadros number
    c=299792458 #speed of light
    e=1.60217662*10**(-19) #electron charge
    me=9.10938*10**(-31) #electron mass
    pi=math.pi
    epsvac=8.8541878176*10**(-12)
    sigmacm=0.4*8065.544
    k=(NA*e**2)/(np.log(10)*2*me*c**2*epsvac)*np.sqrt(np.log(2)/pi)*10**(-1)
    inv_cm_to_au = 1/219474.63068
    Save_Dict = dict()
        # PLOT SETUP DONE
    for file in input_file:
        filename = file.replace('.out','')
        plotname = f'{filename}-uvvis.{UVVIS}'
        title = filename.replace("_"," ")

        excitations = np.array([x for x in Extracted_Values[file]['exc_energies'] if x != 'NaN' and x != 'Not implemented'])
        oscillations = np.array([x for x in Extracted_Values[file]['osc_strengths'] if x != 'NaN' and x != 'Not implemented'])

        if len(excitations) == 0 and len(oscillations) == 0:
            if not(suppressed):
                print(f'{Fore.RED}Excitation energies and oscillator strengths{Style.RESET_ALL} have either not been implemented for this {Style.BRIGHT}{Fore.CYAN}output type{Style.RESET_ALL}, or none were found in the file {Style.BRIGHT}{Fore.YELLOW}{filename}')
            continue
        elif len(excitations) == 0:
            if not(suppressed):
                print(f'{Fore.RED}Excitation energies{Style.RESET_ALL} have either not been implemented for this {Style.BRIGHT}{Fore.CYAN}output type{Style.RESET_ALL}, or none were found in the file {Style.BRIGHT}{Fore.YELLOW}{filename}')
            continue
        elif len(oscillations) == 0:
            if not(suppressed):
                print(f'{Fore.RED}Oscillator strengths{Style.RESET_ALL} have either not been implemented for this {Style.BRIGHT}{Fore.CYAN}output type{Style.RESET_ALL}, or none were found in the file {Style.BRIGHT}{Fore.YELLOW}{filename}')
            continue
        elif len(excitations) > len(oscillations):
            excitations = excitations[0:len(oscillations)]
        elif len(oscillations) > len(excitations):
            oscillations = oscillations[0:len(excitations)]

        excitations = 1E7/(excitations/ inv_cm_to_au)   # From a.u. to cm^-1 to nm
        span = np.linspace(min(excitations)-20, max(excitations)+20, N, endpoint=True) # exctinction coefficient (wavelength range)

        graph = UVVIS_Spectrum(span, excitations, oscillations, k, sigmacm)
        plt.title(title)
        plt.plot(span, graph)
        plt.ylim((0,max(graph)*1.2))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel('Wavelength $(nm)$', fontsize=xlabel_font)
        plt.ylabel('Extinction coefficient', fontsize=ylabel_font)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.savefig(plotname, format=f'{UVVIS}', dpi=600)
        plt.close()

        if SAVE:
            Save_Dict[file] = [span,graph]
    if SAVE:
        np.savez('UVVIS.npz', **Save_Dict)
        print(f'{Style.BRIGHT}{Fore.LIGHTGREEN_EX}UVVIS spectrum data has been saved in UVVIS.npz')

    del Save_Dict # Deletes dictionary from memory since it is no longer needed


class Output_type:
    def __init__(self, filename: str, NeededArguments: dict = None, Quiet: bool = False, Temperature: float = 298.15):
        self.filename = filename

        with open(self.filename,'r') as read:
            lines = read.readlines()[:10]

        if '* O   R   C   A *' in lines[4]: #File type = ORCA
            self.filetype = orca(self.filename, NeededArguments, Quiet, Temperature)
            self.input = 'ORCA'

        if '*************** Dalton - An Electronic Structure Program ***************' in lines[3]:  #File type = DALTON
            self.filetype = dal(self.filename, NeededArguments, Quiet, Temperature)
            self.input = 'DALTON'

        if 'Gaussian, Inc.  All Rights Reserved.' in lines[6]:  #File type = GAUSSIAN
            self.filetype = gaus(self.filename, NeededArguments, Quiet, Temperature)
            self.input = 'GAUSSIAN'

        if '**********  LSDalton - An electronic structure program  **********' in lines[2]:    #File type = LSDALTON
            self.filetype = lsdal(self.filename, NeededArguments, Quiet, Temperature)
            self.input = 'LSDALTON'

        del lines


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


class gaus:
    def __init__(self, filename: str, NeededArguments: dict = None, Quiet: bool = False, Temperature: float = 298.15):
        self.filename = filename
        self.NeededArguments = NeededArguments
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self):
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _Energy(self):
        linenumber = Forward_search_last(self.filename, 'Sum of electronic and zero-point Energies=', 'final energy', err=False)
        if type(linenumber) == int:
            self.tot_energy = float(self.lines[linenumber].split()[-1]) - float(self.lines[linenumber-4].split()[-2])
            return
        linenumber = Forward_search_last(self.filename, 'SCF Done:', 'final energy')
        if type(linenumber) == int:
            self.tot_energy = float(self.lines[linenumber].split()[4])
            return
        self.tot_energy = 'NaN'

    def _ZPV(self):
        linenumber = Forward_search_last(self.filename, 'Sum of electronic and zero-point Energies=', 'ZPV energy')
        if type(linenumber) == int:
            self.zpv = float(self.lines[linenumber].split()[-1])
            return
        self.zpv = 'NaN'

    def _Dipole_moments(self):
        linenumber = Forward_search_last(self.filename, 'Electric dipole moment (input orientation):', 'dipole moments')
        if type(linenumber) == int:
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[linenumber+4].split()[1].replace('D','E')), float(self.lines[linenumber+5].split()[1].replace('D','E')), float(self.lines[linenumber+6].split()[1].replace('D','E')), float(self.lines[linenumber+3].split()[1].replace('D','E'))
            return
        self.dipolex = self.dipoley = self.dipolez = self.total_dipole = 'NaN'

    def _Polarizabilities(self):
        linenumber = ['NaN', 'NaN', 'NaN', 'NaN']
        searchwords = [' xx ', ' yy ', ' zz ', ' iso ']
        for i in range(len(searchwords)):
            linenumber[i] = Forward_search_after_last(self.filename, 'Dipole polarizability, Alpha (input orientation).', searchwords[i], 15, 'polarizabilities')
        if linenumber != ['NaN', 'NaN', 'NaN', 'NaN']:
            self.polx, self.poly, self.polz, self.iso_polar = float(self.lines[linenumber[0]].split()[1].replace('D','E')), float(self.lines[linenumber[1]].split()[1].replace('D','E')), float(self.lines[linenumber[2]].split()[1].replace('D','E')), float(self.lines[linenumber[3]].split()[1].replace('D','E'))
            return
        self.polx = self.poly = self.polz = self.iso_polar = 'NaN'

    def _Frequencies(self):
        self.freq = []
        linenumbers = Forward_search_all(self.filename, 'Frequencies --', 'frequencies')
        if type(linenumbers) == list:
            for i in linenumbers:
                for j in self.lines[i].split()[2:]:
                    self.freq.append(float(j)* self.constants.inv_cm_to_au)
        if len(self.freq) == 0:
            self.freq = ['NaN'] * abs(self.NeededArguments['_Frequencies'])
        if len(self.freq) < self.NeededArguments['_Frequencies']:
            self.freq += ['NaN'] * (self.NeededArguments['_Frequencies'] - len(self.freq))

    def _Excitation_energies(self):
        self.exc_energies = []
        linenumber = Forward_search_last(self.filename, 'Excitation energies and oscillator strengths:', 'excitation energies', err=False)
        linenumbers = Forward_search_all(self.filename, 'Excited State', 'excitation energies')
        if type(linenumber) == int:
            linenumbers = [i for i in linenumbers if i > linenumber]
            for i in linenumbers:
                self.exc_energies.append(float(self.lines[i].split()[4])* self.constants.ev_to_au)
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN'] * abs(self.NeededArguments['_Excitation_energies'])
        if len(self.exc_energies) < self.NeededArguments['_Excitation_energies']:
            self.exc_energies += ['NaN'] * (self.NeededArguments['_Excitation_energies'] - len(self.exc_energies))

    def _Oscillator_strengths(self):
        self.osc_strengths = []
        linenumber = Forward_search_last(self.filename, 'Excitation energies and oscillator strengths:', 'oscillator strengths', err=False)
        linenumbers = Forward_search_all(self.filename, 'Excited State', 'oscillator strengths')
        if type(linenumber) == int:
            linenumbers = [i for i in linenumbers if i > linenumber]
            for i in linenumbers:
                self.osc_strengths.append(float(self.lines[i].split()[-1].replace('f=','')))
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN'] * abs(self.NeededArguments['_Excitation_energies'])
        if len(self.osc_strengths) < self.NeededArguments['_Excitation_energies']:
            self.osc_strengths += ['NaN'] * (self.NeededArguments['_Excitation_energies'] - len(self.osc_strengths))

    def _RotationalConsts(self):
        self.rots = []
        linenumbers = Forward_search_all(self.filename, 'Rotational constants (GHZ):', 'rotational constants')
        for i in self.lines[linenumbers[-2]].split()[3:]:
            self.rots.append(float(i))
        self.rots = np.array(self.rots)
        self.rots = self.rots[self.rots != 0.0]

    def _Mass(self):
        self.mass = 0.0
        linenumber = Forward_search_last(self.filename, 'Molecular mass', 'molecular mass')
        if type(linenumber) == int:
            self.mass = float(self.lines[linenumber].split()[2])

    def _SymmetryNumber(self):
        self.symnum = 0
        linenumber = Forward_search_last(self.filename, 'Rotational symmetry number', 'rotational symmetry number')
        if type(linenumber) == int:
            self.symnum = int(self.lines[linenumber].split()[-1].replace('.',''))

    def _Multiplicity(self):
        self.multi = 0
        linenumber = Forward_search_first(self.filename, 'Multiplicity', 'multiplicity')
        if type(linenumber) == int:
            self.multi = int(self.lines[linenumber].split()[-1])

    def _PartitionFunctions(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}partition function{Style.RESET_ALL} calculation")
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

    def _Enthalpy(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}partition function{Style.RESET_ALL} calculation")
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

    def _Entropy(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}partition function{Style.RESET_ALL} calculation")
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

    def _Gibbs(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}free energy{Style.RESET_ALL} energy calculation")
            self.gibbs = 'NaN'
            return
        self.gibbs = self.enthalpy - self.T*self.entropy / self.constants.au_to_kJmol


class orca:
    def __init__(self, filename: str, NeededArguments: dict = None, Quiet: bool = False, Temperature: float = 298.15):
        self.filename = filename
        self.NeededArguments = NeededArguments
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self):
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _Energy(self):
        linenumber = Forward_search_last(self.filename, 'Electronic energy', 'Final energy', err = False)
        if type(linenumber) == int:
            self.tot_energy = float(self.lines[linenumber].split()[-2])
            return
        linenumber = Forward_search_last(self.filename, 'FINAL SINGLE POINT ENERGY', 'Final energy')
        if type(linenumber) == int:
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        self.tot_energy = 'NaN'

    def _ZPV(self):
        linenumber = Forward_search_last(self.filename, 'Electronic energy', 'ZPV energy')
        if type(linenumber) == int:
            self.zpv = float(self.lines[linenumber].split()[-2]) + float(self.lines[linenumber+1].split()[-4])
            return
        self.zpv = 'NaN'

    def _Enthalpy(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}partition function{Style.RESET_ALL} calculation")
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

    def _Gibbs(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}free energy{Style.RESET_ALL} energy calculation")
            self.gibbs = 'NaN'
            return
        self.gibbs = self.enthalpy - self.T*self.entropy / self.constants.au_to_kJmol

    def _Dipole_moments(self):
        linenumber = Forward_search_last(self.filename, 'Total Dipole Moment', 'dipole moment')
        if type(linenumber) == int:
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[linenumber].split()[-3]), float(self.lines[linenumber].split()[-2]), float(self.lines[linenumber].split()[-1]), float(self.lines[linenumber+2].split()[-1])
            return
        self.dipolex, self.dipoley, self.dipolez, self.total_dipole = 'NaN'

    def _Polarizabilities(self):
        linenumber = Forward_search_after_last(self.filename, 'THE POLARIZABILITY TENSOR', "'diagonalized tensor:'", 10, 'polarizability')
        if type(linenumber) == int:
            self.polx, self.poly, self.polz, self.iso_polar = float(self.lines[linenumber+1].split()[0]), float(self.lines[linenumber+1].split()[1]), float(self.lines[linenumber+1].split()[2]), float(self.lines[linenumber+7].split()[-1])
            return
        self.polx = self.poly = self.polz = self.iso_polar = 'NaN'

    def _Excitation_energies(self):
        self.exc_energies = []
        linenumbers = Forward_search_all(self.filename, 'STATE ', 'excitation energies')
        if type(linenumbers) == list:
            for i in linenumbers:
                self.exc_energies.append(float(self.lines[i].split()[3]))
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN'] * abs(self.NeededArguments['_Excitation_energies']  )
        if len(self.exc_energies) < self.NeededArguments['_Excitation_energies']:
            self.exc_energies += ['NaN'] * (self.NeededArguments['_Excitation_energies'] - len(self.exc_energies))

    def _Oscillator_strengths(self):
        self.osc_strengths = []
        linenumber = Forward_search_last(self.filename, 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS', 'oscillator strengths')
        if type(linenumber) == int:
            for i in range(len(self.exc_energies)):
                if len(self.lines[linenumber+5+i].split()) > 6:
                    self.osc_strengths.append(float(self.lines[linenumber+5+i].split()[3]))
                else:
                    self.osc_strengths.append('NaN')
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN'] * abs(self.NeededArguments['_Excitation_energies'])
        if len(self.osc_strengths) < self.NeededArguments['_Excitation_energies']:
            self.osc_strengths += ['NaN'] * (self.NeededArguments['_Excitation_energies'] - len(self.osc_strengths))

    def _Frequencies(self):
        self.freq = []
        linenumber = Forward_search_last(self.filename, "VIBRATIONAL FREQUENCIES", 'frequencies')
        if type(linenumber) == int:
            for j in self.lines[linenumber+7: self.end]:
                if ": " and " 0.00 " in j:
                    pass
                elif ": " in j and not " 0.00 " in j:
                    self.freq.append(float(j.split()[1])* self.constants.inv_cm_to_au)
                else:
                    break
        if len(self.freq) == 0:
            self.freq = ['NaN'] * abs(self.NeededArguments['_Frequencies'])
        if len(self.freq) < self.NeededArguments['_Frequencies']:
            self.freq += ['NaN'] * (self.NeededArguments['_Frequencies'] - len(self.freq))

    def _RotationalConsts(self):
        self.rots = []
        linenumbers = Forward_search_first(self.filename, 'Rotational constants in MHz', 'rotational constants')
        for i in self.lines[linenumbers].split()[-3:]:
            self.rots.append(float(i))
        self.rots = np.array(self.rots) * 1E-3
        self.rots = self.rots[self.rots != 0.0]

    def _Mass(self):
        self.mass = 0.0
        linenumber = Forward_search_last(self.filename, 'Total Mass', 'molecular mass')
        if type(linenumber) == int:
            self.mass = float(self.lines[linenumber].split()[-2])

    def _SymmetryNumber(self):
        self.symnum = 0
        linenumber = Forward_search_last(self.filename, 'Symmetry Number', 'rotational symmetry number')
        if type(linenumber) == int:
            self.symnum = int(self.lines[linenumber].split()[-1])

    def _Multiplicity(self):
        self.multi = 0
        linenumber = Forward_search_first(self.filename, 'Multiplicity', 'multiplicity')
        if type(linenumber) == int:
            self.multi = int(self.lines[linenumber].split()[-1])

    def _PartitionFunctions(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}partition function{Style.RESET_ALL} calculation")
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

    def _Entropy(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}partition function{Style.RESET_ALL} calculation")
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


class dal:
    def __init__(self, filename: str, NeededArguments: dict = None, Quiet: bool = False, Temperature: float = 298.15):
        self.filename = filename
        self.NeededArguments = NeededArguments
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self):
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _Energy(self):
        linenumber = Forward_search_last(self.filename, '@    Final .* energy:', 'final energy', err=False)
        if type(linenumber) == int:
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        linenumber = Forward_search_last(self.filename, '@ Energy at final geometry is', 'final energy')
        if type(linenumber) == int:
            self.tot_energy = float(self.lines[linenumber].split()[-2])
            return
        self.tot_energy = 'NaN'

    def _ZPV(self):
        linenumber = Forward_search_last(self.filename, 'Total Molecular Energy', 'zero-point energy')
        if type(linenumber) == int:
            self.zpv = float(self.lines[linenumber+5].split()[1])
            return
        self.zpv = 'NaN'

    def _Dipole_moments(self):
        linenumber = Forward_search_last(self.filename, 'Dipole moment components', 'dipole moment')
        if type(linenumber) == int:
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[linenumber+5].split()[1]), float(self.lines[linenumber+6].split()[1]), float(self.lines[linenumber+7].split()[1]), float(self.lines[linenumber-3].split()[0])
            return
        self.dipolex = self.dipoley = self.dipolez = self.total_dipole = 'NaN'

    def _Polarizabilities(self):
        linenumber = Forward_search_last(self.filename, 'SECOND ORDER PROPERTIES', 'polarizabilities')
        if type(linenumber) == int:
            self.polx, self.poly, self.polz = float(self.lines[linenumber+2].split()[-1]), float(self.lines[linenumber+5].split()[-1]), float(self.lines[linenumber+7].split()[-1])
            self.iso_polar = (self.polx + self.poly + self.polz)/3.
            return
        self.polx = self.poly = self.polz = self.iso_polar = 'NaN'

    def _Excitation_energies(self):
        self.exc_energies = []
        self.exc_type = None
        linenumber = Forward_search_last(self.filename, '@  Oscillator strengths are dimensionless.', 'excitation energies', err=False)
        if type(linenumber) == int:
            self.exc_type = '.EXCITA'
            for i in self.lines[linenumber+5: self.end]:
                if "@ "in i:
                    self.exc_energies.append(float(i.split()[3])* self.constants.ev_to_au)
                else:
                    break
        linenumbers = Forward_search_all(self.filename, '@ Excitation energy', 'excitation energies')
        if type(linenumbers) == list:
            self.exc_type = 'MCTDHF'
            for i in linenumbers:
                self.exc_energies.append(float(self.lines[i].split()[-2]))
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN'] * abs(self.NeededArguments['_Excitation_energies'])
        if len(self.exc_energies) < self.NeededArguments['_Excitation_energies']:
            self.exc_energies += ['NaN'] * (self.NeededArguments['_Excitation_energies'] - len(self.exc_energies))

    def _Oscillator_strengths(self):
        self.osc_strengths = []
        if self.exc_type == '.EXCITA':
            linenumber = Forward_search_last(self.filename, '@  Oscillator strengths are dimensionless.', 'oscillator strengths')
            if type(linenumber) == int:
                for i in self.lines[linenumber+5: self.end]:
                    if "@ " in i:
                        self.osc_strengths.append(float(i.split()[-1]))
                    else:
                        break
        if self.exc_type == 'MCTDHF':
            linenumbers = Forward_search_all(self.filename, '@ Excitation energy', 'oscillator strengths')
            if type(linenumbers) == list:
                for i in linenumbers:
                    osc = 0
                    for j in self.lines[i:i+15]:
                        if '@ Oscillator strength' in j:
                            osc += float(j.split()[5])**2
                    self.osc_strengths.append(osc**0.5)
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN'] * len(self.exc_energies)
            return
        if len(self.osc_strengths) < self.NeededArguments['_Excitation_energies']:
            self.osc_strengths += ['NaN'] * (self.NeededArguments['_Excitation_energies'] - len(self.osc_strengths))

    def _Frequencies(self):
        self.freq = []
        linenumber = Forward_search_last(self.filename, 'Vibrational Frequencies and IR Intensities', 'frequencies')
        if type(linenumber) == int:
            for i in self.lines[linenumber+7: self.end]:
                if len(i.split()) < 1:
                    break
                self.freq.append(float(i.split()[3]))
        if len(self.freq) == 0:
            self.freq = ['NaN'] * abs(self.NeededArguments['_Frequencies'])
        if len(self.freq) < self.NeededArguments['_Frequencies']:
            self.freq += ['NaN'] * (self.NeededArguments['_Frequencies'] - len(self.freq))

    def _RotationalConsts(self):
        self.rots = []
        linenumbers = Forward_search_last(self.filename, 'Rotational constants', 'rotational constants')
        for i in self.lines[linenumbers+7].split()[:-1]:
            self.rots.append(float(i))
        self.rots = np.array(self.rots) * 1E-3
        self.rots = self.rots[self.rots != 0.0]

    def _Mass(self):
        self.mass = 0.0
        linenumber = Forward_search_last(self.filename, 'Total mass:', 'molecular mass')
        if type(linenumber) == int:
            self.mass = float(self.lines[linenumber].split()[-2])

    #Symmetry checking not implemented by default in Dalton
    #def _SymmetryNumber(self):
    #
    #    self.symnum = 0
    #    linenumber = Forward_search_last(self.file, 'Symmetry Number', 'rotational symmetry number')
    #    if type(linenumber) == int:
    #        self.symnum = int(self.lines[linenumber].split()[-1])

    def _Multiplicity(self):
        self.multi = 0
        linenumber = Forward_search_last(self.filename, 'Spatial symmetry', 'multiplicity')
        if type(linenumber) == int:
            self.multi = int(self.lines[linenumber].split()[2])

    def _PartitionFunctions(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}partition function{Style.RESET_ALL} calculation")
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

    def _Entropy(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}partition function{Style.RESET_ALL} calculation")
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

    def _Enthalpy(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}partition function{Style.RESET_ALL} calculation")
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

    def _Gibbs(self):
        if CheckForOnlyNans(np.array(self.freq)):
            if not(self.quiet):
                print(f"{Fore.RED}No frequencies{Style.RESET_ALL} found in {Style.BRIGHT}{Fore.YELLOW}{self.filename}{Style.RESET_ALL}, skipping {Style.BRIGHT}{Fore.RED}free energy{Style.RESET_ALL} energy calculation")
            self.gibbs = 'NaN'
            return
        self.gibbs = self.enthalpy - self.T*self.entropy / self.constants.au_to_kJmol


class lsdal:
    def __init__(self, filename: str, NeededArguments: dict = None, Quiet: bool = False, Temperature: float = 298.15):
        self.filename = filename
        self.NeededArguments = NeededArguments
        self.quiet = Quiet
        self.T = Temperature
        self.constants = Constants()

        self.ReadFile()

        self.end = len(self.lines)

    def ReadFile(self):
        with open(self.filename, "r") as file:
            self.lines = file.readlines()

    def _Energy(self):
        linenumber = Forward_search_last(self.filename, 'Final .* energy:', 'final energy', err=False)
        if type(linenumber) == int:
            self.tot_energy = float(self.lines[linenumber].split()[-1])
            return
        linenumber = Forward_search_last(self.filename, 'ENERGY SUMMARY', 'final energy')
        if type(linenumber) == int:
            for i in self.lines[linenumber+3:self.end]:
                if 'E: ' in i:
                    self.tot_energy = float(i.split()[-1])
                else:
                    return
        self.tot_energy = 'NaN'

    def _Dipole_moments(self):
        linenumber = Forward_search_last(self.filename, 'Permanent dipole moment', 'dipole moment')
        if type(linenumber) == int:
            self.dipolex, self.dipoley, self.dipolez, self.total_dipole = self.lines[linenumber+9].split()[1], self.lines[linenumber+10].split()[1], self.lines[linenumber+11].split()[1], self.lines[linenumber+3].split()[0]
            return
        self.dipolex = self.dipoley = self.dipolez = self.total_dipole = 'NaN'

    def _Polarizabilities(self):
        linenumber = Forward_search_last(self.filename, '*          POLARIZABILITY TENSOR RESULTS (in a.u.)          *', 'polarizability')
        if type(linenumber) == int:
            self.polx, self.poly, self.polz, self.iso_polar = self.lines[linenumber+10].split()[-3], self.lines[linenumber+11].split()[-2], self.lines[linenumber+12].split()[-1], self.lines[linenumber+14].split()[-1]
            return
        self.polx = self.poly = self.polz = self.iso_polar = 'NaN'

    def _Excitation_energies(self):
        self.exc_energies = []
        linenumber = Forward_search_last(self.filename, '*                   ONE-PHOTON ABSORPTION RESULTS (in a.u.)                  *', 'excitation energies')
        if type(linenumber) == int:
            for i in range(linenumber+8,self.end):
                if len(self.lines[i].split()) < 1:
                    break
                self.exc_energies.append(float(self.lines[i].split()[0]))
        if len(self.exc_energies) == 0:
            self.exc_energies = ['NaN'] * abs(self.NeededArguments['_Excitation_energies']  )
        if len(self.exc_energies) < self.NeededArguments['_Excitation_energies']:
            self.exc_energies += ['NaN'] * (self.NeededArguments['_Excitation_energies'] - len(self.exc_energies))

    def _Oscillator_strengths(self):
        self.osc_strengths = []
        linenumber = Forward_search_last(self.filename, '*                   ONE-PHOTON ABSORPTION RESULTS (in a.u.)                  *', 'oscillator strengths')
        if type(linenumber) == int:
            for i in range(len(self.exc_energies)):
                self.osc_strengths.append(float(self.lines[linenumber+8+i].split()[-1]))
        if len(self.osc_strengths) == 0:
            self.osc_strengths = ['NaN'] * abs(self.NeededArguments['_Excitation_energies'])
        if len(self.osc_strengths) < self.NeededArguments['_Excitation_energies']:
            self.osc_strengths += ['NaN'] * (self.NeededArguments['_Excitation_energies'] - len(self.osc_strengths))
