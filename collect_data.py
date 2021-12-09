
import sys
from typing import Set #For argv
import numpy as np
import argparse
from csv import writer
import os


#*************************** INPUT PARSING ****************************


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
A script designed to make it easier to exctract data from output files

          Currently the output formats supported are
          ------------------------------------------
            -  ORCA
            -  DALTON
            -  GAUSSIAN
            -  LSDALTON

          This script can currently extract the data
          ------------------------------------------
            -  Total energies
            -  Zero-Point Vibrational energies
            -  Enthalpies
            -  Entropies
            -  Gibbs Free energies
            -  Dipole moments
            -  Polarizability
            -  Excitation energies
            -  Oscillator strengths
            -  Frequencies
            -  Partition functions at a given temperature

Though not all data types have been implemented for all of the output formats

All values are extracted as Atomic Units where applicable'''
,epilog='''
The following is not implemented for ORCA
 -  Entropies

The following is not implemented for DALTON
 -  Non DFT total energies
 -  Zero-Point Vibrational energies
 -  Enthalpies
 -  Entropies
 -  Gibbs Free energies
 -  Frequencies
 -  Partition functions with symmetry

The following is not implemented for GAUSSIAN
 -  Entropies
 -  Excitation energies
 -  Oscillator strengths

The following is not implemented for LSDALTON
 -  Probably some total energies
 -  Zero-Point Vibrational energies
 -  Enthalpies
 -  Gibbs Free energies
 -  Entropies
 -  Dipole moments
 -  Polarizabilities
 -  Excitation energies
 -  Oscillator strengths
 -  Frequencies
 -  Partition functions

For help contact
  Theo Juncker von Buchwald
  fnc970@alumni.ku.dk

  Magnus Bukhave Johansen
  qhw298@alumni.ku.dk''')

parser.add_argument('infile', type=str, nargs='+', help='The file(s) to exctract data from', metavar='File')

parser.add_argument('-q', '--quiet', action='store_true', help='Include for the script to stay silent - This will not remove error messages or the printing of data')
parser.add_argument('-s', '--suppress', action='store_true', help='Include to suppress all print statements (including errors) except the data output')
parser.add_argument('-csv', action='store_true', help='Include to write the found values in csv file(s)')
parser.add_argument('-E', '--energy', action='store_true', help='Include to extract the Total Energy')
parser.add_argument('-Z', '--zpv', action='store_true', help='Include to extract the Zero-Point Vibrational Energy')
parser.add_argument('-H', '--enthalpy', action='store_true', help='Include to extract the Enthalpy')
parser.add_argument('-S', '--entropy', action='store_true', help='Include to extract the Entropy')
parser.add_argument('-G', '--gibbs', action='store_true', help='Include to extract the Gibbs Free Energy')
parser.add_argument('-D', '--dipole', action='store_true', help='Include to extract the Dipole Moment')
parser.add_argument('-P', '--polar', action='store_true', help='Include to extract the Polarizability')
parser.add_argument('-X', '--exc', const=-1, type=int, help='Include to exctract the Excitation Energies. Add a number to exctract that amount of Excitation Energies. It will extract all Excitation energies as default',nargs='?')
parser.add_argument('-O', '--osc', action='store_true', help='Include to extract the Oscillator Strengths')
parser.add_argument('-F', '--freq', const=-1, type=int, help='Include to exctract the Frequencies. Add a number to exctract that amount of Frequencies. It will extract all Frequencies as default', nargs='?')
parser.add_argument('-Q', '--partfunc', const=298.15, type=float, help='Include to calculate partition functions. Add a temperature to calculate at.',nargs='?')


#******************************* SETUP ********************************


if __name__ == "__main__":
    args = parser.parse_args()
    input_file = args.infile
    CSV = args.csv

    Arguments = {
        '_Energy' : args.energy,
        '_ZPV' : args.zpv,
        '_Enthalpy' : args.enthalpy,
        '_Entropy' : args.entropy,
        '_Gibbs' : args.gibbs,
        '_Dipole_moments' : args.dipole,
        '_Polarizabilities' : args.polar,
        '_Excitation_energies' : args.exc,
        '_Oscillator_strengths' : args.osc,
        '_Frequencies' : args.freq,
        '_PartitionFunctions' : args.partfunc
    }

    Outputs = {
        '_Energy' : ['tot_energy'],
        '_ZPV' : ['zpv'],
        '_Enthalpy' : ['enthalpy'],
        '_Entropy' : ['entropy'],
        '_Gibbs' : ['gibbs'],
        '_Dipole_moments' : ['dipolex', 'dipoley', 'dipolez', 'total_dipole'],
        '_Polarizabilities' : ['polx', 'poly', 'polz', 'iso_polar'],
        '_Excitation_energies' : ['exc_energies'],
        '_Oscillator_strengths' : ['osc_strengths'],
        '_Frequencies' : ['freq'],
        '_PartitionFunctions' : ['qTotal']
    }

    Header_text = {
        'tot_energy' : 'Energy',
        'zpv' : 'ZPV Energy',
        'enthalpy' : 'Enthalpy',
        'entropy' : 'Entropy',
        'gibbs' : 'Gibbs Free Energy',
        'dipolex' : 'Dipole x',
        'dipoley' : 'Dipole y',
        'dipolez' : 'Dipole z',
        'total_dipole' : 'Total Dipole',
        'polx' : 'Polarizability xx',
        'poly' : 'Polarizability yy',
        'polz' : 'Polarizability zz',
        'iso_polar' : 'Isotropic Polarizability',
        'exc_energies' : 'Exc. energy',
        'osc_strengths' : 'Osc. strength',
        'freq' : 'Frequency',
        'qTotal' : 'Total molar partition function'
    }

    silence = args.quiet
    suppressed = args.suppress

    if Arguments['_Oscillator_strengths'] != False and Arguments['_Excitation_energies'] == None:
        if suppressed == False:
            print('Excitation energies will be found as well, since you are trying to extract oscillator strengths')
        Arguments['_Excitation_energies'] = -1

    if Arguments['_Oscillator_strengths'] == True:
        Arguments['_Oscillator_strengths'] = Arguments['_Excitation_energies']

    if Arguments['_PartitionFunctions'] != None and Arguments['_Frequencies'] == None:
        if suppressed == False:
            print('Frequencies will be found as well, since you are trying to extract partition functions')
        Arguments['_Frequencies'] = -1

    Variable_arrays = dict([item for item in Arguments.items() if type(item[1]) == int])


#******************************* CLASSES ******************************


class gaus:
    def __init__(self,input):
        file = open(input, "r")
        self.lines = file.readlines()
        file.close

    def _Energy(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Sum of electronic and zero-point Energies=" in self.lines[i]:
                self.tot_energy = float(self.lines[i].split()[-1]) - float(self.lines[i-4].split()[-2])
                energy_test = True
                break
        try:
            energy_test
        except NameError:
            if suppressed == False:
                print(f"No Total energy in {infile}")
            self.tot_energy = 'NaN'
            return

    def _ZPV(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Sum of electronic and zero-point Energies=" in self.lines[i]:
                self.zpv_energy = float(self.lines[i].split()[-1])
                zpv_test = True
                break
        try:
            zpv_test
        except NameError:
            if suppressed == False:
                print(f"No ZPV energy in {infile}")
            self.zpv_energy = 'NaN'
            return

    def _Enthalpy(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Sum of electronic and thermal Enthalpies=" in self.lines[i]:
                self.enthalpy = float(self.lines[i].split()[-1])
                enthalpy_test = True
                break
        try:
            enthalpy_test
        except NameError:
            print(f"No Enthalpy in {infile}")
            self.enthalpy = 'NaN'
            return

    def _Gibbs(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Sum of electronic and thermal Free Energies=" in self.lines[i]:
                self.gibbs = float(self.lines[i].split()[-1])
                gibbs_test = True
                break
        try:
            gibbs_test
        except NameError:
            if suppressed == False:
                print(f"No Gibbs energy in {infile}")
            self.gibbs = 'NaN'
            return

    def _Dipole_moments(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Electric dipole moment (input orientation):" in self.lines[i]:
                dipole_test = True
                self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[i+4].split()[1].replace('D','E')), float(self.lines[i+5].split()[1].replace('D','E')), float(self.lines[i+6].split()[1].replace('D','E')), float(self.lines[i+3].split()[1].replace('D','E'))
                break
        try:
            dipole_test
        except NameError:
            if suppressed == False:
                print(f"No dipole moments in {infile}")
            self.dipolex = self.dipoley = self.dipolez = self.total_dipole = 'NaN'

    def _Polarizabilities(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Dipole polarizability, Alpha (input orientation)" in self.lines[i]:
                start_polar = i+2
                end = len(self.lines)
                break
        try:
            start_polar
        except NameError:
            if suppressed == False:
                print(f"No polarizabilities in {infile}")
            self.polx = self.poly = self.polz = self.iso_polar = 'NaN'
            return
        else:
            for i in self.lines[start_polar:end]:
                if " iso " in i:
                    self.iso_polar = float(i.split()[1].replace('D','E'))
                if " xx " in i:
                    self.polx = float(i.split()[1].replace('D','E'))
                if " yy " in i:
                    self.poly = float(i.split()[1].replace('D','E'))
                if " zz " in i:
                    self.polz = float(i.split()[1].replace('D','E'))
                    break

    def _Frequencies(self):
        self.freq = []
        for i in range(len(self.lines)):
            if "Frequencies --" in self.lines[i]:
                for j in self.lines[i].split()[2:]:
                    self.freq.append(float(j)*inv_cm_to_au)
            elif "- Thermochemistry -" in self.lines[i]:
                break
        if len(self.freq) == 0:
            if suppressed == False:
                print(f"No frequencies in {infile}")
            if Arguments['_Frequencies'] == -1:
                self.freq = ['NaN']
            else:
                self.freq = ['NaN'] * Arguments['_Frequencies']

    def _RotationalConsts(self):
        self.rots = []
        for i in range(len(self.lines)):
            if "Rotational constant" in self.lines[i]:
                self.rots = np.array(list(map(float,set(self.lines[i].split()[3:]))))
                self.rots = self.rots[self.rots != 0.0]
            elif "- Thermochemistry -" in self.lines[i]:
                 break

    def _Mass(self):
        self.mass = 0.0
        for i in range(len(self.lines)):
            if "Molecular mass" in self.lines[i]:
                self.mass = float(self.lines[i].split()[2])
                break
        if self.mass == 0.0 and suppressed == False:
            print(f"No molecular mass found in {infile}")

    def _SymmetryNumber(self):
        self.symnum = 0
        for i in range(len(self.lines)):
            if "Rotational symmetry number" in self.lines[i]:
                self.symnum = int(self.lines[i].split()[-1].replace('.',''))
                break
        if self.symnum == 0 and suppressed == False:
            print(f"No rotational symmetry number found in {infile}")

     
    def _PartitionFunctions(self):
        if checkForOnlyNans(np.array(self.freq)) and suppressed == False:
            print(f"No frequencies found in {infile}, skipping partition function calculation")
            self.qTotal = 'NaN'
            return
        self._RotationalConsts()
        self._Mass()
        self._SymmetryNumber()       
        self.qT = trans_const_fac * self.mass ** (1.5) * Arguments['_PartitionFunctions'] ** (2.5)
        if len(self.rots) == 1:
            self.qR = rot_lin_const * Arguments['_PartitionFunctions'] / (self.symnum * self.rots[0])
        else:
            self.qR = rot_poly_const * Arguments['_PartitionFunctions'] ** (1.5) / ( self.symnum * np.prod(np.array(self.rots)) ** (0.5))
        realfreq = np.array([x for x in self.freq if x != 'NaN'])
        realfreq = realfreq[realfreq > 0.0]
        self.qV = np.prod(1 / (1 - np.exp( - vib_const * realfreq / Arguments['_PartitionFunctions'])))
        self.qE = 1 #Good approximation for most closed-shell molecules 
        self.qTotal = self.qT*self.qR*self.qV*self.qE





class orca:
    def __init__(self,input):
        file = open(input, "r")
        self.lines = file.readlines()
        file.close

    def _Energy(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Electronic energy" in self.lines[i]:
                self.tot_energy = float(self.lines[i].split()[-2])
                energy_test = True
                break
        try:
            energy_test
        except NameError:
            if suppressed == False:
                print(f"No final energy in {infile}")
            self.tot_energy = 'NaN'
            return

    def _ZPV(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Electronic energy" in self.lines[i]:
                self.zpv_energy = float(self.lines[i].split()[-2]) + float(self.lines[i+1].split()[-4])
                zpv_test = True
                break
        try:
            zpv_test
        except NameError:
            if suppressed == False:
                print(f"No ZPV energy in {infile}")
            self.zpv_energy = 'NaN'
            return

    def _Enthalpy(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Total Enthalpy" in self.lines[i]:
                self.enthalpy = float(self.lines[i].split()[-2])
                enthalpy_test = True
                break
        try:
            enthalpy_test
        except NameError:
            if suppressed == False:
                print(f"No Enthalpy in {infile}")
            self.enthalpy = 'NaN'
            return

    def _Gibbs(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Final Gibbs free energy" in self.lines[i]:
                self.gibbs = float(self.lines[i].split()[-2])
                gibbs_test = True
                break
        try:
            gibbs_test
        except NameError:
            if suppressed == False:
                print(f"No Gibbs Free energy in {infile}")
            self.gibbs = 'NaN'
            return

    def _Dipole_moments(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Total Dipole Moment" in self.lines[i]:
                dipole_test = True
                self.dipolex, self.dipoley, self.dipolez, self.total_dipole = float(self.lines[i].split()[-3]), float(self.lines[i].split()[-2]), float(self.lines[i].split()[-1]), float(self.lines[i+2].split()[-1])
                break
        try:
            dipole_test
        except NameError:
            if suppressed == False:
                print(f"No dipole moments in {infile}")
            self.dipolex = self.dipoley = self.dipolez = self.total_dipole = 'NaN'

    def _Polarizabilities(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "THE POLARIZABILITY TENSOR" in self.lines[i]:
                start_polar = i
                end = len(self.lines)
                break
        try:
            start_polar
        except NameError:
            if suppressed == False:
                print(f"No polarizabilities in {infile}")
            self.polx = self.poly = self.polz = self.iso_polar = 'NaN'
            return
        for i in range(start_polar,end):
            if "diagonalized tensor" in self.lines[i]:
                self.polx, self.poly, self.polz, self.iso_polar = float(self.lines[i+1].split()[0]), float(self.lines[i+1].split()[1]), float(self.lines[i+1].split()[2]), float(self.lines[i+7].split()[-1])
                break

    def _Excitation_energies(self):
        self.exc_energies = []
        for i in self.lines:
            if "STATE " in i:
                self.exc_energies.append(float(i.split()[3]))
        if len(self.exc_energies) == 0:
            if suppressed == False:
                print(f"No Excitation energies in {infile}")
            if Arguments['_Excitation_energies'] == -1:
                self.exc_energies = ['NaN']
            else:
                self.exc_energies = ['NaN'] * Arguments['_Excitation_energies']

    def _Oscillator_strengths(self):
        self.osc_strengths = []
        for i in range(len(self.lines)):
            if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in self.lines[i]:
                for j in range(len(self.exc_energies)):
                    self.osc_strengths.append(float(self.lines[i+5+j].split()[3]))
        if len(self.osc_strengths) == 0:
            if suppressed == False:
                print(f"No Oscillator strengths in {infile}")
            self.osc_strengths = ['NaN'] * len(self.exc_energies)


    def _Frequencies(self):
        for i in range(len(self.lines)):
            if "VIBRATIONAL FREQUENCIES" in self.lines[i]:
                start_freq = i
                end = len(self.lines)
                break
        try:
            start_freq
        except NameError:
            if suppressed == False:
                print(f"No frequencies in {infile}")
            if Arguments['_Frequencies'] == -1:
                self.freq = ['NaN']
            else:
                self.freq = ['NaN'] * Arguments['_Frequencies']
            return
        else:
            self.freq = []
            for i in range(start_freq,end):
                if "VIBRATIONAL FREQUENCIES" in self.lines[i]:
                    for j in self.lines[i+7: end]:
                        if ": " and " 0.00 " in j:
                            pass
                        elif ": " in j and not " 0.00 " in j:
                            self.freq.append(float(j.split()[1])*inv_cm_to_au)
                        else:
                            break




class dal:
    def __init__(self,input):
        file = open(input, "r")
        self.lines = file.readlines()
        file.close

    def _Energy(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "@    Final DFT energy:" in self.lines[i]:
                self.tot_energy = float(self.lines[i].split()[-1])
                energy_test = True
                break
            if "@ Energy at final geometry is" in self.lines[i]:
                self.tot_energy = float(self.lines[i].split()[-2])
                energy_test = True
                break
        try:
            energy_test
        except NameError:
            if suppressed == False:
                print(f"No final energy in {infile}")
            self.tot_energy = 'NaN'

    def _Dipole_moments(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "Dipole moment\n" in self.lines[i]:
                self.total_dipole = float(self.lines[i+4].split()[0])
                start_dipole = i+4
                end = len(self.lines)
                break
        try:
            start_dipole
        except NameError:
            if suppressed == False:
                print(f"No dipole moments in {infile}")
            self.dipolex = self.dipoley = self.dipolez = self.total_dipole = 'NaN'
        else:
            for i in self.lines[start_dipole: end]:
                if " x " in i:
                    self.dipolex = float(i.split()[1])
                if " y " in i:
                    self.dipoley = float(i.split()[1])
                if " z " in i:
                    self.dipolez = float(i.split()[1])
                    break
            for i in (self.dipolex, self.dipoley, self.dipolez):
                try:
                    i
                except NameError:
                    i = 'NaN'

    def _Polarizabilities(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "SECOND ORDER PROPERTIES" in self.lines[i]:
                start_polar = i
                end = len(self.lines)
                break
        try:
            start_polar
        except NameError:
            if suppressed == False:
                print(f"No polarizabilities in {infile}")
            self.polx = self.poly = self.polz = self.iso_polar = 'NaN'
            return
        else:
            for i in self.lines[start_polar: end]:
                if "XDIPLEN  ; XDIPLEN" in i:
                    self.polx = float(i.split()[-1])
                if "YDIPLEN  ; YDIPLEN" in i:
                    self.poly = float(i.split()[-1])
                if "ZDIPLEN  ; ZDIPLEN" in i:
                    self.polz = float(i.split()[-1])
                    break
            self.iso_polar = (self.polx+self.poly+self.polz)/3.

    def _Excitation_energies(self):
        for i in range(len(self.lines)):
            if "*** ABACUS - Excitation energies (.EXCITA) ***" in self.lines[i]:
                start_exc = i+2
                end = len(self.lines)
                self.exc_type = '.EXCITA'
                break
            if "--- EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) ---" in self.lines[i]:
                start_exc = i+2
                end = len(self.lines)
                self.exc_type = 'MCTDHF'
                break
        try:
            start_exc
        except NameError:
            if suppressed == False:
                print(f"No excitation energies in {infile}")
            if Arguments['_Excitation_energies'] == -1:
                self.exc_energies = ['NaN']
            else:
                self.exc_energies = ['NaN'] * Arguments['_Excitation_energies']
            return
        else:  
            self.exc_energies = []
            if self.exc_type == 'MCTDHF':
                for i in self.lines[start_exc: end]:
                    if "@ Excitation energy" in i:
                        self.exc_energies.append(float(i.split()[-2]))
            elif self.exc_type == '.EXCITA':
                for i in range(start_exc,end):
                    if "@  Oscillator strengths are dimensionless." in self.lines[i]:
                        start_exc = i+5
                        for j in self.lines[i+5: end]:
                            if "@ "in j:
                                self.exc_energies.append(float(j.split()[3])*ev_to_au)
                            else:
                                break

    def _Oscillator_strengths(self):
        self.osc_strengths = []
        if self.exc_type == 'MCTDHF':
            for i in range(len(self.lines)):
                if "@ Oscillator strength (LENGTH)" in self.lines[i]:
                    self.osc_strengths.append(float(self.lines[i].split()[5]))
        if self.exc_type == '.EXCITA':
            for i in range(len(self.lines)):
                if "@  Oscillator strengths are dimensionless." in self.lines[i]:
                    start_osc = i+5
                    end = len(self.lines)
                    for j in self.lines[start_osc: end]:
                        if "@ " in j:
                            self.osc_strengths.append(float(j.split()[-1]))
                        else:
                            break
        if len(self.osc_strengths) == 0:
            if suppressed == False:
                 print(f"No oscillator strengths in {infile}")
            self.osc_strengths = ['NaN'] * len(self.exc_energies)




class lsdal:
    def __init__(self,input):
        file = open(input, "r")
        self.lines = file.readlines()
        file.close

    def _Energy(self):
        for i in range(len(self.lines)-1,-1,-1):
            if "ENERGY SUMMARY" in self.lines[i]:
                start_energy = i+3
                end = len(self.lines)
                break
        try:
            start_energy
        except NameError:
            if suppressed == False:
                print(f"No final energy could be found in {infile}")
            self.tot_energy = 'NaN'
            return
        for i in self.lines[start_energy:end]:
            if "E:" in i:
                self.tot_energy = float(i.split()[-1])
            else:
                break


#**************************** FUNCTIONS *******************************


def resize(array):
    max_size = 0
    for i in range(len(array[:])):
        if len(array[i]) > max_size:
            max_size = len(array[i])

    for i in range(len(array[:])):
        if array[i] == ['Not implemented']:
            array[i] *= max_size
        else:
            array[i] += ['NaN'] * (max_size - len(array[i]))


def checkForOnlyNans(array):
    for i in array:
        if i != 'NaN':
           return False
    return True


#******************************* CODE *********************************

if __name__ == "__main__":
    Wanted_Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]
    Set_of_values = dict(item for item in Outputs.items() if item[0] in Wanted_Values)
    Extracted_values = dict()
    array_input = list()

    ev_to_au = 0.036749405469679
    inv_cm_to_au = 1/219474.63068
    trans_const_fac = 1.5625517342018425307E+22 #Molar value assuming 1 bar standard pressure
    rot_lin_const = 20.83661793 #Assuming rigid, linear rotor and T>>Rotational temperature and rotational constant in GHz
    rot_poly_const = 168.5837766 #Assuming rigid, polyatomic rotor and T>>Rotational temperature and rotational constant in GHz
    vib_const = 3.157750419E+05 #Assuming harmonic oscillator and frequency in au 
    Barrier = '\n**********************************************\n'
    
    if silence == True:
        print('')

    count = len(input_file)

    for infile in input_file:
        array_input.append([infile])
        if not(silence == False or suppressed == False) == True:
            print(Barrier)
            print(f'Collecting data from {infile}')

#   --------- DEFINITION OF WHICH PROGRAM THE OUTPUT IS FROM ----------

        with open(infile,'r') as read:
            lines = read.readlines()[:10]
        if '* O   R   C   A *' in lines[4]: 							#File type = ORCA
            file_text = orca(infile)
            input_type = 'ORCA'       
    
        if '*************** Dalton - An Electronic Structure Program ***************' in lines[3]:	#File type = DALTON
            file_text = dal(infile)
            input_type = 'DALTON'

        if 'Gaussian, Inc.  All Rights Reserved.' in lines[6]:					#File type = GAUSSIAN
            file_text = gaus(infile)
            input_type = 'GAUSSIAN'

        if '**********  LSDalton - An electronic structure program  **********' in lines[2]:	#File type = LSDALTON
            file_text = lsdal(infile)
            input_type = 'LSDALTON'

        del lines

        input_no_ext = infile.replace('.out','')

#   ---------------- EXTRACTING THE APPROPIATE VALUES -----------------

        for i in Wanted_Values:
            try:
                method = getattr(type(file_text),i)
                method(file_text)
            except AttributeError:
                if suppressed == False:
                    print(f'{infile}: {i} has not been implemented for {input_type}')
        
#   -------------- COLLETING THE VALUES IN DICTIONARIES ---------------

        dict_keys = [*file_text.__dict__.keys()]
        collection_dict = dict()
        
        for i in dict_keys[1:]:
            collection_dict[i] = file_text.__dict__[i]

        Extracted_values[infile] = collection_dict

    if not(silence == False or suppressed == False) == True:
        print(Barrier)

#   ---------------- FINDING FUNCTIONS NOT IMPLEMENTED ----------------

    for infile in input_file:
        for key in Set_of_values:
            for val in Set_of_values[key]:
                try:
                    Extracted_values[infile][val]
                except KeyError:
                    Extracted_values[infile][val] = ['Not implemented']

#   ------------------ TURNING EVERYTHING INTO LISTS ------------------

    for key_outer in Extracted_values.keys():
        for key_inner in Extracted_values[key_outer].keys():
            if type(Extracted_values[key_outer][key_inner]) != list:
                Extracted_values[key_outer][key_inner] = [Extracted_values[key_outer][key_inner]]

#   --------- COLLECTING ALL VALUES IN ARRAYS IN A DICTIONARY ---------

    Final_arrays = dict()

    for key in Set_of_values:
        for val in Set_of_values[key]:
            Final_arrays[val] = []
            for infile in input_file:
                Final_arrays[val].append(Extracted_values[infile][val])

#   -------------------------- RESIZES ARRAYS -------------------------

    for key in Final_arrays.keys():
        if type(Final_arrays[key]) == list:
            resize(Final_arrays[key])

#   ------------------- FIXING VARIABLE SIZE ARRAYS -------------------

    for item in Variable_arrays.items():
        if item[1] > 0:
            for val in Outputs[item[0]]:
                for file in range(0,count):
                    Final_arrays[val][file] = Final_arrays[val][file][0:item[1]]

#   --------------------- CREATES THE HEADER ROW ----------------------
    
    header = ['File']
    temp = [val[1] for val in Header_text.items() if val[0] in Set_of_values]
    for key in Set_of_values.keys():
        for val in Set_of_values[key]:
            if len(Final_arrays[val][0]) > 1:
                for i in range(len(Final_arrays[val][0])):
                    header.append(f'{Header_text[val]} {i+1}')
            else:
                header.append(Header_text[val])

#  --- CREATES AN ARRAY OF THE CORRECT SIZE FILLED WITH THE HEADER ---
#  ----------------- USEFUL FOR TROUBLESHOOTING ----------------------

    output_array = np.array([header] * (count + 1), dtype=object)

#   ------------ FILLS THE ARRAY WITH THE CORRECT VALUES -------------

    if count == 1:
        output_array[1,0] = array_input[0][0]
    else:
        output_array[1:,0] = np.array(np.concatenate(array_input))

    col = 1
    for key in Set_of_values.keys():
        for val in Set_of_values[key]:
            output_array[1:,col:col+len(np.array(Final_arrays[val][0]))] = np.array(Final_arrays[val])
            col += len(np.array(Final_arrays[val][0]))

#   ------------ IF CHOSEN PRINTS THE OUTPUT IN A CSV FILE ------------
#   ---------- ELSE THE RESULTS ARE DUMPED INTO THE TERMINAL ----------

    if CSV == True:
        np.savetxt('data.csv', output_array, delimiter=',', fmt='%s')
        print('Data has been saved in data.csv')
    else:
        print(output_array)

    #os.system("sed -i s/$(printf '\r')\$// data.csv") #Removes 
    #(Carriage return) from the data.csv file

#EOF