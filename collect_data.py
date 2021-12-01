
import sys #For argv
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
            -  Total energy
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
parser.add_argument('-s', '--suppress', action='store_true', help='Include to supress all print statements (including errors) except the data output')
parser.add_argument('-csv', action='store_true', help='Include to write the found values in csv file(s)')
parser.add_argument('-E', '--energy', action='store_true', help='Include to extract the Total Energy')
parser.add_argument('-Z', '--zpv', action='store_true', help='Include to extract the Zero-Point Vibrational Energy')
parser.add_argument('-H', '--enthalpy', action='store_true', help='Include to extract the Enthalpy')
parser.add_argument('-S', '--entropy', action='store_true', help='Include to extract the Entropy')
parser.add_argument('-G', '--gibbs', action='store_true', help='Include to extract the Gibbs Free Energy')
parser.add_argument('-D', '--dipole', action='store_true', help='Include to extract the Dipole Moment')
parser.add_argument('-P', '--polar', action='store_true', help='Include to extract the Polarizability')
parser.add_argument('-X', '--exc', const=0, help='Include to exctract the Excitation Energies. Add a number to exctract that amount of Excitation Energies. It will extract all Excitation energies as default',nargs='?')
parser.add_argument('-O', '--osc', action='store_true', help='Include to extract the Oscillator Strengths')
parser.add_argument('-F', '--freq', const=0, help='Include to exctract the Frequencies. Add a number to exctract that amount of Frequencies. It will extract all Frequencies as default', nargs='?')
parser.add_argument('-Q', '--partfunc', const=298.15, help='Include to calculate partition functions. Add a temperature to calculate at.',nargs='?')


args = parser.parse_args()
input_file = args.infile
CSV = args.csv
TOTAL_ENERGY = args.energy
ZPV_ENERGY = args.zpv
ENTHALPY = args.enthalpy
ENTROPY = args.entropy
GIBBS = args.gibbs

DIPOLE_MOMENT = args.dipole
POLARIZABILITY = args.polar

if args.exc == None:
    EXCITATION_ENERGIES = False
    AMOUNT_EX_ENERGIES = 0
else:
   EXCITATION_ENERGIES = True
   AMOUNT_EX_ENERGIES = int(args.exc)

OSCILLATOR_STRENGTHS = args.osc

if args.freq == None:
    FREQUENCIES = False
    AMOUNT_FREQ = 0
else:
    FREQUENCIES = True
    AMOUNT_FREQ = int(args.freq)

if args.partfunc == None:
   PARTITION_FUNC = False 
else:
   PARTITION_FUNC = True
   TEMPERATURE = float(args.partfunc)

silence = args.quiet
suppressed = args.suppress

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
                self.dipolex, self.dipoley, self.dipolez, self.total_dipole = np.float(self.lines[i+4].split()[1].replace('D','E')), np.float(self.lines[i+5].split()[1].replace('D','E')), np.float(self.lines[i+6].split()[1].replace('D','E')), np.float(self.lines[i+3].split()[1].replace('D','E'))
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

    def _Frequencies(self,num):
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
            if num == 0:
                self.freq = ['NaN']
            else:
                self.freq = ['NaN'] * num
     
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

    def _PartitionFunctions(self,T):
        if checkForOnlyNans(np.array(self.freq)) and suppressed == False:
            print(f"No frequencies found in {infile}, skipping partition function calculation")
            self.qTotal = 'NaN'
        else:        
            self.qT = trans_const_fac * self.mass ** (1.5) * T ** (2.5)
            if len(self.rots) == 1:
                self.qR = rot_lin_const * T / (self.symnum * self.rots[0])
            else:
                self.qR = rot_poly_const * T ** (1.5) / ( self.symnum * np.prod(np.array(self.rots)) ** (0.5))
            self.realfreq = np.array([x for x in self.freq if x != 'NaN'])
            self.realfreq = self.realfreq[self.realfreq > 0.0]
            self.qV = np.prod(1 / (1 - np.exp( - vib_const * self.realfreq / T)))
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

    def _Excitation_energies(self,num):
        self.exc_energies = []
        for i in self.lines:
            if "STATE " in i:
                self.exc_energies.append(float(i.split()[3]))
        if len(self.exc_energies) == 0:
            if suppressed == False:
                print(f"No Excitation energies in {infile}")
            if num == 0:
                self.exc_energies = ['NaN']
            else:
                self.exc_energies = ['NaN'] * num

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


    def _Frequencies(self,num):
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
            if num == 0:
                self.freq = ['NaN']
            else:
                self.freq = ['NaN'] * num
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
            if "Dipole moment" in self.lines[i]:
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
            self.dipolex = self.dipoley = self.dipolez = 0.
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
                    i = 'Nan'

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

    def _Excitation_energies(self,num):
        for i in range(len(self.lines)):
            if "*** ABACUS - Excitation energies (.EXCITA) ***" in self.lines[i]:
                start_exc = i+2
                end = len(self.lines)
                self._exc_type = '.EXCITA'
                break
            if "--- EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) ---" in self.lines[i]:
                start_exc = i+2
                end = len(self.lines)
                self._exc_type = 'MCTDHF'
                break
        try:
            start_exc
        except NameError:
            if suppressed == False:
                print(f"No excitation energies in {infile}")
            if num == 0:
                self.exc_energies = ['NaN']
            else:
                self.exc_energies = ['Nan'] * num
            return
        else:  
            self.exc_energies = []
            if self._exc_type == 'MCTDHF':
                for i in self.lines[start_exc: end]:
                    if "@ Excitation energy" in i:
                        self.exc_energies.append(float(i.split()[-2]))
            elif self._exc_type == '.EXCITA':
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
        if self._exc_type == 'MCTDHF':
            for i in range(len(self.lines)):
                if "@ Oscillator strength (LENGTH)" in self.lines[i]:
                    self.osc_strengths.append(float(self.lines[i].split()[5]))
        if self._exc_type == '.EXCITA':
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
        array[i] += ['NaN'] * (max_size - len(array[i]))


def checkForOnlyNans(array):
    for i in array:
        if i != 'NaN':
           return False
    return True

#******************************* CODE *********************************


ev_to_au = 0.036749405469679
inv_cm_to_au = 1/219474.63068
trans_const_fac = 1.5625517342018425307E+22 #Molar value assuming 1 bar standard pressure
rot_lin_const = 20.83661793 #Assuming rigid, linear rotor and T>>Rotational temperature and rotational constant in GHz
rot_poly_const = 168.5837766 #Assuming rigid, polyatomic rotor and T>>Rotational temperature and rotational constant in GHz
vib_const = 3.157750419E+05 #Assuming harmonic oscillator and frequency in au 
Barrier = '\n**********************************************\n'
count = 0
for infile in input_file:
    count += 1

if silence == True:
    print('')

count = 0
for infile in input_file:
    if not(silence == False or suppressed == False) == True:
        print(Barrier)
        print(f'Collecting data from {infile}')
    count += 1

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

#   - COLLECTING THE APPROPIATE VALUES AND STORING THEM IN THE OUTPUT -

    if TOTAL_ENERGY == True:
        try:
            file_text._Energy()
        except AttributeError:
            if suppressed == False:
                print(f"Energy extraction not implemented for {input_type}")
            temp_energy = ['Not implemented']
        else:
            temp_energy = [file_text.tot_energy]

    if ZPV_ENERGY == True:
        try:
            file_text._ZPV()
        except AttributeError:
            if suppressed == False:
                print(f"ZPV Energy extraction not implemented for {input_type}")
            temp_zpv = ['Not implemented']
        else:
            temp_zpv = [file_text.zpv_energy]

    if ENTHALPY == True:
        try:
            file_text._Enthalpy()
        except AttributeError:
            if suppressed == False:
                print(f"Enthalpy extraction not implemented for {input_type}")
            temp_enthalpy = ['Not implemented']
        else:
            temp_enthalpy = [file_text.enthalpy]

    if ENTROPY == True:
        try:
            file_text._Entropy()
        except AttributeError:
            if suppressed == False:
                print(f"Entropy extraction not implemented for {input_type}")
            temp_entropy = ['Not implemented']
        else:
            temp_entropy = [file_text.entropy]

    if GIBBS == True:
        try:
            file_text._Gibbs()
        except AttributeError:
            if suppressed == False:
                print(f"Gibbs Free Energy extraction not implemented for {input_type}")
            temp_gibbs = ['Not implemented']
        else:
            temp_gibbs = [file_text.gibbs]

    if DIPOLE_MOMENT == True:
        try:
            file_text._Dipole_moments()
        except AttributeError:
            if suppressed == False:
                print(f"Dipole moment extraction not implemented for {input_type}")
            temp_dipole = ['Not implemented'] * 4
        else:
            temp_dipole = [file_text.dipolex, file_text.dipoley, file_text.dipolez, file_text.total_dipole]

    if POLARIZABILITY == True:
        try:
            file_text._Polarizabilities()
        except AttributeError:
            if suppressed == False:
                print(f"Polarizability extraction not implemented for {input_type}")
            temp_polar = ['Not implemented'] * 4
        else:
            temp_polar = [file_text.polx, file_text.poly, file_text.polz, file_text.iso_polar]

    if EXCITATION_ENERGIES == True:
        try:
            file_text._Excitation_energies(AMOUNT_EX_ENERGIES)
        except AttributeError:
            if suppressed == False: 
                print(f"Excitation energy extraction not implemented for {input_type}")
            if AMOUNT_EX_ENERGIES == 0:
                temp_exc = ['Not implemented']
            else:
                temp_exc = ['Not implemented'] * AMOUNT_EX_ENERGIES
        else:
            if AMOUNT_EX_ENERGIES == 0:
                temp_exc = file_text.exc_energies
            elif abs(AMOUNT_EX_ENERGIES) > len(file_text.exc_energies):
                excess = AMOUNT_EX_ENERGIES - len(file_text.exc_energies)
                if suppressed == False:
                    print(f"You are trying to extract more Excitation energies than are calculated for {infile}")
                temp_exc = file_text.exc_energies + ['NaN']*excess
            else:
                temp_exc = file_text.exc_energies[0:AMOUNT_EX_ENERGIES]
        if OSCILLATOR_STRENGTHS == True:
            try:
                file_text._Oscillator_strengths()
            except AttributeError:
                if suppressed == False:
                    print(f"Oscillator strength extraction not implemented for {input_type}")
                temp_osc = ['Not implemented'] * AMOUNT_EX_ENERGIES
            else:
                if AMOUNT_EX_ENERGIES == 0:
                    temp_osc = file_text.exc_energies
                elif abs(AMOUNT_EX_ENERGIES) > len(file_text.exc_energies):
                    temp_osc = file_text.osc_strengths + ['NaN']*excess
                else:
                    temp_osc = file_text.exc_energies[0:AMOUNT_EX_ENERGIES]

    if EXCITATION_ENERGIES == False and OSCILLATOR_STRENGTHS == True and suppressed == False:
        print('You cannot exctract Oscillator Strengths without also extracting the Excitation Energies')

    if FREQUENCIES == True:
        try:
            file_text._Frequencies(AMOUNT_FREQ)
        except AttributeError:
            if suppressed == False:
                print(f"Frequency extraction not implemented for {input_type}")
            if AMOUNT_FREQ == 0:
                temp_freq = ['Not implemented']
            else:
                temp_freq = ['Not implemented'] * AMOUNT_FREQ
        else:
            if AMOUNT_FREQ == 0:
                temp_freq = file_text.freq
            elif abs(AMOUNT_FREQ) > len(file_text.freq):
                excess = AMOUNT_FREQ - len(file_text.freq)
                if suppressed == False:
                    print(f"You are trying to extract more Frequencies than are calculated for {infile}")
                temp_freq = file_text.freq + ['NaN']*excess
            else:
                temp_freq = file_text.freq[0:AMOUNT_FREQ]
        if PARTITION_FUNC == True:
            try:
                file_text._RotationalConsts()
                file_text._Mass()
                file_text._SymmetryNumber()
                file_text._PartitionFunctions(TEMPERATURE)
            except AttributeError:
                if suppressed == False:
                    print(f"Partition function calculation not implemented for {input_type}")
                temp_partfunc = ['Not implemented']
            else:
                temp_partfunc = [file_text.qTotal]


    if FREQUENCIES == False and PARTITION_FUNC == True and suppressed == False:
        print('You cannot calculate partition functions without also extracting the vibrational frequencies')


#   ----- CREATES ARRAY CONSISTING OF ALL THE VALUES -----


    if count == 1:
        array_input = [input_no_ext]
        if TOTAL_ENERGY == True:
            array_energy = temp_energy
        if ZPV_ENERGY == True:
            array_zpv = temp_zpv
        if ENTHALPY == True:
            array_enthalpy = temp_enthalpy
        if ENTROPY == True:
            array_entropy = temp_entropy
        if GIBBS == True:
            array_gibbs = temp_gibbs
        if DIPOLE_MOMENT == True:
            array_dipole = temp_dipole
        if POLARIZABILITY == True:
            array_polar = temp_polar
        if EXCITATION_ENERGIES == True:
            array_exc = temp_exc
            if OSCILLATOR_STRENGTHS == True:
                array_osc = temp_osc
        if FREQUENCIES == True:
            array_freq = temp_freq
            if PARTITION_FUNC == True:
                array_part_func = temp_partfunc
    elif count == 2:
        array_input = [array_input, [input_no_ext]]
        if TOTAL_ENERGY == True:
            array_energy = [array_energy, temp_energy]
        if ZPV_ENERGY == True:
            array_zpv = [array_zpv, temp_zpv]
        if ENTHALPY == True:
            array_enthalpy = [array_enthalpy, temp_enthalpy]
        if ENTROPY == True:
            array_entropy = [array_entropy, temp_entropy]
        if GIBBS == True:
            array_gibbs = [array_gibbs, temp_gibbs]
        if DIPOLE_MOMENT == True:
            array_dipole = [array_dipole, temp_dipole]
        if POLARIZABILITY == True:
            array_polar = [array_polar, temp_polar]
        if EXCITATION_ENERGIES == True:
            array_exc = [array_exc, temp_exc]
            if OSCILLATOR_STRENGTHS == True:
                array_osc = [array_osc, temp_osc]
        if FREQUENCIES == True:
            array_freq = [array_freq, temp_freq]
            if PARTITION_FUNC == True:
                array_part_func = [array_part_func,temp_partfunc]
    else:
        array_input = [*array_input, [input_no_ext]]
        if TOTAL_ENERGY == True:
            array_energy = [*array_energy, temp_energy]
        if ZPV_ENERGY == True:
            array_zpv = [*array_zpv, temp_zpv]
        if ENTHALPY == True:
            array_enthalpy = [*array_enthalpy, temp_enthalpy]
        if ENTROPY == True:
            array_entropy = [*array_entropy, temp_entropy]
        if GIBBS == True:
            array_gibbs = [*array_gibbs, temp_gibbs]
        if DIPOLE_MOMENT == True:
            array_dipole = [*array_dipole, temp_dipole]
        if POLARIZABILITY == True:
            array_polar = [*array_polar, temp_polar]
        if EXCITATION_ENERGIES == True:
            array_exc = [*array_exc, temp_exc]
            if OSCILLATOR_STRENGTHS == True:
                array_osc = [*array_osc, temp_osc]
        if FREQUENCIES == True:
            array_freq = [*array_freq, temp_freq]
            if PARTITION_FUNC == True:
                array_part_func = [*array_part_func,temp_partfunc]

    del file_text.lines

if not(silence == False or suppressed == False) == True:
    print(Barrier)


#   ------- EXCTEND ARRAYS THAT MAY VARY IN SIZE --------

if not(silence == False or suppressed == False) == True:
    print('Consolidating data - will soon print')

if suppressed == False:
    print('')

if count > 1:
    if EXCITATION_ENERGIES == True:
        resize(array_exc)
        if OSCILLATOR_STRENGTHS == True:
            resize(array_osc)
    if FREQUENCIES == True:
        resize(array_freq)

#   -------------- CREATES THE HEADER ROW ---------------

header = ['File']
if TOTAL_ENERGY == True:
    header += ['Energy']
if ZPV_ENERGY == True:
    header += ['ZPV Energy']
if ENTHALPY == True:
    header += ['Enthalpy']
if ENTROPY == True:
    header += ['Entropy']
if GIBBS == True:
    header += ['Gibbs Free energy']
if DIPOLE_MOMENT == True:
    header += ['Dipole x', 'Dipole y', 'Dipole z', 'Total dipole']
if POLARIZABILITY == True:
    header += ['Polarizablity xx', 'Polarizability yy', 'Polarizability zz', 'Isotropic polarizability']
if count == 1:
    if EXCITATION_ENERGIES == True:
        header += [f'Exc. energy {i+1}' for i in range(len(array_exc))]
        if OSCILLATOR_STRENGTHS == True:
            header += [f'Osc. strength {i+1}' for i in range(len(array_osc))]
    if FREQUENCIES == True:
        header += [f'Frequency {i+1}' for i in range(len(array_freq))]
        if PARTITION_FUNC == True:
           header += ['Total molar partition function']
else:
    if EXCITATION_ENERGIES == True:
        header += [f'Exc. energy {i+1}' for i in range(len(array_exc[0]))]
        if OSCILLATOR_STRENGTHS == True:
            header += [f'Osc. strength {i+1}' for i in range(len(array_osc[0]))]
    if FREQUENCIES == True:
        header += [f'Frequencies {i+1}' for i in range(len(array_freq[0]))]
        if PARTITION_FUNC == True:
            header += ['Total molar partition function']

#  CREATES AN ARRAY OF THE CORRECT SIZE FILLED WITH THE HEADER
#                 USEFUL FOR TROUBLESHOOTING

output_array = np.array([header] * (count + 1), dtype=object)

#   ------ FILLS THE ARRAY WITH THE CORRECT VALUES ------

if count == 1:
    output_array[1,0] = array_input[0]
else:
    output_array[1:,0] = np.array(np.concatenate(array_input))

col = 1
if TOTAL_ENERGY == True:
    output_array[1:,col] = np.array(array_energy).T
    col += 1
if ZPV_ENERGY == True:
    output_array[1:,col] = np.array(array_zpv).T
    col += 1
if ENTHALPY == True:
    output_array[1:,col] = np.array(array_enthalpy).T
    col += 1
if ENTROPY == True:
    output_array[1:,col] = np.array(array_entropy).T
    col += 1
if GIBBS == True:
    output_array[1:,col] = np.array(array_gibbs).T
    col += 1
if DIPOLE_MOMENT == True:
    output_array[1:,col:col+4] = np.array(array_dipole)
    col += 4
if POLARIZABILITY == True:
    output_array[1:,col:col+4] = np.array(array_polar)
    col += 4
if EXCITATION_ENERGIES == True:
    if count == 1:
        output_array[1:,col:col+len(np.array(array_exc))] = np.array(array_exc)
        col += len(np.array(array_exc))
        if OSCILLATOR_STRENGTHS == True:
            output_array[1:,col:col+len(np.array(array_osc))] = np.array(array_osc)
    else:
        output_array[1:,col:col+len(np.array(array_exc)[0])] = np.array(array_exc)
        col += len(np.array(array_exc[0]))
        if OSCILLATOR_STRENGTHS == True:
            output_array[1:,col:col+len(np.array(array_osc)[0])] = np.array(array_osc) 
if FREQUENCIES == True:
    if count == 1:
        output_array[1:,col:col+len(np.array(array_freq))] = np.array(array_freq)
        col += len(np.array(array_freq))
    else:
        output_array[1:,col:col+len(np.array(array_freq)[0])] = np.array(array_freq)
        col += len(np.array(array_freq[0]))
    if PARTITION_FUNC == True:
        output_array[1:,col] = np.array(array_part_func).T
        col+=1


#   ----- IF CHOSEN PRINTS THE OUTPUT IN A CSV FILE -----
#   --- ELSE THE RESULTS ARE DUMPED INTO THE TERMINAL ---

if CSV == True:
    np.savetxt('data.csv', output_array, delimiter=',', fmt='%s')
    print('Data has been saved in data.csv')
else:
    print(output_array)

#os.system("sed -i s/$(printf '\r')\$// data.csv") #Removes  #(Carriage return) from the data.csv file

#EOF

