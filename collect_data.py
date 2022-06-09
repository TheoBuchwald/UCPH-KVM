
import argparse
import os
import time
import numpy as np
from KurtGroup.Kurt import output_processing as op
from functools import partial
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
from matplotlib import rc
from types import FunctionType
import json
import copy


#*************************** INPUT PARSING ****************************


class TerminalInformation():

    def __init__(self, total_nr: int, max_filename_length: int = None):
        self.full_bar = "▇"
        self.empty_bar = " "
        self.total_nr = total_nr
        self.max_filename_length = max_filename_length
        self.terminal_width, _ = os.get_terminal_size()

        self.UP = lambda x : f"\x1B[{x}A"
        self.CLR = "\x1B[0K"
        self.CLRnl = "\x1B[0K\n"
        self.additionalLines = 0

        self.setProgressbarWidth(True)
        self.updateProgressbar(0, False, False)

    def setProgressbarWidth(self, timing: bool):
        self.progressbar_scale_factor = self.terminal_width
        self.progressbar_scale_factor -= 2 # Corresponds to removing the sidebars
        self.progressbar_scale_factor -= len(str(self.total_nr)) * 2 - 2 # Removes the counter
        if timing:
            self.progressbar_scale_factor -= 23 # Removes the timing parts except the average files pr. second time
            self.progressbar_scale_factor -= 5 # Estimate for average files pr. second time
        self.progressbar_scale_factor -= 10 # Additional just in case

    def updateProgressbar(self, nr: int, print_filename: bool, timing: bool, *, filename: str = None):
        self.progress = nr/self.total_nr*self.progressbar_scale_factor
        self.additionalLines = 0

        self.bar = int(self.progress)*self.full_bar + (self.progressbar_scale_factor-int(self.progress))*self.empty_bar

        self.progressbar = f"|{self.bar}| {nr}/{self.total_nr}"

        if timing:
            self.progressbar += self.timer(nr)

        if print_filename:
            if nr == 1:
                print("")
            self.additionalLines += 1
            self.progressbar = f"<<< {filename:{self.max_filename_length}} >>>{self.CLRnl}{self.progressbar}"

        lines = self.progressbar.split('\n')
        if isinstance(lines, list) and len(lines) > 1:
            self.progressbar = f"{self.UP(len(lines)-1)}{self.progressbar}"

        print(self.progressbar + self.CLR, end="\r")

    def start_timer(self):
        start_time = time.perf_counter()
        self.times = [start_time]
        self.time_differences = []

    def timer(self, nr):
        self.times.append(time.perf_counter())

        self.time_differences.append(self.times[-1] - self.times[-2])

        minutes = lambda x : x // 60
        seconds = lambda x : x - x // 60

        if len(self.times) > 1:
            time_average = np.mean(self.time_differences)
            time_spent = sum(self.time_differences)
            time_left = time_average*(self.total_nr - nr)
            return f"; {minutes(time_spent):02.0f}:{seconds(time_spent):02.0f}<{minutes(time_left):02.0f}:{seconds(time_left):02.0f}; {1/time_average:.2f} files/s"
        else:
            return "; 00:00<00:00; 0 files/s"

# Stolen from https://geekflare.com/flatten-list-python/ and modified slightly
def flatten_list(data, flat_list):
    # iterating over the data
    for element in data:
        # checking for list
        if type(element) == list:
            # calling the same function with current element as new argument
            flatten_list(element, flat_list)
        else:
            flat_list.append(element)

def Resize(array: list) -> None:
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

def Data_Extraction(infile, Needed_Values: dict, quiet: bool = False, Temperature: float = 298.15) -> dict:
    Extracted_values = dict()

    infile = op.OutputType(str(infile), Quiet=quiet, Temperature=Temperature)

    # Extracting data
    Extract_data(quiet, Needed_Values, infile.filename, infile.extract, infile.input)

    # Removing the file text from memory as this is no longer needed (and fills in memory)
    infile.extract.__delattr__('lines')

    # List of all dictionary keys for infile.extract
    dict_keys = [*infile.extract.__dict__.keys()]

    # Dictionary
    collection_dict = dict()

    # Collecting the data in dictionaries
    for i in dict_keys[1:]:
        collection_dict[i] =  infile.extract.__dict__[i]

    # Assigning to the Extracted_values dictionary with the filename as key so all data can be easily found in the future
    Extracted_values[infile.filename] = collection_dict

    return Extracted_values

def Extract_data(suppressed: bool, Wanted_Values: dict, infile: str, file_text: dict, input_type: str) -> None:
    # Loops over all requested values and runs the corresponding function
    # If the function has not been implemented it will print an error message
    for i in Wanted_Values:
        try:
            method = getattr(type(file_text),i)
            method(file_text)
        except AttributeError:
            if not(suppressed):
                with open("collect_data.log", "a") as logfile:
                    logfile.write(f'{infile}: {i} has not been implemented for {input_type}\n')

def Check_if_Implemented(input_file: dict, Set_of_values: dict, Extracted_values: dict) -> None:
    # Checks to see if the keys of a double dictionary exists
    # If they don't it is assumed that the function related to the data hasn't been implemented
    for infile in input_file:
        for values in Set_of_values.values():
            for val in values:
                try:
                    Extracted_values[infile][val]
                except KeyError:
                    Extracted_values[infile][val] = ['Not implemented']

def Collect_and_sort_data(input_file: str, Set_of_values: dict, Extracted_values: dict) -> dict:
    Final_arrays = dict()

    # All data from Extracted_values is sorted by the Set_of_values
    # This way excess values in Extracged_values are ignored
    for key in Set_of_values:
        for val in Set_of_values[key]:
            Final_arrays[val] = []
            for infile in input_file:
                Final_arrays[val].append(Extracted_values[infile][val])
    return Final_arrays

def Downsizing_variable_arrays(Outputs: dict, Variable_arrays: dict, count: int, Final_arrays: dict) -> None:
    # Downsizes arrays if requested
    # If not requested nothing happens
    for item in Variable_arrays.items():
        if item[1] > 0:
            for val in Outputs[item[0]]:
                for file in range(0,count):
                    Final_arrays[val][file] = Final_arrays[val][file][0:item[1]]

def Upsizing_variable_arrays(Outputs: dict, Variable_arrays: dict, count: int, Final_arrays: dict, Arguments: dict) -> None:
    # Upsizes arrays if requested
    # If not requested nothing happens
    for key, arg in Variable_arrays.items():
        for val in Outputs[key]:
            for file in range(0,count):
                Final_arrays[val][file] += ['NaN'] * (arg-len(Final_arrays[val][file]))


def Create_Header(Header_text: dict, Set_of_values: dict, Final_arrays: dict) -> list:
    header = ['File']
    # Adds to the header row all relevant headers for the data-points requested
    for key in Set_of_values.keys():
        for val in Set_of_values[key]:
            if len(Final_arrays[val][0]) > 1:
                for i in range(len(Final_arrays[val][0])):
                    header.append(f'{Header_text[val]} {i+1}')
            else:
                header.append(Header_text[val])
    return header

def Fill_output_array(Set_of_values: dict, array_input: dict, count: int, Final_arrays: dict, output_array: list) -> None:
    # Fills the output array with the header row
    if count == 1:
        output_array[1,0] = array_input[0][0]
    else:
        output_array[1:,0] = np.array(np.concatenate(array_input))

    # Fills the output array with the requested data
    # As it is not known the length of all lists in the Final_arrays the col is used to slowly move through the array
    col = 1
    for key in Set_of_values.keys():
        for val in Set_of_values[key]:
            output_array[1:,col:col+len(np.array(Final_arrays[val][0]))] = np.array(Final_arrays[val])
            col += len(np.array(Final_arrays[val][0]))

def Make_complex_propagator_spectrum(input_file: list, suppressed: bool, Format: str, Extracted_Values: dict, SAVE: bool = True) -> None:
    # A LOT OF PLOT SETUP
    rc('text', usetex=True)
    xlabel_font = ylabel_font = title_font = 16
    plt.rc('font', size=12) # x/y axis font size
    NA=6.02214199E23 #a vogadros number
    c=299792458 # speed of light
    eps0=8.8541878176E-12
    aufreq=4.134137334E16
    hartreetohz=6.579683920502E15
    Save_Dict = dict()
        # PLOT SETUP DONE
    for file in input_file:
        filename = file.replace('.out','')
        plotname = f'{filename}-complex-propagator.{Format}'
        title = filename.replace("_"," ")

        frequencies = np.array([x[0] for x in Extracted_Values[file]['complex_propagator'] if x != 'NaN' and x != 'Not implemented'])
        imaginary_part = np.array([x[2] for x in Extracted_Values[file]['complex_propagator'] if x != 'NaN' and x != 'Not implemented'])

        if len(frequencies) == 0:
            if not(suppressed):
                print(f'Polarizability dampening frequencies have either not been implemented for this output type, or none were found in the file {filename}')
            continue

        span = c/(frequencies*hartreetohz)*10**9
        graph = (frequencies*aufreq*10000)/(c*eps0) * (imaginary_part*1.64877727436*10**(-41)*NA)/(2.3*1000)

        plt.title(title)
        plt.plot(span, graph)
        plt.ylim((0,max(graph)*1.2))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel('Wavelength $(nm)$', fontsize=xlabel_font)
        plt.ylabel('Extinction coefficient', fontsize=ylabel_font)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.savefig(plotname, format=f'{Format}', dpi=600)
        plt.close()

        if SAVE:
            Save_Dict[file] = [span,graph]
    if SAVE:
        np.savez('complex_propagator.npz', **Save_Dict)
        print(f'Complex propagator spectrum data has been saved in complex_propagator.npz')

def UVVIS_Spectrum(t: list, l: list, f: list, k: float, sigmacm: float):
    lambda_tot = np.zeros(len(t))
    for x in range(1,len(t)):
        lambda_tot[x] = sum((k/sigmacm)*f*np.exp(-4*np.log(2)*((1/t[x]-1/l)/(1E-7*sigmacm))**2))
    return lambda_tot

def Make_uvvis_spectrum(input_file: list, suppressed: bool, UVVIS_Spectrum: FunctionType, Format: str, Extracted_Values: dict, SAVE: bool = True) -> None:
    # A LOT OF PLOT SETUP
    rc('text', usetex=True)
    xlabel_font = ylabel_font = title_font = 16
    plt.rc('font', size=12) # x/y axis font size
    N=1000 # number of calculated points in curve
    NA=6.02214199*10**23 #avogadros number
    c=299792458 #speed of light
    e=1.60217662*10**(-19) #electron charge
    me=9.10938*10**(-31) #electron mass
    pi=np.pi
    epsvac=8.8541878176*10**(-12)
    sigmacm=0.4*8065.544
    k=(NA*e**2)/(np.log(10)*2*me*c**2*epsvac)*np.sqrt(np.log(2)/pi)*10**(-1)
    inv_cm_to_au = 1/219474.63068
    Save_Dict = dict()
        # PLOT SETUP DONE
    for file in input_file:
        filename = file.replace('.out','')
        plotname = f'{filename}-uvvis.{Format}'
        title = filename.replace("_"," ")

        excitations = np.array([x for x in Extracted_Values[file]['exc_energies'] if x != 'NaN' and x != 'Not implemented'])
        oscillations = np.array([x for x in Extracted_Values[file]['osc_strengths'] if x != 'NaN' and x != 'Not implemented'])

        if len(excitations) == 0 and len(oscillations) == 0:
            if not(suppressed):
                print(f'Excitation energies and oscillator strengths have either not been implemented for this output type, or none were found in the file {filename}')
            continue
        elif len(excitations) == 0:
            if not(suppressed):
                print(f'Excitation energies have either not been implemented for this output type, or none were found in the file {filename}')
            continue
        elif len(oscillations) == 0:
            if not(suppressed):
                print(f'Oscillator strengths have either not been implemented for this output type, or none were found in the file {filename}')
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
        plt.savefig(plotname, format=f'{Format}', dpi=600)
        plt.close()

        if SAVE:
            Save_Dict[file] = [span,graph]
    if SAVE:
        np.savez('UVVIS.npz', **Save_Dict)
        print(f'UVVIS spectrum data has been saved in UVVIS.npz')

    del Save_Dict # Deletes dictionary from memory since it is no longer needed

def Spectra(args):
    """
    This function is used to do any methods related to the Spectra keyword
    """
    InputFiles = args.infile
    UVVis = args.uvvis
    ComplexPropagator = args.complex_propagator
    Format = args.format
    Save = args.save
    Quiet = args.quiet
    Multiprocessing = args.multiprocessing

    if UVVis:
        NeededArguments = {'_Excitation_energies': -1, '_Oscillator_strengths': -1}
        ArgumentsToValues = {'_Excitation_energies': ['exc_energies'], '_Oscillator_strengths': ['osc_strengths']}

        if Multiprocessing:
            with Pool(int(cpu_count()/2)) as pool:
                ExtractedValues = pool.map(partial(Data_Extraction, Needed_Values=NeededArguments, quiet=Quiet), InputFiles)
                ExtractedValues = {key: value for dictionary in ExtractedValues for key, value in dictionary.items()} # Reformatting Extracted_values
        else:
            ExtractedValues = dict()
            for infile in InputFiles:
                ExtractedValues[infile] = Data_Extraction(infile, NeededArguments, Quiet)[infile]

        Check_if_Implemented(InputFiles, ArgumentsToValues, ExtractedValues)   # Finding functions not implemented

        for OuterKey in ExtractedValues:   # Turn everything into lists
            for InnerKey, Value in ExtractedValues[OuterKey].items():
                if not isinstance(Value, list):
                    ExtractedValues[OuterKey][InnerKey] = [Value]

        Make_uvvis_spectrum(InputFiles, Quiet, UVVIS_Spectrum, Format, ExtractedValues, Save)
        return

    elif ComplexPropagator:
        NeededArguments = {'_Complex_propagator': True}
        ArgumentsToValues = {'_Complex_propagator': ['complex_propagator']}

        if Multiprocessing:
            with Pool(int(cpu_count()/2)) as pool:
                ExtractedValues = pool.map(partial(Data_Extraction, Needed_Values=NeededArguments, quiet=Quiet), InputFiles)
                ExtractedValues = {key: value for dictionary in ExtractedValues for key, value in dictionary.items()} # Reformatting Extracted_values
        else:
            ExtractedValues = dict()
            for infile in InputFiles:
                ExtractedValues[infile] = Data_Extraction(infile, NeededArguments, Quiet)[infile]

        Check_if_Implemented(InputFiles, ArgumentsToValues, ExtractedValues)   #Finding functions not implemented

        for OuterKey, Dict in ExtractedValues.items():   # Turn everything into lists
            for InnerKey, Value in Dict.items():
                if not isinstance(Value, list):
                    ExtractedValues[OuterKey][InnerKey] = [Value]

        Make_complex_propagator_spectrum(InputFiles, Quiet, Format, ExtractedValues, Save)
        return

    print('Doing nothing - please use --uvvis or --complex-propagator')



def Extract(args):
    """
    This function is used for any methods related to the extract keyword
    """

    # Setting the input_files variable to all files
    InputFiles = args.infile

    # These are all the possible arguments that extract data
    # They are here set to correspond the argument passed
    RequestedArguments = {
        '_Energy': args.energy,
        '_ZPV': args.zpv,
        '_Dipole_moments': args.dipole,
        '_Polarizabilities': args.polar,
        '_Excitation_energies': args.exc,
        '_Oscillator_strengths': args.osc,
        '_Frequencies': args.freq,
        '_Enthalpy': args.enthalpy,
        '_Entropy': args.entropy,
        '_Gibbs': args.gibbs,
        '_PartitionFunctions': args.partfunc,
        '_CPUS': args.cpu_time,
        '_Optimized_Geometry': args.optgeom
    }

    # These are the datapoints that will be extracted per argument
    # The value for each key is how that value is saved in the Classes
    Outputs = {
        '_Energy': ['tot_energy'],
        '_ZPV': ['zpv'],
        '_Enthalpy': ['enthalpy'],
        '_Entropy': ['entropy'],
        '_Gibbs': ['gibbs'],
        '_Dipole_moments': ['dipolex', 'dipoley', 'dipolez', 'total_dipole'],
        '_Polarizabilities': ['polx', 'poly', 'polz', 'iso_polar'],
        '_Excitation_energies': ['exc_energies'],
        '_Oscillator_strengths': ['osc_strengths'],
        '_Frequencies': ['freq'],
        '_PartitionFunctions': ['qTotal'],
        '_CPUS': ['total_cpu_time', 'wall_cpu_time'],
    }

    # These are what will be written in the header for each data-point
    # The keys are what the data-points are saved as in the classes
    # The values are what the data-points should be written under in the output header
    HeaderText = {
        'tot_energy': 'Energy',
        'zpv': 'ZPV EnergTotaly',
        'enthalpy': 'Enthalpy',
        'entropy': 'Entropy (kJ/(mol*K)',
        'gibbs': 'Gibbs Free Energy',
        'dipolex': 'Dipole x',
        'dipoley': 'Dipole y',
        'dipolez': 'Dipole z',
        'total_dipole': 'Total Dipole',
        'polx': 'Polarizability xx',
        'poly': 'Polarizability yy',
        'polz': 'Polarizability zz',
        'iso_polar': 'Isotropic Polarizability',
        'exc_energies': 'Exc. energy',
        'osc_strengths': 'Osc. strength',
        'freq': 'Frequency',
        'qTotal': 'Total molar partition function',
        'total_cpu_time': f'Total CPU time ({args.cpu_time})',
        'wall_cpu_time': f'Wall CPU time ({args.cpu_time})'
    }

    # Remainder of arguments
    T = args.temp
    Save = args.save
    Quiet = args.quiet
    Multiprocessing = args.multiprocessing
    ProgressBar = args.progressbar
    UnitTesting = args.unittest
    SaveName = args.savename

    # Making a copy of RequestedArguments
    # This is so arguments that are dependent on others can be called independently
    NeededArguments = RequestedArguments.copy()

    # Ensuring that excitation energies are calculated when oscillator strengths are
    if NeededArguments['_Oscillator_strengths'] and NeededArguments['_Excitation_energies'] == None:
        NeededArguments['_Oscillator_strengths'] = NeededArguments['_Excitation_energies'] = -1
    elif NeededArguments['_Oscillator_strengths']:
        NeededArguments['_Oscillator_strengths'] = NeededArguments['_Excitation_energies']

    # Ensuring that enthalpies and entropies are calculated as these are needed to calculate the Gibbs free energy
    if NeededArguments['_Gibbs']:
        NeededArguments['_Enthalpy'] = True
        NeededArguments['_Entropy'] = True

    # Ensuring that frequencies are calculated as these are needed to calculate the Partition Functions
    if NeededArguments['_PartitionFunctions']:
        NeededArguments['_Frequencies'] = -1

    # Ensuring that frequencies are calculated as these are needed to calculate the enthalpy
    if NeededArguments['_Enthalpy']:
        NeededArguments['_Frequencies'] = -1
        NeededArguments['_Energy'] = True

    # Ensuring that frequencies are calculated as these are needed to calculate the entropy
    if NeededArguments['_Entropy']:
        NeededArguments['_Frequencies'] = -1

    # Dictionary of data where the amount of values printed can be changed
    # Examples of this are the Excitation energies and the Frequencies
    VariableArrays = dict([item for item in NeededArguments.items() if type(item[1]) == int])

    # List of arguments that have been requested
    WantedValues = [key for key, val in RequestedArguments.items() if not(val == None or val == False)]

    # List of arguments that are needed for those requested
    # This is done so you can request as an example the Gibbs free energies withou also having the enthalpies and entropies printed
    NeededValues = [key for key, val in NeededArguments.items() if not(val == None or val == False)]

    # Data-points that will be written to the terminal or save file
    # These are found from the Outputs dictionary by comparing with the Wanted_Values list
    ArgumentsToValues = {key: val for key, val in Outputs.items() if key in WantedValues}

    Values = []
    flatten_list([val for val in ArgumentsToValues.values()], Values)

    # How many files to run the script on
    Count = len(InputFiles)

    # If multiprocessing is enabled it will be run using half of the available CPUS
    # Else they will be run in a linear fashion
    if ProgressBar:
        max_filename_length = len(max(InputFiles, key=len))
        TerminalOutput = TerminalInformation(Count, max_filename_length)
        TerminalOutput.start_timer()
    if Multiprocessing:
        with Pool(int(cpu_count()/2)) as pool:
            ExtractedValues = []
            for i, result in enumerate(pool.imap(partial(Data_Extraction, Needed_Values=NeededValues, quiet=Quiet, Temperature=T), InputFiles), start=1):
                if ProgressBar:
                    TerminalOutput.updateProgressbar(i, False, True)
                ExtractedValues.append(result)
            ExtractedValues = {key: value for dictionary in ExtractedValues for key, value in dictionary.items()} # Reformatting Extracted_values
    else:
        ExtractedValues = dict()
        for i, file in enumerate(InputFiles, start=1):
            if ProgressBar:
                TerminalOutput.updateProgressbar(i, True, True, filename=file)
            ExtractedValues[file] = Data_Extraction(file, NeededValues, Quiet, T)[file]

    # Creating Input_Array where all values are put in lists
    InputArray = [[i] for i in ExtractedValues]

    # Checking if some functions have not been implemented for the relevant extraction types
    Check_if_Implemented(InputFiles, ArgumentsToValues, ExtractedValues)

    # If the CPU time has been requested in second or hours instead of minutes
    # then the calculation from minutes to the requested unit is done here
    if RequestedArguments['_CPUS'] == 's':
        for file in InputFiles:
            if ExtractedValues[file]['total_cpu_time'] != ['Not Implemented'] or ExtractedValues[file]['total_cpu_time'] != ['Nan']:
                continue
            ExtractedValues[file]['total_cpu_time'] *= 60
            ExtractedValues[file]['wall_cpu_time'] *= 60
    elif NeededArguments['_CPUS'] == 'h':
        for file in InputFiles:
            if ExtractedValues[file]['total_cpu_time'] != ['Not Implemented'] or ExtractedValues[file]['total_cpu_time'] != ['Nan']:
                continue
            ExtractedValues[file]['total_cpu_time'] /= 60
            ExtractedValues[file]['wall_cpu_time'] /= 60

    # This is done purely for the unittest script to work correctly
    if UnitTesting:
        SaveDict = copy.deepcopy(ExtractedValues)
        for key_outer, dictionary in ExtractedValues.items():
            for key_inner in dictionary:
                if key_inner not in Values:
                    SaveDict[key_outer].pop(key_inner)
        return SaveDict

    # If something is at this point not in a list somehow they will be after this
    for OuterKey, Dict in ExtractedValues.items():
        for InnerKey, Value in Dict.items():
            if not isinstance(Value, list):
                ExtractedValues[OuterKey][InnerKey] = [Value]

    # Collecting all values in arrays in a dictionary
    # Some values in the Extracted_Values dictionary may not have been requested, so these are removed here
    FinalArrays = Collect_and_sort_data(InputFiles, ArgumentsToValues, ExtractedValues)

    # Resizing arrays
    # An example is Excitation energies where there may be more of them in one output file than another
    # By doing this it fits properly in what is printed to the terminal
    for key in FinalArrays:
        if isinstance(FinalArrays[key], list):
            Resize(FinalArrays[key])

    # Fixing the size of variable size arrays so that they match what was requested
    Downsizing_variable_arrays(Outputs, VariableArrays, Count, FinalArrays)
    Upsizing_variable_arrays(Outputs, VariableArrays, Count, FinalArrays, RequestedArguments)

    # Creation of header row
    Header = Create_Header(HeaderText, ArgumentsToValues, FinalArrays)

    # Create output array from header row
    OutputArray = np.array([Header] * (Count + 1), dtype=object)

    # Filling the output array with the extracted data
    Fill_output_array(ArgumentsToValues, InputArray, Count, FinalArrays, OutputArray)

#   ------------ IF CHOSEN PRINTS THE OUTPUT IN A CSV FILE ------------
#   ---------- ELSE THE RESULTS ARE DUMPED INTO THE TERMINAL ----------

    if ProgressBar:
        print("")

    # If this statement is true, then only the filenames have been written to the Output_Array
    if len(OutputArray) == OutputArray.size:
       print("No data was extracted, therefore nothing more will be printed")
       return

    elif Save == 'return':
        SaveDict = copy.deepcopy(ExtractedValues)
        for key_outer, dictionary in ExtractedValues.items():
            for key_inner in dictionary:
                if key_inner not in Values:
                    SaveDict[key_outer].pop(key_inner)
        return SaveDict

    elif Save == 'csv':
        np.savetxt(f'{SaveName}.csv', OutputArray, delimiter=',', fmt='%s')
        print(f'Data has been saved in {SaveName}.csv')
        return

    elif Save == 'npz':
        SaveDict = {i[0]: i[1:] for i in OutputArray}
        np.savez(f'{SaveName}.npz', **SaveDict)
        print(f'Data has been saved in {SaveName}.npz')
        return

    elif Save == 'json':
        SaveDict = copy.deepcopy(ExtractedValues)
        for key_outer, dictionary in ExtractedValues.items():
            for key_inner in dictionary:
                if key_inner not in Values:
                    SaveDict[key_outer].pop(key_inner)
        json_object = json.dumps(SaveDict, indent=4)

        with open(f"{SaveName}.json", "w") as outfile:
            outfile.write(json_object)
        print(f"Data has been saved in {SaveName}.json")

        return

    print(OutputArray)


def main():
    #---------------------------
    # Creating main parser
    #---------------------------
    Parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=f'''
    A script designed to make it easier to extract data from output files

    To extract data you have to use the keyword extract before any other keywords
    To make spectra you will have to use the keyword spectra before any other keywords

            Currently the output formats supported are
            ------------------------------------------
                -  ORCA
                -  DALTON
                -  GAUSSIAN
                -  LSDALTON
                -  VELOXCHEM
                -  AMSTERDAM MODELING SUITE
''', epilog=f'''
For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk

    Magnus Bukhave Johansen
    qhw298@alumni.ku.dk
''')

    #---------------------------
    # Creating spectra subparser
    #---------------------------
    subparser = Parser.add_subparsers(dest='pars')
    subparser.required = False

    SpectraSubparser = subparser.add_parser('spectra', formatter_class=argparse.RawDescriptionHelpFormatter, description=f'''
    This part of the script is used for making spectra from data collected from output files

            Currently the output formats supported are
            ------------------------------------------
                -  ORCA
                -  DALTON
                -  GAUSSIAN
                -  LSDALTON
                -  VELOXCHEM
                -  AMSTERDAM MODELING SUITE

    It is currently possible to make UVVIS spectra using excitation energies and complex propagator theory

    The following is not implemented for ORCA
    -  UVVIS based on complex propagator theory

    The following is not implemented for DALTON
    -  UVVIS based on complex propagator theory

    The following is not implemented for GAUSSIAN
    -  UVVIS based on complex propagator theory

    The following is not implemented for VELOXCHEM
    -  UVVIS based on excitation energies
    -  UVVIS based on complex propagator theory

    The following is not implemented for AMSTERDAM MODELING SUITE
    -  UVVIS based on excitation energies
    -  UVVIS based on complex propagator theory
''', help='Use to make spectra such as UVVIS from excitation energies or complex propagator theory')

    # Setting the Spectra function to be run if spectra is used
    SpectraSubparser.set_defaults(func=Spectra)

    # Adding arguments
    SpectraSubparser.add_argument('infile', type=str, nargs='+', help='The file(s) to make spectra from', metavar='File')

    SpectraGroup = SpectraSubparser.add_mutually_exclusive_group()
    SpectraGroup.add_argument('--uvvis', action='store_true', help='Include to make UVVIS spectra')
    SpectraGroup.add_argument('--complex-propagator', action='store_true', help='Include to use complex propagator theory to make the spectra')

    SpectraDataProcessingGroup = SpectraSubparser.add_argument_group('Data processing commands')
    SpectraDataProcessingGroup.add_argument('--format', default='png', const='png', type=str, help='Include this to change the picture format. Will use png as default. If \'--save\' is used together with this, the processed data will be saved in a .npz file', nargs='?', choices=['png', 'eps', 'pdf', 'svg', 'ps'])
    SpectraDataProcessingGroup.add_argument('-s', '--save', action='store_true', help='Saves extracted and processed data in a npz file')

    SpectraAdditionalCommandsGroup = SpectraSubparser.add_argument_group('Additional commands')
    SpectraAdditionalCommandsGroup.add_argument('-q', '--quiet', action='store_true', help='Include for the script to stay silent - This will not remove error messages or the printing of data')
    SpectraAdditionalCommandsGroup.add_argument('-mp','--multiprocessing', action='store_true', help='Include to use the multiprocessing library for data extraction')

    #---------------------------
    # Creating extract subparser
    #---------------------------
    ExtractionSubparser = subparser.add_parser('extract', formatter_class=argparse.RawDescriptionHelpFormatter, description=f'''
    This part of the script is for extracting data from output files and either printing it in the terminal or saving it to a file

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
                -  CPU time used
                -  Optimized geometries (or last geometry in file)

    Though not all data types have been implemented for all of the output formats

    All values are extracted as Atomic Units where applicable

    The following is not implemented for LSDALTON
    -  Zero-Point Vibrational energies
    -  Enthalpies
    -  Entropies
    -  Gibbs Free energies
    -  Frequencies
    -  Partition functions

    The following is not implemented for VELOXCHEM
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
    -  CPU time used

    The following is not implemented for AMSTERDAM MODELING SUITE
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
    -  CPU time used
'''
    , help='Use to extract data from output files')

    # Setting the Extract function to be run if extract is used
    ExtractionSubparser.set_defaults(func=Extract)

    # Adding arguments
    ExtractionSubparser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='File')

    ExtractionGroup = ExtractionSubparser.add_argument_group('Data extraction commands')
    ExtractionGroup.add_argument('-E', '--energy', action='store_true', help='Include to extract the Total Energy')
    ExtractionGroup.add_argument('-Z', '--zpv', action='store_true', help='Include to extract the Zero-Point Vibrational Energy')
    ExtractionGroup.add_argument('-H', '--enthalpy', action='store_true', help='Include to calculate molar enthalpies (kJ/mol)')
    ExtractionGroup.add_argument('-S', '--entropy', action='store_true', help='Include to calculate molar entropes (kJ/(mol*K)')
    ExtractionGroup.add_argument('-G', '--gibbs', action='store_true', help='Include to extract the Gibbs Free Energy')
    ExtractionGroup.add_argument('-D', '--dipole', action='store_true', help='Include to extract the Dipole Moment')
    ExtractionGroup.add_argument('-P', '--polar', action='store_true', help='Include to extract the Polarizability')
    ExtractionGroup.add_argument('-X', '--exc', const=-1, type=int, help='Include to extract the Excitation Energies. Add a number to extract that amount of Excitation Energies. It will extract all Excitation energies as default',nargs='?')
    ExtractionGroup.add_argument('-O', '--osc', action='store_true', help='Include to extract the Oscillator Strengths')
    ExtractionGroup.add_argument('-F', '--freq', const=-1, type=int, help='Include to extract the Frequencies. Add a number to extract that amount of Frequencies. It will extract all Frequencies as default', nargs='?')
    ExtractionGroup.add_argument('-Q', '--partfunc', action='store_true', help='Include to calculate molar partition functions.')
    ExtractionGroup.add_argument('-T', '--temp', const=298.15, default=298.15, type=float, help='Include to calculate at a different temperature. Default is 298.15 K', nargs='?')
    ExtractionGroup.add_argument('-C', '--cpu_time', const=['m'], help='Include to extract total cpu time and pr. cpu time. You can change the output from being in seconds, minutes and hours, where the default is minutes', nargs='?', choices=['s', 'm', 'h'])
    ExtractionGroup.add_argument('-geom', '--optgeom', action='store_true',help='Include to extract optimized geometries and save to \'filename_opt.xyz\'.')

    ExtractionDataProcessingGroup = ExtractionSubparser.add_argument_group('Data processing commands')
    ExtractionDataProcessingGroup.add_argument('-s', '--save', const='csv', type=str, help='Saves extracted and processed data. The extracted data is by default saved in a csv file', nargs='?', choices=['csv', 'npz', 'json', 'return'])
    ExtractionDataProcessingGroup.add_argument('--name', default='data', const='data', type=str, help='Define the name of the datafile where the extracted data is stored', nargs='?', dest='savename')

    ExtractionAdditionalCommandsGroup = ExtractionSubparser.add_argument_group('Additional commands')
    ExtractionAdditionalCommandsGroup.add_argument('-q', '--quiet', '--no-log', action='store_true', help="Include to not print error messages to the 'collect_data.log' file", dest='quiet')
    ExtractionAdditionalCommandsGroup.add_argument('-mp','--multiprocessing', action='store_true', help='Include to use the multiprocessing library for data extraction')
    ExtractionAdditionalCommandsGroup.add_argument('--no-progressbar', action='store_false', help='Include to deactivate progress bar', dest='progressbar')
    ExtractionAdditionalCommandsGroup.add_argument('--unittest', action='store_true', help=argparse.SUPPRESS)

    # Parses the arguments
    args = Parser.parse_args()

    # The arguments are sent to the correct function
    # The function may be one of Spectra, Extract, ...
    args.func(args)


if __name__ == "__main__":
    main()
