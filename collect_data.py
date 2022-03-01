
import argparse
import numpy as np
import dependencies.output_processing as op
from functools import partial
from multiprocessing import Pool, cpu_count
from colorama import Fore, Style, init

init(autoreset=True)


#*************************** INPUT PARSING ****************************


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=f'''
    A script designed to make it easier to extract data from output files

            Currently the output formats supported are
            ------------------------------------------
                -  {Style.BRIGHT}{Fore.CYAN}ORCA{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.CYAN}DALTON{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.CYAN}GAUSSIAN{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.CYAN}LSDALTON{Style.RESET_ALL}

            This script can currently extract the data
            ------------------------------------------
                -  {Style.BRIGHT}{Fore.GREEN}Total energies{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Zero-Point Vibrational energies{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Enthalpies{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Entropies{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Gibbs Free energies{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Dipole moments{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Polarizability{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Excitation energies{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Oscillator strengths{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Frequencies{Style.RESET_ALL}
                -  {Style.BRIGHT}{Fore.GREEN}Partition functions at a given temperature{Style.RESET_ALL}

    Though not all data types have been implemented for all of the output formats

    All values are extracted as {Style.BRIGHT}{Fore.MAGENTA}Atomic Units{Style.RESET_ALL} where applicable'''
    ,epilog=f'''
    The following is not implemented for {Style.BRIGHT}{Fore.CYAN}ORCA{Style.RESET_ALL}
    -  {Fore.RED}Entropies{Style.RESET_ALL}

    The following is not implemented for {Style.BRIGHT}{Fore.CYAN}DALTON{Style.RESET_ALL}
    -  {Fore.RED}Non DFT total energies{Style.RESET_ALL}
    -  {Fore.RED}Zero-Point Vibrational energies{Style.RESET_ALL}
    -  {Fore.RED}Enthalpies{Style.RESET_ALL}
    -  {Fore.RED}Entropies{Style.RESET_ALL}
    -  {Fore.RED}Gibbs Free energies{Style.RESET_ALL}
    -  {Fore.RED}Frequencies{Style.RESET_ALL}
    -  {Fore.RED}Partition functions with symmetry{Style.RESET_ALL}

    The following is not implemented for {Style.BRIGHT}{Fore.CYAN}GAUSSIAN{Style.RESET_ALL}
    -  {Fore.RED}Entropies{Style.RESET_ALL}
    -  {Fore.RED}Excitation energies{Style.RESET_ALL}
    -  {Fore.RED}Oscillator strengths{Style.RESET_ALL}

    The following is not implemented for {Style.BRIGHT}{Fore.CYAN}LSDALTON{Style.RESET_ALL}
    -  {Fore.RED}Probably some total energies{Style.RESET_ALL}
    -  {Fore.RED}Zero-Point Vibrational energies{Style.RESET_ALL}
    -  {Fore.RED}Enthalpies{Style.RESET_ALL}
    -  {Fore.RED}Gibbs Free energies{Style.RESET_ALL}
    -  {Fore.RED}Entropies{Style.RESET_ALL}
    -  {Fore.RED}Dipole moments{Style.RESET_ALL}
    -  {Fore.RED}Polarizabilities{Style.RESET_ALL}
    -  {Fore.RED}Excitation energies{Style.RESET_ALL}
    -  {Fore.RED}Oscillator strengths{Style.RESET_ALL}
    -  {Fore.RED}Frequencies{Style.RESET_ALL}
    -  {Fore.RED}Partition functions{Style.RESET_ALL}

    For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk

    Magnus Bukhave Johansen
    qhw298@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='File')

    ExtractionGroup = parser.add_argument_group('Data extraction commands')
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

    DataProcessingGroup = parser.add_argument_group('Data processing commands')
    DataProcessingGroup.add_argument('-s', '--save', const='csv', type=str, help='Saves extracted and processed data. The extracted data is by default saved in a csv file', nargs='?', choices=['csv','npz'])
    DataProcessingGroup.add_argument('--uvvis', const='png', type=str, help='Include to calculate a UV-VIS spectrum. Tou can give the picture format as an optional argument. Will use png as default. If \'--save\' is used togrther with this, the processed data will be saved in a .npz file', nargs='?', choices=['png', 'eps', 'pdf', 'svg', 'ps'])

    AdditionalCommandsGroup = parser.add_argument_group('Additional commands')
    AdditionalCommandsGroup.add_argument('-q', '--quiet', action='store_true', help='Include for the script to stay silent - This will not remove error messages or the printing of data')


#******************************* SETUP ********************************


    args = parser.parse_args()
    input_file = args.infile
    SAVE = args.save
    UVVIS = args.uvvis
    T = args.temp

    RequestedArguments = {   #These are all the possible arguments that extract data
        '_Energy' : args.energy,
        '_ZPV' : args.zpv,
        '_Dipole_moments' : args.dipole,
        '_Polarizabilities' : args.polar,
        '_Excitation_energies' : args.exc,
        '_Oscillator_strengths' : args.osc,
        '_Frequencies' : args.freq,
        '_Enthalpy' : args.enthalpy,
        '_Entropy' : args.entropy,
        '_Gibbs' : args.gibbs,
        '_PartitionFunctions' : args.partfunc,
    }

    Outputs = {     #These are the datapoints that will be extracted per argument
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
        '_PartitionFunctions' : ['qTotal'],
    }

    Header_text = { #These are what will be written in the header for each datapoint
        'tot_energy' : 'Energy',
        'zpv' : 'ZPV Energy',
        'enthalpy' : 'Enthalpy',
        'entropy' : 'Entropy (kJ/(mol*K)',
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
        'qTotal' : 'Total molar partition function',
    }

    quiet = args.quiet

    NeededArguments = RequestedArguments.copy()

    if UVVIS:
        NeededArguments['_Oscillator_strengths'] = NeededArguments['_Excitation_energies'] = -1

    if NeededArguments['_Oscillator_strengths'] and NeededArguments['_Excitation_energies'] == None:   #Ensuring that excitation energies are calculated when oscillator strengths are
        NeededArguments['_Oscillator_strengths'] = NeededArguments['_Excitation_energies'] = -1
    elif NeededArguments['_Oscillator_strengths']:
        NeededArguments['_Oscillator_strengths'] = NeededArguments['_Excitation_energies']

    if NeededArguments['_Gibbs']:
        NeededArguments['_Enthalpy'] = True
        NeededArguments['_Entropy'] = True

    if NeededArguments['_PartitionFunctions']:  #Ensuring that frequencies are calculated as these are needed to calculate the Partition Functions
        NeededArguments['_Frequencies'] = -1

    if NeededArguments['_Enthalpy']:  #Ensuring that frequencies are calculated as these are needed to calculate the Partition Functions
        NeededArguments['_Frequencies'] = -1
        NeededArguments['_Energy'] = True

    if NeededArguments['_Entropy']:  #Ensuring that frequencies are calculated as these are needed to calculate the Partition Functions
        NeededArguments['_Frequencies'] = -1

    Variable_arrays = dict([item for item in NeededArguments.items() if type(item[1]) == int])    #Dictionary of data where the amount of values printed can be changed
    Wanted_Values = [item[0] for item in RequestedArguments.items() if not(item[1] == None or item[1] == False)]
    Needed_Values = [item[0] for item in NeededArguments.items() if not(item[1] == None or item[1] == False)]
    Set_of_Values = dict(item for item in Outputs.items() if item[0] in Wanted_Values)

    count = len(input_file)

    with Pool(int(cpu_count()/2)) as pool:
        Extracted_Values = pool.map(partial(op.Data_Extraction, Needed_Values=Needed_Values, NeededArguments=NeededArguments, quiet=quiet, Temperature=T), input_file)

    Extracted_Values = {key: value for dictionary in Extracted_Values for key, value in dictionary.items()} # Reformatting Extracted_values
    Input_Array = [[i] for i in Extracted_Values] # Creating array_input

    op.Check_if_Implemented(input_file, Set_of_Values, Extracted_Values)   #Finding functions not implemented

    for key_outer in Extracted_Values.keys():   # Turn everything into lists
        for key_inner in Extracted_Values[key_outer].keys():
            if type(Extracted_Values[key_outer][key_inner]) != list:
                Extracted_Values[key_outer][key_inner] = [Extracted_Values[key_outer][key_inner]]

    Final_arrays = op.Collect_and_sort_data(input_file, Set_of_Values, Extracted_Values)   #Collecting all values in arrays in a dictionary

    for key in Final_arrays.keys(): #Resizing arrays
        if type(Final_arrays[key]) == list:
            op.Resize(Final_arrays[key])

    if UVVIS:
        op.Make_uvvis_spectrum(input_file, quiet, op.UVVIS_Spectrum, count, Final_arrays, UVVIS, Extracted_Values, SAVE)

    op.Downsizing_variable_arrays(Outputs, Variable_arrays, count, Final_arrays)   #Fixing the size of variable size arrays

    Header = op.Create_Header(Header_text, Set_of_Values, Final_arrays)    #Creation of header row

    Output_Array = np.array([Header] * (count + 1), dtype=object)   #Create output array from header row

    op.Fill_output_array(Set_of_Values, Input_Array, count, Final_arrays, Output_Array)    #Filling the output array with the extracted data

#   ------------ IF CHOSEN PRINTS THE OUTPUT IN A CSV FILE ------------
#   ---------- ELSE THE RESULTS ARE DUMPED INTO THE TERMINAL ----------

    if SAVE == 'csv':
        np.savetxt('data.csv', Output_Array, delimiter=',', fmt='%s')
        print(f'{Style.BRIGHT}{Fore.LIGHTGREEN_EX}Data has been saved in data.csv')
    elif SAVE == 'npz':
        Save_Dict = {i[0]: i[1:] for i in Output_Array}
        np.savez('data.npz', **Save_Dict)
        print(f'{Style.BRIGHT}{Fore.LIGHTGREEN_EX}Data has been saved in data.npz')
    else:
        print(Output_Array)


#EOF
