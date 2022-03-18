
import argparse
import numpy as np
from KurtGroup.Kurt import output_processing as op
from functools import partial
from multiprocessing import Pool, cpu_count


#*************************** INPUT PARSING ****************************

def Spectra(args):
    input_file = args.infile
    UVVIS = args.uvvis
    complex_propagator = args.complex_propagator
    format = args.format
    SAVE = args.save
    quiet = args.quiet
    MULTIPROCESSING = args.multiprocessing

    if UVVIS:
        NeededArguments = {'_Excitation_energies': -1, '_Oscillator_strengths': -1}
        Set_of_Values = {'_Excitation_energies': ['exc_energies'], '_Oscillator_strengths': ['osc_strengths']}

        if MULTIPROCESSING:
            with Pool(int(cpu_count()/2)) as pool:
                Extracted_Values = pool.map(partial(op.Data_Extraction, Needed_Values=NeededArguments, NeededArguments=NeededArguments, quiet=quiet), input_file)
                Extracted_Values = {key: value for dictionary in Extracted_Values for key, value in dictionary.items()} # Reformatting Extracted_values
        else:
            Extracted_Values = dict()
            for infile in input_file:
                Extracted_Values[infile] = op.Data_Extraction(infile, NeededArguments, NeededArguments, quiet)[infile]

        op.Check_if_Implemented(input_file, Set_of_Values, Extracted_Values)   #Finding functions not implemented

        for key_outer in Extracted_Values.keys():   # Turn everything into lists
            for key_inner in Extracted_Values[key_outer].keys():
                if type(Extracted_Values[key_outer][key_inner]) != list:
                    Extracted_Values[key_outer][key_inner] = [Extracted_Values[key_outer][key_inner]]

        op.Make_uvvis_spectrum(input_file, quiet, op.UVVIS_Spectrum, format, Extracted_Values, SAVE)

    elif complex_propagator:
        NeededArguments = {'_Complex_propagator': True}
        Set_of_Values = {'_Complex_propagator': ['complex_propagator']}

        if MULTIPROCESSING:
            with Pool(int(cpu_count()/2)) as pool:
                Extracted_Values = pool.map(partial(op.Data_Extraction, Needed_Values=NeededArguments, NeededArguments=NeededArguments, quiet=quiet), input_file)
                Extracted_Values = {key: value for dictionary in Extracted_Values for key, value in dictionary.items()} # Reformatting Extracted_values
        else:
            Extracted_Values = dict()
            for infile in input_file:
                Extracted_Values[infile] = op.Data_Extraction(infile, NeededArguments, NeededArguments, quiet)[infile]

        op.Check_if_Implemented(input_file, Set_of_Values, Extracted_Values)   #Finding functions not implemented

        for key_outer in Extracted_Values.keys():   # Turn everything into lists
            for key_inner in Extracted_Values[key_outer].keys():
                if type(Extracted_Values[key_outer][key_inner]) != list:
                    Extracted_Values[key_outer][key_inner] = [Extracted_Values[key_outer][key_inner]]

        op.Make_complex_propagator_spectrum(input_file, quiet, format, Extracted_Values, SAVE)

    else:
        print('Doing nothing - please use --uvvis or --complex-propagator')



def Extract(args):
    input_file = args.infile
    SAVE = args.save
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
        '_CPUS' : args.cpu_time,
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
        '_CPUS' : ['total_cpu_time', 'wall_cpu_time'],
    }

    Header_text = { #These are what will be written in the header for each datapoint
        'tot_energy' : 'Energy',
        'zpv' : 'ZPV EnergTotaly',
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
        'total_cpu_time': f'Total CPU time ({args.cpu_time})',
        'wall_cpu_time' : f'Wall CPU time ({args.cpu_time})'
    }

    quiet = args.quiet
    MULTIPROCESSING = args.multiprocessing

    NeededArguments = RequestedArguments.copy()

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

    if MULTIPROCESSING:
        with Pool(int(cpu_count()/2)) as pool:
            Extracted_Values = pool.map(partial(op.Data_Extraction, Needed_Values=Needed_Values, NeededArguments=NeededArguments, quiet=quiet, Temperature=T), input_file)
            Extracted_Values = {key: value for dictionary in Extracted_Values for key, value in dictionary.items()} # Reformatting Extracted_values
    else:
        Extracted_Values = dict()
        for infile in input_file:
            Extracted_Values[infile] = op.Data_Extraction(infile, Needed_Values, NeededArguments, quiet, T)[infile]

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

    op.Downsizing_variable_arrays(Outputs, Variable_arrays, count, Final_arrays)   #Fixing the size of variable size arrays

    Header = op.Create_Header(Header_text, Set_of_Values, Final_arrays)    #Creation of header row

    Output_Array = np.array([Header] * (count + 1), dtype=object)   #Create output array from header row

    op.Fill_output_array(Set_of_Values, Input_Array, count, Final_arrays, Output_Array)    #Filling the output array with the extracted data

#   ------------ IF CHOSEN PRINTS THE OUTPUT IN A CSV FILE ------------
#   ---------- ELSE THE RESULTS ARE DUMPED INTO THE TERMINAL ----------

    if SAVE == 'csv':
        np.savetxt('data.csv', Output_Array, delimiter=',', fmt='%s')
        print(f'Data has been saved in data.csv')
    elif SAVE == 'npz':
        Save_Dict = {i[0]: i[1:] for i in Output_Array}
        np.savez('data.npz', **Save_Dict)
        print(f'Data has been saved in data.npz')
    else:
        print(Output_Array)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=f'''
    A script designed to make it easier to extract data from output files

    To extract data you have to use the keyword extract before any other keywords
    To make spectra you will have to use the keyword spectra before any other keywords

            Currently the output formats supported are
            ------------------------------------------
                -  ORCA
                -  DALTON
                -  GAUSSIAN
                -  LSDALTON
                -  VELOXCHEM'''
    , epilog=f'''
For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk

    Magnus Bukhave Johansen
    qhw298@alumni.ku.dk''')
    subparser = parser.add_subparsers(dest='pars')
    subparser.required = False

    Spectra_subparser = subparser.add_parser('spectra', formatter_class=argparse.RawDescriptionHelpFormatter, description=f'''
    This part of the script is used for making spectra from data collected from output files

            Currently the output formats supported are
            ------------------------------------------
                -  ORCA
                -  DALTON
                -  GAUSSIAN
                -  LSDALTON
                -  VELOXCHEM

    It is currently possible to make UVVIS spectra using excitation energies and complex propagator theory

    The following is not implemented for ORCA
    -  UVVIS based on complex propagator theory

    The following is not implemented for DALTON
    -  UVVIS based on complex propagator theory

    The following is not implemented for GAUSSIAN
    -  UVVIS based on complex propagator theory

    The following is not implemented for VELOXCHEM
    -  UVVIS based on excitation energies
    -  UVVIS based on complex propagator theory'''
    , help='Use to make spectra such as UVVIS from excitation energies or complex propagator theory')
    Spectra_subparser.set_defaults(func=Spectra)
    Spectra_subparser.add_argument('infile', type=str, nargs='+', help='The file(s) to make spectra from', metavar='File')

    Spectra_group = Spectra_subparser.add_mutually_exclusive_group()
    Spectra_group.add_argument('--uvvis', action='store_true', help='Include to make UVVIS spectra')
    Spectra_group.add_argument('--complex-propagator', action='store_true', help='Include to use complex propagator theory to make the spectra')

    SpectraDataProcessingGroup = Spectra_subparser.add_argument_group('Data processing commands')
    SpectraDataProcessingGroup.add_argument('--format', default='png', const='png', type=str, help='Include this to change the picture format. Will use png as default. If \'--save\' is used together with this, the processed data will be saved in a .npz file', nargs='?', choices=['png', 'eps', 'pdf', 'svg', 'ps'])
    SpectraDataProcessingGroup.add_argument('-s', '--save', action='store_true', help='Saves extracted and processed data in a npz file')

    SpectraAdditionalCommandsGroup = Spectra_subparser.add_argument_group('Additional commands')
    SpectraAdditionalCommandsGroup.add_argument('-q', '--quiet', action='store_true', help='Include for the script to stay silent - This will not remove error messages or the printing of data')
    SpectraAdditionalCommandsGroup.add_argument('-mp','--multiprocessing', action='store_true', help='Include to use the multiprocessing library for data extraction')

    Extraction_subparser = subparser.add_parser('extract', formatter_class=argparse.RawDescriptionHelpFormatter, description=f'''
    This part of the script is for extracting data from output files and either printing it in the terminal or saveing it to a file

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
'''
    , help='Use to extract data from output files')
    Extraction_subparser.set_defaults(func=Extract)
    Extraction_subparser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='File')

    ExtractionGroup = Extraction_subparser.add_argument_group('Data extraction commands')
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

    ExtractionDataProcessingGroup = Extraction_subparser.add_argument_group('Data processing commands')
    ExtractionDataProcessingGroup.add_argument('-s', '--save', const='csv', type=str, help='Saves extracted and processed data. The extracted data is by default saved in a csv file', nargs='?', choices=['csv','npz'])

    ExtractionAdditionalCommandsGroup = Extraction_subparser.add_argument_group('Additional commands')
    ExtractionAdditionalCommandsGroup.add_argument('-q', '--quiet', action='store_true', help='Include for the script to stay silent - This will not remove error messages or the printing of data')
    ExtractionAdditionalCommandsGroup.add_argument('-mp','--multiprocessing', action='store_true', help='Include to use the multiprocessing library for data extraction')

    args = parser.parse_args()
    args.func(args)


#EOF
