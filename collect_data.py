
import argparse
import numpy as np
from KurtGroup.Kurt import output_processing as op
from functools import partial
from multiprocessing import Pool, cpu_count


#*************************** INPUT PARSING ****************************

def Spectra(args):
    """
    This function is used to do any methods related to the Spectra keyword
    """
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
    """
    This function is used for any methods related to the extract keyword
    """

    # Setting the input_files variable to all files
    input_files = args.infile

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
    Header_text = {
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
    SAVE = args.save
    quiet = args.quiet
    MULTIPROCESSING = args.multiprocessing

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
    Variable_arrays = dict([item for item in NeededArguments.items() if type(item[1]) == int])

    # List of arguments that have been requested
    Wanted_Values = [item[0] for item in RequestedArguments.items() if not(item[1] == None or item[1] == False)]

    # List of arguments that are needed for those requested
    # This is done so you can request as an example the Gibbs free energies withou also having the enthalpies and entropies printed
    Needed_Values = [item[0] for item in NeededArguments.items() if not(item[1] == None or item[1] == False)]

    # Data-points that will be written to the terminal or save file
    # These are found from the Outputs dictionary by comparing with the Wanted_Values list
    Set_of_Values = dict(item for item in Outputs.items() if item[0] in Wanted_Values)

    # How many files to run the script on
    count = len(input_files)

    # If multiprocessing is eneabled it will be run using half of the available CPUS
    # Else they will be run in a linear fashion
    if MULTIPROCESSING:
        with Pool(int(cpu_count()/2)) as pool:
            Extracted_Values = pool.map(partial(op.Data_Extraction, Needed_Values=Needed_Values, NeededArguments=NeededArguments, quiet=quiet, Temperature=T), input_files)
            Extracted_Values = {key: value for dictionary in Extracted_Values for key, value in dictionary.items()} # Reformatting Extracted_values
    else:
        Extracted_Values = dict()
        for infile in input_files:
            Extracted_Values[infile] = op.Data_Extraction(infile, Needed_Values, NeededArguments, quiet, T)[infile]

    # Creating Input_Array where all values are put in lists
    Input_Array = [[i] for i in Extracted_Values]

    # Checking if some functions have not been implemented for the relevant extraction types
    op.Check_if_Implemented(input_files, Set_of_Values, Extracted_Values)

    # If something is at this point not in a list somehow they will be after this
    for key_outer in Extracted_Values.keys():
        for key_inner in Extracted_Values[key_outer].keys():
            if type(Extracted_Values[key_outer][key_inner]) != list:
                Extracted_Values[key_outer][key_inner] = [Extracted_Values[key_outer][key_inner]]

    # Collecting all values in arrays in a dictionary
    # Some values in the Extracted_Values dictionary may not have been requested, so these are removed here
    Final_arrays = op.Collect_and_sort_data(input_files, Set_of_Values, Extracted_Values)

    # Resizing arrays
    # An example is Excitation energies where there may be more of them in one output file than another
    # By doing this it fits properly in what is printed to the terminal
    for key in Final_arrays.keys():
        if type(Final_arrays[key]) == list:
            op.Resize(Final_arrays[key])

    # Fixing the size of variable size arrays so that they match what was requested
    op.Downsizing_variable_arrays(Outputs, Variable_arrays, count, Final_arrays)

    # Creation of header row
    Header = op.Create_Header(Header_text, Set_of_Values, Final_arrays)

    # Create output array from header row
    Output_Array = np.array([Header] * (count + 1), dtype=object)

    # Filling the output array with the extracted data
    op.Fill_output_array(Set_of_Values, Input_Array, count, Final_arrays, Output_Array)

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


def main():
    #---------------------------
    # Creating main parser
    #---------------------------
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

    #---------------------------
    # Creating spectra subparser
    #---------------------------
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

    # Setting the Spectra function to be run if spectra is used
    Spectra_subparser.set_defaults(func=Spectra)

    # Adding arguments
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

    #---------------------------
    # Creating extract subparser
    #---------------------------
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

    # Setting the Extract function to be run if extract is used
    Extraction_subparser.set_defaults(func=Extract)

    # Adding arguments
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

    # Parses the arguments
    args = parser.parse_args()

    # The arguments are sent to the correct function
    # The function may be one of Spectra, Extract, ...
    args.func(args)


if __name__ == "__main__":
    main()

#EOF
