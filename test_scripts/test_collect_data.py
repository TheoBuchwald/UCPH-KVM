
from ...KurtGroup.Kurt import output_processing as op
import unittest
import json
import os

def ReadJSONFile():
    BASE_DIR = os.path.dirname(__file__)
    with open(f'{BASE_DIR}/test_data.json', 'r') as json_file:
        dictionary = json.load(json_file)
    return dictionary


class Test_output_processing(unittest.TestCase):

    def test_Energy_Extraction(self):

        Arguments = {'_Energy': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_Energy': ['tot_energy']}, Extracted_Values)

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['tot_energy'], data_file[infile]['tot_energy'])

    def test_ZPV_Extraction(self):

        Arguments = {'_ZPV': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_ZPV': ['zpv']}, Extracted_Values)

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['zpv'], data_file[infile]['zpv'])


    def test_Dipole_Extraction(self):

        Arguments = {'_Dipole_moments': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_Dipole_moments': ['total_dipole']}, Extracted_Values)

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['dipolex'], data_file[infile]['dipolex'])
            self.assertAlmostEqual(Extracted_Values[infile]['dipoley'], data_file[infile]['dipoley'])
            self.assertAlmostEqual(Extracted_Values[infile]['dipolez'], data_file[infile]['dipolez'])
            self.assertAlmostEqual(Extracted_Values[infile]['total_dipole'], data_file[infile]['total_dipole'])

    def test_Polarizability_Extraction(self):

        Arguments = {'_Polarizabilities': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_Polarizabilities': ['iso_polar']}, Extracted_Values)

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['polx'], data_file[infile]['polx'])
            self.assertAlmostEqual(Extracted_Values[infile]['poly'], data_file[infile]['poly'])
            self.assertAlmostEqual(Extracted_Values[infile]['polz'], data_file[infile]['polz'])
            self.assertAlmostEqual(Extracted_Values[infile]['iso_polar'], data_file[infile]['iso_polar'])

    def test_Excitation_Extraction(self):

        Arguments = {'_Excitation_energies': -1}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_Excitation_energies': ['exc_energies']}, Extracted_Values)

        for infile in data_file:
            self.assertListEqual(Extracted_Values[infile]['exc_energies'], data_file[infile]['exc_energies'])

    def test_Oscillator_Extraction(self):

        Arguments = {'_Excitation_energies': -1, '_Oscillator_strengths': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_Oscillator_strengths': ['osc_strengths']}, Extracted_Values)

        for infile in data_file:
            self.assertListEqual(Extracted_Values[infile]['osc_strengths'], data_file[infile]['osc_strengths'])

    def test_Frequency_Extraction(self):

        Arguments = {'_Frequencies': -1}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_Frequencies': ['freq']}, Extracted_Values)

        for infile in data_file:
            self.assertListEqual(Extracted_Values[infile]['freq'], data_file[infile]['freq'])

    def test_Enthalpy(self):

        Arguments = {'_Energy': True, '_Frequencies': -1, '_Enthalpy': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_Enthalpy': ['enthalpy']}, Extracted_Values)

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['enthalpy'], data_file[infile]['enthalpy'])

    def test_Entropy_Extraction(self):

        Arguments = {'_Frequencies': -1, '_Entropy': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_Entropy': ['entropy']}, Extracted_Values)

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['entropy'], data_file[infile]['entropy'])

    def test_Gibbs_Extraction(self):

        Arguments = {'_Energy': True, '_Frequencies': -1, '_Enthalpy': True, '_Entropy': True, '_Gibbs': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_Gibbs': ['gibbs']}, Extracted_Values)

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['gibbs'], data_file[infile]['gibbs'])

    def test_Partitionfunction_Extraction(self):

        Arguments = {'_Frequencies': -1, '_PartitionFunctions': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_PartitionFunctions': ['qTotal']}, Extracted_Values)

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['qTotal'], data_file[infile]['qTotal'])

    def test_CPUtime_Extraction(self):

        Arguments = {'_CPUS': True}
        Values = [item[0] for item in Arguments.items() if not(item[1] == None or item[1] == False)]

        data_file = ReadJSONFile()

        Extracted_Values = dict()
        for infile in data_file:
            file = f'test_systems/{infile}'
            Extracted_Values[infile] = op.Data_Extraction(file, Needed_Values=Values, NeededArguments=Arguments, quiet=True)[file]

        op.Check_if_Implemented(data_file, {'_CPUS' : 'wall_cpu_time'}, Extracted_Values)

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['total_cpu_time'], data_file[infile]['total_cpu_time'])
            self.assertAlmostEqual(Extracted_Values[infile]['wall_cpu_time'], data_file[infile]['wall_cpu_time'])


if __name__ == '__main__':
    unittest.main()
