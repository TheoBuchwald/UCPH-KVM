
import unittest
import json
import os
import sys

sys.path.append('')

import KurtGroup.Kurt.output_processing as op
from collect_data import Check_if_Implemented, Data_Extraction

def ReadJSONFile():
    BASE_DIR = os.path.dirname(__file__)
    with open(f'{BASE_DIR}/test_data.json', 'r') as json_file:
        dictionary = json.load(json_file)
    return dictionary

def Extraction(data_file, func: str):
    Extracted_Values = dict()
    for infile in data_file:
        Extracted_Values[infile] = {}
        file = f'test_systems/{infile}'
        outfile = op.OutputType(file, Quiet=True)

        Values = getattr(outfile, func)()
        # Enthalpy = outfile.getEnthalpy()
        if Values:
            Extracted_Values[infile] = {'test': Values}

    Check_if_Implemented(data_file, {'_': ['test']}, Extracted_Values)

    return Extracted_Values

# def Check_if_Implemented(input_file: str, Set_of_values: dict, Extracted_values: dict) -> None:
#     # Checks to see if the keys of a double dictionary exists
#     # If they don't it is assumed that the function related to the data hasn't been implemented
#     for infile in input_file:
#         for key in Set_of_values:
#             for val in Set_of_values[key]:
#                 try:
#                     Extracted_values[infile][val]
#                 except KeyError:
#                     Extracted_values[infile][val] = ['Not implemented']


class Test_output_processing(unittest.TestCase):

    def test_Energy_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getEnergy")

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], data_file[infile]['tot_energy'])

    def test_ZPV_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getZeroPointVibrationalEnergy")

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], data_file[infile]['zpv'])


    def test_Dipole_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getDipoleMoment")

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['test'][0], data_file[infile]['dipolex'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][1], data_file[infile]['dipoley'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][2], data_file[infile]['dipolez'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][3], data_file[infile]['total_dipole'])

    def test_Polarizability_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getPolarizability")

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['test'][0], data_file[infile]['polx'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][1], data_file[infile]['poly'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][2], data_file[infile]['polz'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][3], data_file[infile]['iso_polar'])

    def test_Excitation_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getExcitationEnergies")

        for infile in data_file:
            self.assertListEqual(Extracted_Values[infile]['test'], data_file[infile]['exc_energies'])

    def test_Oscillator_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getOscillatorStrengths")

        for infile in data_file:
            self.assertListEqual(Extracted_Values[infile]['test'], data_file[infile]['osc_strengths'])

    def test_Frequency_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getFrequencies")

        for infile in data_file:
            self.assertListEqual(Extracted_Values[infile]['test'], data_file[infile]['freq'])

    def test_Enthalpy(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getEnthalpy")

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], data_file[infile]['enthalpy'])

    def test_Entropy_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getEntropy")

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], data_file[infile]['entropy'])

    def test_Gibbs_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getGibbsFreeEnergy")

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], data_file[infile]['gibbs'])

    def test_Partitionfunction_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getPartitionFunction")

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], data_file[infile]['qTotal'])

    def test_CPUtime_Extraction(self):

        data_file = ReadJSONFile()

        Extracted_Values = Extraction(data_file, "getCPUTime")

        for infile in data_file:
            self.assertAlmostEqual(Extracted_Values[infile]['test'][0], data_file[infile]['total_cpu_time'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][1], data_file[infile]['wall_cpu_time'])


if __name__ == '__main__':
    unittest.main()
