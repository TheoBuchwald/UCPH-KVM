
import unittest
import json
import os
import sys

sys.path.append('')

import KurtGroup.Kurt.output_processing as op
from collect_data import Check_if_Implemented

def ReadJSONFile(filename: str):
    BASE_DIR = os.path.dirname(__file__)
    with open(f'{BASE_DIR}/{filename}', 'r') as json_file:
        dictionary = json.load(json_file)
    return dictionary

def Extraction(data_file, func: str):
    Extracted_Values = dict()
    for infile in data_file:
        Extracted_Values[infile] = {}
        file = f'test_systems/{infile}'
        outfile = op.OutputType(file, Quiet=True)

        Values = getattr(outfile, func)()

        if Values:
            Extracted_Values[infile] = {'test': Values}

    Check_if_Implemented(data_file, {'_': ['test']}, Extracted_Values)

    return Extracted_Values


class Test_output_processing(unittest.TestCase):

    def test_Energy_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getEnergy")

        for infile in DATA_FILE:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['tot_energy'])

    def test_ZPV_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getZeroPointVibrationalEnergy")

        for infile in DATA_FILE:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['zpv'])


    def test_Dipole_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getDipoleMoment")

        for infile in DATA_FILE:
            self.assertAlmostEqual(Extracted_Values[infile]['test'][0], DATA_FILE[infile]['dipolex'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][1], DATA_FILE[infile]['dipoley'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][2], DATA_FILE[infile]['dipolez'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][3], DATA_FILE[infile]['total_dipole'])

    def test_Polarizability_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getPolarizability")

        for infile in DATA_FILE:
            self.assertAlmostEqual(Extracted_Values[infile]['test'][0], DATA_FILE[infile]['polx'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][1], DATA_FILE[infile]['poly'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][2], DATA_FILE[infile]['polz'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][3], DATA_FILE[infile]['iso_polar'])

    def test_Excitation_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getExcitationEnergies")

        for infile in DATA_FILE:
            self.assertListEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['exc_energies'])

    def test_Oscillator_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getOscillatorStrengths")

        for infile in DATA_FILE:
            self.assertListEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['osc_strengths'])

    def test_Frequency_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getFrequencies")

        for infile in DATA_FILE:
            self.assertListEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['freq'])

    def test_Enthalpy(self):

        Extracted_Values = Extraction(DATA_FILE, "getEnthalpy")

        for infile in DATA_FILE:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['enthalpy'])

    def test_Entropy_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getEntropy")

        for infile in DATA_FILE:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['entropy'])

    def test_Gibbs_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getGibbsFreeEnergy")

        for infile in DATA_FILE:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['gibbs'])

    def test_Partitionfunction_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getPartitionFunction")

        for infile in DATA_FILE:
            self.assertAlmostEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['qTotal'])

    def test_CPUtime_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getCPUTime")

        for infile in DATA_FILE:
            self.assertAlmostEqual(Extracted_Values[infile]['test'][0], DATA_FILE[infile]['total_cpu_time'])
            self.assertAlmostEqual(Extracted_Values[infile]['test'][1], DATA_FILE[infile]['wall_cpu_time'])

TEST_DATA = "test_data.json"
DATA_FILE = ReadJSONFile(TEST_DATA)

if __name__ == '__main__':
    unittest.main()
