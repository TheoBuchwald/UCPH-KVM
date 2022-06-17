
from argparse import Namespace
import unittest
import json
import os
import sys

sys.path.append('')

import KurtGroup.Kurt.output_processing as op
import collect_data as cd

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

    cd.Check_if_Implemented(data_file, {'_': ['test']}, Extracted_Values)

    return Extracted_Values


class Test_output_processing(unittest.TestCase):

    def test_Energy_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getEnergy")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['tot_energy'])

    def test_ZPV_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getZeroPointVibrationalEnergy")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['zpv'])


    def test_Dipole_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getDipoleMoment")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'][0], DATA_FILE[infile]['dipolex'])
            self.assertEqual(Extracted_Values[infile]['test'][1], DATA_FILE[infile]['dipoley'])
            self.assertEqual(Extracted_Values[infile]['test'][2], DATA_FILE[infile]['dipolez'])
            self.assertEqual(Extracted_Values[infile]['test'][3], DATA_FILE[infile]['total_dipole'])

    def test_Polarizability_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getPolarizability")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'][0], DATA_FILE[infile]['polx'])
            self.assertEqual(Extracted_Values[infile]['test'][1], DATA_FILE[infile]['poly'])
            self.assertEqual(Extracted_Values[infile]['test'][2], DATA_FILE[infile]['polz'])
            self.assertEqual(Extracted_Values[infile]['test'][3], DATA_FILE[infile]['iso_polar'])

    def test_Excitation_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getExcitationEnergies")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['exc_energies'])

    def test_Oscillator_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getOscillatorStrengths")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['osc_strengths'])

    def test_Frequency_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getFrequencies")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['freq'])

    def test_Enthalpy(self):

        Extracted_Values = Extraction(DATA_FILE, "getEnthalpy")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['enthalpy'])

    def test_Entropy_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getEntropy")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['entropy'])

    def test_Gibbs_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getGibbsFreeEnergy")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['gibbs'])

    def test_Partitionfunction_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getPartitionFunction")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'], DATA_FILE[infile]['qTotal'])

    def test_CPUtime_Extraction(self):

        Extracted_Values = Extraction(DATA_FILE, "getCPUTime")

        for infile in DATA_FILE:
            self.assertEqual(Extracted_Values[infile]['test'][0], DATA_FILE[infile]['total_cpu_time'])
            self.assertEqual(Extracted_Values[infile]['test'][1], DATA_FILE[infile]['wall_cpu_time'])


class Test_collect_data(unittest.TestCase):

    def test_Extract(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=['m'],
            dipole=True,
            energy=True,
            enthalpy=True,
            entropy=True,
            exc=-1,
            freq=-1,
            gibbs=True,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=True,
            osc=True,
            partfunc=True,
            polar=True,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=True,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        self.assertTrue(Values)

    def test_Extract_multiprocessing(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=['m'],
            dipole=True,
            energy=True,
            enthalpy=True,
            entropy=True,
            exc=-1,
            freq=-1,
            gibbs=True,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=True,
            optgeom=True,
            osc=True,
            partfunc=True,
            polar=True,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=True,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        self.assertTrue(Values)

    def test_Extract_energy(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=True,
            enthalpy=False,
            entropy=False,
            exc=None,
            freq=None,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['tot_energy'], DATA_FILE[infile]['tot_energy'])

    def test_Extract_ZPV(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=False,
            enthalpy=False,
            entropy=False,
            exc=None,
            freq=None,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=True,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['zpv'], DATA_FILE[infile]['zpv'])

    def test_Extract_dipole(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=True,
            energy=False,
            enthalpy=False,
            entropy=False,
            exc=None,
            freq=None,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['dipolex'], DATA_FILE[infile]['dipolex'])
            self.assertEqual(Values[f'test_systems/{infile}']['dipoley'], DATA_FILE[infile]['dipoley'])
            self.assertEqual(Values[f'test_systems/{infile}']['dipolez'], DATA_FILE[infile]['dipolez'])
            self.assertEqual(Values[f'test_systems/{infile}']['total_dipole'], DATA_FILE[infile]['total_dipole'])

    def test_Extract_polarizability(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=False,
            enthalpy=False,
            entropy=False,
            exc=None,
            freq=None,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=True,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['polx'], DATA_FILE[infile]['polx'])
            self.assertEqual(Values[f'test_systems/{infile}']['poly'], DATA_FILE[infile]['poly'])
            self.assertEqual(Values[f'test_systems/{infile}']['polz'], DATA_FILE[infile]['polz'])
            self.assertEqual(Values[f'test_systems/{infile}']['iso_polar'], DATA_FILE[infile]['iso_polar'])

    def test_Extract_excitation(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=False,
            enthalpy=False,
            entropy=False,
            exc=-1,
            freq=None,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['exc_energies'], DATA_FILE[infile]['exc_energies'])

    def test_Extract_oscillator(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=False,
            enthalpy=False,
            entropy=False,
            exc=-1,
            freq=None,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=True,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['osc_strengths'], DATA_FILE[infile]['osc_strengths'])

    def test_Extract_frequencies(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=False,
            enthalpy=False,
            entropy=False,
            exc=None,
            freq=-1,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['freq'], DATA_FILE[infile]['freq'])

    def test_Extract_enthalpy(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=True,
            enthalpy=True,
            entropy=False,
            exc=None,
            freq=-1,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['enthalpy'], DATA_FILE[infile]['enthalpy'])

    def test_Extract_entropy(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=False,
            enthalpy=False,
            entropy=True,
            exc=None,
            freq=-1,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['entropy'], DATA_FILE[infile]['entropy'])

    def test_Extract_gibbs(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=True,
            enthalpy=True,
            entropy=True,
            exc=None,
            freq=-1,
            gibbs=True,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['gibbs'], DATA_FILE[infile]['gibbs'])

    def test_Extract_enthalpy(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=None,
            dipole=False,
            energy=False,
            enthalpy=False,
            entropy=False,
            exc=None,
            freq=-1,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=True,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['qTotal'], DATA_FILE[infile]['qTotal'])

    def test_Extract_CPUtime(self):
        files = ['CCSD_Ethanol_dal.out', 'CCSD_Ethanol_exci_gaus.out', 'CCSD_Ethanol_gaus.out', 'CCSD_Ethanol_lsdal.out', 'CCSD_Ethanol_orca.out', 'CCSD_Methane_dal.out', 'CCSD_Methane_exci_gaus.out', 'CCSD_Methane_gaus.out', 'CCSD_Methane_lsdal.out', 'CCSD_Methane_orca.out', 'CCSD_Water_dal.out', 'CCSD_Water_exci_gaus.out', 'CCSD_Water_gaus.out', 'CCSD_Water_lsdal.out', 'CCSD_Water_orca.out', 'DFT_Ethanol_exci_dal.out', 'DFT_Ethanol_exci_gaus.out', 'DFT_Ethanol_exci_lsdal.out', 'DFT_Ethanol_exci_orca.out', 'DFT_Ethanol_gaus.out', 'DFT_Ethanol_lsdal.out', 'DFT_Ethanol_opt_lsdal.out', 'DFT_Ethanol_opt_velox.out', 'DFT_Ethanol_orca.out', 'DFT_Ethanol_pol_lsdal.out', 'DFT_Ethanol_pol_velox.out', 'DFT_Ethanol_vib_dal.out', 'DFT_Methane_exci_dal.out', 'DFT_Methane_exci_gaus.out', 'DFT_Methane_exci_lsdal.out', 'DFT_Methane_exci_orca.out', 'DFT_Methane_gaus.out', 'DFT_Methane_lsdal.out', 'DFT_Methane_opt_lsdal.out', 'DFT_Methane_opt_velox.out', 'DFT_Methane_orca.out', 'DFT_Methane_pol_lsdal.out', 'DFT_Methane_pol_velox.out', 'DFT_Methane_vib_dal.out', 'DFT_Water_exci_dal.out', 'DFT_Water_exci_gaus.out', 'DFT_Water_exci_lsdal.out', 'DFT_Water_exci_orca.out', 'DFT_Water_gaus.out', 'DFT_Water_lsdal.out', 'DFT_Water_opt_lsdal.out', 'DFT_Water_opt_velox.out', 'DFT_Water_orca.out', 'DFT_Water_pol_lsdal.out', 'DFT_Water_pol_velox.out', 'DFT_Water_vib_dal.out', 'HF_Ethanol_dal.out', 'HF_Ethanol_gaus.out', 'HF_Ethanol_lsdal.out', 'HF_Ethanol_opt_dal.out', 'HF_Methane_dal.out', 'HF_Methane_gaus.out', 'HF_Methane_lsdal.out', 'HF_Methane_opt_dal.out', 'HF_Water_dal.out', 'HF_Water_gaus.out', 'HF_Water_lsdal.out', 'HF_Water_opt_dal.out', 'MP2_Ethanol_dal.out', 'MP2_Ethanol_gaus.out', 'MP2_Ethanol_lsdal.out', 'MP2_Methane_dal.out', 'MP2_Methane_gaus.out', 'MP2_Methane_lsdal.out', 'MP2_Water_dal.out', 'MP2_Water_gaus.out', 'MP2_Water_lsdal.out', 'RIMP2_Ethanol_lsdal.out', 'RIMP2_Methane_lsdal.out', 'RIMP2_Water_lsdal.out']

        args = Namespace(cpu_time=['m'],
            dipole=False,
            energy=False,
            enthalpy=False,
            entropy=False,
            exc=None,
            freq=None,
            gibbs=False,
            infile=[f'test_systems/{file}' for file in files],
            multiprocessing=False,
            optgeom=False,
            osc=False,
            partfunc=False,
            polar=False,
            quiet=True,
            save='return',
            temp=298.15,
            zpv=False,
            progressbar=False,
            unittest=True)

        Values = cd.Extract(args)

        for infile in DATA_FILE:
            self.assertEqual(Values[f'test_systems/{infile}']['total_cpu_time'], DATA_FILE[infile]['total_cpu_time'])
            self.assertEqual(Values[f'test_systems/{infile}']['wall_cpu_time'], DATA_FILE[infile]['wall_cpu_time'])


TEST_DATA = "test_data.json"
DATA_FILE = ReadJSONFile(TEST_DATA)

if __name__ == '__main__':
    unittest.main()
