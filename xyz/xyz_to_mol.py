import argparse
from collections import Counter #For number of unique elements
from Kurt import chemical_information as ci

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''A script to convert xyz files to mol files for DALTON''', epilog='''For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='.xyz file')
    parser.add_argument('--charge', default=[0], nargs=1, type=int, help='Include to specify charge - 0 if not included')
    parser.add_argument('--basis', default=['pc-1'], nargs=1, type=str, help='Include to specify basis set of the molecular atoms - pc-1 if not included')
    parser.add_argument('--RIbasis', nargs=1, type=str, help='Include to specify basis set of the molecular atoms - RI-BASIS if not included')

    args = parser.parse_args()

    input_files = args.infile
    basis = args.basis[0]
    charge = args.charge[0]

    basis_sets_dal = ['3-21G', '3-21G*', '3-21++G', '3-21++G*', '4-31G', '6-311+G*', '6-311G', '6-311G*', '6-311G**', '6-311++G**', '6-311++G(2d,2p)', '6-311G(2df,2pd)', '6-311++G(3df,3pd)', '6-31+G', '6-31+G*', '6-31G', '6-31G*', '6-31G**', '6-31++G', '6-31++G*', '6-31++G**', '6-31G(3df,3pd)', 'Ahlrichs-Coulom-Fit', 'Ahlrichs-pVDZ', 'Ahlrichs-TZV', 'Ahlrichs-VDZ', 'Ahlrichs-VTZ', 'Almlof-Taylor-ANO', 'ano-1', 'ano-2', 'ano-3', 'ano-4', 'ANO-DK3', 'ANO-RCC', 'apc-0', 'apc-1', 'apc-2', 'apc-3', 'apc-4', 'aug-ccJ-pVTZ', 'aug-cc-pCV5Z', 'aug-cc-pCVDZ', 'aug-cc-pCVQZ', 'aug-cc-pCVTZ', 'aug-cc-pCVTZ-CTOCD-uc', 'aug-cc-pV(5+d)Z', 'aug-cc-pV5Z', 'aug-cc-pV5Z-DK', 'aug-cc-pV5Z_OPTRI', 'aug_cc_pv5z_pp', 'aug-cc-pV5Z-RI', 'aug-cc-pV(6+d)Z', 'aug-cc-pV6Z', 'aug-cc-pV6Z-RI', 'aug-cc-pV(D+d)Z', 'aug-cc-pVDZ', 'aug-cc-pVDZ-DK', 'aug-cc-pVDZ_OPTRI', 'aug_cc_pvdz_pp', 'aug-cc-pVDZ-RI', 'aug-cc-pV(Q+d)Z', 'aug-cc-pVQZ', 'aug-cc-pVQZ-DK', 'aug-cc-pVQZ_OPTRI', 'aug_cc_pvqz_pp', 'aug-cc-pVQZ-RI', 'aug-cc-pV(T+d)Z', 'aug-cc-pVTZ', 'aug-cc-pVTZ-DK', 'aug-cc-pVTZ-J', 'aug-cc-pVTZ-lresc', 'aug-cc-pVTZ_OPTRI', 'aug_cc_pvtz_pp', 'aug-cc-pVTZ-RI', 'aug-cc-pwCV5Z-RI', 'aug-cc-pwCVQZ-RI', 'aug-cc-pwCVTZ-RI', 'aug-pc-0_emsl', 'aug-pc-1_emsl', 'aug-pc-2_emsl', 'aug-pc-3_emsl', 'aug-pc-4_emsl', 'aug-pcH-1', 'aug-pcH-2', 'aug-pcH-3', 'aug-pcH-4', 'aug-pcJ-0', 'aug-pcJ-1', 'aug-pcJ-2', 'aug-pcJ-3', 'aug-pcJ-4', 'aug-pcseg-0', 'aug-pcseg-1', 'aug-pcseg-2', 'aug-pcseg-3', 'aug-pcseg-4', 'aug-pV7Z', 'aug-Test', 'ccJ-pV5Z', 'ccJ-pVDZ', 'ccJ-pVQZ', 'ccJ-pVTZ', 'cc-pCV5Z', 'cc-pCVDZ', 'cc-pCVQZ', 'cc-pCVTZ', 'cc-pV(5+d)Z', 'cc-pV5Z', 'cc-pV5Zdenfit', 'cc-pV5Z-DK', 'cc_pv5z_pp', 'cc-pV5Z-RI', 'cc-pV(6+d)Z', 'cc-pV6Z', 'cc-pV6Z-RI', 'cc-pV(D+d)Z', 'cc-pVDZ', 'cc-pVDZ-DK', 'cc-pVDZ_emsl', 'cc-pVDZ-F12', 'cc-pVDZ-F12_OPTRI', 'cc_pvdz_pp', 'cc-pVDZ-RI', 'cc-pVDZ-uncontracted', 'cc-pV(Q+d)Z', 'cc-pVQZ', 'cc-pVQZdenfit', 'cc-pVQZ-DK', 'cc-pVQZ-F12', 'cc-pVQZ-F12_OPTRI', 'cc_pvqz_pp', 'cc-pVQZ-RI', 'cc-pV(T+d)Z', 'cc-pVTZ', 'cc-pVTZ_apra', 'cc-pVTZdenfit', 'cc-pVTZ-DK', 'cc-pVTZ_emsl', 'cc-pVTZ-F12', 'cc-pVTZ-F12_OPTRI', 'cc_pvtz_pp', 'cc-pVTZ-RI', 'cc-pwCV5Z', 'cc-pwCV5Z-DK', 'cc_pwcv5z_pp', 'cc-pwCVDZ', 'cc_pwcvdz_pp', 'cc-pwCVQZ', 'cc-pwCVQZ-DK', 'cc_pwcvqz_pp', 'cc-pwCVTZ', 'cc-pwCVTZ-DK', 'cc_pwcvtz_pp', 'crenl_ecp', 'crens_ecp', 'def2_qzvp', 'def2_qzvpp', 'def2_sv_p', 'def2_svp', 'def2_tzvp', 'def2_tzvpd', 'def2_tzvpp', 'df-def2', 'diffs', 'DZ(Dunning)', 'DZP+Diffuse(Dunning)', 'DZP(Dunning)', 'DZP+Ryderg(Dunning)', 'dzq', 'DZ+Ryderg(Dunning)', 'ecp-sdd-DZ', 'GAMESS-PVTZ', 'GAMESS-VTZ', 'hay_wadt_m_n1_ecp', 'hay_wadt_vdz_n1_ecp', 'Huckel', 'Huckel_old', 'Huz-II', 'Huz-III', 'Huz-IIIsu3', 'Huz-IIsu2', 'Huz-IV', 'Huz-IVsu4', 'lanl08', 'lanl08d', 'lanl08_f', 'lanl08_p', 'lanl2dzdp_ecp', 'lanl2dz_ecp', 'lanl2tz', 'lanl2tz_f', 'lanl2tz_p', 'loprop-6-31+G*', 'loprop-6-31G*', 'loprop-aug-cc-pVDZ', 'loprop-aug-cc-pVQZ', 'loprop-aug-cc-pVTZ', 'loprop-augx-cc-pVDZ', 'loprop-cc-pVDZ', 'loprop-cc-pVTZ', 'McLean-Chandler-VTZ', 'MINI(Huzinaga)', 'MINI(Scaled)', 'modified_lanl2dz', 'N07D:B3LYP', 'NASA-Ames-ANO', 'NQvD', 'pc-0', 'pc-0_emsl', 'pc-1', 'pc-1_emsl', 'pc-2', 'pc-2_emsl', 'pc-3', 'pc-3_emsl', 'pc-4', 'pc-4_emsl', 'pcH-1', 'pcH-2', 'pcH-3', 'pcH-4', 'pcJ-0', 'pcJ-1', 'pcJ-2', 'pcJ-3', 'pcJ-4', 'pcS-0', 'pcS-1', 'pcS-2', 'pcS-3', 'pcS-4', 'pcseg-0', 'pcseg-1', 'pcseg-2', 'pcseg-3', 'pcseg-4', 'pV6Z', 'pV7Z', 'raf-r', 'Sadlej-pVTZ', 'Sadlej-pVTZ-J', 'skjc_polarized_p_2d_lfk', 'skjc_vdz_ecp', 'sd_aug_cc_pvqz', 'sd_aug_cc_pvtz', 'sd_cc_pvqz', 'sd_cc_pvtz', 'STO-2G', 'STO-3G', 'STO-6G', 'stuttgart_rlc_ecp', 'stuttgart_rsc_1997_ecp', 'stuttgart_rsc_ano_ecp', 'stuttgart_rsc_segmented_ecp', 'SV+DouleRyderg(Dunning-Hay)', 'SV(Dunning-Hay)', 'SVP+Diffuse(Dunning-Hay)', 'SVP(Dunning-Hay)', 'SV+Ryderg(Dunning-Hay)', 'Turomole-DZ', 'Turomole-DZP', 'Turomole-QZV', 'Turomole-QZVP', 'Turomole-QZVPP', 'Turomole-SV', 'Turomole-SVP', 'Turomole-TZ', 'Turomole-TZP', 'Turomole-TZV', 'Turomole-TZV2P', 'Turomole-TZVP', 'Turomole-TZVPP', 'Turomole-TZVPPP', 'TZ(Dunning)', 'UnitTest_genD', 'UnitTest_genF', 'UnitTest_genP', 'UnitTest_genS', 'UnitTest_segD', 'UnitTest_segD1p', 'UnitTest_segF', 'UnitTest_segF1p', 'UnitTest_segP', 'UnitTest_segP1p', 'UnitTest_segS', 'UnitTest_segS1p', 'UnitTest_segSP', 'Wachtersa+f']
    BSE = True

    if basis in basis_sets_dal:
        BSE = False

    if BSE:
        #Use the BSE API to get the basis set
        print("Need to get basis set from BSE!")

    if args.RIbasis:
        RIbasis = args.RIbasis[0]
    else:
        RIbasis = f'{basis}-RI'

    for molfile in input_files:
        namesmol = []
        molx, moly, molz = [], [], []

        #Convert label into atomic number

        #Read in xyz coordinates
        with open(molfile, "r") as f:
            lines = f.readlines()
        for i in range(2, len(lines)):
            x = lines[i].split()
            namesmol.append(x[0])
            molx.append(float(x[1]))
            moly.append(float(x[2]))
            molz.append(float(x[3]))

        #Lines to be inserted in .mol file
        lines_mol =[]

        #Starting lines
        lines_mol.append('ATOMBASIS\n')
        lines_mol.append('./'+molfile+'\n')
        lines_mol.append('Hej Theo\n')
        lines_mol.append('Atomtypes='+str(len(set(namesmol)))+f' Charge={charge} NoSymmetry Angstrom\n')

        unique_indices = list(Counter(namesmol).keys())

        #Write coordinates and basis set
        for i in unique_indices:
            ind = unique_indices.index(i)
            count = list(Counter(namesmol).values())[ind]
            if not BSE:
                lines_mol.append(f'  {ci.AtomicInformation(i).atomnr():.4f}     {count} Bas={basis} Aux={RIbasis}\n')
            else:
                BasisSet = ci.BasisSet()
                try:
                    basis_mol = BasisSet.AtomBasisSet('dalton', basis, i, SupressHeader=True)
                except RuntimeError:
                    print("Failed to get basis set from BSE. Please check the spelling, upper-/lowercase IS important")
                    print('The problem may also be that the basis set does not exist for the given atom')
                    exit()
                BlockTypes = ['s functions', 'p functions', 'd functions', 'f functions', 'g functions', 'h functions', 'i functions', 'j functions', 'k functions']
                blocks = 0
                for j in BlockTypes:
                    if j in basis_mol:
                        blocks += 1
                Block = f'{blocks}'
                for j in BlockTypes[:blocks]:
                    Block += f' {basis_mol.count(j)}'
                lines_mol.append(f'  Charge={ci.AtomicInformation(i).atomnr():.4f}     Atoms={count}     Blocks={Block}\n')
                basis_mol = basis_mol.replace('H','').split('\n')[5:-2]
                basis_mol = [i for i in basis_mol if 'functions' not in i]
            for j, atom in enumerate(namesmol):
                if atom == i:
                    lines_mol.append(''.join([atom.ljust(2),' ',f"{molx[j]:.9f}".rjust(14),' ', f"{moly[j]:.9f}".rjust(19), ' ',f"{molz[j]:.9f}".rjust(19) ,'\n']))
            if BSE:
                for j in basis_mol:
                    lines_mol.append(f'{j}\n')

        #Output is .mol
        with open(molfile[:-4] + '_' + basis + '.mol','w') as f:
            f.writelines(lines_mol)
