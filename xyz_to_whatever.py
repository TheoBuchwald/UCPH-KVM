import numpy as np
import sys
from collections import Counter #For number of unique elements
import os
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='''A script to make inputfiles for orca, dalton and molpro using a xyz file''',epilog='''For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to convert to inputs', metavar='.xyz file')

    parser.add_argument('-dal', '--dalton', action='store_true', help='Include to create .mol files for DALTON')
    # parser.add_argument('-molpro', '--molpro', action='store_true', help='Include to create MOLPRO files')
    parser.add_argument('-orca', action='store_true', help='Include to create .inp files for ORCA')

    parser.add_argument('-hf','--hartree', action='store_true', help='Include to do a RHF calculation')
    parser.add_argument('-mp2', action='store_true', help='Include to do a MP2 calculation')
    parser.add_argument('-ccsd', action='store_true', help='Include to do a CCSD calculation')
    parser.add_argument('-ccsdpt', action='store_true', help='Include to do a CCSD(T) calculation')
    parser.add_argument('-cc2', action='store_true', help='Include to do a CC2 calculation')
    parser.add_argument('-cc3', action='store_true', help='Include to do a CC3 calculation')
    parser.add_argument('-dft', '--functional', type=str, nargs='*', help='Include to do a DFT calculation. Write the functionals wished after')
    parser.add_argument('-basis','--basisset', type=str, nargs=1, help='The basis set used', default='pc-1')

    parser.add_argument('--unrestricted', action='store_true', help='Include to change to unrestricted methods')

    parser.add_argument('-P', action='store_true', help='Include to calculate polarizabilities')
    parser.add_argument('-X', action='store_true', help='Include to calculate excitation energies')
    parser.add_argument('-F', action='store_true', help='Include to calculate frequencies')
    parser.add_argument('-g', action='store_true', help='Include to do a geometry calculation')
    parser.add_argument('-add', '--additional', type = str, nargs='*', help='Include to store additional arguments')

    args = parser.parse_args()

    input_file = args.infile

    Wavefuntion_Methods = {
        'RHF' : args.hartree,
        'MP' : args.mp2,
        'CCSD' : args.ccsd,
        'CCSD(T)' : args.ccsdpt,
        'CC2' : args.cc2,
        'CC3' : args.cc3,
    }

    Basis = args.basisset
    DFT = args.functional
    Additional = args.additional
    Unrestricted = args.unrestricted

    Calc = {
        'geo' : args.g,
        'freq' : args.F,
        'pol' : args.P,
        'exc' : args.X
    }

    if DFT == None:
        Methods = [i[0] for i in Wavefuntion_Methods.items() if i[1] != False]
    else:
        Methods = [i[0] for i in Wavefuntion_Methods.items() if i[1] != False] + DFT


class orca:
    def __init__(self, input):
        file = open(input, "r")
        self.lines = file.readlines()
        file.close
        self.text = '! '
        if Unrestricted == True:
            self.text += 'U'
        self.text
    
    def options(self, method):
        if method in Wavefuntion_Methods:
            self.text += method + ' ' + Basis
            if Additional != None:
                self.text += ' ' + ' ' .join(Additional)
                return
        self.text += method + ' ' + Basis
        if Additional != None:
            self.text += 'KS ' + ' ' .join(Additional)

    

    

class dal:
    def __init__(self,input):
        file = open(input, "r")
        self.lines = file.readlines()
        file.close


if __name__ == '__main__':

    Outputs = {
        dal : args.dalton,
        # 'molpro' : args.molpro,
        orca : args.orca
    }
    
    Out = [i[0] for i in Outputs.items() if i[1] != False]

    for prog in Out:
        for infile in input_file:
            for method in Methods:
                file = prog(infile)
                file.options(method)

                print(file.text)