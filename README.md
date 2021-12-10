# Quantum-chemistry-output-parser
This repository was created to collect all the scripts and programs developed and usen in Proffessor Kurt V. Mikkelsens group at the University of Copenhagen.
## collect_data.py
<details><summary> Program information </summary>
<p>
  
  A script designed to make it easier to extract data from .out files
  
  Currently the following has been implemented:
  
  | Data types | ORCA | GAUSSIAN | DALTON | LSDALTON |
  |:---|:---:|:---:|:---:|:---:|
  | Total energies | x | x | x | x|
  | Zero-Point Vibrational energies | x | x |   |   |
  | Enthalpies | x | x |   |   |
  | Entropies |   |   |   |   |
  | Gibbs Free energies | x | x |   |   |
  | Dipole moments | x | x | x |   |
  | Polarizabilities | x | x | x |   |
  | Excitation energies | x |   | x |   |
  | Oscillator strengths | x |   | x |   |
  | Frequencies | x | x |   |   |
  | Partition functions |   | x |   |   |
  
  Some more advanced functions are:
  - UV/VIS Spectra
    - Requires excitation energies and oscillator strengths in the .out file
  
  For more information about this script you can contact
  
  Theo Juncker von Buchwald @ fnc970@alumni.ku.dk / theo.buchmail@gmail.com
  
  Magnus Bukhave Johansen @ qhw298@alumni.ku.dk
</p>
</details>
