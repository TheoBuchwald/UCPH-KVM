# Quantum-chemistry-output-parser
This repository was created to collect all the scripts and programs developed and usen in Proffessor Kurt V. Mikkelsens group at the University of Copenhagen.
## [collect_data.py](./collect_data.py)
<details><summary> Program information </summary>
<p>
  A script designed to make it easier to extract data from .out files

  Currently the following has been implemented:<br/>
  | Data types                      |       ORCA       |     GAUSSIAN     |      DALTON      |     LSDALTON     |
  |:--------------------------------|:----------------:|:----------------:|:----------------:|:----------------:|
  | Total energies                  |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|
  | Zero-Point Vibrational energies |:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |
  | Enthalpies                      |:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |
  | Entropies                       |        :x:       |        :x:       |        :x:       |        :x:       |
  | Gibbs Free energies             |:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |
  | Dipole moments                  |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        :x:       |
  | Polarizabilities                |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        :x:       |
  | Excitation energies             |:heavy_check_mark:|        :x:       |:heavy_check_mark:|        :x:       |
  | Oscillator strengths            |:heavy_check_mark:|        :x:       |:heavy_check_mark:|        :x:       |
  | Frequencies                     |:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |
  | Partition functions             |        :x:       |:heavy_check_mark:|        :x:       |        :x:       |

  Some more advanced functions are:
  - UV/VIS Spectra
    - Requires excitation energies and oscillator strengths in the .out file
</p>
</details>

## [sandwich.py](./sandwich.py)
<details><summary> Program information </summary>
<p>
  A script designed to make nanoparticles on either side of a molecule
  
  Takes the molecule as a xyz file, the two atoms the nanoparticles will be aligned with and the diameter of the particles (in that order).
  
  #### Keywords
  
  By default the atomnumbers used to choose alignment is those shown in molden. If instead you wish to choose by the linenumbers as they are in the xyz file you can use the *-l* or *--linenumber* keywords. <br/>
  As default the basis set pc-1 and the RI basis set pc-1-RI will be used. This can be changed with the keywords *--basis and* --ribasis accordingly. <br/>
  An xyz file containing all the information about the junvtion will also be saved, this can be turned off by supplying the keyword *--returnxyz*. <br/>
  If the nanoparticles are spherical in nature (such as Au, Ag & Cu contrary to TiO<sub>2</sub> which is a slab) they will by default turn inwards towards the molecule. For the   nanoparticles to turn outwards the keyword *--outwards* can be supplied. <br/>
  Furthermore the charge of the molecule in the junction is by default 0, this can be changed using the *--charge* keyword <\br>
</p>
</details>
