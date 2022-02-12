# Quantum-chemistry-output-parser
This repository was created to collect all the scripts and programs developed and used in Professor Kurt V. Mikkelsens group at the University of Copenhagen.
## [collect_data.py](./collect_data.py)
<details><summary> Program information </summary>
<p>
  A script designed to make it easier to extract data from output files

  Currently the following has been implemented:<br/>
  | Data types                      |       ORCA       |     GAUSSIAN     |      DALTON      |     LSDALTON     |
  |:--------------------------------|:----------------:|:----------------:|:----------------:|:----------------:|
  | Total energies                  |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|
  | Zero-Point Vibrational energies |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |
  | Enthalpies                      |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |
  | Entropies                       |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |
  | Gibbs Free energies             |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |
  | Dipole moments                  |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|
  | Polarizabilities                |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|
  | Excitation energies             |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|
  | Oscillator strengths            |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|
  | Frequencies                     |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |
  | Partition functions             |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |

  **N/A means not applicable*

  Some more advanced functions are:
  - UV/VIS Spectra
    - Requires excitation energies and oscillator strengths in the output file
    - It is possible to choose between different formats for the figure

  The data you want extracted is done using keywords when calling the script. The keywords you call will be printed either in the terminal or written to a csv file.
</p>
</details>

## [sandwich.py](./sandwich.py)
<details><summary> Program information </summary>
<p>
  A script designed to place nanoparticles on either side of a molecule

  Takes the molecule as a xyz file, the two atoms the nanoparticles will be aligned with and the diameter of the particles (in that order).

  #### Keywords

  By default the atomnumbers used to choose alignment is those shown in molden. If instead you wish to choose by the linenumbers as they are in the xyz file you can use the *-l* or *--linenumber* keywords. <br/>
  As default the basis set pc-1 and the RI basis set pc-1-RI will be used. This can be changed with the keywords *--basis and* --ribasis accordingly. <br/>
  An xyz file containing all the information about the junvtion will also be saved, this can be turned off by supplying the keyword *--returnxyz*. <br/>
  If the nanoparticles are spherical in nature (such as Au, Ag & Cu contrary to TiO<sub>2</sub> which is a slab) they will by default turn inwards towards the molecule. For the   nanoparticles to turn outwards the keyword *--outwards* can be supplied. <br/>
  Furthermore the charge of the molecule in the junction is by default 0, this can be changed using the *--charge* keyword <\br>
</p>
</details>

## [xyz2povray.py](./xyz2povray.py)
<details><summary> Program information </summary>
<p>
  A script designed to convert a xyz file to a pov file for the programme POV-Ray which can be used ot ake visually pretty graphics

  The only argument you have to provide is the xyz file

  Apart from this the script will also automatically start generating the figures requested using some antialiasing settings applied in the script
</p>
</details>

## [pov-editor.py](./pov-editor.py)
<details><summary> Program information </summary>
<p>
  A script designed to take the camera position of an existing pov file and update the graphics arguments of said file

  You need to supply two arguments. The pov file wherein the camera position is located and the xyz file so the script can generate the updated graphics.

  This script is especially useful in conjunction with either imol (which only exist for Mac) or Avogadro. In both programes you can export a certain view as a pov file. This is where the camera position is located.

  Apart from this the script will also automatically start generating the figures requested using some antialiasing settings applied in the script
</p>
</details>