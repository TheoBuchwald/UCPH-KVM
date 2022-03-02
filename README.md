# Quantum-chemistry-output-parser
This repository was created to collect all the scripts and programs developed and used in Professor Kurt V. Mikkelsens group at the University of Copenhagen.

Once you have cloned the repository you will have to run ./build.sh to initialize the repository.

You should also do this when pulling new updates, as they may have changed the build package.

When you wish to push an update, to ensure that it works, you may want to use ./test.sh in the root directory of the repository and check for any errors or failures in the different tests.

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

  The data you want extracted is done using keywords when calling the script. The keywords you call will be printed either in the terminal or written to a csv or npz file.
</p>
</details>

## [sandwich.py](./junctions/sandwich.py)
<details><summary> Program information </summary>
<p>
  A script designed to place nanoparticles on either side of a molecule

  Takes the molecule as a xyz file, the two atoms the nanoparticles will be aligned with and the diameter of the particles (in that order).

  #### Keywords

  By default the atomnumbers used to choose alignment is those shown in molden. If instead you wish to choose by the linenumbers as they are in the xyz file you can use the *-l* or *--linenumber* keywords. <br/>
  As default the basis set pc-1 will be used. This can be changed with the keyword *--basis*. <br/>
  An xyz file containing all the information about the junction will also be saved, this can be turned off by supplying the keyword *--returnxyz*. <br/>
  If the nanoparticles are spherical in nature (such as Au, Ag & Cu contrary to TiO<sub>2</sub> which is a slab) they will by default turn inwards towards the molecule. For the nanoparticles to turn outwards the keyword *--outwards* can be supplied. <br/>
  Furthermore the charge of the molecule in the junction is by default 0, this can be changed using the *--charge* keyword <\br>
</p>
</details>

## [leftright.py](./junctions/leftright.py)
<details><summary> Program information </summary>
<p>
  A script designed to place nanoparticles on either side of a molecule in two separate files

  Takes the molecule as a xyz file, the two atoms the nanoparticles will be aligned with and the diameter of the particles (in that order).

  #### Keywords

  By default the atomnumbers used to choose alignment is those shown in molden. If instead you wish to choose by the linenumbers as they are in the xyz file you can use the *-l* or *--linenumber* keywords. <br/>
  As default the basis set pc-1 will be used on the atoms in the molecule while the LANL2DZ and LANL-ECP basis sets will be used on the atoms in the nanoparticles. This can be changed with the keywords *--basis*, *--NPbasis*, and *--ECPbasis* accordingly. <br/>
  The CPU and memory options can be changed from the default of 16 CPU and 16 GB memory with the keywords *--cpu* and *--mem*. <br/>
  An xyz file containing all the information about the junction will also be saved, this can be turned off by supplying the keyword *--returnxyz*. <br/>
  If the nanoparticles are spherical in nature (such as Au, Ag & Cu contrary to TiO<sub>2</sub> which is a slab) they will by default turn inwards towards the molecule. For the nanoparticles to turn outwards the keyword *--outwards* can be supplied. <br/>
  Furthermore the charge of the molecule in the junction is by default 0, this can be changed using the *--charge* keyword <\br>
</p>
</details>

## [xyz_to_gauss.py](./xyz/xyz_to_gauss.py)
<details><summary> Program information </summary>
<p>
  A script designed to convert a xyz file to a com input file for the Gaussian suite of programs

  You will need to supply the xyz file and keywords. Other options can be added via the command line. Use -h on the script to see the available options.

  You can also supply basis sets not implemented in Gaussian, in which case an API to BSE (https://www.basissetexchange.org/) is used.
</p>
</details>

## [xyz_to_mol.py](./xyz/xyz_to_mol.py)
<details><summary> Program information </summary>
<p>
  A script designed to convert a xyz file to a mol file for the programme DALTON

  You will need to supply the xyz file

  Apart from this, you can also supply a basis set, RI-basis set and the charge with the keywords *--basis*, *--RIbasis*, and *--charge*
</p>
</details>

## [xyz_to_molpro.py](./xyz/xyz_to_molpro.py)
<details><summary> Program information </summary>
<p>
  A script designed to convert a xyz file to a molpro file

  You will need to supply the xyz file as well as a keywords nr. to determine the options for the programme
</p>
</details>

## [xyz_to_orca.py](./xyz/xyz_to_orca.py)
<details><summary> Program information </summary>
<p>
  A script designed to convert a xyz file to a inp file for the programme ORCA

  You will need to supply the xyz file as well as a keywords nr. to determine the options for the programme

  Apart from this, you can also supply a charge and memory limits with the keywords *--charge* and *--mem*

  If you want extra calculations you can supply either of the keywords *--extra1* and *--extra2*
</p>
</details>

## [xyz_to_povray.py](./visualization/xyz_to_povray.py)
<details><summary> Program information </summary>
<p>
  A script designed to convert a xyz file to a pov file for the programme POV-Ray which can be used ot ake visually pretty graphics

  The only argument you have to provide is the xyz file(s)

  Apart from this the script will also automatically start generating the figures requested using some antialiasing settings applied in the script. Those settings are:

    +A0.1 +AM2 +AG0 +R5 -J

  +A0.1: Antialliasing set to 0.1 threshold<br/>
  +AM2: Antialiasing method 2<br/>
  +AG0: Gamma set to 0<br/>
  +R5: Depth set to 5<br/>
  -J: Jitter set to off
</p>
</details>

## [pov_editor.py](./visualization/pov_editor.py)
<details><summary> Program information </summary>
<p>
  A script designed to take the camera position of an existing pov file and update the graphics arguments of said file

  You need to supply two arguments. The pov file wherein the camera position is located and the xyz file so the script can generate the updated graphics.

  This script is especially useful in conjunction with either imol (which only exist for Mac) or Avogadro. In both programes you can export a certain view as a pov file. This is where the camera position is located.

  Apart from this the script will also automatically start generating the figures requested using some antialiasing settings applied in the script Those settings are:

    +A0.1 +AM2 +AG0 +R5 -J

  +A0.1: Antialliasing set to 0.1 threshold<br/>
  +AM2: Antialiasing method 2<br/>
  +AG0: Gamma set to 0<br/>
  +R5: Depth set to 5<br/>
  -J: Jitter set to off
</p>
</details>
