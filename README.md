# Quantum-chemistry-output-parser

This repository was created to collect all the scripts and programs developed and used in Professor, Ph.D., Dr. Scient. Kurt V. Mikkelsens group at the University of Copenhagen.

Many of the scripts in this repository are dependent on functions also found in the [KurtGroup](https://pypi.org/project/KurtGroup/) Python package. It is not necessary to install this Python package to use these scripts, however if you wish to use some of the functionalities this package supplies for your own personal scripts it might be an idea to install it using either of

```
pip install KurtGroup
pip install KurtGroup --user
```

To install a specific version you can use

```
pip install KurtGroup==version
```

It can be updated using one of the following commands

```
pip install KurtGroup -U
pip install KurtGroup --user -U
```

The KurtGroup package information and source code is located in the folder [KurtGroup](./KurtGroup/). More information about this package can be found on the [PyPI web page](https://pypi.org/project/KurtGroup) or in the [README](./KurtGroup/README.md) in said folder.

When you wish to push an update, to ensure that it works, you may want to use ./[test.sh](./test.sh) in the root directory of the repository and check for any errors or failures in the different tests.

## [collect_data.py](./collect_data.py)
<details><summary> Program information </summary>
<p>

  A script designed to make it easier to extract data from output files

  Currently the following has been implemented:<br/>
  | Data types                      |       ORCA       |     GAUSSIAN     |      DALTON      |     LSDALTON     |     VeloxChem    |        AMS       |      Q-Chem      |      DIRAC       |
  |:--------------------------------|:----------------:|:----------------:|:----------------:|:----------------:|:----------------:|:----------------:|:----------------:|:----------------:|
  | Total energies                  |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|
  | Zero-Point Vibrational energies |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |        :x:       |        :x:       |        :x:       |        :x:       |
  | Enthalpies                      |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |        :x:       |        :x:       |        :x:       |        :x:       |
  | Entropies                       |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |        :x:       |        :x:       |        :x:       |        :x:       |
  | Gibbs Free energies             |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |        :x:       |        :x:       |        :x:       |        :x:       |
  | Dipole moments                  |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |
  | Polarizabilities                |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |        :x:       |
  | Excitation energies             |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |:heavy_check_mark:|:heavy_check_mark:|
  | Oscillator strengths            |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |:heavy_check_mark:|:heavy_check_mark:|
  | Frequencies                     |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |        :x:       |        :x:       |        :x:       |        :x:       |
  | Partition functions             |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        N/A       |        :x:       |        :x:       |        :x:       |        :x:       |
  | CPU time used                   |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |:heavy_check_mark:|        :x:       |
  | Optimized geometries            |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |

  **N/A means not applicable*

  When it comes to more advanced features the following has been implemented:

  | Data processing                 |       ORCA       |     GAUSSIAN     |      DALTON      |     LSDALTON     |     VeloxChem    |        AMS       |      Q-Chem      |
  |:--------------------------------|:----------------:|:----------------:|:----------------:|:----------------:|:----------------:|:----------------:|:----------------:|
  | UVVIS using excitation energies |:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|:heavy_check_mark:|        :x:       |        :x:       |:heavy_check_mark:|
  | UVVIS using complex propagators |        :x:       |        :x:       |:heavy_check_mark:|        :x:       |        :x:       |        :x:       |        :x:       |

  It is possible to choose between multiple formats for the spectra (png, eps,...)

  The graph data can also be saved in a npz file using the *-s* or *--save* keyword. Here it will be saved as the [wavelength span, extinction coefficient] for the UVVIS spectra

  Additionally the unit of the x-axis can be changed between nm, eV, and cm^-1 using the *-u* or *--unit* keyword

  The data you want extracted is done using keywords when calling the script. The keywords you call will be printed either in the terminal or written to a csv or npz file.

</p>
</details>

## [permutation_checker.py](./permutation_checker.py)
<details><summary> Program information </summary>
<p>

  A script designed to check and compare the indicies of equations derived from Box 13.2 in *Molecular Electronic Structure Theory*

  #### Keywords

  The keywords -P, -bra, and -E are required arguments and must be given as in the examples:
    -P cde klm or -P cd kl ...
    -bra ai bj or -bra -ai ...
    -E dn or -E dn cl ...

  The keywords -F, -L, -g, -t, -LV, -RV, and -sum are optional and must be given as in the examples:
    -F ci
    -L cile
    -g cile
    -t cile or -t cile dlem ...
    -LV ci
    -RV ck
    -sum cdeklm or -sum -clmedk ...

  If -sum is not provided the unique permutations will not be found


</p>
</details>

## [sandwich.py](./sandwich.py)
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

## [leftright.py](./leftright.py)
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

## [LevelSystem5.py](./LevelSystem5.py)
<details><summary> Program information </summary>
<p>

  For this script you need to manually edit the coupling elements and energies. After having done this you can run the script and generate the Coulomb stairs for your five level system.

</p>
</details>

## [xyz_align.py](./xyz_align.py)
<details><summary> Program information </summary>
<p>

  This script is designed to align two molecules so that the RMSD between them is as small as possible.

  If you find this script useful for any publishable work, please cite the corresponding paper:

  - Berhane Temelso, Joel M. Mabey, Toshiro Kubota, Nana Appiah-padi, George C. Shields
  J. Chem. Info. Model. 2017, 57(5), 1045-1054

</p>
</details>

## [xyz_to_gauss.py](./xyz_to_gauss.py)
<details><summary> Program information </summary>
<p>

  A script designed to convert a xyz file to a com input file for the Gaussian suite of programs

  You will need to supply the xyz file and keywords. Other options can be added via the command line. Use -h on the script to see the available options.

  You can also supply basis sets not implemented in Gaussian, in which case an API to the [Basis Set Exchange](https://www.basissetexchange.org/) is used.

</p>
</details>

## [xyz_to_mol.py](./xyz_to_mol.py)
<details><summary> Program information </summary>
<p>

  A script designed to convert a xyz file to a mol file for the program DALTON

  You will need to supply the xyz file

  Apart from this, you can also supply a basis set and the charge with the keywords *--basis* and *--charge*

  You can also supply basis sets not implemented in DALTON, in which case an API to the [Basis Set Exchange](https://www.basissetexchange.org/) is used.

</p>
</details>

## [xyz_to_molpro.py](./xyz_to_molpro.py)
<details><summary> Program information </summary>
<p>

  A script designed to convert a xyz file to a molpro file

  You will need to supply the xyz file as well as a keywords nr. to determine the options for the program

</p>
</details>

## [xyz_to_orca.py](./xyz_to_orca.py)
<details><summary> Program information </summary>
<p>

  A script designed to convert a xyz file to a inp file for the program ORCA

  You will need to supply the xyz file as well as a keywords nr. to determine the options for the program

  Apart from this, you can also supply a charge and memory limits with the keywords *--charge* and *--mem*

  If you want extra calculations you can supply either of the keywords *--extra1* and *--extra2*

</p>
</details>

## [xyz_to_povray.py](./xyz_to_povray.py)
<details><summary> Program information </summary>
<p>

  A script designed to convert a xyz file to a pov file for the program POV-Ray which can be used to make visually pretty graphics

  The only argument you have to provide is the xyz file(s)

  Apart from this the script will also automatically start generating the figures requested using some antialiasing settings applied in the script. Those settings are:

```
+A0.1 +AM2 +AG0 +R5 -J
```

  +A0.1: Antialliasing set to 0.1 threshold<br/>
  +AM2: Antialiasing method 2<br/>
  +AG0: Gamma set to 0<br/>
  +R5: Depth set to 5<br/>
  -J: Jitter set to off

</p>
</details>

## [pov_editor.py](./pov_editor.py)
<details><summary> Program information </summary>
<p>

  A script designed to take the camera position of an existing pov file and update the graphics arguments of said file

  You need to supply two arguments. The pov file wherein the camera position is located and the xyz file so the script can generate the updated graphics.

  This script is especially useful in conjunction with either imol (which only exist for Mac) or Avogadro. In both programes you can export a certain view as a pov file. This is where the camera position is located.

  Apart from this the script will also automatically start generating the figures requested using some antialiasing settings applied in the script Those settings are:

```
+A0.1 +AM2 +AG0 +R5 -J
```

  +A0.1: Antialliasing set to 0.1 threshold<br/>
  +AM2: Antialiasing method 2<br/>
  +AG0: Gamma set to 0<br/>
  +R5: Depth set to 5<br/>
  -J: Jitter set to off

</p>
</details>

# Deprecated scripts

## [getopt_adf.py](./getopt_adf.py)
<details><summary> Program information </summary>
<p>

  This script is designed to extract the optimized geometry from a geometry
  optimization run in ADF. The same functionality can be found in the [collect_data.py](./collect_data.py) script.

</p>
</details>

## [getopt_gauss.py](./getopt_gauss.py)
<details><summary> Program information </summary>
<p>

  This script is designed to extract the optimized geometry from a geometry
  optimization run in Gaussian. The same functionality can be found in the [collect_data.py](./collect_data.py) script.

</p>
</details>

## [getopt_orca.py](./getopt_orca.py)
<details><summary> Program information </summary>
<p>

  This script is designed to extract the optimized geometry from a geometry
  optimization run in Orca. The same functionality can be found in the [collect_data.py](./collect_data.py) script.

</p>
</details>
