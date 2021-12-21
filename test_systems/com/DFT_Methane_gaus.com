%chk=DFT_Methane_gaus.chk
%nprocshared=4
%mem=4GB
#p opt freq b3lyp/cc-pvdz geom=connectivity polar

Methane.xyz

0 1
 C                  1.01825000    0.09315000    0.03147000
 H                  2.11045000    0.09315000    0.03147000
 H                  0.65419000   -0.60280000    0.79043000
 H                  0.65418000   -0.21615000   -0.95071000
 H                  0.65418000    1.09840000    0.25470000

 1 2 1.0 3 1.0 4 1.0 5 1.0
 2
 3
 4
 5

