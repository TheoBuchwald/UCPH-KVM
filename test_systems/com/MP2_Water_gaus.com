%chk=MP2_Water_gaus.chk
%nprocshared=4
%mem=4GB
#p opt freq mp2/cc-pvdz geom=connectivity polar

Water.xyz

0 1
 O                  0.90644000   -0.06978000    0.02183000
 H                  1.87433000   -0.04844000   -0.01941000
 H                  0.62759000    0.35671000   -0.80237000

 1 2 1.0 3 1.0
 2
 3

