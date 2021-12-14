%chk=Methane.chk
%nprocshared=4
%mem=4GB
#p opt freq hf/cc-pvdz geom=connectivity polar

Methane.xyz

0  1
C           1.01825         0.09315         0.03147
H           2.11045         0.09315         0.03147
H           0.65419        -0.60280         0.79043
H           0.65418        -0.21615        -0.95071
H           0.65418         1.09840         0.25470

