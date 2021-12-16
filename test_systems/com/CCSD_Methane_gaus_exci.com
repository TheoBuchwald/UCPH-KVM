%chk=CCSD_Methane_gaus.chk
%nprocshared=4
%mem=4GB
#p eomccsd=(nstate=10)/cc-pvdz guess=read geom=allcheck polar


