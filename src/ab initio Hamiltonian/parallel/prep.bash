#!/bin/bash


echo "********** The Quantum Espresso calculations have been successfully completed. ***********"
echo " "
echo "                Now generating all the executables"
echo " "
make
echo " All done with make"
echo " "
echo " 1.) ******  Generating a bunch of input files ******** "

## Get number of kpoints from the DFT system.save
nkpt=($(ls system.save/|wc -l))
nkpts=$(($nkpt - 3))
echo "nkpts = ",$nkpt
echo
echo " actual nkpts =",$nkpts
echo "$nkpts" >> nkpt.ipt

./datagen.x
echo " "
echo "   : successfully generated all the input files "
echo " " 
cp efermiinrydberg.ipt efryd.ipt

echo " 2.) ******  Generating the transition matrix elements file - tmels ******  "
./tmelnew.x
echo " "
echo "   : successfully generated the tmels file "
echo " " 
echo " 3.) ******  Generating more input files and reordering eigen values and tmels ****** "
./setup2.x
echo " "
echo "   : successfully completed setup "
echo " " 
echo " 4.) ****** Converting u1.dat from u of g to u of x  ******"
./conugtoux.x
echo " "
echo "   : successfully converted u of g to u of x and saved in u1.dat "
echo " " 
echo " 5.) ****** Orthonormalizing the wavefunctions and storing in u2.dat ******"
./orthog.x
echo " "
echo "   : successfully orthonormalized wavefunciton coefficients "
echo " " 
echo " 6.) ****** Converting the density from rhoofr to rho : x, y, z ******"
./condens2.x
echo " "
echo "   : successfully converted the density file "
echo " " 
echo " 6.) ****** Calculating the HLL dielectric function : generating W.els file ******"
./ll.x
echo " "
echo "   : successfully generated the W.els file "
echo " " 
echo " 7.) ****** Calculating the HLL dielectric function : generating W0.els file ******"
./whom0.x
echo " "
echo "   : successfully generated the W0.els file "
echo " " 
echo " 8.) ****** Calculating the GW correction to DFT bandgap ******"
./bridgegw.x
echo " "
echo "   : successfully generated the eb.dat, bs.dat files "
echo " " 
echo " " 
echo " 9.) ****** Calculating the ladder.dat file ******"
./bridgelad.x
echo " "
echo "   : successfully generated the ladder.dat file "
echo " " 
echo " " 
