#! /bin/bash

# Simple script to do clean up of files

make clean
rm -r system.save
rm *ipt  *.txt *.xpts *.els fort.* *.dat
rm bvecs cpbd decut ebdat enkfile epsilon err gvectors gwipt ladcap ldaclips pdadat prog rhoofr tmels vecy
rm kandb.h niter.h omega.h rsval.h eps1 eps2 inds loss loss_SE opcons refl gap  
