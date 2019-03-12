#!/bin/csh -f
set echo
cd /Users/markh/Desktop/usf/mapman
setenv MAPSIZE 10000000
./mapman << EOF 
read m1 mapman_test.xE ccp4
! normalise m1
drange m1 -2.0 5.0
mappage m1 mapman_test.omap
quit
