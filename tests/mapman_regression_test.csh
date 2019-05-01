#! /bin/tcsh -f
#
#    Quick test to make sure the installed version of mapman can actually read 
#    and spline interpolate CCP4 map files
#
#
# allow user to specify executable
set mapman = "$1"
if("$mapman" == "") set mapman = lx_mapman


if(! $?CCP4) then
    echo "Sorry, need CCP4 set up to run this test."
    exit 9
endif


# create a simple HKL file
echo "0 0 2 1 0" >! F.hkl
# turn it into an MTZ
f2mtz hklin F.hkl hklout F.mtz << EOF >! f2mtz.log
CELL 10 10 10 90 90 90
SYMM 1
labout H K L F PHI
ctypou H H H F P
EOF
# calculate a map
fft hklin F.mtz mapout ffted.map << EOF >! fft.log
labin F1=F PHI=PHI
RESO 2
EOF
# expand the map a bit to avoid edge effects in the interpolation
mapmask mapin ffted.map mapout F.map << EOF >! mapmask.log
xyzlim -0.2 1.2 -0.2 1.2 -0.2 1.2
axis X Y Z
EOF
# invent a probe pdb with a line of atoms to do the interpolation
awk 'BEGIN{for(x=0;x<10;x+=0.1){\
  printf("ATOM     10  OW  WAT X  10       0.000   0.000%8.3f  1.00  5.00\n",x)}}' |\
cat >! probe.pdb

# create a simple input file for mapman
cat << EOF >! mapman.in
read map1 F.map ccp4
peek value map1 probe.pdb probed.pdb spline ;
quit
EOF

# make sure mapman can handle the map
setenv MAPSIZE `ls -l F.map | awk '{printf "%d", $5}'`
cat mapman.in | $mapman -b mapsize $MAPSIZE > mappeek.log

# now examine the expected results
cat mappeek.log |\
awk '/PEEK-A-BOO/{++n;\
  x=(n-1)*0.1;\
  pi=4*atan2(1,1);\
  print x,$NF,0.002*cos(x*4*pi/10)}' |\
cat >! plotme.txt

# check deviation from expectation
set rmsd = `awk '{++n;sum+=($2-$3)**2} END{if(n)print sqrt(sum/n)}' plotme.txt`
set rmsa = `awk '{++n;sum+=($2+$3)**2} END{if(n)print sqrt(sum/n)}' plotme.txt`

# now test that we actually got a result, that it wasn't exactly zero, and it was smaller than 1% error
set test = `echo $rmsd $rmsa | awk '{print ( NF==2 && ( $1 > 0) && ( $1 < $2/100 ) )}'`

# short and sweet result
if("$test" == "1") then
    echo "$mapman passed"
else
    echo "ERROR: $mapman failed the test"
    exit 9
endif

