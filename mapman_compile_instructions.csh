#! /bin/tcsh -f
#
#   get mapman working using git repositories
#

# needed after default RedHat 8.4 install (run as root)
if ( `whoami` == "root" ) then
yum install tcsh 
yum install make git m4 patch
yum install gcc gcc-gfortran libgfortran
# hopefully your sysadmin already did this
endif

# regular user start here

# re-compile ccp4 libraries
git clone https://github.com/dials/ccp4io
cd ccp4io/libccp4/
chmod a+x configure 
./configure
make
cd ../..


# now get the version of USF being maintained by Martyn Winn
git clone https://github.com/martynwinn/Uppsala-Software-Factory
cd Uppsala-Software-Factory

# a few corrections...
mkdir bin
mkdir ccp4libs_latest_m64_linux

# copy in CCP4 libs
cp ../ccp4io/libccp4/fortran/.libs/libccp4f.a ccp4libs_latest_m64_linux/
cp ../ccp4io/libccp4/ccp4/.libs/libccp4c.a ccp4libs_latest_m64_linux/

# do the build
./make_all.csh mapman -64

# test it
wget https://raw.githubusercontent.com/fraser-lab/holton_scripts/master/map_bender/mapman_regression_test.csh
chmod a+x mapman_regression_test.csh
./mapman_regression_test.csh ./bin/mapman


