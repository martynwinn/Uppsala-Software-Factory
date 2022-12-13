#!/bin/csh -f

# Build script for gklib USF library. MRH Feb 2009

set echo 
set int_size = "32"
set OS = "OSX"

if ($#argv > 0) then
    set OS = $argv[1]
endif

if ($#argv > 1) then
    set int_size = $argv[2]
endif

set OS = `echo $OS | tr '[a-z]' '[A-Z]'`

set sys=`uname -s`
if ("$OS" == "OSX" && $sys == Darwin) then
echo "Building OS X gklib library"
else
goto LINUX
endif


# OS X build  ###################################################################### OS X Build

set lib=osx_kleylib
##set subdir=os

set fort="gfortran -D$OS -m$int_size -O -u -check_bounds -std=legacy -ffixed-line-length-132"

set ccom="gcc -c -g -D$OS -m$int_size"

# backup existing library
if (-e $lib) mv $lib {$lib}.ckp

#if (! -d $subdir) then
#  echo WARNING ... creating subdirectory $subdir
#  mkdir $subdir
#endif

##cd $subdir

echo Compiling ...

$fort -c *.f |& tee compilation.log | grep -i error
grep -i warn compilation.log
grep -i info compilation.log

if ("$OS" == "OSX") then
  echo Compiling fmalloc.c ...
  $ccom fmalloc.c
endif

ar rsv $lib *.o
ranlib $lib

cd ..

exit $status

## LINUX build

LINUX:

set sys=`uname -s`
if ($sys != Linux) then
  echo ERROR ... this machine is not a Linux box nor OS X
  exit -1
endif
echo "Building LINUX gklib library"

set lib=linux_kleylib

set fort="gfortran -DLINUX -O -u -g -m"$int_size" -std=legacy -ffixed-line-length-132"


#set ccom='cc -c -g -DOSX'
set ccom="gcc -c -g -DLINUX -m"$int_size

# backup existing library
if (-e $lib) mv $lib {$lib}.ckp

echo Compiling gklib ...

$fort -c *.f |& tee compilation.log | grep -i error
grep -i warn compilation.log
grep -i info compilation.log

echo Compiling fmalloc.c ...
$ccom fmalloc.c

ar rsv $lib *.o
ranlib $lib

cd ..

exit 0
