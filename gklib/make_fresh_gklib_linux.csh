#!/bin/csh -f

# make a new subroutine library for linux

set echo

set sys=`uname -s`
if ($sys != Linux) then
  echo ERROR ... this is not a Linux box
  exit -1
endif

set lib=linux_kleylib

set fort='gfortran -DLINUX -O -u -g -check_bounds -fbacktrace -fbounds-check -m32 -std=legacy -ffixed-line-length-132'


#set ccom='cc -c -g -DOSX'
set ccom='gcc -c -g -DLINUX -m32 '

# backup existing library
if (-e $lib) mv $lib {$lib}.ckp

echo Compiling ...

$fort -c *.f |& tee compilation.log | grep -i error
grep -i warn compilation.log
grep -i info compilation.log

echo Compiling fmalloc.c ...
$ccom fmalloc.c

ar rsv $lib *.o
ranlib $lib

cd ..
##touch make_fresh_osxlib

# NOTE: if the creation of the library fails try this:
#
# cd osx
# ar rsv /Volumes/GERARD/$lib *.o
# ranlib /Volumes/GERARD/$lib
# ls /Volumes/GERARD/$lib
# cd ..
# cp /Volumes/GERARD/$lib .

exit 0

