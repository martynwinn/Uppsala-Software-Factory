#!/bin/tcsh -f

# make_all.csh - build file for USF programs. MRH Jan 2010
# Programname as parameter builds only that program

# Need to set home dir below.

######################################################
# Check param list

set int_size = "32"
set one_prog = ""
set skiplib = ""
set static_flag = ""

set iopt = 0
if ($#argv == 0) goto start

argloop:
    @ iopt++
    if ($iopt> $#argv) goto lastarg # Loop through args

    set alen = `echo $argv[$iopt] | grep -v "-" | wc -m` # Look for arg with no minus sign
    if ($alen > 0) set one_prog = $argv[$iopt]

help:
    if ("$argv[$iopt]" == "help" || "$argv[$iopt]" == "-h" || "$argv[$iopt]" == "-help" || "$argv[$iopt]" == "--" || $#argv == 0) then # -help
	echo "make_all.csh [progname] [-32] [-64] [-skiplib] [-static]"

	exit 0

	else if ("$argv[$iopt]" == "-64") then # -m64
	set int_size = "64"

	else if ("$argv[$iopt]" == "-32") then # -m32
	set int_size = "32"

	else if ("$argv[$iopt]" == "-skiplib") then # -skiplib
	set skiplib = "Yes"

	else if ("$argv[$iopt]" == "-static") then # -static
	set static_flag = "-static"

    endif

goto argloop

lastarg:

# Check OS
######################################################
start:

#set echo
date

set sys=`uname -s`
if ($sys == Darwin) then
  set OS = osx

  #set ROOT_DIR = /Users/markh/Desktop/usf

else if ($sys == Linux) then
  set OS = linux

  #set ROOT_DIR = /home/markh/usf

else 
  echo "Unknown OS : $sys"
   exit 0
endif

set ROOT_DIR = `echo $cwd`

set makefile = Makefile_$OS

set BIN_DIR = $ROOT_DIR/bin
set LIB_DIR = $ROOT_DIR/gklib

echo ... MAKE-ing $OS $int_size bit executables ...

# Go through each prog directory
######################################################

foreach prog  (*)

    if ("$one_prog" != "" && "$one_prog" != "$prog") continue

 if ("$prog" == "gklib" || "$prog" == "bin" || "$prog" == "rave_package") continue;
  if (-d $prog) then
    cd $prog

    # copy include file. Needed by some progs. Later on should change them to look upstairs for it.
    cp $ROOT_DIR/maxdim.incl .

    # Copy makefile template
    cat $ROOT_DIR/Makefile_${OS}_template | sed -e "s/PROGRAMNAME/$prog/g" > Makefile_${OS}
    # See if it has sub_dirs
    if (-e ${prog}_subs.f) then 
	cp Makefile_${OS} /tmp/make_tmp ; cat /tmp/make_tmp  | sed -e 's/##SUBS## //g' > Makefile_${OS}; 
    endif
    # See if it has incls
    if (-e ${prog}.incl) then 
	cp Makefile_${OS} /tmp/make_tmp ; cat /tmp/make_tmp  | sed -e 's/##INCL## //g' > Makefile_${OS}; 
    endif

    # Set 32 or 64 bit
    if ("$int_size" == "64") then 
	cp Makefile_${OS} /tmp/make_tmp ; cat /tmp/make_tmp  | sed -e 's/m32/m64/g' > Makefile_${OS}; 
    endif


    # Set static flag
    if ("$static_flag" != "") then 
	cp Makefile_${OS} /tmp/make_tmp ; cat /tmp/make_tmp  | sed -e 's/##STATIC##//g' > Makefile_${OS}; 
    endif


    # Make fpp links for OS X if not yet made
    if ("$OS" == "osx") then 
	foreach file  (*)
	    if ($file:e == "f") then
	       set file2 = `echo $file | sed -e 's/\.f/.fpp/'`
		##echo "File2 : $file2"
            if (! -e $file2) ln -s $file $file2
	    endif
	end
    endif

    cd ..
  endif
end


# Make gklib
######################################################

set echo

cd $LIB_DIR

if ("$skiplib" == "") make_fresh_gklib.csh $OS $int_size

if ($status < 0) then
    echo "Error - gklib build failed"
    exit -1
endif

# Make others
######################################################

cd ..

foreach prog  (*)
    if ("$one_prog" != "" && "$one_prog" != "$prog") continue

  if (-d $prog && -e $prog/$makefile) then
    echo "Building $prog"
    cd $prog
    rm $prog $prog.o
    cat $makefile
    make -f $makefile
    cp $prog $BIN_DIR
    cd ..
  endif
end

goto exit

######################################################

exit:
# Go for coffee
exit 0
