Uppsala Software Factory. March 2019
---------------------------------------------

These are the famous Uppsala Software Factory programs, rescued into GitHub.

The following is the original README, now out-of-date!!


Introduction
----------

The USF program suite is no longer maintained or supported, but the source code is being made available along with a build system for Linux and OS X on Intel Macs.
Most users should be able to simply download the relevant executables and run these without any need for compilation, but if you use Linux or wish to develop the code, you will need the full distribution kit.

Download 
--------

The distribution kit can be downloaded from : http://xray.bmc.uu.se/markh/usf

Organisation of the distribution kit
----------------------------

35 programs are available, and each has its own directory with the same name as the program, under the parent directory 'usf'. In the top-level 'usf' directory there is a shell script called 'make_all.csh', and big surprise, if you execute that file it will will build all the programs, leaving executables in both the programs' own subdirectories and also in the collective directory 'usf/bin'. 
Building under OS X

The package has been built under Mac OS X v10.6 using the gfortran compiler.
In order to build yourself you will need to install XCODE which can be downloaded from http://developer.apple.com/technology/xcode.html, and gfortran which can be found here http://www.macresearch.org/gfortran-leopard. You will also find a gfortran dmg file in the distribution. 
Follow the instructions that come with these installation packages. 
To be able to download these packages, you will have to register as an Apple Developer.

USF Program list
-------------

aconio 
ave 
cello 
coma 
comap 
comdem 
crave 
dataman 
essens 
findncs 
imp 
lsqman 
mama
mapfix 
mapman 
mappage
maskit 
mave
moleman 
moleman2 
ncs6d
o2d
odbman 
oops 
oops2 
prof
seaman 
site2rt 
sod 
solex
spancsi
ssencs 
voidoo
xpand 
xplo2d 


Programs in the Uppsala Software Factory were written by Gerard Kleywegt, 
while the archive is currently administered by Mark Harris.
