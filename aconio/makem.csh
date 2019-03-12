#!/bin/csh -f

set echo

alias f77_alpha_comp 'f77 -Olimit 2000 -C -O -u -check bounds -align dcommons -warn argument_checking -c'
alias f77_alpha_link 'f77 -non_shared -om -O -math_library accurate'
alias f77_6d_comp  'f77 -Olimit 2000 -C -O -u -check_bounds -w0 -c'
alias f77_6d_link  'f77 -O'
alias f77_lx_comp  'f77 -C -O -u -check_bounds -c'
alias f77_lx_link  'f77 -O'
alias f77_osx_comp  'g77 -C -O -u -check_bounds -c'
alias f77_osx_link  'g77 -O'

set sys=`uname -s`

switch ($sys)

# dec alpha/osf1
  case OSF1:
f77_alpha_comp aconio.f
f77_alpha_link -o DEC_ACONIO aconio.o ../gklib/alpha_kleylib
strip DEC_ACONIO
    breaksw

# silicon graphics
  case IRIX64:
f77_6d_comp aconio.f
f77_6d_link -o SGI_ACONIO aconio.o ../gklib/6d_kleylib
strip SGI_ACONIO
../man2html/MAN2HTML << EOF
aconio.txt
ACONIO
aconio_man.html
EOF
    breaksw

# linux
  case Linux:
f77_lx_comp aconio.f
f77_lx_link -o LX_ACONIO aconio.o ../gklib/lx_kleylib
strip LX_ACONIO
    breaksw

# mac osx
  case Darwin:
f77_osx_comp aconio.f
f77_osx_link -o OSX_ACONIO aconio.o ../gklib/osx_kleylib
strip OSX_ACONIO
    breaksw

  default:
    echo ERROR ... Not implemented on this type of machine $sys
    exit 1
endsw

unset echo

exit 0

