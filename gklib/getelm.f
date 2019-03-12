c
c
c
      subroutine getelm (atmnam,el,i,mass)
c
c ... GETELM - try to deduce the element nr and get the mass
c     ATMNAM *4 = PDB atom name
c     EL     *2 = PDB element type
c     I         = element number or < 1 if not found
c     MASS      = mass of element
c
      implicit none
c
      real mass,radius
c
      integer i
c
      logical lhydro
c
      character atmnam*(*),el*(*),try*6,line*80
c
code ...
c
c ... (1) is it hydrogen ?
c
      if (lhydro(atmnam)) then
        try = ' H'
        call elinfo (try,line,i,mass,radius,.false.)
        if (i .gt. 0) return
      end if
c
c ... (2) try intact atom name
c
      try = atmnam(1:2)
      call elinfo (try,line,i,mass,radius,.false.)
      if (i .gt. 0) return
c
c ... (3) try PDB element symbol
c
      try = el(1:2)
      call elinfo (try,line,i,mass,radius,.false.)
      if (i .gt. 0) return
c
c ... (4) try space + second char of atom name
c
      try = ' '//atmnam(2:2)
      call elinfo (try,line,i,mass,radius,.false.)
      if (i .gt. 0) return
c
c ... (5) try space + first char of atom name
c
      try = ' '//atmnam(1:1)
      call elinfo (try,line,i,mass,radius,.false.)
      if (i .gt. 0) return
c
c ... alas
c
      i = -999
      mass = 0.0
c
      return
      end
