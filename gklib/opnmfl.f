c
c
c
      subroutine opnmfl(fname,lform,iunit,ierr)
C     ====================================
C Find out if input file FNAME is formatted or not
C Leave file closed
c
      implicit none
c
      character*(*) fname
C LFORM is true if input file in formatted
      logical lform
c
      integer index, kdummy, kfail, iunit, ierr
      character*80 line
c
code ...
c
      lform = .false.
      kdummy = 0
      kfail = 1
      call ccpdpn(iunit,fname,'READONLY','F',kdummy,kfail)
      if (kfail .eq. -1) then
        ierr = -1
        return
      end if
c
C Read 1st line as is formatted, should contain string MAPEXCHANGE if so
      read (iunit,6001,err=100) line
 6001 format(A)
c
      if (index(line,'MAPEXCHANGE') .gt. 0) then
         lform = .true.
      else
         lform = .false.
      endif
c
C Close file
 100  close (unit=iunit)
      ierr = 0
c
      return 
      end
