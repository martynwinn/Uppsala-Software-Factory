c
c
c
      subroutine maskqq (lun,mask,ext1,ext2,ext3,space,how,ierr)
c
c --- Read a mask
c
      implicit none
c
      integer lun, ext1, ext2, ext3, space
      integer mask (ext1, ext2, ext3) 
      integer i, j, k, ierr
c
      character how*(*)
c
code ...
c
      ierr = -1
c
      if (ext1*ext2*ext3 .gt. space) then
        call errcon (' Space exceeded in MASKQQ')
        call jvalut (' Required :',1,(ext1*ext2*ext3))
        call jvalut (' Maximum  :',1,space)
        return
      end if
c
      if (how(1:3) .eq. 'NEW') then
        read (lun, 40,err=9000,end=9000)
     +    (((mask(i,j,k), i=1,ext1),j=1,ext2),k=1,ext3)
      else
        read (lun, 20,err=9000,end=9000)
     +    (((mask(i,j,k), i=1,ext1),j=1,ext2),k=1,ext3)
      end if
c
      call prompt (' Mask read OK')
      ierr = 0
c
      return
c
 9000 continue
      return
c
20    format (40i2)
40    format (80i1)
c
      end
