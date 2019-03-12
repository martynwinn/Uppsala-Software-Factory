c
c
c
      subroutine priref (nr,str1,str2,str3,str4,str5,str6)
c
      implicit none
c
      integer nr,length
c
      character str1*(*),str2*(*),str3*(*)
      character str4*(*),str5*(*),str6*(*)
c
code ...
c
      if (nr .le. 0) then
        call prompt (' Reference(s) for this program:')
        nr = 0
      end if
c
      if (length(str1) .le. 0) return
c
      nr = nr + 1
c
      write (*,6000) nr,str1(1:length(str1))
      if (length(str2).gt.0) write(*,6010) str2(1:length(str2))
      if (length(str3).gt.0) write(*,6010) str3(1:length(str3))
      if (length(str4).gt.0) write(*,6010) str4(1:length(str4))
      if (length(str5).gt.0) write(*,6010) str5(1:length(str5))
      if (length(str6).gt.0) write(*,6010) str6(1:length(str6))
c
 6000 format (/' * ',i2,' * ',a)
 6010 format (8x,a)
c
      return
      end
