c
c
c
      subroutine xhp_stamp (isiz,x,y,prog)
c
      implicit none
c
      integer isiz,length,ll
c
      real x,y
c
      character stamp*24,prog*(*),mytext*120,user*40
c
code ...
c
      call xhp_txtsiz (isiz)
c
      call xhp_moveto (x,y)
c
      call gkdate (stamp)
c
      call gkuser (user)
      ll = length (user)
      if (ll .lt. 1) user = 'an UNKNOWN user'
c
      write (mytext,*) 'Produced by ',prog(1:length(prog)),
     +  ' at ',stamp(1:length(stamp)),' for ',
     +  user(1:length(user))
c
      call pretty (mytext)
      ll = length(mytext)
      call xhp_text (mytext,ll)
c
      return
      end
