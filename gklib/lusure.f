c
c ========================================================================
c
      logical function lusure(prompt)
c
c ... ask the user whether or he/or she is sure (of whatever)
c
c ... LUSURE = "Are you sure" in Japanese
c
c === G J Kleywegt @ 910926/920512
c
      implicit none
c
      integer mylen,length
c
      character answer*1, prompt*(*), myprom*80
c
code ...
c
      answer = 'N'
      mylen  = length(prompt)
c
      if (mylen.eq.0) then
        myprom = ' Are you sure (Y/N) ?'
      else
        mylen  = min (mylen,70)
        myprom = ' '//prompt(1:mylen)//' (Y/N) ?'
      endif
c
      call textin (myprom,answer)
      lusure = (answer.eq.'Y' .or. answer.eq.'y')
c
      return
      end
