c
c
c
      subroutine dosymb (nopt,optpar,ldone)
c
      integer maxsym,maxval,maxnam
      parameter (maxsym=100, maxval=256, maxnam=20)
c
      integer nopt,nsym,i,j,k,ll,length
      integer symlen(maxsym)
c
      logical ldone
c
      character optpar(nopt)*(*)
      character symnam(maxsym)*(maxnam),symval(maxsym)*(maxval)
      character line*(maxval)
c
      common /cbsymc/ symnam,symval
      common /cbsymi/ symlen,nsym
c
      data nsym /0/
c
code ...
c
      ldone = .false.
c
      if (nopt .lt. 1) return
c
c ... substitute symbols (but not second parameter if & command)
c
      if (nsym .ge. 1) then
        k = 2
        if (optpar(1) .eq. '&') k = 3
        do i=k,nopt
          if (optpar(i)(1:1) .eq. '$') then
            line = optpar(i)(2:)
            call upcase (line)
            call remspa (line)
            ll = length(line)
            do j=1,nsym
              if (ll .eq. symlen(j)) then
                if (line(1:ll) .eq. symnam(j)(1:ll)) then
                  optpar (i) = symval (j)
                  goto 10
                end if
              end if
            end do
            call errcon ('Unrecognised symbol')
            call textut (' Name :',optpar(i))
   10       continue
          end if
        end do
      end if
c
      if (optpar(1) .eq. '&') then
        ldone = .true.
c
c ... list symbols ?
c
        if (nopt .ge. 2) then
          if (optpar(2)(1:1) .eq. '?') then
            call ivalut (' Nr of defined symbols :',1,nsym)
            if (nsym .gt. 0) then
              do i=1,nsym
                line = ' Symbol '//
     +                 symnam(i)(1:symlen(i))//' :'
                call textut (line,symval(i))
              end do
            end if
            return
          end if
        end if
c
        if (nopt .lt. 2 .or. length(optpar(2)) .lt. 1) then
          call errcon ('No argument(s) for & command')
          return
        end if
c
c ... existing symbol ?
c
        if (nsym .ge. 1) then
          line = optpar(2)
          call upcase (line)
          call remspa (line)
          ll = length (line)
          do j=1,nsym
            if (ll .eq. symlen(j)) then
              if (line(1:ll) .eq. symnam(j)(1:ll)) then
                k = j
                goto 1000
              end if
            end if
          end do
        end if
c
c ... if here, new symbol
c
        if (nsym .ge. maxsym) then
          call errcon ('Too many symbol definitions')
          call ivalut (' Max :',1,maxsym)
          return
        end if
c
        nsym = nsym + 1
        symnam (nsym) = optpar(2)
        call upcase (symnam(nsym))
        call remspa (symnam(nsym))
        symlen (nsym) = length(symnam(nsym))
        symval (nsym) = ' '
        k = nsym
c
c ... get (new) value
c
 1000   continue
        if (nopt .ge. 3) then
          symval (k) = optpar(3)
        else
          line = ' Value for symbol '//
     +           symnam(k)(1:symlen(k))//' ?'
          call textin (line,symval(k))
        end if
        line = ' Symbol '//symnam(k)(1:symlen(k))//' :'
        call textut (line,symval(k))
      end if
c
      return
      end
