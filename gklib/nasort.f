c
c
c
      subroutine nasort (na,nam,ind)
c
c ... NASORT - simple bubble-sort by atom name (hydrogens to the bottom)
c
      implicit none
c
      integer na,ind(*),i,j,k
c
      logical swaped,lhydro,l1,l2
c
      character*4 nam(*),dumnam*4
c
code ...
c
c ... replace G(amma) by C etc.
c
      do i=1,na
        if (nam(i)(3:3) .eq. 'G') nam(i)(3:3) = 'C'
        if (nam(i)(3:3) .eq. 'Z') nam(i)(3:3) = 'F'
        if (nam(i)(3:3) .eq. 'H') nam(i)(3:3) = 'G'
      end do
c
      i = 1
c
   10 continue
      swaped = .false.
      l1 = lhydro(nam(i))
c
      do j=i+1,na
c
        l2 = lhydro(nam(j))
        if (l1 .and. l2 .or. (.not. l1 .and. .not. l2) ) then
          if (nam(i)(3:4) .gt. nam(j)(3:4)) goto 20
          if (nam(i)(3:4) .lt. nam(j)(3:4)) goto 30
          if (nam(i)(1:2) .gt. nam(j)(1:2)) goto 20
        else if (.not. l2 .and. l1) then
           goto 20
        end if
        goto 30
c
   20   continue
        dumnam = nam(i)
        nam(i) = nam(j)
        nam(j) = dumnam
        k = ind(i)
        ind(i) = ind(j)
        ind(j) = k
        goto 10
c
   30   continue
c
      end do
c
      i = i + 1
      if (i .lt. na) goto 10
c
c ... replace C by G(amma) etc.
c
      do i=1,na
        if (nam(i)(3:3) .eq. 'G') nam(i)(3:3) = 'H'
        if (nam(i)(3:3) .eq. 'F') nam(i)(3:3) = 'Z'
        if (nam(i)(3:3) .eq. 'C') nam(i)(3:3) = 'G'
      end do
c
      return
      end
