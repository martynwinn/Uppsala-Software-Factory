c
c ... seq_subs.f - subroutines for dealing with amino-acid sequences
c
c ... G J Kleywegt - 1997-2005
c
c
      subroutine inicmp (cmpmat,maxtyp)
c
c ... INICMP () - initialise BLOSUM-45 matrix; assumes:
c
c      data typ3lc /'ALA','ARG','ASN','ASP','CYS','GLU','GLN',
c     +             'GLY','HIS','ILE','LEU','LYS','MET','PHE',
c     +             'PRO','SER','THR','TRP','TYR','VAL'/
c
c      data typ1lc /'A',  'R',  'N',  'D',  'C',  'E',  'Q',
c     +             'G',  'H',  'I',  'L',  'K',  'M',  'F',
c     +             'P',  'S',  'T',  'W',  'Y',  'V'/
c
      implicit none
c
      integer maxtyp,i,j
c
      real cmpmat(maxtyp,maxtyp)
c
code ...
c
      do i=1,maxtyp
        do j=1,maxtyp
          cmpmat (i,j) = -9999.999
        end do
      end do
c
      if (maxtyp .lt. 20) then
        call errstp (
     +    'INICMP - Need at least 20 types - E-mail Gerard !')
        return
      end if
c
      CMPMAT ( 1, 1) =   5.0
      CMPMAT ( 1, 2) =  -2.0
      CMPMAT ( 1, 3) =  -1.0
      CMPMAT ( 1, 4) =  -2.0
      CMPMAT ( 1, 5) =  -1.0
      CMPMAT ( 1, 6) =  -1.0
      CMPMAT ( 1, 7) =  -1.0
      CMPMAT ( 1, 8) =   0.0
      CMPMAT ( 1, 9) =  -2.0
      CMPMAT ( 1,10) =  -1.0
      CMPMAT ( 1,11) =  -1.0
      CMPMAT ( 1,12) =  -1.0
      CMPMAT ( 1,13) =  -1.0
      CMPMAT ( 1,14) =  -2.0
      CMPMAT ( 1,15) =  -1.0
      CMPMAT ( 1,16) =   1.0
      CMPMAT ( 1,17) =   0.0
      CMPMAT ( 1,18) =  -2.0
      CMPMAT ( 1,19) =  -2.0
      CMPMAT ( 1,20) =   0.0
      CMPMAT ( 2, 2) =   7.0
      CMPMAT ( 2, 3) =   0.0
      CMPMAT ( 2, 4) =  -1.0
      CMPMAT ( 2, 5) =  -3.0
      CMPMAT ( 2, 6) =   1.0
      CMPMAT ( 2, 7) =   0.0
      CMPMAT ( 2, 8) =  -2.0
      CMPMAT ( 2, 9) =   0.0
      CMPMAT ( 2,10) =  -3.0
      CMPMAT ( 2,11) =  -2.0
      CMPMAT ( 2,12) =   3.0
      CMPMAT ( 2,13) =  -1.0
      CMPMAT ( 2,14) =  -2.0
      CMPMAT ( 2,15) =  -2.0
      CMPMAT ( 2,16) =  -1.0
      CMPMAT ( 2,17) =  -1.0
      CMPMAT ( 2,18) =  -2.0
      CMPMAT ( 2,19) =  -1.0
      CMPMAT ( 2,20) =  -2.0
      CMPMAT ( 3, 3) =   6.0
      CMPMAT ( 3, 4) =   2.0
      CMPMAT ( 3, 5) =  -2.0
      CMPMAT ( 3, 6) =   0.0
      CMPMAT ( 3, 7) =   0.0
      CMPMAT ( 3, 8) =   0.0
      CMPMAT ( 3, 9) =   1.0
      CMPMAT ( 3,10) =  -2.0
      CMPMAT ( 3,11) =  -3.0
      CMPMAT ( 3,12) =   0.0
      CMPMAT ( 3,13) =  -2.0
      CMPMAT ( 3,14) =  -2.0
      CMPMAT ( 3,15) =  -2.0
      CMPMAT ( 3,16) =   1.0
      CMPMAT ( 3,17) =   0.0
      CMPMAT ( 3,18) =  -4.0
      CMPMAT ( 3,19) =  -2.0
      CMPMAT ( 3,20) =  -3.0
      CMPMAT ( 4, 4) =   7.0
      CMPMAT ( 4, 5) =  -3.0
      CMPMAT ( 4, 6) =   0.0
      CMPMAT ( 4, 7) =   2.0
      CMPMAT ( 4, 8) =  -1.0
      CMPMAT ( 4, 9) =   0.0
      CMPMAT ( 4,10) =  -4.0
      CMPMAT ( 4,11) =  -3.0
      CMPMAT ( 4,12) =   0.0
      CMPMAT ( 4,13) =  -3.0
      CMPMAT ( 4,14) =  -4.0
      CMPMAT ( 4,15) =  -1.0
      CMPMAT ( 4,16) =   0.0
      CMPMAT ( 4,17) =  -1.0
      CMPMAT ( 4,18) =  -4.0
      CMPMAT ( 4,19) =  -2.0
      CMPMAT ( 4,20) =  -3.0
      CMPMAT ( 5, 5) =  12.0
      CMPMAT ( 5, 6) =  -3.0
      CMPMAT ( 5, 7) =  -3.0
      CMPMAT ( 5, 8) =  -3.0
      CMPMAT ( 5, 9) =  -3.0
      CMPMAT ( 5,10) =  -3.0
      CMPMAT ( 5,11) =  -2.0
      CMPMAT ( 5,12) =  -3.0
      CMPMAT ( 5,13) =  -2.0
      CMPMAT ( 5,14) =  -2.0
      CMPMAT ( 5,15) =  -4.0
      CMPMAT ( 5,16) =  -1.0
      CMPMAT ( 5,17) =  -1.0
      CMPMAT ( 5,18) =  -5.0
      CMPMAT ( 5,19) =  -3.0
      CMPMAT ( 5,20) =  -1.0
      CMPMAT ( 6, 6) =   6.0
      CMPMAT ( 6, 7) =   2.0
      CMPMAT ( 6, 8) =  -2.0
      CMPMAT ( 6, 9) =   1.0
      CMPMAT ( 6,10) =  -2.0
      CMPMAT ( 6,11) =  -2.0
      CMPMAT ( 6,12) =   1.0
      CMPMAT ( 6,13) =   0.0
      CMPMAT ( 6,14) =  -4.0
      CMPMAT ( 6,15) =  -1.0
      CMPMAT ( 6,16) =   0.0
      CMPMAT ( 6,17) =  -1.0
      CMPMAT ( 6,18) =  -2.0
      CMPMAT ( 6,19) =  -1.0
      CMPMAT ( 6,20) =  -3.0
      CMPMAT ( 7, 7) =   6.0
      CMPMAT ( 7, 8) =  -2.0
      CMPMAT ( 7, 9) =   0.0
      CMPMAT ( 7,10) =  -3.0
      CMPMAT ( 7,11) =  -2.0
      CMPMAT ( 7,12) =   1.0
      CMPMAT ( 7,13) =  -2.0
      CMPMAT ( 7,14) =  -3.0
      CMPMAT ( 7,15) =   0.0
      CMPMAT ( 7,16) =   0.0
      CMPMAT ( 7,17) =  -1.0
      CMPMAT ( 7,18) =  -3.0
      CMPMAT ( 7,19) =  -2.0
      CMPMAT ( 7,20) =  -3.0
      CMPMAT ( 8, 8) =   7.0
      CMPMAT ( 8, 9) =  -2.0
      CMPMAT ( 8,10) =  -4.0
      CMPMAT ( 8,11) =  -3.0
      CMPMAT ( 8,12) =  -2.0
      CMPMAT ( 8,13) =  -2.0
      CMPMAT ( 8,14) =  -3.0
      CMPMAT ( 8,15) =  -2.0
      CMPMAT ( 8,16) =   0.0
      CMPMAT ( 8,17) =  -2.0
      CMPMAT ( 8,18) =  -2.0
      CMPMAT ( 8,19) =  -3.0
      CMPMAT ( 8,20) =  -3.0
      CMPMAT ( 9, 9) =  10.0
      CMPMAT ( 9,10) =  -3.0
      CMPMAT ( 9,11) =  -2.0
      CMPMAT ( 9,12) =  -1.0
      CMPMAT ( 9,13) =   0.0
      CMPMAT ( 9,14) =  -2.0
      CMPMAT ( 9,15) =  -2.0
      CMPMAT ( 9,16) =  -1.0
      CMPMAT ( 9,17) =  -2.0
      CMPMAT ( 9,18) =  -3.0
      CMPMAT ( 9,19) =   2.0
      CMPMAT ( 9,20) =  -3.0
      CMPMAT (10,10) =   5.0
      CMPMAT (10,11) =   2.0
      CMPMAT (10,12) =  -3.0
      CMPMAT (10,13) =   2.0
      CMPMAT (10,14) =   0.0
      CMPMAT (10,15) =  -2.0
      CMPMAT (10,16) =  -2.0
      CMPMAT (10,17) =  -1.0
      CMPMAT (10,18) =  -2.0
      CMPMAT (10,19) =   0.0
      CMPMAT (10,20) =   3.0
      CMPMAT (11,11) =   5.0
      CMPMAT (11,12) =  -3.0
      CMPMAT (11,13) =   2.0
      CMPMAT (11,14) =   1.0
      CMPMAT (11,15) =  -3.0
      CMPMAT (11,16) =  -3.0
      CMPMAT (11,17) =  -1.0
      CMPMAT (11,18) =  -2.0
      CMPMAT (11,19) =   0.0
      CMPMAT (11,20) =   1.0
      CMPMAT (12,12) =   5.0
      CMPMAT (12,13) =  -1.0
      CMPMAT (12,14) =  -3.0
      CMPMAT (12,15) =  -1.0
      CMPMAT (12,16) =  -1.0
      CMPMAT (12,17) =  -1.0
      CMPMAT (12,18) =  -2.0
      CMPMAT (12,19) =  -1.0
      CMPMAT (12,20) =  -2.0
      CMPMAT (13,13) =   6.0
      CMPMAT (13,14) =   0.0
      CMPMAT (13,15) =  -2.0
      CMPMAT (13,16) =  -2.0
      CMPMAT (13,17) =  -1.0
      CMPMAT (13,18) =  -2.0
      CMPMAT (13,19) =   0.0
      CMPMAT (13,20) =   1.0
      CMPMAT (14,14) =   8.0
      CMPMAT (14,15) =  -3.0
      CMPMAT (14,16) =  -2.0
      CMPMAT (14,17) =  -1.0
      CMPMAT (14,18) =   1.0
      CMPMAT (14,19) =   3.0
      CMPMAT (14,20) =   0.0
      CMPMAT (15,15) =   9.0
      CMPMAT (15,16) =  -1.0
      CMPMAT (15,17) =  -1.0
      CMPMAT (15,18) =  -3.0
      CMPMAT (15,19) =  -3.0
      CMPMAT (15,20) =  -3.0
      CMPMAT (16,16) =   4.0
      CMPMAT (16,17) =   2.0
      CMPMAT (16,18) =  -4.0
      CMPMAT (16,19) =  -2.0
      CMPMAT (16,20) =  -1.0
      CMPMAT (17,17) =   5.0
      CMPMAT (17,18) =  -3.0
      CMPMAT (17,19) =  -1.0
      CMPMAT (17,20) =   0.0
      CMPMAT (18,18) =  15.0
      CMPMAT (18,19) =   3.0
      CMPMAT (18,20) =  -3.0
      CMPMAT (19,19) =   8.0
      CMPMAT (19,20) =  -1.0
      CMPMAT (20,20) =   5.0
c
      do i=1,maxtyp
        do j=1,maxtyp
          if (cmpmat (i,j) .lt. -1000.0) then
            cmpmat (i,j) = cmpmat (j,i)
          end if
        end do
      end do
c
      call prompt (' *** BLOSUM-45 substitution matrix loaded ***')
c
      return
      end
