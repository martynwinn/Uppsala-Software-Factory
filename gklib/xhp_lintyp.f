c
      SUBROUTINE xhp_lintyp(PENTYP,LINE)
c
C ---	DEFINE LINE TYPE AS A CERTAIN PEN AND TYPE OF LINE
C
C ---	'PENTYP' MUST HAVE VALUE 1-4
C ---	LINE = 1 ->SOLID LINE
C	       2-6 -> VARIOUS DASHED LINES
C	       7 ->DOT AT POINTS
C	       8 ->DOTTED
C
      include 'xhp_graf.incl'
c
!mrh
      BYTE HPLINE(8)
      character*1 HPLINE_c(8)
c
      INTEGER PENTYP,line
c
      BYTE P(3),PL(2)
      character*1 P_c(3),PL_c(2)
      equivalence(HPLINE, HPLINE_c),(P,P_c),(PL,PL_c)
c
      DATA HPLINE_c/' ','2','3','4','5','6','0','1'/ , P_c/'L','T',';'/
      DATA PL_c/',','1'/
c
      IF(LINE .LE. 0  .OR.  LINE .GT. 8)LINE = 1
      IF(LINE .NE. 1)WRITE(LUN,10)P(1),P(2),HPLINE(LINE),PL,P(3)
      IF(LINE .EQ. 1)WRITE(LUN,10)P
  10  FORMAT(1X,6A1)
C
      IF(PENTYP .LT. 1  .OR.  PENTYP .GT. 4)RETURN
      CALL xhp_PEN(3,PENTYP)
c
      RETURN
      END
