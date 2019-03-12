c
c
c
	SUBROUTINE getr1 (LUN, CELL, NSYM, IERROR)
c
        implicit none
C
	REAL*4    CELL(6)     
	INTEGER*2 IP	    ,NB
	INTEGER*4 NSYM	  ,IBLID	 ,NWR	   ,LUN
	INTEGER*4 IERROR	,I
c
code ...
c
C
C   Read the type 1 record
C
	READ (LUN,ERR=900) IBLID,IP,NB,NWR,CELL,NSYM
	IF(IBLID.NE.1) GOTO 910
C
C  Convert COS(ANGLE) to ANGLE
C
	DO 100 I  = 4,6
	IF(ABS(CELL(I)).GT.1.000) GO TO 920
	CELL(I) = (180./3.141593)*ACOS(CELL(I))
  100 CONTINUE
	RETURN
C
C   Errors
C
  900 IERROR	 = 1
	WRITE (6,1000)
	RETURN
C
  910 IERROR	 = 2
	WRITE(6,1100) IBLID
	RETURN
C
  920 IERROR	 = 4
	WRITE (6,1200)  (CELL(I),I=4,6)
	RETURN
C
 1000 FORMAT(' MAP> Error reading R1 record')
 1100 FORMAT(' MAP> Record R1 block-id = ',I5,' (should be 1) ')
 1200 FORMAT(' MAP> Angle cosine(s) out of range:', 3F15.6)
c
	END
