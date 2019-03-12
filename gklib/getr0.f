c
c
c
	SUBROUTINE getr0 ( LUN, TITLE, IERROR)
c
        implicit none
c
	INTEGER*4  IBLID	 ,NWR	   ,IFT	   ,TITLE(20) 
	INTEGER*4  LUN	   ,IERROR	,I
	INTEGER*2  IP	    ,NB
c
code ...
c
C
C   Blank title
C
	DO 100 I  = 1,20
	TITLE(I)    = 0
  100 CONTINUE
C
C   Read the type zero record
C
	READ (LUN,ERR=900) IBLID,IP,NB,NWR,
     1    (TITLE(I),I=1,4),IFT,(TITLE(I),I=6,7)
	IF (IBLID.NE.0)    GOTO 910
	IF (IFT.NE.3)	GOTO 920
	RETURN
C 
C   Errors
C
  900 IERROR	  = 1
	WRITE (6,1000)
	RETURN
C
  910 IERROR	  = 2
	WRITE(6,1100) IBLID
	RETURN
C
  920 IERROR	  = 3
	WRITE(6,1200) IFT
	RETURN
C
 1000 FORMAT (' MAP> Error reading R0 record')
 1100 FORMAT (' MAP> Record R0 block-id  = ',I5,' (should be 0)')
 1200 FORMAT (' MAP> Record R0 file type = ',I5,' (should be 3)')
c
	END
