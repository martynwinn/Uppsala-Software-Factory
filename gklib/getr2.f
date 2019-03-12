c
c
c
	SUBROUTINE getr2 (LUN, IPERM, IGRID, MIFD, IERROR)
c
        implicit none
C
	INTEGER*4 IPERM(3)    ,IGRID(3)    ,MIFD(9)     ,IBLID
	INTEGER*4 NWR	   ,IPATT	 ,LUN	   ,IERROR
	INTEGER*2 IP	 ,NB
c
code ...
c
C
C   Read the type 2 record
C
	READ (LUN,ERR=900) IBLID,IP,NB,NWR,IPATT,IPERM,IGRID,MIFD
	IF(IBLID.NE.2) GO TO 910
	IF (IPATT.LT.0 .OR. IPATT.GT.1) GO TO 920
	RETURN
C
C   Errors
C
  900 IERROR	 = 1
	WRITE(6,1000)
	RETURN
C
  910 IERROR	 = 2
	WRITE (6,1100) IBLID
	RETURN
C
  920 IERROR	 = 5
	WRITE (6,1200)  IPATT
	RETURN
C
 1000 FORMAT (' MAP> Error reading R2 record')
 1100 FORMAT (' MAP> Record R2 block-id = ',I5,' (should be 2)')
 1200 FORMAT (' MAP> Map IPATT type = ',I5,' (should be 0 or 1)')
c
	END
