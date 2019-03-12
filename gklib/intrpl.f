c
c
c
      subroutine intrpl (ed, nxa, nya, nza, x, densty, errcod)
c
c ---	Wayne's interpolation
c
      implicit none
c
      integer nxa, nya, nza, errcod
      real ed (nxa, nya, nza), x(3), densty
c
      real rx(4,4), rxy(4), spline,dx,dy,dz,zk0,zk1
      integer jx,jy,jz,ifx,ify,ifz,ii,jj,i,j
c
code ...
c
C
C    DETERMINE OFFSETS
c
      JX=int(X(1))
      JY=int(X(2))
      JZ=int(X(3))
c
      errcod = 1
ccc
c      if (jx.ge.(nxa).or.jy.ge.(nya).or.jz.ge.(nza)) then
c        write (*,*) 'Border ',jx,jy,jz
c      end if
ccc
	if(jx .le. 0) return
	if(jy .le. 0) return
	if(jz .le. 0) return
c
	if (jx .ge. nxa) return
	if (jy .ge. nya) return
	if (jz .ge. nza) return
c
	if (abs(x(1)-float(jx)) .lt. 0.001) then
          if (abs(x(2)-float(jy)) .lt. 0.001) then
            if (abs(x(3)-float(jz)) .lt. 0.001) then
	      densty = ed(jx, jy, jz)
	      errcod = 0
	      return
            end if
          end if
	end if
c
	errcod = 0
c
      DX=X(1)-FLOAT(JX)
      DY=X(2)-FLOAT(JY)
      DZ=X(3)-FLOAT(JZ)
C
C    IS THIS A BORDER OR INTERIOR POINT
c
      IFX=2
      IFY=2
      IFZ=2
      IF(JX.EQ.1) IFX=1
      IF(JX.EQ.NXA-1) IFX=3
      IF(JY.EQ.1) IFY=1
      IF(JY.EQ.NYA-1) IFY=3
      IF(JZ.EQ.1) IFZ=1
      IF(JZ.EQ.NZA-1) IFZ=3
C
C    X INTERPOLATION
c
      II=JY-IFY
      JJ=JZ-IFZ
      DO 120 I=1,4
      DO 120 J=1,4
      GO TO(112,115,114),IFX
  112 ZK0 = 2.*ED(JX,  II+I,JJ+J) -
     +      5.*ED(JX+1,II+I,JJ+J) +
     +      4.*ED(JX+2,II+I,JJ+J) -
     +         ED(JX+3,II+I,JJ+J)
      ZK1 =    ED(JX,  II+I,JJ+J) -
     +      2.*ED(JX+1,II+I,JJ+J) +
     +         ED(JX+2,II+I,JJ+J)
      GO TO 120
  114 ZK0 =    ED(JX-1,II+I,JJ+J) -
     +      2.*ED(JX,  II+I,JJ+J) +
     +         ED(JX+1,II+I,JJ+J)
      ZK1 = 2.*ED(JX+1,II+I,JJ+J) -
     +      5.*ED(JX,  II+I,JJ+J) +
     +      4.*ED(JX-1,II+I,JJ+J) -
     +         ED(JX-2,II+I,JJ+J)
      GO TO 120
  115 ZK0 =    ED(JX-1,II+I,JJ+J) -
     +      2.*ED(JX,  II+I,JJ+J) +
     +         ED(JX+1,II+I,JJ+J)
      ZK1 =    ED(JX,  II+I,JJ+J) -
     +      2.*ED(JX+1,II+I,JJ+J) +
     +         ED(JX+2,II+I,JJ+J)
  120 RX(I,J) = SPLINE (DX,ED(JX,II+I,JJ+J),
     +                  ED(JX+1,II+I,JJ+J),ZK0,ZK1)
C
C    Y INTERPOLATION
c
      DO 130 J=1,4
      GO TO(122,125,124),IFY
  122 ZK0 = 2.*RX(1,J) - 5.*RX(2,J) + 4.*RX(3,J) - RX(4,J)
      ZK1 =    RX(1,J) - 2.*RX(2,J) +    RX(3,J)
      GO TO 130
  124 ZK0 =    RX(2,J) - 2.*RX(3,J) +    RX(4,J)
      ZK1 = 2.*RX(4,J) - 5.*RX(3,J) + 4.*RX(2,J) - RX(1,J)
      GO TO 130
  125 ZK0 =    RX(1,J) - 2.*RX(2,J) +    RX(3,J)
      ZK1 =    RX(2,J) - 2.*RX(3,J) +    RX(4,J)
  130 RXY(J) = SPLINE (DY,RX(IFY,J),RX(IFY+1,J),ZK0,ZK1)
C
C    Z INTERPOLATION
c
      GO TO(132,135,134),IFZ
  132 ZK0 = 2.*RXY(1) - 5.*RXY(2) + 4.*RXY(3) - RXY(4)
      ZK1 =    RXY(1) - 2.*RXY(2) +    RXY(3)
      GO TO 140
  134 ZK0 =    RXY(2) - 2.*RXY(3) +    RXY(4)
      ZK1 = 2.*RXY(4) - 5.*RXY(3) + 4.*RXY(2) - RXY(1)
      GO TO 140
  135 ZK0 =    RXY(1) - 2.*RXY(2) +    RXY(3)
      ZK1 =    RXY(2) - 2.*RXY(3) +    RXY(4)
  140 DENSTY = SPLINE (DZ,RXY(IFZ),RXY(IFZ+1),ZK0,ZK1)
c
      RETURN
      END
