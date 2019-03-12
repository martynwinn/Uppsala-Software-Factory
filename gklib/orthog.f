c
c ===========================================================================
c
      subroutine orthog (cell,a,isw)
c
c ... Alwyn's orthogonalisation routine
c
c ---	Orthogonalise_matrix_from_unit_cell
c ---	Used to be called Ortho
c ---	convert fractional coords to cartesian , isw<=0
c ---	convert cartesian coords to fractional , isw>0
c ---	cell is a,b,c,alpha,beta,gamma in degrees
c ---	returns a      11   12   13
c ---	               21   22   23
c ---	               31   32   33
c ---	direction of x-axis is kept
c
      implicit none
c
      real degtor
      parameter (degtor = 6.2831853071796/360.0)
c
      real cell(6),a(3,3)
      integer isw
c
      real cabg(3),sabg(3),cabgs(3),abcs(3),vol,sabgs1
c
      integer i,j
c
code ...
c
      do i=1,3
        cabg(i)=cos(cell(i+3)*degtor)
        sabg(i)=sin(cell(i+3)*degtor)
      end do
c
      cabgs(1)=(cabg(2)*cabg(3)-cabg(1))/(sabg(2)*sabg(3))
      cabgs(2)=(cabg(3)*cabg(1)-cabg(2))/(sabg(3)*sabg(1))
      cabgs(3)=(cabg(1)*cabg(2)-cabg(3))/(sabg(1)*sabg(2))
c
      vol=cell(1)*cell(2)*cell(3)*sqrt(1.0+2.0*cabg(1)*cabg(2)*cabg(3)
     *   -cabg(1)**2-cabg(2)**2-cabg(3)**2)
c
      abcs(1)=cell(2)*cell(3)*sabg(1)/vol
      abcs(2)=cell(1)*cell(3)*sabg(2)/vol
      abcs(3)=cell(1)*cell(2)*sabg(3)/vol
c
      sabgs1=sqrt(1.0-cabgs(1)**2)
c
      do i=1,3
        do j=1,3
          a(i,j)=0.0
        end do
      end do
c
      if(isw .le. 0) then
        a(1,1)=cell(1)
        a(1,2)=cabg(3)*cell(2)
        a(1,3)=cabg(2)*cell(3)
        a(2,2)=sabg(3)*cell(2)
        a(2,3)=-sabg(2)*cabgs(1)*cell(3)
        a(3,3)=sabg(2)*sabgs1*cell(3)
      else
        a(1,1)=1.0/cell(1)
        a(1,2)=-cabg(3)/(sabg(3)*cell(1))
        a(1,3)=-(cabg(3)*sabg(2)*cabgs(1)+cabg(2)*sabg(3))/
     *          (sabg(2)*sabgs1*sabg(3)*cell(1))
        a(2,2)=1.0/(sabg(3)*cell(2))
        a(2,3)=cabgs(1)/(sabgs1*sabg(3)*cell(2))
        a(3,3)=1.0/(sabg(2)*sabgs1*cell(3))
      end if
c
      return
      end
