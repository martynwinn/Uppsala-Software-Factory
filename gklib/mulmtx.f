c
c ===========================================================================
c
	subroutine mulmtx (a, b, c, dim1, dim2, dim3)
c
c ---	MATrix_MULtiply
c	Multiply 2 matrices to gether and store in a 3rd matrix.
c	The matrices are C(dim1,dim3) = a(dim1,dim2)*b(dim2,dim3)
c	They must be different storeage locations.
c ---	Written by 
	implicit none
	integer dim1, dim2, dim3
	real a(dim1,dim2), b(dim2, dim3), c(dim1,dim3)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	integer i, j, k
	do 100 i=1,dim1
	do 100 j=1,dim3
	  c(i,j) = 0.
	  do 100 k=1,dim2
100	    c(i,j) = c(i,j)+ a(i,k)* b(k,j)
	return
	end
