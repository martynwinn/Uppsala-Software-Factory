c
c ===========================================================================
c
	subroutine vecrtv (xa, xb, ct, r, t)                                   
c ---	VECtor_Rotate_Translate_Vector                                    
c ---	Calculate transformed coords.                                     
c ----	XA and XB can be same array.                                     
c 	xa = untransformed coords                                            
c	xb = the transformed coords                                           
c	ct = number points                                                    
c	r = rotation matrix, t=translation vector                                  
	implicit none                                                          
	integer ct                                                             
	real xa(3,ct), xb(3,ct), r(3,3), t(3)                                  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	integer i,j                                                            
	real x(3)                                                              
	do 100 i=1,ct                                                          
	  do 110 j=1,3                                                         
110	    x(j) = xa(1,i)*r(j,1)+xa(2,i)*r(j,2)+xa(3,i)*r(j,3)+t(j)        
	  do 100 j=1,3                                                         
100	  xb(j,i) = x(j)                                                    
	return                                                                 
	end                                                                    
