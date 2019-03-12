c
c=======================================================================
c
	real function distce (x1, x2)                                          
c ---	Calculate distance between 2 points                               
	implicit none                                                          
	real x1(3), x2(3)                                                      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	integer i                                                              
	do 110 i=1,3                                                           
	  if(x1(i) .ne. x2(i))goto 100                                         
110	continue                                                            
	distce = 0.                                                            
	return                                                                 
100	distce =    sqrt((x1(1)-x2(1))**2                                   
     $	                +(x1(2)-x2(2))**2                                
     $	                +(x1(3)-x2(3))**2)                               
	return                                                                 
	end                                                                    
