All references to the standard are to IEEE754 2008 

all functions are standard compliant with the exception of signaling errors

not implemented parts of standard:  
	different roundings                              ; the only one available is TOWARD_ZERO  
	remainder(x, y)                                  ;  
	convert to/from decimal character sequence       ;  
	anything related to error handling and signaling ;  
	fused multiply add                               ; I am too tired to do it in the nearest future  
  
no reason to implement   
	radix                                            ; everything is in base 2   
	is canonical                                     ; every number accessible to user thru functions other than pack/unpack is canonical  
	minmag, maxmag, totalordermag                    ; not needed since user can use usual functions with abs() on arguments  
	copysign(x)                                      ; x = x   


internal representation is binary float with extended precision  
this way there is no need to handle subnormals in different way than normals, pack() and unpack() take care of that   

