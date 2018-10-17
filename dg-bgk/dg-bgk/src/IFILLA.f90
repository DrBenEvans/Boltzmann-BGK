      SUBROUTINE IFILLA(MA,NA,LA,KA,val) 
! 
      INTEGER MA(NA,LA,KA),val 
! 
	  DO 1002 K=1,KA 
      DO 1000 J=1,LA 
      DO 1001 I=1,NA 
      MA(I,J,K)=val 
 1001 CONTINUE 
 1000 CONTINUE 
 1002 CONTINUE 
! 
       RETURN 
       END 
