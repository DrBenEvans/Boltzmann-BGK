      SUBROUTINE IFILLM(MA,NA,LA,K) 
! 
      INTEGER MA(NA,LA) 
! 
      DO 1000 J=1,LA 
      DO 1001 I=1,NA 
      MA(I,J)=K 
 1001 CONTINUE 
 1000 CONTINUE 
! 
       RETURN 
       END 
