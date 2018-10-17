      SUBROUTINE IDENTM(M,N,A,B) 
      REAL A(M,N),B(M,N) 
      DO 2000 J=1,N 
      DO 1000 I=1,M 
      A(I,J)=B(I,J) 
 1000 CONTINUE 
 2000 CONTINUE 
      RETURN 
      END 
