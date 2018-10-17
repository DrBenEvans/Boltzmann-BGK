      SUBROUTINE IDENTA(L,M,N,A,B) 
      REAL A(L,M,N),B(L,M,N) 
      DO 3000 K=1,N 
      DO 2000 J=1,M 
      DO 1000 I=1,L 
      A(I,J,K)=B(I,J,K) 
 1000 CONTINUE 
 2000 CONTINUE 
 3000 CONTINUE 
      RETURN 
      END 
