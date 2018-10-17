      SUBROUTINE ADDVEC(A,B,C,M) 
! 
      REAL A(M),B(M),C(M) 
! 
! *** VECTOR A + VECTOR B ==> VECTOR C 
! 
      DO 1000 I=1,M 
      C(I)=A(I)+B(I) 
 1000 CONTINUE 
! 
      RETURN 
      END 
