       SUBROUTINE GETFOR(LIFT,DRAG,NBOUN,NPOIN,NBNOI,PS,COORD,BSIDO)
! 
! *** THIS SUBROUTINE COMPUTES THE LIFT AND DRAG COEFFICIENT OF THE BODY
! BEING STUDIED 
! 
      IMPLICIT NONE 
      INCLUDE 'mpif.h' 
! 
       INTEGER IB,IP1,IP2
       INTEGER NBNOI,NBOUN,NPOIN
       INTEGER ROTMAT(2,2)
       INTEGER BSIDO(NBNOI,NBOUN)
       
       REAL LIFT,DRAG,X1,X2,Y1,Y2,VEC(2),LEN,PERP(2)
       REAL PS_AV,MAG
       REAL PS(NPOIN),COORD(2,NPOIN)
!
! 90deg anticlckwise ROTATION MATRIX
!
     ROTMAT(1,1)=0
     ROTMAT(1,2)=-1
     ROTMAT(2,1)=1
     ROTMAT(2,2)=0 
!
! INITIALISE LIFT AND DRAG 
!
      LIFT = 0.0
      DRAG = 0.0
      DO 1000 IB=1,NBOUN
         IF(BSIDO(4,IB).EQ.4)THEN
            IP1=BSIDO(1,IB)
            IP2=BSIDO(2,IB)
            X1=COORD(1,IP1)
            Y1=COORD(2,IP1)
            X2=COORD(1,IP2)
            Y2=COORD(2,IP2)
            LEN=SQRT((X2-X1)**2+(Y2-Y1)**2)
            VEC(1)=X2-X1
            VEC(2)=Y2-Y1  
            PERP(1)=ROTMAT(1,1)*VEC(1)+ROTMAT(1,2)*VEC(2)
            PERP(2)=ROTMAT(2,1)*VEC(1)+ROTMAT(2,2)*VEC(2)
            MAG=SQRT(PERP(1)*PERP(1)+PERP(2)*PERP(2))
            PERP(1)=PERP(1)/MAG
            PERP(2)=PERP(2)/MAG
            PS_AV=0.5*(PS(IP1)+PS(IP2))
            LIFT=LIFT-PS_AV*LEN*PERP(2) 
            DRAG=DRAG-PS_AV*LEN*PERP(1)           
         ENDIF 
 1000 CONTINUE 
! 
      RETURN 
      END 
                 
             
