       SUBROUTINE FLUCON2(NELEM_PP,VNPNT,DISNF_PP,NBNOI,NBOUN_PP,& 
     &                  BSIDO_PP,NBNOR,RSIDO_PP,VCORD,rv,& 
     &                 ETA,NSIDE_PP,ISIDE_PP,RORDER,TORDER,SUMWEIGHT,R,&
     &           MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &           VSPACE_FIRST,VSPACE_LAST,MPI_RANK_P) 
! 
! *** THIS SUBROUTINE CALCULATE THE MASS CONSERVATION PARAMETER (ETA) AT 
! *** EACH BOUNDARY NODE FOR USE IN GETBOU.f WHEN CALCULATING THE FLUX OF 
! *** MOLECULES THAT HAVE BOUNCED OFF THE WALL 
!
      IMPLICIT NONE
      include 'mpif.h' 
! 
      INTEGER NELEM_PP,VNPNT,NBNOI,NBOUN_PP,IP1,IP2 
      INTEGER NBNOR,IB,IV,IS,IV1,IV2,IETA,IZETA 
      INTEGER NSIDE_PP,RORDER,TORDER,IV3,IV4
      INTEGER BSIDO_PP(NBNOI,NBOUN_PP),IEL,ISIDE_PP(8,NSIDE_PP) 
      INTEGER IN1,IN2,IP1T,IP2T,MPI_IERR,ITYPE 

      INTEGER MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V
      INTEGER MPI_RANK_P
      INTEGER VNPNT_PART,VSPACE_FIRST,VSPACE_LAST

      REAL RSIDO_PP(NBNOR,NBOUN_PP),ETA1,ZETA1,RT,THETA,THETA2 
      REAL VCORD(3,VNPNT) 
      REAL rv,ANX,ANY,TEST,VNORM,DIFF1,DIFF2,DIFF3,DIFF4 
      REAL CX,CY,PI,JAC,WEIGHT,SUMWEIGHT,FRAC 
      REAL Mw,nf,ETATEST1,ETATEST2,ZETATEST1,ZETATEST2
      REAL nf1,nf2,THETA1,ETAf,ZETAf,nf11,nf12 
      REAL SPEED,R,C1,C2,C3,C4,Co4
      REAL SUMBOTTOM,ETA(NBOUN_PP)
      REAL DISNF_PP(3,VSPACE_FIRST:VSPACE_LAST,NELEM_PP)
      REAL SUMTOP_ARR(NBOUN_PP),SUMTOPG_ARR(NBOUN_PP)
      REAL CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4
      CHARACTER filename*80
! 
! *** SET UP GAS CONSTANT AS A PARAMETER 
! 
      PARAMETER(PI=3.1416,Co4=4.0)
      FRAC=Co4/SUMWEIGHT
! 
! *** LOOP OVER THE BOUNDARY SIDES 
!
      DO IB=1, NBOUN_PP ! djahdfasdfjhalkhdfa
        IP1=BSIDO_PP(1,IB) 
        IP2=BSIDO_PP(2,IB) 
        IEL=BSIDO_PP(3,IB)
        ITYPE=BSIDO_PP(4,IB)
        IF(ITYPE.EQ.4)THEN !fgsdfgskjfgkljsfglkjsdfg
          ANX=RSIDO_PP(1,IB) 
          ANY=RSIDO_PP(2,IB) 
          IF((ANX.GT.0).AND.(ANY.GE.0))THEN
            THETA1=ATAN(ANY/ANX)-(PI/2)
          ELSEIF((ANX.LE.0).AND.(ANY.GT.0))THEN
            IF(ANX.EQ.0.0)THEN
              THETA1=PI/2
              GOTO 159
            ENDIF
            THETA1=(PI/2)-ATAN(ABS(ANY/ANX))
          ELSEIF((ANX.LT.0).AND.(ANY.LE.0))THEN
            THETA1=ATAN(ABS(ANY/ANX))+(PI/2)
          ELSEIF((ANX.GE.0).AND.(ANY.LT.0))THEN
            IF(ANX.EQ.0)THEN
              THETA1=-PI/2
              GOTO 159
            ENDIF  
            THETA1=-ATAN(ABS(ANY/ANX))-(PI/2)
          ENDIF
 159      CONTINUE
!
! ***     COMPUTE THETA2
!
          THETA2=THETA1+PI        !WE WILL THEN INTEGRATE BETWEEN THESE TWO ANGLES
          IF(THETA2.GT.(PI+0.01))THETA2=THETA2-2*PI
! 
! ***     FIND BOUNDARY SIDE EDGE IN ISIDE_PP TO GET LOCAL NODE NUMBERS 
! 
          DO IS=1,NSIDE_PP !dkljsdfaadfda
            IP1T=ISIDE_PP(1,IS) 
            IP2T=ISIDE_PP(2,IS) 
            IF((IP1.EQ.IP1T).AND.(IP2.EQ.IP2T))THEN 
              IN1=ISIDE_PP(5,IS) 
              IN2=ISIDE_PP(6,IS) 
              EXIT
            ELSEIF((IP2.EQ.IP1T).AND.(IP1.EQ.IP2T))THEN 
              IN1=ISIDE_PP(6,IS) 
              IN2=ISIDE_PP(5,IS) 
              EXIT
            ENDIF 
          ENDDO ! DO IS=1,NSIDE_PP !dkljsdfaadfda
          IF(IS.EQ.(NSIDE_PP+1))THEN
            WRITE(*,*) "MPI",MPI_RANK_P,"SIDE",IP1,IP2,IB,"NOT FOUND"
            STOP
          ENDIF  
! 
! ***     LOOP OVER ALL VELOCITY SPACE NODES IN THE TRANSFORMED HEMISPHERE THAT 
! ***     BEGINS AT THETA1
! 
          SUMTOP_ARR(IB)=0.0
          DO IV=VSPACE_FIRST,VSPACE_LAST !fdfgdvsdvjkkljnffgddd
            ETA1=VCORD(1,IV)                
            ZETA1=VCORD(2,IV)              
            RT=ETA1*(rv/2)+(rv/2)       !MAP BACK TO POLAR COORDINATES 
            THETA=ZETA1*PI
            IF(ANX.GE.0)THEN
              IF((THETA.GE.THETA1).AND.(THETA.LE.THETA2))THEN
                CX=RT*COS(THETA)            !CONVERT TO CARTESIANS 
                CY=RT*SIN(THETA)
! 
! ***           IS THE FLUX INTO OR AWAY FROM THE WALL? 
! 
                TEST=CX*ANX+CY*ANY
                VNORM=ABS(TEST) 
                JAC=PI*rv*rv*0.25*(ETA1+1) 
                WEIGHT=VCORD(3,IV)
                nf1=DISNF_PP(IN1,IV,IEL)
                nf2=DISNF_PP(IN2,IV,IEL)
                nf=0.5*(nf1+nf2)
!
                SUMTOP_ARR(IB)=SUMTOP_ARR(IB)+VNORM*nf*JAC*WEIGHT*FRAC
              ELSE
                CONTINUE
              ENDIF 
            ELSEIF(ANX.LE.0)THEN
              IF((THETA.GE.THETA1).OR.(THETA.LE.THETA2))THEN
                CX=RT*COS(THETA)            !CONVERT TO CARTESIANS 
                CY=RT*SIN(THETA)
!   
!   ***         IS THE FLUX INTO OR AWAY FROM THE WALL? 
!   
                TEST=CX*ANX+CY*ANY
                VNORM=ABS(TEST) 
                JAC=PI*rv*rv*0.25*(ETA1+1) 
                WEIGHT=VCORD(3,IV)
                nf1=DISNF_PP(IN1,IV,IEL)
                nf2=DISNF_PP(IN2,IV,IEL)
                nf=0.5*(nf1+nf2)
                SUMTOP_ARR(IB)=SUMTOP_ARR(IB)+VNORM*nf*JAC*WEIGHT*FRAC
              ELSE
                CONTINUE
              ENDIF 
            ENDIF       

! ***       END THE LOOP OVER VELOCITY SPACE 
          ENDDO ! DO IV=VSPACE_FIRST,VSPACE_LAST !fdfgdvsdvjkkljnffgddd
        ENDIF! IF(ITYPE.EQ.4)THEN !fgsdfgskjfgkljsfglkjsdfg
! 
! ***   END LOOP OVER THE BOUNDARY SIDES 
! 
      ENDDO ! DO IB=1, NBOUN_PP !djahdfasdfjhalkhdfa
!
      CALL MPI_ALLREDUCE(SUMTOP_ARR,SUMTOPG_ARR,NBOUN_PP,MPI_REAL,&
     &      MPI_SUM,MPI_COMM_V,MPI_IERR)

      DO IB=1, NBOUN_PP ! djahdfasdfjhalkhdfa2
        ITYPE=BSIDO_PP(4,IB)
        IF(ITYPE.EQ.4)THEN !fgsdfgskjfgkljsfglkjsdfg2

! ***     USE AN ANALYTICAL DERIVATION FOR THE DENOMINATOR
          SUMBOTTOM=SQRT(2*PI*(R**3)*(BSIDO_PP(5,IB))**3)
! 
! ***     CALCULATE ETA 
! 
          IF(SUMBOTTOM.LT.0.1)THEN 
            WRITE(*,500) 
            WRITE(*,501) 
            STOP 
          ENDIF
!	
          ETA(IB)=SUMTOPG_ARR(IB)/SUMBOTTOM
!
        ELSE ! IF(ITYPE.EQ.4)THEN !fgsdfgskjfgkljsfglkjsdfg2
          ETA(IB)=0.0
        ENDIF! IF(ITYPE.EQ.4)THEN !fgsdfgskjfgkljsfglkjsdfg2
! 
! ***   END LOOP OVER THE BOUNDARY SIDES 
! 
      ENDDO ! DO IB=1, NBOUN_PP !djahdfasdfjhalkhdfa2
! 
! *** FORMAT STATEMENTS 
! 
  500 FORMAT('ERROR CALCULATING FLUX CONSERVATION PARAMETER IN FLUCON') 
  501 FORMAT('PROGRAM STOPPED!') 
      RETURN 
      END 
         
 
         
       
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
       
