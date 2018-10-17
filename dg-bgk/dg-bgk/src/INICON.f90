       SUBROUTINE INICON(NPOIN,NNODE,NELEM,NELEM_PP,INTMA_PP,ELGRP,& 
     &     NPOIN_PP,COORD_PP,VNPNT,VCORD,DISNF_PP,IPCOM_PP,rv,IVD,&
     &                 MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
     &                 MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &                 VSPACE_FIRST,VSPACE_LAST) 
!
! 
! *** THIS SUBROUTINE ESTABLISHES THE INITIAL CONDITIONS OF THE PROBLEM 
! 
       IMPLICIT NONE 
      INCLUDE 'mpif.h' 
      TYPE InputVariables
          INTEGER :: TORDER !Quadrature order in Theta for V-SPACE
          INTEGER :: NTIME !Number of timesteps
          INTEGER :: IMMAT !Lumped mass matrices (no/yes)
          INTEGER :: RS !Using Restart data (no/yes)
          INTEGER :: INF !Is there an inflow (no/yes)
          INTEGER :: NVSPACEPART ! number of VSPACE partitions
          REAL :: CSAFM !Safety factor applied to timestep (Courant)
          REAL :: rv !Radial extent of the V-SPACE
          REAL :: T1 !Initial condition temp
          REAL :: P1 !Initial condition pressure
          REAL :: U0 !Initial condition X-vel
          REAL :: V0 !Initial condition Y-vel
          REAL :: CINF(4)
          REAL :: W
          REAL :: ALPHA
      END TYPE

! *** DECLARE VARIABLES 
! 
      TYPE(InputVariables) :: IVD
      INTEGER IN,INV,RS,IV,IP,NELEM,NPOIN,NNODE,maxNELEM_PP 
      INTEGER IE,TAG,IP_PP,I,IE_PP,NELEM_PP 
      INTEGER VNPNT,NPOIN_pp,IPT,INTMA_PP(NNODE,NELEM_PP) 
      INTEGER TEST1,TEST2,NPART,IPCOM_PP(NPOIN_PP) 
      INTEGER ELGRP(NELEM,2),MPI_IERR,MPI_STATUS 
! 
      REAL VCORD(3,VNPNT),COORD_PP(2,NPOIN_PP) 
      REAL T1,P1,GC,RHO1,n1,NA,C0,ETA,ZETA,R,THETA 
      REAL M,BETA1,TMP1,MAGX,rv,C 
      REAL UX,UY,TMP,SPEED,CO,PI,F0,XPOS 
      REAL SUM,W,L 

      ! mpi-stuff for position space partitioning
      INTEGER MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,GROUP_P
      ! mpi-stuff for velocity space partitioning
      INTEGER MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,GROUP_V
      INTEGER VSPACE_FIRST,VSPACE_LAST

      REAL DISNF_PP(NNODE,VSPACE_FIRST:VSPACE_LAST,NELEM_PP)
      REAL DISNFPARCEL(NNODE) 

! 
      CHARACTER filename*80,TEXT*80 
! 
! *** SET PI, GAS CONSTANT (GC),MOLAR MASS AND AVAGADRO'S NUMBER AS PARAMETERS 
! 
      PARAMETER(GC=0.287E+03,NA=6.022E+026,M=29,PI=3.1416)
! 
! *** INITIALISE DISNF_PP 
! 
      IF(MPI_RANK_P.NE.0)THEN 
       DO IE_PP=1,NNODE
          DO IV=VSPACE_FIRST,VSPACE_LAST
            DO IN=1,NNODE
              DISNF_PP(IN,IV,IE_PP) = 0.0
            ENDDO
          ENDDO
        ENDDO
      ENDIF
! 
! *** ASK THE USER IF THEY ARE USING A RESTART FILE AS INITIAL CONDITIONS 
! 
      RS = IVD%RS 
! 
! *** IF USING A RESTART FILE ASK FOR THE NAME OF THE FILE 
! 
      IF(RS.EQ.1)THEN ! sdflgjsdfglkjsfdgs
        IF(MPI_RANK_P.EQ.0)THEN !cdifjaoncaoda
151       READ(21,*) filename 
          OPEN  (11, file=filename, status='old', err=151) 
          WRITE(*,*) 'READING RESTART FILE =',filename 
          WRITE(*,*) 
! 
! ***   READ IN DATA FROM RESTRT FILE 
! 
          READ(11,*) TEXT 
          READ(11,*) TEST1,TEST2 
! 
! ***   CHECK DATA 
! 
          IF((TEST1.NE.NPOIN).OR.(TEST2.NE.VNPNT))THEN 
            WRITE(*,208) 
            WRITE(*,209) 
            STOP 
          ENDIF 
! 
        ENDIF ! IF(MPI_RANK_P.EQ.0)THEN !cdifjaoncaoda
! 
! ***   READ DISNF_PP 
! 
        ! THIS IS SUPER SLOW BUT FAST TO IMPLEMENT
        DO 1010 IV=1,VNPNT 
          DO 1011 IE=1,NELEM 
            IF(MPI_RANK_P.EQ.0)THEN !qqwerqwerqw
              READ(11,*) (DISNFPARCEL(I),I=1,NNODE) 
! ***     CHECK THE SLAVE ALLOCATION OF THIS ELEMENT 
              NPART=ELGRP(IE,1) 
              IE_PP=ELGRP(IE,2) 
            ENDIF   !qqwerqwerqw
! ***     BROADCAST 
            IF((IV.GE.VSPACE_FIRST).AND.(IV.LE.VSPACE_LAST))THEN !dsfkjshd
              CALL MPI_BCAST(NPART,1,MPI_INTEGER,0,& 
     &          MPI_COMM_P,MPI_IERR) 
              CALL MPI_BCAST(IE_PP,1,MPI_INTEGER,0,& 
     &          MPI_COMM_P,MPI_IERR) 
              CALL MPI_BCAST(DISNFPARCEL,NNODE,MPI_REAL,0,& 
     &          MPI_COMM_P,MPI_IERR) 
              IF((MPI_RANK_P+1).EQ.NPART)THEN !ljvvkjkddfgd
                DO IN=1,NNODE 
                  DISNF_PP(IN,IV,IE_PP)=DISNFPARCEL(IN) 
                ENDDO 
              ENDIF   !ljvvkjkddfgd  
            ENDIF  !dsfkjshd
 1011     CONTINUE 
 1010   CONTINUE 
! 
! ***   CLOSE RESTART FILE 
! 
        IF(MPI_RANK_P.EQ.0)CLOSE(11) 
! 
      ELSE ! IF(RS.EQ.1)THEN ! sdflgjsdfglkjsfdgs
! 
! ***   INITIAL CONDITIONS FOR A GAS EXPANSION PROBLEM: 
! 
! ***   ASK THE USER FOR THE PRE-EXPANSION GAS CONDITIONS: 
! 
          IF(MPI_RANK_P.EQ.0)THEN    ! dsfkjccsdfsdfaaa 
! 
            T1 = IVD%T1 
            P1 = IVD%P1 
            W = IVD%W 
            WRITE(*,*) 
            L=(W*1E-09)/2 
! 
! ***    CALCULATE GAS PROPERTIES USED FOR MAXWELL DISTRIBUTION 
!        
            RHO1=(P1*(10**05))/(GC*T1)     !DENSITIES 
            n1=RHO1*NA/M                   !NUMBER DENSITIES 
            TMP1=RHO1/(2*P1*(10**05)) 
            BETA1=SQRT(TMP1)               !BETA VARIABLE 
          ENDIF   !IF(MPI_RANK_P.EQ.0)THEN    ! dsfkjccsdfsdfaaa 
! 
! ***   BROADCAST VARIABLES 
! 
        CALL MPI_BCAST(L,1,MPI_REAL,0,MPI_COMM_P,MPI_IERR) 
        CALL MPI_BCAST(BETA1,1,MPI_REAL,0,MPI_COMM_P,MPI_IERR) 
        CALL MPI_BCAST(n1,1,MPI_REAL,0,MPI_COMM_P,MPI_IERR)
! 
! ***   BEGIN LOOP OVER THE PHYSICAL SPACE DISCONTINUOUS NODES 
!
          DO 1005 IE=1,NELEM_PP
            DO 1006 IN=1,NNODE 
              IP_PP=INTMA_PP(IN,IE) 
! 
! ***         X-COORDINATE OF NODE 
! 
              XPOS=COORD_PP(1,IP_PP)
! 
! ***   INSIDE THE GAS CLOUD 
! 
              MAGX=ABS(XPOS)
              IF(MAGX.LE.L)THEN ! sdfkacakjnaa
! ***    LOOP OVER ALL VELOCITY SPACE NODES 
                C=(BETA1**2)/(PI)         !COEFFICIENT OF THE MAXWELL DISTRIBUTION FUNCTION
                DO 1002 INV=VSPACE_FIRST,VSPACE_LAST
                  ETA=VCORD(1,INV) 
                  ZETA=VCORD(2,INV) 
! ***    HERE WE TRANSFORM THE (ETA,ZETA) COORDINATES BACK TO (CX,CY) COORDINATES 
                  R=ETA*(rv/2)+(rv/2)     !FIRST TO POLAR COORDINATES 
                  THETA=ZETA*PI 
                  UX=R*COS(THETA)          !THEN TO CARTESIANS 
                  UY=R*SIN(THETA) 
                  TMP=UX*UX+UY*UY         
                  SPEED=SQRT(TMP)  !MOLECULAR SPEED AT THIS COORDINATE IN VELSPACE MESH 
                  TMP=-((BETA1**2)*(SPEED**2)) 
                  F0=C*EXP(TMP)       !DISTRIBUTION FUNCTION 
                  DISNF_PP(IN,INV,IE)=n1*F0     !INITIAL CONDITION (nf) MATRIX 
! ***    END LOOP OVER VELOCITY SPACE NODES 
 1002           CONTINUE
! ***    IF OUTSIDE THE GAS CLOUD, SET nf TO ZERO 
              ELSE ! IF(MAGX.LE.L)THEN ! sdfkacakjnaa
                DO 1003 INV=VSPACE_FIRST,VSPACE_LAST
                  DISNF_PP(IN,INV,IE)=0.0
 1003           CONTINUE
              ENDIF ! IF(MAGX.LE.L)THEN ! sdfkacakjnaa
! 
! ***    END LOOP OVER THE PHYSICAL SPACE DISCONTINUOUS NODES 
! 
 1006       CONTINUE
 1005     CONTINUE
!
      ENDIF !  IF(RS.EQ.1)THEN ! sdflgjsdfglkjsfdgs
! 
! *** FORMAT STATEMENTS 
! 
 201  FORMAT('WHAT IS THE TEMPERATURE(IN K) OF THE GAS?', $) 
 202  FORMAT('WHAT IS THE PRESSURE (IN bar) OF THE GAS?', $) 
 203  FORMAT('WHAT IS THE WIDTH OF THE GAS CLOUD (in nanometers) ?', $) 
 204  FORMAT('DO WISH TO START THE RUN USING A RESTART FILE? ') 
 205  FORMAT('(YES = 1, NO = 0)',$) 
 206  FORMAT('WHAT IS THE NAME OF THE RESTART FILE?  ',$) 
 208  FORMAT('THIS RESTART FILE IS NOT COMPATIBLE WITH YOUR MESHES!') 
 209  FORMAT('PROGRAM STOPPED') 
! 
      RETURN 
      END 
