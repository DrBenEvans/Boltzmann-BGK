        SUBROUTINE GETLOC(NSIDE,NGEOM,NELEM,maxNELEM_PP,NPOIN,NBNOI,&
     &      NBOUN,NBNOR,NNODE,MXCOM_PP,NGRPS,ELGRP,NEGRP,&
     &      GCOMM,NX,NY,EL,GEOME,MMAT,CMMAT,INTMA,&
     &      COORD,BSIDO,RSIDO,ISIDE,NELEM_PP,&
     &      maxNBOUN_PP,NBOUN_PP,maxNPOIN_PP,NPOIN_PP,&
     &      maxNSIDE_PP,NSIDE_PP,ISIDE_PP,&
     &      NX_PP,NY_PP,EL_PP,GEOME_PP,MMAT_PP,CMMAT_PP,&
     &      INTMA_PP,COORD_PP,BSIDO_PP,RSIDO_PP,IMMAT,&
     &      MPI_RANK_P,LCOMM_PP,IPCOM_PP,IBCOM_PP,ISCOM_PP,&
     &      NPGRP,NCOMM_PP,SDCOM_PP,MPI_COMM_P)
!
! *** THIS SUBROUTINE CONVERTS GLOBAL VARIABLES/ARRAYS INTO LOCAL PROCESSOR VARIABLES/ARRAYS
!
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!
      INTEGER FLAG,NNODE,NSIDE_PP,maxNSIDE_PP,NPOIN_PP,NCOMM_PP
      INTEGER NBOUN_PP,maxNBOUN_PP,NELEM_PP,IP,IP_PP,maxNELEM_PP
      INTEGER NSIDE,NGEOM,NELEM,NPOIN,NBNOI,NBOUN,NBNOR
      INTEGER NGRPS,RANK,TAG1,TAG2,TAG3,SIZE1,SIZE2,SIZE3,TAG
      INTEGER TAG4,TAG5,SIZE4,SIZE5,TAG6,SIZE6,TAG7,SIZE7,TAG8
      INTEGER SIZE8,TAG9,SIZE9,TAG10,SIZE10,TAG11,SIZE11,TAG12
      INTEGER SIZE12,TAG13,SIZE13,TAG14,SIZE14,TAG15,SIZE15
      INTEGER TAG16,SIZE16,TAG17,SIZE17,MPI_RANK_P,MPI_IERR
      INTEGER TAG18,SIZE18
      INTEGER IE1,IE2,IS,RL,RR,IR,IE,NELEM_G,CO,IEG,IG,IGT,I
      INTEGER IEL,IER,RANKL,RANKR,IPT,J,IB,maxNPOIN_PP
      INTEGER IP1,IP2,IP3,IMMAT,MPI_STATUS(MPI_STATUS_SIZE)
      INTEGER MXCOM_PP,NPROC,GCOMM(NGRPS,NGRPS)
      INTEGER ELGRP(NELEM,2),ISIDE(8,NSIDE)
      INTEGER INTMA(NNODE,NELEM),BSIDO(NBNOI,NBOUN)
      INTEGER INTMA_PP(NNODE,maxNELEM_PP),SDCOM_PP(3,maxNSIDE_PP)
      INTEGER BSIDO_PP(NBNOI,maxNBOUN_PP),NPGRP(NGRPS)
      INTEGER NEGRP(NGRPS),ISIDE_PP(8,maxNSIDE_PP)
      INTEGER IPCOM_PP(maxNPOIN_pp),IBCOM_PP(maxNBOUN_pp) 
      INTEGER ISCOM_PP(maxNSIDE_PP),LCOMM_PP(maxNSIDE_PP)
      INTEGER MPI_COMM_P
!
      REAL NX(NSIDE),NY(NSIDE),EL(NSIDE),GEOME(NGEOM,NELEM)
      REAL MMAT(NNODE,NELEM),CMMAT(3,3,NELEM),COORD(2,NPOIN)
      REAL RSIDO(NBNOR,NBOUN),NX_PP(maxNSIDE_PP),NY_PP(maxNSIDE_PP)
      REAL EL_PP(maxNSIDE_PP),GEOME_PP(NGEOM,maxNELEM_PP)
      REAL MMAT_PP(NNODE,maxNELEM_PP),CMMAT_PP(3,3,maxNELEM_PP)
      REAL COORD_PP(2,maxNPOIN_PP),RSIDO_PP(NBNOR,maxNBOUN_PP) 
      REAL COO,CX,CY

      INTEGER NPOIN_PP_CP,NBOUN_PP_CP,NCOMM_PP_CP
      INTEGER NSIDE_PP_CP,NELEM_PP_CP
      INTEGER,ALLOCATABLE :: INTMA_PP_CP(:,:),IPCOM_PP_CP(:)
      INTEGER,ALLOCATABLE ::BSIDO_PP_CP(:,:),ISIDE_PP_CP(:,:)
      REAL,ALLOCATABLE :: RSIDO_PP_CP(:,:),COORD_PP_CP(:,:)
      REAL,ALLOCATABLE :: GEOME_PP_CP(:,:),MMAT_PP_CP(:,:)
      REAL,ALLOCATABLE :: CMMAT_PP_CP(:,:,:),NX_PP_CP(:),NY_PP_CP(:)
      REAL,ALLOCATABLE :: EL_PP_CP(:)
      INTEGER,ALLOCATABLE :: LCOMM_PP_CP(:),IBCOM_PP_CP(:)
      INTEGER,ALLOCATABLE :: ISCOM_PP_CP(:)
!
      IF(MPI_RANK_P.EQ.0)THEN ! vnsdsjkdkja
        ALLOCATE(INTMA_PP_CP(NNODE,maxNELEM_PP))
        ALLOCATE(IPCOM_PP_CP(maxNPOIN_PP))
        ALLOCATE(COORD_PP_CP(2,maxNPOIN_PP))

        DO 10001 IG=1,NGRPS
! ***     NELEM_G IS THE NUMBER OF ELEMENTS IN GROUP OF ELEMENTS IG
          NELEM_G=NEGRP(IG)
! ***     INTITIALISE INTMA_PP,IPCOM_PP AND IBCOM_PP
          CO=0
          CALL IFILLM(INTMA_PP_CP,NNODE,NELEM_G,CO)
          CALL IFILLV(IPCOM_PP_CP,maxNPOIN_pp,CO)
!
          FLAG=1 ! local point index
          I=0
!
          DO 1001 IE=1,NELEM_G    !LOOP OVER THE PROCESSOR'S ELEMENTS          
 1002       CONTINUE
            I=I+1
            IGT=ELGRP(I,1)
            IF(IG.EQ.IGT)THEN
              IEG=ELGRP(I,2)
              IP1=INTMA(1,I)
              IP2=INTMA(2,I)
              IP3=INTMA(3,I)
! ***         FILL IN IPCOM_PP_CP AND INTMA_PP_CP
              ! scan the local point array to find IP1 in IPCOM_PP_CP
              DO 1003 IP_PP=1,maxNPOIN_pp
                IP=IPCOM_PP_CP(IP_PP)
                IF(IP.EQ.IP1)THEN     ! IP1 found
                  INTMA_PP_CP(1,IEG)=IP_PP 
                  GOTO 2000
                ENDIF
 1003         CONTINUE              ! IP1 not found
              INTMA_PP_CP(1,IEG)=FLAG
              IPCOM_PP_CP(FLAG)=IP1    ! add it to IPCOM_PP_CP
              FLAG=FLAG+1
              IF(FLAG.GT.maxNPOIN_pp)THEN
                WRITE(*,101);WRITE(*,102);STOP
              ENDIF
 2000         CONTINUE
              ! scan the local point array to find IP2 in IPCOM_PP_CP
              DO 1004 IP_PP=1,maxNPOIN_pp
                IP=IPCOM_PP_CP(IP_PP)
                IF(IP.EQ.IP2)THEN     ! IP2 found
                  INTMA_PP_CP(2,IEG)=IP_PP
                  GOTO 2001
                ENDIF
 1004         CONTINUE              ! IP2 not found
              INTMA_PP_CP(2,IEG)=FLAG
              IPCOM_PP_CP(FLAG)=IP2    ! add it to IPCOM_PP_CP
              FLAG=FLAG+1
              IF(FLAG.GT.maxNPOIN_pp)THEN
                WRITE(*,101);WRITE(*,102);STOP
              ENDIF
 2001         CONTINUE
!
              ! scan the local point array to find IP3 in IPCOM_PP_CP
              DO 1005 IP_PP=1,maxNPOIN_pp
                IP=IPCOM_PP_CP(IP_PP)
                IF(IP.EQ.IP3)THEN      ! IP3 found
                  INTMA_PP_CP(3,IEG)=IP_PP
                  GOTO 2002
                ENDIF
 1005         CONTINUE               ! IP3 not found
              INTMA_PP_CP(3,IEG)=FLAG
              IPCOM_PP_CP(FLAG)=IP3     ! add it to IPCOM_PP_CP
              FLAG=FLAG+1
              IF(FLAG.GT.maxNPOIN_pp)THEN
                WRITE(*,101);WRITE(*,102);STOP
              ENDIF
 2002         CONTINUE
            ELSE ! RELATED TO IF(IG.EQ.IGT)
              GOTO 1002 ! BASICALLY "CYCLE", USELESS
            ENDIF! RELATED TO IF(IG.EQ.IGT)

! ***       END LOOP OVER PROCESSOR ELEMENTS
 1001     CONTINUE
          NELEM_PP_CP=NEGRP(IG)
          NPOIN_PP_CP=FLAG-1
          NPGRP(IG)=NPOIN_PP_CP            
! ***     MPI COMMAND TO SEND NPOIN_PP_CP,INTMA_PP_CP AND IPCOM_PP_CP TO PROCESSOR RANK IG !!!!!!!!!!!!!!!
          TAG1=1
          SIZE1=1
          TAG2=2
          SIZE2=NNODE*NELEM_PP_CP
          TAG3=3
          SIZE3=maxNPOIN_PP
          IF(IG.NE.1)THEN
            CALL MPI_SEND(NPOIN_PP_CP,SIZE1,MPI_INTEGER,IG-1,TAG1,&
     &             MPI_COMM_P,MPI_IERR)
            CALL MPI_SEND(INTMA_PP_CP,SIZE2,MPI_INTEGER,IG-1,TAG2,&
     &              MPI_COMM_P,MPI_IERR)
            CALL MPI_SEND(IPCOM_PP_CP,SIZE3,MPI_INTEGER,IG-1,TAG3,&
     &                MPI_COMM_P,MPI_IERR)
          ELSE
            NPOIN_PP = NPOIN_PP_CP
            INTMA_PP(:,1:NELEM_PP) = INTMA_PP_CP(:,1:NELEM_PP)
            IPCOM_PP(1:NPOIN_PP_CP)   = IPCOM_PP_CP(1:NPOIN_PP_CP)  
          ENDIF
10001   CONTINUE
        DEALLOCATE(INTMA_PP_CP)
      ENDIF ! RELATED TO IF(MPI_RANK_P.EQ.0)  vnsdsjkdkja
!      CALL MPI_BARRIER(MPI_COMM_P,MPI_IERR)
!
      IF(MPI_RANK_P.NE.0)THEN
        TAG1=1
        SIZE1=1
        TAG2=2
        SIZE2=NNODE*NELEM_PP
        CALL MPI_RECV(NPOIN_PP,SIZE1,MPI_INTEGER,0,& 
     &           TAG1,MPI_COMM_P,MPI_STATUS,MPI_IERR )
        CALL MPI_RECV(INTMA_PP,SIZE2,MPI_INTEGER,0,& 
     &           TAG2,MPI_COMM_P,MPI_STATUS,MPI_IERR )
        TAG3=3
        SIZE3=maxNPOIN_PP
        CALL MPI_RECV(IPCOM_PP,SIZE3,MPI_INTEGER,0,& 
     &           TAG3,MPI_COMM_P,MPI_STATUS,MPI_IERR )
      ENDIF
!
! *** GET COORD_PP FOR THIS PROCESSOR
!
      print *,'Here we go'
      CALL MPI_BARRIER(MPI_COMM_P,MPI_IERR)

      COO=0.0
      IF(MPI_RANK_P.EQ.0)THEN !ccmfcjskha
        DO 10003 IG=1,NGRPS
          NELEM_PP_CP=NEGRP(IG)
              TAG=IG*10
          IF(IG.NE.1) THEN
          CALL MPI_RECV(IPCOM_PP_CP,maxNPOIN_PP,MPI_INTEGER,IG-1,&
     &               TAG,MPI_COMM_P,MPI_STATUS,MPI_IERR)
          ELSE
            IPCOM_PP_CP(1:maxNPOIN_PP) = IPCOM_PP(1:maxNPOIN_PP)
          ENDIF
          CALL RFILLM(COORD_PP_CP,2,maxNPOIN_PP,COO)
!
          NPOIN_PP_CP=NPGRP(IG)
          DO 1006 IP_PP=1,NPOIN_PP_CP
            IP=IPCOM_PP_CP(IP_PP)
            CX=COORD(1,IP)
            CY=COORD(2,IP)
            COORD_PP_CP(1,IP_PP)=CX
            COORD_PP_CP(2,IP_PP)=CY
 1006     CONTINUE
!
! ***     MPI COMMAND TO SEND COORD_PP TO PROCESSOR RANK IG !!!!!!!!!!!!!!!!!!!!!!!!!
!
          TAG4=4
          SIZE4=2*NPOIN_PP_CP
          IF(IG.NE.1) THEN
            CALL MPI_SEND(COORD_PP_CP,SIZE4,MPI_REAL,IG-1,& 
     &             TAG4,MPI_COMM_P,MPI_IERR )
          ELSE
              COORD_PP = COORD_PP_CP
          ENDIF  
10003   CONTINUE
        DEALLOCATE(COORD_PP_CP)
      ENDIF  ! IF(MPI_RANK_P.EQ.0)THEN !ccmfcjskha
!
      IF(MPI_RANK_P.NE.0)THEN ! ccmdxkaljfa
        TAG=(MPI_RANK_P+1)*10
        CALL MPI_SEND(IPCOM_PP,maxNPOIN_PP,MPI_INTEGER,&
     &              0,TAG,MPI_COMM_P,MPI_IERR)
        TAG4=4
        SIZE4=2*NPOIN_PP
        CALL MPI_RECV(COORD_PP,SIZE4,MPI_REAL,0,& 
     &           TAG4,MPI_COMM_P,MPI_STATUS,MPI_IERR )
      ENDIF !IF(MPI_RANK_P.NE.0)THEN ! ccmdxkaljfa
!
! *** NEXT DETERMINE BSIDO_PP AND RSIDO_PP
!
      IF(MPI_RANK_P.EQ.0)THEN  ! diescaa
        ALLOCATE(BSIDO_PP_CP(NBNOI,maxNBOUN_PP))
        ALLOCATE(RSIDO_PP_CP(NBNOR,maxNBOUN_PP))
        ALLOCATE(IBCOM_PP_CP(maxNBOUN_pp))
        DO 10005 IG=1,NGRPS
          NELEM_PP_CP=NEGRP(IG)
          CALL IFILLM(BSIDO_PP_CP,NBNOI,maxNBOUN_PP,CO)
          CALL RFILLM(RSIDO_PP_CP,NBNOR,maxNBOUN_PP,COO)
          CALL IFILLV(IBCOM_PP_CP,maxNBOUN_PP,CO)
          TAG=IG*10
          IF(IG.NE.1)THEN
            CALL MPI_RECV(IPCOM_PP_CP,maxNPOIN_PP,MPI_INTEGER,IG-1,&
     &               TAG,MPI_COMM_P,MPI_STATUS,MPI_IERR)
          ELSE
            IPCOM_PP_CP(1:maxNPOIN_PP) = IPCOM_PP(1:maxNPOIN_PP)
          ENDIF
!
! *** L  OOP OVER ALL THE BOUNDARY SIDES
!
          FLAG=0
          DO 1007 IB=1,NBOUN
            IE=BSIDO(3,IB)
            RANK=ELGRP(IE,1)
            IF(RANK.EQ.IG)THEN
              FLAG=FLAG+1
              IF(FLAG.GT.maxNBOUN_PP)THEN
                WRITE(*,103);WRITE(*,104);STOP
              ENDIF
              NPOIN_PP_CP=NPGRP(IG)
              DO I=1,NPOIN_PP_CP ! scan for BSIDO(1,IB) om IPCOM_PP
                IPT=IPCOM_PP_CP(I)
                IP=BSIDO(1,IB)
                IF(IP.EQ.IPT)THEN ! found 
                  BSIDO_PP_CP(1,FLAG)=I
                  GOTO 1030
                ENDIF
              ENDDO
 1030         CONTINUE
              DO I=1,NPOIN_PP_CP  ! scan for BSIDO(1,IB) om IPCOM_PP
                IPT=IPCOM_PP_CP(I)
                IP=BSIDO(2,IB)
                IF(IP.EQ.IPT)THEN ! found
                  BSIDO_PP_CP(2,FLAG)=I
                  GOTO 1031
                ENDIF
              ENDDO
 1031         CONTINUE
              IEL=BSIDO(3,IB)
              BSIDO_PP_CP(3,FLAG)=ELGRP(IEL,2)
              BSIDO_PP_CP(4,FLAG)=BSIDO(4,IB)
              BSIDO_PP_CP(5,FLAG)=BSIDO(5,IB)
              BSIDO_PP_CP(6,FLAG)=BSIDO(6,IB)
              DO I=1,NBNOR
                RSIDO_PP_CP(I,FLAG)=RSIDO(I,IB)
              ENDDO
              IBCOM_PP_CP(FLAG)=IB
            ENDIF
 1007     CONTINUE
          NBOUN_PP_CP=FLAG
!
! *** M  PI COMMAND TO SEND NBOUN_PP_CP,BSIDO_PP_CP,RSIDO_PP_CP AND IBCOM_PP_CP TO PROCESSOR RANK IR !!!!!!!
!
          TAG5=5
          TAG6=6
          SIZE6=NBNOI*maxNBOUN_PP
          TAG7=7
          SIZE7=NBNOR*maxNBOUN_PP
          TAG8=8
          SIZE8=maxNBOUN_PP
         
          IF(IG.NE.1)THEN 
            CALL MPI_SEND(NBOUN_PP_CP,1    ,MPI_INTEGER,IG-1,& 
     &             TAG5,MPI_COMM_P,MPI_IERR)
            CALL MPI_SEND(BSIDO_PP_CP,SIZE6,MPI_INTEGER,IG-1,& 
     &             TAG6,MPI_COMM_P,MPI_IERR)
            CALL MPI_SEND(RSIDO_PP_CP,SIZE7,MPI_REAL,IG-1,&
     &             TAG7,MPI_COMM_P,MPI_IERR) 
            CALL MPI_SEND(IBCOM_PP_CP,SIZE8,MPI_INTEGER,IG-1,&
     &             TAG8,MPI_COMM_P,MPI_IERR)
          ELSE 
            NBOUN_PP = NBOUN_PP_CP
            BSIDO_PP(:,1:NBOUN_PP_CP) = BSIDO_PP_CP(:,1:NBOUN_PP_CP)
            RSIDO_PP(:,1:NBOUN_PP_CP) = RSIDO_PP_CP(:,1:NBOUN_PP_CP)
            IBCOM_PP(1:NBOUN_PP_CP) = IBCOM_PP_CP(1:NBOUN_PP_CP)
          ENDIF
10005   CONTINUE
        DEALLOCATE(BSIDO_PP_CP)
        DEALLOCATE(RSIDO_PP_CP)
        DEALLOCATE(IBCOM_PP_CP)
      ENDIF    ! IF(MPI_RANK_P.EQ.0)THEN  ! diescaa
!
      IF(MPI_RANK_P.NE.0)THEN  !dsjkhadfxx
        TAG=(MPI_RANK_P+1)*10
        CALL MPI_SEND(IPCOM_PP,maxNPOIN_PP,MPI_INTEGER,&
     &    0,TAG,MPI_COMM_P,MPI_IERR)
        TAG5=5
        CALL MPI_RECV(NBOUN_PP,1    ,MPI_INTEGER,0,& 
     &   TAG5,MPI_COMM_P,MPI_STATUS,MPI_IERR )
        TAG6=6
        SIZE6=NBNOI*maxNBOUN_PP
        CALL MPI_RECV(BSIDO_PP,SIZE6,MPI_INTEGER,0,&
     &    TAG6,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        TAG7=7
        SIZE7=NBNOR*maxNBOUN_PP
        CALL MPI_RECV(RSIDO_PP,SIZE7,MPI_REAL,0,&
     &    TAG7,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        TAG8=8
        SIZE8=maxNBOUN_PP
        CALL MPI_RECV(IBCOM_PP,SIZE8,MPI_INTEGER,0,&
     &     TAG8,MPI_COMM_P,MPI_STATUS,MPI_IERR)
      ENDIF  ! IF(MPI_RANK_P.NE.0)THEN  !dsjkhadfxx
!
! *** NEXT DETERMINE GEOME_PP,MMAT_PP,CMMAT_PP
!
      CALL MPI_BARRIER(MPI_COMM_P,MPI_IERR)

      IF(MPI_RANK_P.EQ.0)THEN ! jdhdcaaaaaa
        ALLOCATE(GEOME_PP_CP(NGEOM,maxNELEM_PP))
        ALLOCATE(MMAT_PP_CP(NNODE,maxNELEM_PP))
        ALLOCATE(CMMAT_PP_CP(3,3,maxNELEM_PP))
        DO 10007 IG=1,NGRPS
          NELEM_PP_CP=NEGRP(IG)
!
          CALL RFILLM(GEOME_PP_CP,NGEOM,maxNELEM_PP,COO)
          CALL RFILLM(MMAT_PP_CP,NNODE,maxNELEM_PP,COO)
          CALL RFILLA(CMMAT_PP_CP,3,3,maxNELEM_PP,COO)
!
          DO 1008 IE=1,NELEM 
            RANK=ELGRP(IE,1)
            IEG=ELGRP(IE,2)
            IF(RANK.EQ.IG)THEN
              DO I=1,NGEOM
                GEOME_PP_CP(I,IEG)=GEOME(I,IE)
              ENDDO
              IF(IMMAT.EQ.1)THEN
                DO I=1,NNODE
                  MMAT_PP_CP(I,IEG)=MMAT(I,IE)
                ENDDO
              ELSE
                DO I=1,3
                  DO J=1,3
                   CMMAT_PP_CP(I,J,IEG)=CMMAT(I,J,IE)
                  ENDDO
                ENDDO
              ENDIF      
            ENDIF
 1008     CONTINUE
!
! ***     MPI COMMAND TO SEND GEOME_PP_CP,MMAT_PP_CP AND CMMAT_PP_CP TO PROCESSOR RANK IG !!!
!
          TAG9=9
          NELEM_PP_CP=NEGRP(IG)
          SIZE9=NGEOM*maxNELEM_PP
          IF(IG.NE.1)THEN
            CALL MPI_SEND(GEOME_PP_CP,SIZE9,MPI_REAL,IG-1,& 
     &               TAG9,MPI_COMM_P,MPI_IERR)
            IF(IMMAT.EQ.1)THEN
              TAG10=10
              SIZE10=maxNELEM_PP*NNODE
              CALL MPI_SEND(MMAT_PP_CP,SIZE10,MPI_REAL,IG-1,&
     &              TAG10,MPI_COMM_P,MPI_IERR)
            ELSE
              TAG10=10
              SIZE10=3*3*maxNELEM_PP
              CALL MPI_SEND(CMMAT_PP_CP,SIZE10,MPI_REAL,IG-1,&
     &          TAG10,MPI_COMM_P,MPI_IERR)
            ENDIF
          ELSE
            GEOME_PP(:,1:maxNELEM_PP) = GEOME_PP_CP(:,1:maxNELEM_PP)

            IF(IMMAT.EQ.1)THEN
              MMAT_PP(:,1:maxNELEM_PP) = MMAT_PP_CP(:,1:maxNELEM_PP)
            ELSE
              CMMAT_PP(:,:,1:maxNELEM_PP) =&
     &               CMMAT_PP_CP(:,:,1:maxNELEM_PP)
            ENDIF
          ENDIF
10007   CONTINUE
        DEALLOCATE(GEOME_PP_CP)
        DEALLOCATE(MMAT_PP_CP)
        DEALLOCATE(CMMAT_PP_CP)
      ENDIF       ! IF(MPI_RANK_P.EQ.0)THEN ! jdhdcaaaaaa
!
      IF(MPI_RANK_P.NE.0)THEN !cncsskdjhfsdf
        CALL RFILLM(GEOME_PP,NGEOM,maxNELEM_PP,COO)
        CALL RFILLM(MMAT_PP,NNODE,maxNELEM_PP,COO)
        CALL RFILLA(CMMAT_PP,3,3,maxNELEM_PP,COO)
        TAG9=9
        SIZE9=NGEOM*maxNELEM_PP
        CALL MPI_RECV(GEOME_PP,SIZE9,MPI_REAL,0,&
     &       TAG9,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        IF(IMMAT.EQ.1)THEN 
          TAG10=10
          SIZE10=maxNELEM_PP*NNODE
          CALL MPI_RECV(MMAT_PP,SIZE10,MPI_REAL,0,&
     &           TAG10,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        ELSE
          TAG10=10
          SIZE10=3*3*maxNELEM_PP
          CALL MPI_RECV(CMMAT_PP,SIZE10,MPI_REAL,0,&
     &        TAG10,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        ENDIF
      ENDIF    ! IF(MPI_RANK_P.NE.0)THEN !cncsskdjhfsdf

!
! *** NEXT DETERMINE ISIDE_PP,NX_PP,NY_PP,EL_PP AND ISCOM_PP
!
      IF(MPI_RANK_P.EQ.0)THEN !xvmnsdnmbfa
        ALLOCATE(ISIDE_PP_CP(8,maxNSIDE_PP))
        ALLOCATE(NX_PP_CP(maxNSIDE_PP))
        ALLOCATE(NY_PP_CP(maxNSIDE_PP))
        ALLOCATE(EL_PP_CP(maxNSIDE_PP))
        ALLOCATE(ISCOM_PP_CP(maxNSIDE_PP))
        ALLOCATE(LCOMM_PP_CP(maxNSIDE_PP))
        DO 10009 IG=1,NGRPS
          CALL IFILLM(ISIDE_PP_CP,8,maxNSIDE_PP,CO)
          CALL RFILLV(NX_PP_CP,maxNSIDE_PP,COO)
          CALL RFILLV(NY_PP_CP,maxNSIDE_PP,COO)
          CALL RFILLV(EL_PP_CP,maxNSIDE_PP,COO)
          CALL IFILLV(ISCOM_PP_CP,maxNSIDE_PP,CO)
          CALL IFILLV(LCOMM_PP_CP,maxNSIDE_PP,CO)
          TAG=IG*10
          IF(IG.NE.1)THEN
          CALL MPI_RECV(IPCOM_PP_CP,maxNPOIN_PP,MPI_INTEGER,IG-1,&
     &           TAG,MPI_COMM_P,MPI_STATUS,MPI_IERR)
          ELSE
            IPCOM_PP_CP(1:maxNPOIN_PP) = IPCOM_PP(1:maxNPOIN_PP)
          ENDIF
!
! ***     LOOP OVER THE ELEMENT EDGES
!
          FLAG=0
          NCOMM_PP_CP=0
          DO 1009 IS=1,NSIDE
            IEL=ISIDE(3,IS)
            IER=ISIDE(4,IS)
            IF(IEL.EQ.0)THEN ! This should not happen...
            RANKL=(-1);GOTO 1100
            ENDIF
            RANKL=ELGRP(IEL,1)
 1100       CONTINUE
            IF(IER.EQ.0)THEN
            RANKR=(-1);GOTO 1101
            ENDIF
            RANKR=ELGRP(IER,1)
 1101       CONTINUE
!
            IF((RANKL.EQ.IG).AND.(RANKR.EQ.IG))THEN!SIDE IS NOT AT A P-SPACE DOMAIN BOUNDARY
              FLAG=FLAG+1
              IF(FLAG.GT.maxNSIDE_PP)THEN
                WRITE(*,105)
                WRITE(*,106)
                STOP
              ENDIF   !OR AT A PROCESSOR DOMAIN BOUNDARY
              ISCOM_PP_CP(FLAG)=IS    !FILL IN ISCOM_PP_CP
              IP1=ISIDE(1,IS)
              DO 1010 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP1)THEN
                  ISIDE_PP_CP(1,FLAG)=I    !FILL IN ISIDE_PP_CP(1,*)
                  GOTO 1011
                ENDIF
 1010         CONTINUE
 1011         CONTINUE
              IP2=ISIDE(2,IS)
              DO 1012 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP2)THEN
                  ISIDE_PP_CP(2,FLAG)=I     !FILL IN ISIDE_PP_CP(2,*)
                  GOTO 1013
                ENDIF
 1012         CONTINUE
 1013         CONTINUE
              ISIDE_PP_CP(3,FLAG)=ELGRP(IEL,2)    !FILL IN ISIDE_PP_CP(3,*)
              ISIDE_PP_CP(4,FLAG)=ELGRP(IER,2)    !FILL IN ISIDE_PP_CP(4,*)
              DO I=5,8
                ISIDE_PP_CP(I,FLAG)=ISIDE(I,IS)!FILL IN ISIDE_PP_CP(5,_) ->> ISIDE(8,*)
              ENDDO       
!
! *** FILL   IN NX_PP_CP,NY_PP_CP,EL_PP_CP
!
              NX_PP_CP(FLAG)=NX(IS)
              NY_PP_CP(FLAG)=NY(IS)
              EL_PP_CP(FLAG)=EL(IS)
            ELSEIF((RANKL.EQ.IG).AND.(RANKR.NE.IG).AND.(RANKR.NE.-1))THEN
              NCOMM_PP_CP=NCOMM_PP_CP+1
              FLAG=FLAG+1
              IF(FLAG.GT.maxNSIDE_PP)THEN
                WRITE(*,105)
                WRITE(*,106)
                STOP
              ENDIF      !SIDE IS AT A PROCESSOR DOMAIN BOUNDARY
              ISCOM_PP_CP(FLAG)=IS
              LCOMM_PP_CP(FLAG)=RANKR
              IP1=ISIDE(1,IS)     !RHS ELEMENT IN OTHER PROCESSOR DOMAIN
              DO 1014 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP1)THEN
                  ISIDE_PP_CP(1,FLAG)=I     !FILL IN ISIDE_PP_CP(1,*)
                  GOTO 1015
                ENDIF
 1014         CONTINUE
 1015         CONTINUE
              IP2=ISIDE(2,IS)
              DO 1016 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP2)THEN
                  ISIDE_PP_CP(2,FLAG)=I      !FILL IN ISIDE_PP_CP(2,*)
                  GOTO 1017
                ENDIF
 1016         CONTINUE
 1017         CONTINUE
              ISIDE_PP_CP(3,FLAG)=ELGRP(IEL,2)    !FILL IN ISIDE_PP_CP(3,*)
              ISIDE_PP_CP(4,FLAG)=-1              !FILL IN ISIDE_PP_CP(4,*)
              DO I=5,6
                ISIDE_PP_CP(I,FLAG)=ISIDE(I,IS)  !FILL IN ISIDE_PP_CP(5,_) ->> ISIDE(8,*)
              ENDDO
              DO I=7,8
                ISIDE_PP_CP(I,FLAG)=-1
              ENDDO
!
! *** FILL   IN NX_PP_CP,NY_PP_CP,EL_PP_CP
!
              NX_PP_CP(FLAG)=NX(IS)
              NY_PP_CP(FLAG)=NY(IS)
              EL_PP_CP(FLAG)=EL(IS)
            ELSEIF((RANKR.EQ.IG).AND.(RANKL.NE.IG).AND.(RANKL.NE.-1))THEN
              FLAG=FLAG+1
              NCOMM_PP_CP=NCOMM_PP_CP+1
              IF(FLAG.GT.maxNSIDE_PP)THEN
                WRITE(*,105)
                WRITE(*,106)
                STOP
              ENDIF      !SIDE IS AT A PROCESSOR DOMAIN BOUNDARY
              ISCOM_PP_CP(FLAG)=IS
              LCOMM_PP_CP(FLAG)=RANKL
              IP1=ISIDE(1,IS)    !LHS ELEMENT IN OTHER PROCESSOR DOMAIN
              DO 1018 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP1)THEN
                  ISIDE_PP_CP(1,FLAG)=I    !FILL IN ISIDE_PP_CP(1,*)
                  GOTO 1019
                ENDIF
 1018         CONTINUE
 1019         CONTINUE
              IP2=ISIDE(2,IS)
              DO 1020 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP2)THEN
                  ISIDE_PP_CP(2,FLAG)=I    !FILL IN ISIDE_PP_CP(2,*)
                  GOTO 1021
                ENDIF
 1020         CONTINUE
 1021         CONTINUE
              ISIDE_PP_CP(3,FLAG)=-1    !FILL IN ISIDE_PP_CP(3,*)
              ISIDE_PP_CP(4,FLAG)=ELGRP(IER,2)    !FILL IN ISIDE_PP_CP(4,*)
              DO I=5,6
                ISIDE_PP_CP(I,FLAG)=-1   !FILL IN ISIDE_PP_CP(5,_) ->> ISIDE(8,*)
              ENDDO
              DO I=7,8
                ISIDE_PP_CP(I,FLAG)=ISIDE(I,IS)
              ENDDO
!
! *** FILL   IN NX_PP_CP,NY_PP_CP,EL_PP_CP
!
              NX_PP_CP(FLAG)=NX(IS)
              NY_PP_CP(FLAG)=NY(IS)
              EL_PP_CP(FLAG)=EL(IS)
            ELSEIF((RANKL.EQ.IG).AND.(RANKR.EQ.-1))THEN
              FLAG=FLAG+1
              IF(FLAG.GT.maxNSIDE_PP)THEN
                WRITE(*,105)
                WRITE(*,106)
                STOP
              ENDIF    !SIDE IS AT A P-SPACE DOMAIN BOUNDARY
              ISCOM_PP_CP(FLAG)=IS
              IP1=ISIDE(1,IS)   !RHS ELEMENT OUTSIDE DOMAIN
              DO 1022 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP1)THEN
                  ISIDE_PP_CP(1,FLAG)=I    !FILL IN ISIDE_PP_CP(1,*)
                  GOTO 1023
                ENDIF
 1022         CONTINUE
 1023         CONTINUE
              IP2=ISIDE(2,IS)
              DO 1024 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP2)THEN
                  ISIDE_PP_CP(2,FLAG)=I      !FILL IN ISIDE_PP_CP(2,*)
                  GOTO 1025
                ENDIF
 1024         CONTINUE
 1025         CONTINUE
              ISIDE_PP_CP(3,FLAG)=ELGRP(IEL,2)    !FILL IN ISIDE_PP_CP(3,*)
              ISIDE_PP_CP(4,FLAG)=0               !FILL IN ISIDE_PP_CP(4,*)
              DO I=5,6
                ISIDE_PP_CP(I,FLAG)=ISIDE(I,IS)     !FILL IN ISIDE_PP_CP(5,_) ->> ISIDE(8,*)
              ENDDO
              DO I=7,8
                ISIDE_PP_CP(I,FLAG)=0
              ENDDO
!
! *** FILL   IN NX_PP_CP,NY_PP_CP,EL_PP_CP
!
              NX_PP_CP(FLAG)=NX(IS)
              NY_PP_CP(FLAG)=NY(IS)
              EL_PP_CP(FLAG)=EL(IS)
            ELSEIF((RANKL.EQ.-1).AND.(RANKR.EQ.IG))THEN
              FLAG=FLAG+1
              IF(FLAG.GT.maxNSIDE_PP)THEN
                WRITE(*,105)
                WRITE(*,106)
                STOP
              ENDIF    !SIDE IS AT A P-SPACE DOMAIN BOUNDARY
              ISCOM_PP_CP(FLAG)=IS
              IP1=ISIDE(1,IS)    !LHS ELEMENT OUTSIDE DOMAIN
              DO 1026 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP1)THEN
                  ISIDE_PP_CP(1,FLAG)=I   !FILL IN ISIDE_PP_CP(1,*)
                  GOTO 1027
                ENDIF
 1026         CONTINUE
 1027         CONTINUE
              IP2=ISIDE(2,IS)
              DO 1028 I=1,maxNPOIN_PP
                IPT=IPCOM_PP_CP(I)
                IF(IPT.EQ.IP2)THEN
                  ISIDE_PP_CP(2,FLAG)=I    !FILL IN ISIDE_PP_CP(2,*)
                  GOTO 1029
                ENDIF
 1028         CONTINUE
 1029         CONTINUE
              ISIDE_PP_CP(3,FLAG)=0    !FILL IN ISIDE_PP_CP(3,*)
              ISIDE_PP_CP(4,FLAG)=ELGRP(IER,2)   !FILL IN ISIDE_PP_CP(4,*)
              DO I=5,6
                ISIDE_PP_CP(I,FLAG)=0    !FILL IN ISIDE_PP_CP(5,_) ->> ISIDE(8,*)
              ENDDO
              DO I=7,8
                ISIDE_PP_CP(I,FLAG)=ISIDE(I,IS)
              ENDDO
! 	          		
! ***         FILL   IN NX_PP_CP,NY_PP_CP,EL_PP_CP
!
              NX_PP_CP(FLAG)=NX(IS)
              NY_PP_CP(FLAG)=NY(IS)
              EL_PP_CP(FLAG)=EL(IS)
            ENDIF
! ***     END LOOP OVER SIDES
 1009     CONTINUE
          NSIDE_PP_CP=FLAG
!
! ***     MPI COMMAND TO SEND NSIDE_PP_CP,NX_PP_CP,NY_PP_CP,EL_PP_CP,ISIDE_PP_CP AND ISCOM_PP_CP TO PROCESSOR RANK IG!!!!!!!!!!!!!
!
          IF(IG.NE.1)THEN
            TAG11=11
            SIZE11=1
            TAG12=12
            SIZE12=maxNSIDE_PP
            TAG13=13
            SIZE13=maxNSIDE_PP
            TAG14=14
            SIZE14=maxNSIDE_PP
            TAG15=15
            SIZE15=8*maxNSIDE_PP
            TAG16=16
            SIZE16=maxNSIDE_PP
            TAG17=17
            SIZE17=maxNSIDE_PP
            TAG18=18
            SIZE18=1
            CALL MPI_SEND(NSIDE_PP_CP,SIZE11,MPI_INTEGER,IG-1,&
     &         TAG11,MPI_COMM_P,MPI_IERR)   
            CALL MPI_SEND(NX_PP_CP,SIZE12,MPI_REAL,IG-1,&
     &         TAG12,MPI_COMM_P,MPI_IERR)  
            CALL MPI_SEND(NY_PP_CP,SIZE13,MPI_REAL,IG-1,&
     &          TAG13,MPI_COMM_P,MPI_IERR)
            CALL MPI_SEND(EL_PP_CP,SIZE14,MPI_REAL,IG-1,&
     &          TAG14,MPI_COMM_P,MPI_IERR)
            CALL MPI_SEND(ISIDE_PP_CP,SIZE15,MPI_INTEGER,IG-1,&
     &       TAG15,MPI_COMM_P,MPI_IERR)
            CALL MPI_SEND(ISCOM_PP_CP,SIZE16,MPI_INTEGER,IG-1,&
     &           TAG16,MPI_COMM_P,MPI_IERR)
            CALL MPI_SEND(LCOMM_PP_CP,SIZE17,MPI_INTEGER,IG-1,&
     &           TAG17,MPI_COMM_P,MPI_IERR)
            CALL MPI_SEND(NCOMM_PP_CP,SIZE18,MPI_INTEGER,IG-1,&
     &           TAG18,MPI_COMM_P,MPI_IERR)
          ELSE
            NSIDE_PP = NSIDE_PP_CP
            NX_PP(1:NSIDE_PP_CP) = NX_PP_CP(1:NSIDE_PP_CP)
            NY_PP(1:NSIDE_PP_CP) = NY_PP_CP(1:NSIDE_PP_CP)
            EL_PP(1:NSIDE_PP_CP) = EL_PP_CP(1:NSIDE_PP_CP)
            ISIDE_PP(:,1:NSIDE_PP_CP) = ISIDE_PP_CP(:,1:NSIDE_PP_CP)
            ISCOM_PP(1:NSIDE_PP_CP) = ISCOM_PP_CP(1:NSIDE_PP_CP)
            LCOMM_PP(1:NSIDE_PP_CP) = LCOMM_PP_CP(1:NSIDE_PP_CP)
            NCOMM_PP = NCOMM_PP_CP
          ENDIF

10009   CONTINUE
        DEALLOCATE(NX_PP_CP)
        DEALLOCATE(NY_PP_CP)
        DEALLOCATE(EL_PP_CP)
        DEALLOCATE(ISIDE_PP_CP)
        DEALLOCATE(ISCOM_PP_CP)
        DEALLOCATE(LCOMM_PP_CP)
      ENDIF   !IF(MPI_RANK_P.EQ.0)THEN !xvmnsdnmbfa

      IF(MPI_RANK_P.NE.0)THEN !ddpddpddpd
        TAG=(MPI_RANK_P+1)*10
        CALL MPI_SEND(IPCOM_PP,maxNPOIN_PP,MPI_INTEGER,&
     &    0,TAG,MPI_COMM_P,MPI_IERR)
        TAG11=11
        SIZE11=1    
        CALL MPI_RECV(NSIDE_PP,SIZE11,MPI_INTEGER,0,&
     &    TAG11,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        TAG12=12
        SIZE12=maxNSIDE_PP
        TAG13=13
        SIZE13=maxNSIDE_PP
        TAG14=14
        SIZE14=maxNSIDE_PP
        CALL MPI_RECV(NX_PP,SIZE12,MPI_REAL,0,&
     &    TAG12,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        CALL MPI_RECV(NY_PP,SIZE13,MPI_REAL,0,&
     &    TAG13,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        CALL MPI_RECV(EL_PP,SIZE14,MPI_REAL,0,&
     &    TAG14,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        TAG15=15
        SIZE15=8*maxNSIDE_PP
        TAG16=16
        SIZE16=maxNSIDE_PP
        TAG17=17
        SIZE17=maxNSIDE_PP
        TAG18=18
        SIZE18=1
        CALL MPI_RECV(ISIDE_PP,SIZE15,MPI_INTEGER,0,&
     &    TAG15,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        CALL MPI_RECV(ISCOM_PP,SIZE16,MPI_INTEGER,0,&
     &    TAG16,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        CALL MPI_RECV(LCOMM_PP,SIZE17,MPI_INTEGER,0,&
     &    TAG17,MPI_COMM_P,MPI_STATUS,MPI_IERR)
        CALL MPI_RECV(NCOMM_PP,SIZE18,MPI_INTEGER,0,&
     &    TAG18,MPI_COMM_P,MPI_STATUS,MPI_IERR)
      ENDIF !  IF(MPI_RANK_P.NE.0)THEN !ddpddpddpd

!
! ***   CALL THE ROUTINE FOR CONSTRUCTING THE EDGFLX COMMUNICATION MATRICES
! 
      CALL EDGCOM(NGRPS,MPI_RANK_P,NCOMM_PP,NSIDE_PP,SDCOM_PP,ISIDE_PP,&
     &               LCOMM_PP,ISCOM_PP,MPI_COMM_P)
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
!
! *** FORMAT STATEMENTS
!
 101  FORMAT('NEED TO INCREASE maxNPOIN_pp IN DG_LIN_CONVEC.f!')
 102  FORMAT('PROGRAM STOPPED')    
 103  FORMAT('NEED TO INCREASE maxNBOUN_pp IN DG_LIN_CONVEC.f!')
 104  FORMAT('PROGRAM STOPPED')  
 105  FORMAT('NEED TO INCREASE maxNSIDE_pp IN DG_LIN_CONVEC.f!')
 106  FORMAT('PROGRAM STOPPED')
!
      RETURN
      END
