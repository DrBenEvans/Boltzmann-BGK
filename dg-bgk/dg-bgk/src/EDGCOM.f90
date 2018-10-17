      SUBROUTINE EDGCOM(NGRPS,MPI_RANK_P,NCOMM_PP,NSIDE_PP,SDCOM_PP,&
     &           ISIDE_PP,LCOMM_PP,ISCOM_PP,MPI_COMM_P)
!
! *** THIS ROUTINE CONSTRUCTS THE SDCOM COMMUNICATION MATRICES TO BE USED 
! *** IN EDGFLX.f90
! 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!
      INTEGER IG,IS,IST,NGRPS,NCOMM_PP,I
      INTEGER MPI_IERR,MPI_RANK_P,MPI_COMM_P
      INTEGER SDCOM_PP(3,NCOMM_PP),NSIDE_PP,TMPNSIDE_PP,ISTEST
      INTEGER OPPPART,ISCOM_PP(NSIDE_PP),TAG,IEL,IER,GLOBIS
      INTEGER ISIDE_PP(8,NSIDE_PP),FLAG,LCOMM_PP(NSIDE_PP)
      INTEGER MPI_STATUS(MPI_STATUS_SIZE)
!
      CALL MPI_BARRIER(MPI_COMM_P,MPI_IERR)
!
! *** LOOP OVER EACH GROUP
!
      DO 1000 IG=1,NGRPS
! ***   INITIALISE SDCOM
        IF((MPI_RANK_P+1).EQ.IG)THEN
          CALL IFILLM(SDCOM_PP,3,NCOMM_PP,0)
        ENDIF
! ***   LOOP OVER EACH SIDE
        IF((MPI_RANK_P+1).EQ.IG)TMPNSIDE_PP=NSIDE_PP
        CALL MPI_BCAST(TMPNSIDE_PP,1,MPI_INTEGER,IG-1,MPI_COMM_P,&
     &                   MPI_IERR)
        IF((MPI_RANK_P+1).EQ.IG)FLAG=0
!
        DO 1100 I=1,TMPNSIDE_PP
          IF((MPI_RANK_P+1).EQ.IG)THEN
            IEL=ISIDE_PP(3,I)
            IER=ISIDE_PP(4,I)
            IF((IEL.EQ.-1).OR.(IER.EQ.-1))THEN ! it's an edge
                                               ! shared between ranks
              FLAG=FLAG+1
              OPPPART=LCOMM_PP(I) ! opposite partition for shared edge
                                  ! that is, rank+1
              GLOBIS=ISCOM_PP(I)  ! global index for shared edge
            ENDIF
          ENDIF
          CALL MPI_BCAST(IEL,1,MPI_INTEGER,IG-1,MPI_COMM_P,MPI_IERR)
          CALL MPI_BCAST(IER,1,MPI_INTEGER,IG-1,MPI_COMM_P,MPI_IERR)
          IF((IEL.EQ.-1).OR.(IER.EQ.-1))THEN
            CALL MPI_BCAST(OPPPART,1,MPI_INTEGER,IG-1,MPI_COMM_P,MPI_IERR)
            CALL MPI_BCAST(GLOBIS,1,MPI_INTEGER,IG-1,MPI_COMM_P,MPI_IERR)
            IF((MPI_RANK_P+1).EQ.OPPPART)THEN
              DO IST=1,NSIDE_PP   ! search ISCOM_PP for 
                                  ! global edge
                ISTEST=ISCOM_PP(IST)
                IF(ISTEST.EQ.GLOBIS)GOTO 999 ! found, "EXIT"
              ENDDO
 999          CONTINUE
              TAG=100
              CALL MPI_SEND(IST,1,MPI_INTEGER,IG-1,&
     &        TAG,MPI_COMM_P,MPI_IERR)
            ENDIF
            IF((MPI_RANK_P+1).EQ.IG)THEN
              TAG=100
              CALL MPI_RECV(IST,1,MPI_INTEGER,OPPPART-1,&
     &            TAG,MPI_COMM_P,MPI_STATUS,MPI_IERR)
              SDCOM_PP(1,FLAG)=I 
              SDCOM_PP(2,FLAG)=IST
              SDCOM_PP(3,FLAG)=OPPPART
            ENDIF
          ENDIF
          CALL MPI_BARRIER(MPI_COMM_P,MPI_IERR)
 1100   CONTINUE ! END LOOP OVER SIDES IN THE GROUP
!
        CALL MPI_BARRIER(MPI_COMM_P,MPI_IERR)
 1000 CONTINUE ! END LOOP OVER GROUPS
!
      RETURN
      END
