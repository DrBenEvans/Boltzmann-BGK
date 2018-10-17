      SUBROUTINE OUTPUT(NNODE,NPOIN_PP,NPOIN,NELEM,NELEM_PP,&
     &              maxNELEM_PP,COORD,INTMA , DISNF_PP, VNPNT,IPCOM_PP,&
     &                     ELGRP,NEGRP,IVD,&
     &                     MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
     &                     MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &                     VSPACE_FIRST,VSPACE_LAST,VCORD,rv)  
! 
		IMPLICIT NONE 
        INCLUDE 'mpif.h' 
!
! CREATE STRUCTURE FOR INPUT VARIABLES
!
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
          REAL :: ALPHA !WALL MOLECULAR REFLECTION PARAMETER
          REAL :: R !GAS CONSTANT
          REAL :: d !MOLECULAR DIAMETER
          REAL :: M !MOLAR MASS
          CHARACTER :: LobattoFile*80
          CHARACTER :: PSpaceFile*80
          CHARACTER :: OutFile*80
          CHARACTER :: PartitionFile*80
          CHARACTER :: RestartInFile*80
          CHARACTER :: ResidualFile*80
          CHARACTER :: ResultsFile1*80
          CHARACTER :: ResultsFile2*80
          CHARACTER :: RestartOutFile*80
          CHARACTER :: GIDMeshFile*80
        END TYPE
! 
	  INTEGER VNPNT,NPOIN,NELEM,NNODE,IP,IN,NPART,NELEM_PP 
      INTEGER maxNELEM_PP  
      INTEGER IPT,J,IE,MPI_IERR,MPI_STATUS(MPI_STATUS_SIZE) 
      INTEGER IE_PP,TAG 
      INTEGER NPOIN_PP,IPCOM_PP(NPOIN_PP),IV,IVTRUE,ELGRP(NELEM,2) 
      REAL COORD(2,NPOIN) 
      INTEGER NEGRP(MPI_SIZE_P)

! 
      INTEGER INTMA(NNODE,NELEM),PRINT_EVERY
      ! mpi-stuff for position space partitioning
      INTEGER MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,GROUP_P
      ! mpi-stuff for velocity space partitioning
      INTEGER MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,GROUP_V,IVP
      INTEGER VNPNT_PART,VSPACE_FIRST,VSPACE_LAST,IG,TSIZE
      REAL DISNF_PP(NNODE,VSPACE_FIRST:VSPACE_LAST,maxNELEM_PP)
      REAL DISNF(NNODE,VSPACE_FIRST:VSPACE_LAST,NELEM)

      REAL ETA,ZETA,WEIGHT,RT,THETA,CX,CY,VCORD(3,VNPNT),rv,PI
      PARAMETER (PI=3.1416)
!
      TYPE(InputVariables) :: IVD
! *** DECLARE THE CHARACTER VARIABLE FILENAME 
      CHARACTER filename*80 
! 
! *** OUTPUT FOR RESTART. 
! 
! *** SET UP CHANNEL TO WRITE TO RESTART FILE 
! 
      PRINT_EVERY = (VSPACE_LAST+1-VSPACE_FIRST)/5
      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN 
        WRITE(*,*) 
        WRITE(*,*) 
        filename = IVD%RestartOutFile
        PRINT*,'OPENING NEW RESTART FILE = ',IVD%RestartOutFile 
        OPEN  (20, file=filename, status='UNKNOWN')  
    
        WRITE(20,185)  
        WRITE(20,*)  NPOIN,'    ',VNPNT  
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
    
      VNPNT_PART = VSPACE_LAST-VSPACE_FIRST+1

      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN 
        WRITE(*,*)"Sending DISNF data to master in MPI_COMM_P..."
      ENDIF
      DO IG=1,MPI_SIZE_P !
        TSIZE = NEGRP(IG)*3*VNPNT_PART
        IF((MPI_RANK_P+1).EQ.IG)THEN
          CALL MPI_SEND(DISNF_PP(:,VSPACE_FIRST:VSPACE_LAST,:),&
     &            TSIZE,MPI_REAL,0,IG,&
     &            MPI_COMM_P,MPI_IERR)
        ENDIF
        IF(MPI_RANK_P.EQ.0)THEN
          IF(IG.NE.1) THEN
            CALL MPI_RECV(DISNF_PP(:,VSPACE_FIRST:VSPACE_LAST,:),&
     &              TSIZE,MPI_REAL,IG-1,IG,&
     &              MPI_COMM_P,MPI_STATUS,MPI_IERR)
          ENDIF
          ! Copying DISNF_PP into DISNF    
          DO IE=1,NELEM ! scanning through DISNF
            ! Checking if we have received the necessary data 
            ! in this iteration of the DO IG,MPI_SIZE_P group
            NPART = ELGRP(IE,1) 
            IF(NPART.EQ.IG)THEN
              IE_PP=ELGRP(IE,2)
              DISNF(:,:,IE) = DISNF_PP(:,:,IE_PP)
            ENDIF
          ENDDO
        ENDIF
      ENDDO


      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN 
        WRITE(*,*)"Sending DISNF data to general master..."
      ENDIF
    
      ! only the master ranks in the MPI_COMM_P communicators
      IF(MPI_RANK_P.EQ.0)THEN !jdfcdknsdasdjkha
        TSIZE =  NELEM*3*VNPNT_PART
        DO IVP=0,MPI_SIZE_V-1 !dvkjdkjdjkshasdjkh
          IF((IVP.EQ.MPI_RANK_V).AND.(IVP.NE.0)) THEN
              CALL MPI_SEND(DISNF,TSIZE,MPI_REAL,0,1,&
     &             MPI_COMM_V,MPI_IERR)
          ENDIF          
          IF(0.EQ.MPI_RANK_V) THEN !dsaaanccmna

            IF(IVP.NE.0) CALL MPI_RECV(DISNF,TSIZE,MPI_REAL,IVP,1,&
     &             MPI_COMM_V,MPI_STATUS,MPI_IERR)

            WRITE(*,"(A41,I2)")&
                 &"Writing to file data from vspace part no.", IVP
            DO IV=VSPACE_FIRST,VSPACE_LAST ! 1 to VNPNT_PART for rank 0
              IVTRUE=IV+IVP*VNPNT_PART
              ETA=VCORD(1,IVTRUE) 
              ZETA=VCORD(2,IVTRUE) 
              WEIGHT=VCORD(3,IVTRUE)
              RT=ETA*(rv/2)+(rv/2)!MAP BACK TO POLAR COORDINATES 
              THETA=ZETA*PI
              CX=RT*COS(THETA)!CONVERT TO CARTESIANS 
              CY=RT*SIN(THETA) 
              IF((MOD(IV,PRINT_EVERY)).EQ.0)THEN 
                WRITE(*,"(A23,I5)"),&
     &                   "V-space iteration no.", IV
              ENDIF

              DO IE=1,NELEM
                WRITE(20,"(I5,5E15.7)")IE,CX,CY,(DISNF(IN,IV,IE),IN=1,NNODE)
              ENDDO
            ENDDO
          ENDIF ! IF(0.EQ.MPI_RANK_V) THEN !dsaaanccmna
        ENDDO ! DO IVP=0,MPI_SIZE_V-1 !dvkjdkjdjkshasdjkh


 
        IF(MPI_RANK_V.EQ.0)THEN !lkjfsdccddna
          WRITE(*,*)"Written DISNF data to file."
          CLOSE(20) 
! ***     OUTPUT FOR GID MESHFILE 
! ***     SET UP CHANNEL 
          filename = IVD%GIDMeshFile 
          PRINT*,'OPENING GIDMESH FILE = ',IVD%GIDMeshFile 
          OPEN  (19, file=filename, status='unknown') 
          WRITE(*,*) 
! 
          WRITE(19,210)  
          WRITE(19,300) 
          DO IP=1,NPOIN 
            WRITE(19,400)IP,(COORD(J,IP),J=1,2) 
          ENDDO
          WRITE(19,220) 
! 
          WRITE(19,230) 
          DO IE=1,NELEM 
            WRITE(19,200)IE,(INTMA(J,IE),J=1,NNODE) 
          ENDDO
          WRITE(19,240) 
! 
          CLOSE(19) 
          WRITE(*,*) 
          WRITE(*,*) 
        ENDIF ! IF(MPI_RANK_V.EQ.0)THEN !lkjfsdccddna
      ENDIF ! IF(MPI_RANK_P.EQ.0)THEN !jdfcdknsdasdjkha
! 
! *** format statements 
!	 
  200 FORMAT(10I8) 
  300 FORMAT('   COORDINATES') 
  400 FORMAT(I10,2(2X,F10.5)) 
  180 FORMAT('WHAT WOULD YOU LIKE THE RESTART FILE TO BE CALLED?& 
     & (____________.RES) :    ', $) 
  185 FORMAT('NPOIN     VNPNT') 
  190 FORMAT('WHAT WOULD YOU LIKE THE GID MESHFILE TO BE CALLED?& 
     & (____________.RES) :    ', $) 
  210 FORMAT('mesh dimension = 2 elemtype triangle nnode = 3') 
  220 FORMAT('end coordinates') 
  230 FORMAT('elements') 
  240 FORMAT('end elements') 
! 
      RETURN 
      END 
