      SUBROUTINE DG_lin_convec(NDIMN ,NNODE ,NPOIN ,NELEM ,& 
     &                     NBOUN ,INTMA ,COORD ,& 
     &                     NTIME ,CSAFM ,& 
     &                     IMMAT ,NBNOI ,BSIDO ,& 
     &                     NBNOR,NGEOM,RORDER,TORDER,& 
     &                     VNPNT,SUMWEIGHT,VCORD,& 
     &                     rv,RHO,UVEL,VVEL,PS,& 
     &                     TEMP,ALPHA,CINF,MXSID,NPROC,IVD,FORCEOUT,&
     &                     d,R,M,&
     &                     MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
     &                     MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &                     VSPACE_FIRST,VSPACE_LAST)

! 
! *** KEYWORDS :    - DISCONTINUOUS GALERKIN   
!				    - SECOND ORDER TAYLOR GALERKIN 
!                   - TWO-STEP SCHEME 
!                   - LINEAR FEM 
!				    - GLOBAL TIMESTEPPING 
      IMPLICIT NONE
      INCLUDE 'mpif.h' 
!
! CREATE STRUCTURE FOR INPUT VARIABLES
!
        TYPE InputVariables
          INTEGER :: TORDER !Quadrature order in Theta for V-SPACE
          INTEGER :: NTIME !Number of timesteps
          INTEGER :: FORCEOUT !output forces at each timestep? 
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
!     DECLARE DYNAMIC ARRAYS
      INTEGER, ALLOCATABLE :: ELGRP(:,:)
      REAL, ALLOCATABLE :: MMAT(:,:)
      REAL, ALLOCATABLE :: DISNF_PP(:,:,:)
      INTEGER, ALLOCATABLE :: NPGRP(:)
      INTEGER, ALLOCATABLE :: NEGRP(:)
      INTEGER, ALLOCATABLE :: GCOMM(:,:)
!
      TYPE(InputVariables) :: IVD
! 
! *** DECLARE INTEGER VARIABLES 
      INTEGER NSIDE, NDIMN,NNODE,NBNOR,NBNOI,NPROC,NCOMM_PP 
      INTEGER NGEOM, NELEM,NPOIN,NBOUN,NTIME,MXCOM_PP 
      INTEGER MXSID,MPI_IERR,IMMAT,VNPNT,RORDER,TORDER
      INTEGER IP1,IP2,IP3,FORCEOUT
      INTEGER maxNELEM_PP,maxNBOUN_PP,maxNSIDE_PP,maxNPOIN_PP 
      INTEGER NELEM_PP,NSIDE_PP,NPOIN_PP,NBOUN_PP 
      ! mpi-stuff for position space partitioning
      INTEGER MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,GROUP_P
      ! mpi-stuff for velocity space partitioning
      INTEGER MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,GROUP_V
      INTEGER VNPNT_PART,VSPACE_FIRST,VSPACE_LAST
      INTEGER IG
      



! *** PARAMETERS 
      PARAMETER (maxNELEM_PP=30000,maxNBOUN_PP=2400,maxNPOIN_PP=16000) 
      PARAMETER (maxNSIDE_PP=81000) 
! *** DECLARE REAL VARIABLES MMAT (MASS MATRIX), UX,UY,CSAFM 
      REAL    UX, UY, CSAFM,ALPHA,DTE 
      REAL    UVEL(NPOIN),VVEL(NPOIN),RHO(NPOIN),PS(NPOIN) 
      REAL    COORD(2,NPOIN),GEOME(NGEOM,NELEM),RSIDO(3,NBOUN) 
      REAL    CMMAT(3,3,NELEM),RELEN(NELEM)      
      REAL    rv,CINF(4),TEMP(NPOIN),VCORD(3,VNPNT) 
      REAL    RHO_PP(maxNPOIN_pp),UVEL_PP(maxNPOIN_PP) 
      REAL    TEMP_PP(maxNPOIN_PP),NX_PP(maxNSIDE_PP) 
      REAL    EL_PP(maxNSIDE_PP) 
      REAL    COORD_PP(2,maxNPOIN_PP),PS_PP(maxNPOIN_PP) 
      REAL    RSIDO_PP(NBNOR,maxNBOUN_PP),MMAT_PP(NNODE,maxNELEM_PP) 
      REAL    CMMAT_PP(3,3,maxNELEM_PP),GEOME_PP(NGEOM,maxNELEM_PP) 
      REAL    DELUN_PP(1,3,maxNELEM_pp),UMEAN_PP(maxNELEM_PP) 
      REAL    NX(mxsid),NY(mxsid) 
      REAL    FLUXP_PP(1,3,maxNELEM_pp),FLUYP_PP(1,3,maxNELEM_pp)  
      REAL    EL(mxsid),VVEL_PP(maxNPOIN_pp),NY_PP(maxNSIDE_PP) 
      REAL    TARRAY(2),TIME,SUMWEIGHT,d,R,M 
! *** DECLARE INTEGER ARRAYS 
      INTEGER ISIDE(8,MXSID)
      INTEGER ISCOM_PP(maxNSIDE_pp)
      INTEGER INTMA(3,NELEM) ,BSIDO(NBNOI,NBOUN)
      INTEGER IBCOM_PP(maxNBOUN_pp) 
      INTEGER LCOMM_PP(maxNSIDE_PP),IPCOM_PP(maxNPOIN_PP) 
      INTEGER INTMA_PP(3,maxNELEM_PP),BSIDO_PP(NBNOI,maxNBOUN_PP) 
      INTEGER ISIDE_PP(8,maxNSIDE_PP),SDCOM_PP(3,maxNSIDE_PP) 
! *** CHARACTER VARIABLES
      CHARACTER filename*80
! *** VARIABLES FROM THE GETSID SUBROUTINE 
      INTEGER LWHER(NPOIN),LHOWM(NPOIN),ICONE(3*NELEM)  
      DOUBLE PRECISION OUTPUT_TIMING



! *** ALLOCATE MEMORY


      VNPNT_PART=VNPNT/MPI_SIZE_V
      ALLOCATE(ELGRP(NELEM,2),MMAT(NNODE,NELEM))
      ALLOCATE(DISNF_PP(NNODE,VSPACE_FIRST:VSPACE_LAST,maxNELEM_PP))
      ALLOCATE(NPGRP(MPI_SIZE_P),NEGRP(MPI_SIZE_P),GCOMM(MPI_SIZE_P,MPI_SIZE_P))
!      
      IF(MPI_RANK_P.EQ.0)THEN   
      ! GET THE ELEMENT SIDES DATA (ROWS 1-4 IN ISIDE) 
          CALL GETSID(NELEM,NPOIN,NSIDE,INTMA,ISIDE,LWHER,& 
     &                      LHOWM,ICONE,MXSID)
      ENDIF
      !BROADCAST NSIDE TO ALL PROCS
      CALL MPI_BCAST(NSIDE,1,MPI_INTEGER,0,MPI_COMM_P,MPI_IERR)


      IF(MPI_RANK_P.EQ.0)THEN ! sfasdasfsda
! ***   GET THE NORMALS AND LENGTHS OF THE ELEMENT SIDES 
        CALL NORMAL(ISIDE,COORD,NSIDE,NPOIN,NX,NY,EL) 
! ***   OBTAIN THE GEOMETRY-ARRAYS 
        CALL GETGEO(NELEM ,NPOIN ,NNODE ,NDIMN ,NGEOM ,INTMA ,COORD ,& 
     &              GEOME ) 
! ***   GET THE REPRESENTATIVE ELEMENT/POINT LENGTHS 
        CALL GETLEN(NNODE ,NELEM ,RELEN ,& 
     &              NGEOM ,GEOME )
! ***   FIND THE ELEMENT LUMPED MASS MATRICES 
        IF(IMMAT.EQ.1)THEN 
               CALL GETMAT(NELEM ,NNODE ,NGEOM ,GEOME ,MMAT ) 
         ELSE 
               CALL CONMAT(NELEM,NGEOM,GEOME,CMMAT) 
         ENDIF
! ***   GET THE NORMALS OF THE MESH 
        CALL GETNOR(NDIMN ,NBNOI ,NBNOR ,NPOIN ,NBOUN ,BSIDO ,RSIDO ,& 
     &              COORD,MPI_RANK_P ) 
      ENDIF !  IF(MPI_RANK_P.EQ.0)THEN ! sfasdasfsda
! 
! *** CONVERT FROM METIS NODAL PARTITION TO ELEMENT PARTITION 
! 
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 

!     The PTMESH routine does not care about velocity space and velocity
!     space partitioning.
      CALL PTMESH(NELEM,INTMA,NPOIN,MPI_SIZE_P,NSIDE,ISIDE,& 
     &        NELEM_PP,NPOIN_PP,NBOUN_PP,NSIDE_PP,NEGRP,& 
     &      LCOMM_PP,GCOMM,MXCOM_PP,ELGRP,& 
     &      NGEOM,NBNOI,NBOUN,NBNOR,NNODE,NX,NY,EL,GEOME,MMAT,& 
     &      CMMAT,COORD,BSIDO,RSIDO,maxNSIDE_PP,maxNBOUN_PP,& 
     &       maxNELEM_PP,maxNPOIN_PP,NX_PP,NY_PP,EL_PP,& 
     &       GEOME_PP,MMAT_PP,CMMAT_PP,INTMA_PP,COORD_PP,BSIDO_PP,& 
     &       RSIDO_PP,IMMAT,MPI_RANK_P,IPCOM_PP,IBCOM_PP,ISCOM_PP,NPGRP,& 
     &       ISIDE_PP,NCOMM_PP,SDCOM_PP,IVD,MPI_COMM_P)

!
! *** SET UP THE INITIAL CONDITIONS OF THE PROBLEM (FOR A GAS EXPANSION) 
!
      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN ! write only on master rank
        WRITE(*,*) "Initialization..."
      ENDIF
!      CALL INICON(NPOIN,NNODE,NELEM,NELEM_PP,INTMA_PP,ELGRP,&
!     &NPOIN_PP,COORD_PP,VNPNT,VCORD,DISNF_PP,IPCOM_PP,rv,&
!     &                 IVD,&
!     &                 MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
!     &                 MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
!     &                 VSPACE_FIRST,VSPACE_LAST) 
!
!      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
! 
! *** FOR A UNIFORM INITIAL CONDITION 
! 

      CALL INICON2(NPOIN,NNODE,NELEM,NELEM_PP,INTMA_PP,&
     &    ELGRP,NPOIN_PP,VNPNT,VCORD,DISNF_PP,IPCOM_PP,CINF,rv,&
     &    IVD,R,M,&
     &                 MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
     &                 MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &                 VSPACE_FIRST,VSPACE_LAST) 

! *** SYNCHRONISE 
! 
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
! 
! *** PERFORM THE TIME STEP ITERATIONS 
!
      CALL ADVNCE(NNODE ,NGEOM ,NBNOI , NBOUN, & 
     &        NBNOR , COORD, BSIDO, NBOUN_PP, NPOIN, NPOIN_PP ,NELEM_PP,NELEM,& 
     &                  NTIME ,BSIDO_PP ,INTMA_PP ,& 
     &                  FLUXP_PP ,FLUYP_PP ,& 
     &                  GEOME_PP ,MMAT_PP  ,CMMAT_PP, RELEN ,DTE ,& 
     &                  DELUN_PP ,& 
     &                  RSIDO_PP ,UMEAN_PP,& 
     &                  IMMAT ,CSAFM ,UX,UY,& 
     &                  NSIDE_PP,ISIDE_PP,NX_PP,NY_PP,EL_PP,DISNF_PP,& 
     &                  VNPNT,SUMWEIGHT,DISNF_PP,rv,VCORD,& 
     &                  RHO_PP,UVEL_PP,VVEL_PP,PS_PP,TEMP_PP,& 
     &                  ALPHA,CINF,LCOMM_PP,GCOMM,NPGRP,IPCOM_PP,& 
     &      ISCOM_PP,MXCOM_PP,NCOMM_PP,maxNELEM_PP,RORDER,TORDER,maxNPOIN_PP,&
     &          SDCOM_PP,IVD,FORCEOUT,d,R,M,&
     &                     MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
     &                     MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &                     VSPACE_FIRST,VSPACE_LAST)
! 
! *** OUTPUT FOR GID AND RESTART FILE 
! 
      OUTPUT_TIMING = MPI_WTIME()
!      CALL OUTPUT(NNODE ,NPOIN_PP ,NPOIN,NELEM ,NELEM_PP,&
!     &      maxNELEM_PP,COORD,INTMA, DISNF_PP, VNPNT,IPCOM_PP,&
!     &                     ELGRP,NEGRP,IVD,&
!     &                     MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
!     &                     MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
!     &                     VSPACE_FIRST,VSPACE_LAST,VCORD,rv)
      OUTPUT_TIMING = MPI_WTIME() - OUTPUT_TIMING
      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN ! write only on master rank
          WRITE(*,"(A20,F6.2,A3)") "OUTPUT written in ", OUTPUT_TIMING,&
     &                             "sec"
      ENDIF
!
       RETURN 
       END 
