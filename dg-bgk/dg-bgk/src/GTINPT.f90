        SUBROUTINE GTINPT(NDIMN ,NNODE ,NPOIN ,NELEM ,& 
     &                  NBOUN ,INTMA ,COORD , BSIDO ,& 
     &                  IELSI , NBNOI,MPI_RANK_V) 
! 
      REAL COORD(NDIMN,NPOIN)   
      REAL TEXT(20) 
! 
      INTEGER INTMA(NNODE,NELEM),NNODE, NBNOI 
      INTEGER NELEM, NPOIN, NBOUN 
      INTEGER BSIDO(NBNOI,NBOUN),IELSI(2,NELEM) 
      INTEGER MPI_RANK_V
! 
! *** READ ALL THE INPUT DATA IN THE FOLLOWING MANNER : 
! 
!     1. CONNECTIVITY MATRIX : INTMA 
!     2. COORDINATES OF THE NODES : COORD 
!     3. BOUNDARY CONDITIONS AND SIDES : BSIDO 
!     4. CONSTANTS OF THE PROBLEM:NIN,NTIME,IMMAT,CSAFM 
! 
! *** 1. CONNECTIVITY MATRIX : INTMA  
!  
      READ(14,1) TEXT
      WRITE(*,"(A3,I2,A31)"),'VSR',MPI_RANK_V,&
     &                ': Reading connectivity matrix' 
      DO 1010 I=1,NELEM 
        READ(14,*) IELEM,(INTMA(J,I),J=1,NNODE) 
 1010 CONTINUE 
! 
! *** 2. COORDINATES 
! 
      READ(14,1) TEXT 
      WRITE(*,"(A3,I2,A31)"),'VSR',MPI_RANK_V,&
     &           ': Reading P-space coordinates'
      DO 1020 I=1,NPOIN 
        READ(14,*) IPOIN,(COORD(J,I),J=1,NDIMN) 
 1020 CONTINUE 
! 
! *** 3. BOUNDARY SIDES : A) NODES ON THE SIDE  
!                         B) ELEMENT CONTAINING THE SIDE 
!                         C) NATURE OF THE BOUNDARY SIDE
!                         D) TEMP IF WALL BOUNDARY
!                         E) ALPHA x 10 IF WALL BOUNDARY 
      READ(14,1)TEXT 
      WRITE(*,"(A3,I2,A31)"),'VSR',MPI_RANK_V,&
     &              ': Reading boundary data'
      DO 1070 I=1,NBOUN 
        READ(14,*)(BSIDO(J,I),J=1,NBNOI) 
 1070 CONTINUE 
! 
!     READ(14,1) TEXT 
! *** NIN = NUMBER OF INITIAL CONDITIONS 
!      READ(14,*) NIN 
!      PRINT*,'Number of initial conditions=',NIN
! *** READ IN NTIME,IMMAT 
! 
!      READ(14,1) TEXT 
!      READ(14, *) NTIME,IMMAT,ALPHA 
!      PRINT*,'NTIME=',NTIME,'IMMAT=',IMMAT,'ALPHA=',ALPHA
! 
! *** 6.2 READ IN :CSAFM (SAFETY FACTOR) 
!      READ(14,1) TEXT 
!      READ(14, *) CSAFM    
!      PRINT*,'Safety factor = ',CSAFM
! 
! *** FILL IN THE ELEMENT-SIDES CONNECTIVITY MATRIX 
! *** FILLS IELSI(2,NELEM) WITH ZEROS 
      CALL IFILLM(IELSI,2,NELEM,0) 
! 
      DO 3070 I=1,NBOUN 
        IE=BSIDO(3,I) 
        IF(IELSI(1,IE).EQ.0) IELSI(1,IE)=I 
        IF(IELSI(1,IE).NE.I) IELSI(2,IE)=I 
 3070 CONTINUE 
! 
! *** CLOSE CHANNEL 5 
! 
      CLOSE(5) 
! 
! ***-----FORMATS 
! 
    1 FORMAT(20A4) 
! 
      RETURN 
      END 
