      SUBROUTINE ALOTIM(NELEM ,CSAFM ,& 
     &                  RELEN ,DTE,MXSPD ) 
! 
      INTEGER NELEM 
      REAL CSAFM 
      REAL DELTE(NELEM) 
      REAL RELEN(NELEM) 
      REAL DTE,MXSPD 
! 
! *** ALLOWABLE TIMESTEP AT ELEMENTS. 
! 
      DO 2000 IE=1,NELEM 
      ALEN=RELEN(IE) 
! 
      DELTE(IE)=(CSAFM*ALEN)/(MXSPD) 
 2000 CONTINUE 
! 
      DTE=MINVAL(DELTE) 
! 
      RETURN 
      END 