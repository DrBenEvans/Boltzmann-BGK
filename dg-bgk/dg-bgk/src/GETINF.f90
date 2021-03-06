       SUBROUTINE GETINF(ITIME,CINF) 
! 
! *** THIS SUBROUTINE CALCULATES THE UNKNOWN VALUE APPLIED AT THE INFLOW BOUNDARY 
!		(SET AT THE MO FOR A UNIFORM INITIAL CONDITION OF 1.0) 
      INTEGER ITIME 
      REAL CINF(3) 
!       
       IF(ITIME.LT.20) THEN 
       CINF(1)=1.0 
       ELSEIF((ITIME.GE.20).AND.(ITIME.LT.30))THEN 
       CINF(1)=1.0+(ITIME-20)*((5.0-1.0)/10.0) 
       ELSEIF((ITIME.GE.30).AND.(ITIME.LT.50))THEN 
       CINF(1)=5.0 
       ELSEIF((ITIME.GE.50).AND.(ITIME.LT.60))THEN 
       CINF(1)=5.0-(ITIME-50)*((5.0-1.0)/10.0) 
       ELSE 
       CINF(1)=1.0 
       ENDIF 
! 
! 
      PRINT*,'ITIME =',ITIME 
      PRINT*,'CINF(1) =',CINF(1) 
      RETURN 
      END 
