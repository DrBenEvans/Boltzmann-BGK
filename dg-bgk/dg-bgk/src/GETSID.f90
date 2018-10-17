      SUBROUTINE GETSID(NELEM,NPOIN,ILOCA,INTMA,ISIDE,LWHER,& 
     &                 LHOWM,ICONE,mxsid) 
! 
      IMPLICIT NONE 
! 
      INTEGER mxsid 
      INTEGER NELEM,ISIDE(8,mxsid),INTMA(3,NELEM) 
      INTEGER IS,IN,INODE,ILOC,IP1L,IP2L,IP1R,IP2R,IEL,IER 
      INTEGER IP1,IP2,IP,IE,JLOCA,ILOC1,IELE,IWHER,IN1 
      INTEGER IPT, J, IN2, NPOIN, ILOCA,LWHER(NPOIN) 
      INTEGER LHOWM(NPOIN), ICONE(3*NELEM) 
! ** FILL IN LHOWM : NR. OF ELEMENTS PER NODE ! How many lines 
! 
      DO 1490 IP=1,NPOIN 
      LHOWM(IP)=0 ! Zero ! thank you for this very informative comment
 1490 CONTINUE 
!                       
      DO 1501 IE=1,NELEM        
      DO 1500 IN=1,3  
 
      IP=INTMA(IN,IE) 
      LHOWM(IP)=LHOWM(IP)+1  
 1500 CONTINUE 
 1501 CONTINUE   
! 
! ** FILL IN LWHER : LOCATION OF EACH NODE INSIDE ICONE ! Where are the lines 
! 
      LWHER(1)=0 
      DO 1600 IP=2,NPOIN 
      LWHER(IP)=LWHER(IP-1)+LHOWM(IP-1) 
 1600 CONTINUE 
! 
! ** FILL IN ICONE : ELEMENTS IN EACH NODE 
! 
      DO 1690 IP=1,NPOIN 
      LHOWM(IP)=0 
 1690 CONTINUE 
!  
      DO 1700 IE=1,NELEM 
      DO 1701 IN=1,3 
      IP=INTMA(IN,IE) 
      LHOWM(IP)=LHOWM(IP)+1 
      JLOCA=LWHER(IP)+LHOWM(IP) 
      ICONE(JLOCA)=IE 
 1701 CONTINUE 
 1700 CONTINUE 
! 
! ** LOOP OVER THE NODES 
! 
      ILOCA=0 
! 
      DO 3000 IP=1,NPOIN 
      ILOC1=ILOCA 
      IELE=LHOWM(IP) 
      IF(IELE.EQ.0) GOTO 3000 
! 
! ** INITIALIZE ISIDE ----> IMPORTANT FOR BOUNDARY SIDES 
! 
      DO 3001 IS=1,IELE+2 
      ISIDE(3,IS+ILOC1)=0 
      ISIDE(4,IS+ILOC1)=0 
 3001 CONTINUE 
! 
      IWHER=LWHER(IP) 
! 
! ** LOOP OVER ELEMENTS SURROUNDING THE POINT IP 
! 
      IP1=IP 
      DO 3090 IEL=1,IELE 
      IE=ICONE(IWHER+IEL) 
! 
! ** FIND OUT POSITION OF IP IN THE CONNECTIVITY MATRIX 
! 
      DO 3091 IN=1,3 
      IN1=IN 
      IPT=INTMA(IN,IE) 
      IF(IPT.EQ.IP) GOTO 3092 ! break, we found the point IP
                              ! in the element, 1 <= IN <= 3
                              ! why not use the "EXIT" statement?
 3091 CONTINUE 
 3092 CONTINUE 
! 
      ! CICLE ON THE OTHER TWO POINTS IN THE ELEMENT
      DO 3100 J=1,2           ! we now find the other two points
      IN2=IN1+J               ! the other two are MOD(IN+1,3)+1
      IF(IN2.GT.3) IN2=IN2-3  ! and  MOD(IN+2,3)+1
      IP2=INTMA(IN2,IE)       !  
      IF(IP2.LT.IP1) GOTO 3100! IF IP2 < IP1, we already checked this
                              ! side, go to the next iteration
                              ! why not use the "CYCLE" statement?   

! ** CHECK THE SIDE ----->  NEW OR OLD 
! 
      IF(ILOCA.EQ.ILOC1) GOTO 7304 ! useless
      DO 5600 IS=ILOC1+1,ILOCA 
      JLOCA=IS 
      IF(ISIDE(2,IS).EQ.IP2) GOTO 7303 
 5600 CONTINUE 
 7304 CONTINUE 
! 
! ** NEW SIDE 
! 
      ILOCA=ILOCA+1 
      ISIDE(1,ILOCA)=IP1 
      ISIDE(2,ILOCA)=IP2 
      ISIDE(2+J,ILOCA)=IE 
      GOTO 3012 
! 
! ** OLD SIDE 
! 
 7303 CONTINUE 
      ISIDE(2+J,JLOCA)=IE 
 3012 CONTINUE 
! 
 3100 CONTINUE ! END OF CICLE ON THE OTHER TWO POINTS IN THE ELEMENT
! 
! ** END LOOP OVER ELEMENTS SURROUNDING POINT IP 
! 
 3090 CONTINUE 
! 
      DO 8000 IS=ILOC1+1,ILOCA 
      IF(ISIDE(3,IS).NE.0) GOTO 8000 ! "CYCLE" ( or "CONTINUE" in C)
                                     ! ""
      ISIDE(3,IS)=ISIDE(4,IS) 
      ISIDE(4,IS)=0 
      ISIDE(1,IS)=ISIDE(2,IS) 
      ISIDE(2,IS)=IP1     
 8000 CONTINUE 
! 
! ** END LOOP OVER POINTS 
!       
 3000 CONTINUE 
! 
! *** CHECK THAT THE NUMBER OF SIDES IS NOT TOO LARGE!
!
        IF(ILOCA.GT.MXSID)THEN
        PRINT*,'MXSID IS NOT BIG ENOUGH YOU MUPPET!'
        STOP
        ENDIF
! *** BEGIN LOOP OVER EACH SIDE 
! 
         DO 1000 IS=1,ILOCA 
! 
                IP1=ISIDE(1,IS) 
                IP2=ISIDE(2,IS) 
                IEL=ISIDE(3,IS) 
                IER=ISIDE(4,IS) 
! 
! *** CHECK WHETHER SIDE IS A BOUNDARY SIDE 
! 
                IF(IER.EQ.0)THEN 
                DO IN=1,3 
                INODE=INTMA(IN,IEL) 
                IF(INODE.EQ.IP1)THEN 
                ILOC=IN   ! FOUND THE POINT  INDEX
                GOTO 500  ! "CYCLE"
                ENDIF 
                ENDDO 
 500            IF(ILOC.EQ.3)THEN    ! that is,IP1L=ILOC and
                IP1L=3               ! IP2L=MOD(ILOC+1,3)+1
                IP2L=1 
                ELSE 
                IP1L=ILOC 
                IP2L=ILOC+1 
                ENDIF 
                ISIDE(5,IS)=IP1L 
                ISIDE(6,IS)=IP2L 
                ISIDE(7,IS)=0 
                ISIDE(8,IS)=0 
                ELSE 
! 
! *** ESTABLISH THE TWO LOCAL LHS NODES ASSOCIATED WITH THE EDGE  
!  
        DO IN=1,3 
      INODE=INTMA(IN,IEL) 
        IF(INODE.EQ.IP1)THEN 
        ILOC=IN 
        GOTO 100 
        ENDIF 
        ENDDO 
 100    IF(ILOC.EQ.3)THEN 
        IP1L=3 
        IP2L=1 
        ELSE 
        IP1L=ILOC 
        IP2L=ILOC+1 
        ENDIF 
         ISIDE(5,IS)=IP1L 
          ISIDE(6,IS)=IP2L 
! 
! *** ESTABLISH THE TWO LOCAL RHS NODES ASSOCIATED WITH THE EDGE  
! 
        DO IN=1,3 
        INODE=INTMA(IN,IER) 
        IF(INODE.EQ.IP2)THEN 
        ILOC=IN 
          GOTO 200 
          ENDIF 
         ENDDO 
 200      IF(ILOC.EQ.3)THEN 
        IP1R=1 
          IP2R=3 
          ELSE 
          IP2R=ILOC 
          IP1R=ILOC+1 
          ENDIF 
        ISIDE(7,IS)=IP1R 
        ISIDE(8,IS)=IP2R 
        ENDIF 
 1000 CONTINUE 
! 

        RETURN 
        END 
