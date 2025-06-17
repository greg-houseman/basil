C*--------------------------------------------------------------------
C*    Basil / Sybil:   stmesh.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE STMESH(IQU,IPLOT,SNTRP,BNTRP,CNTRP,DNTRP,
C    :                  AMESH,BMESH,CMESH,DMESH,EMESH,
     :                  PROFILTMP,PROFILTMP2,EDGEPTS,PTS,
     :                  EX,EY,VHB,UVP,LEM,NOR,IHELP,IBNGH,IBC,
     :                  IHELP2,SSQ,IELFIX,INTV,RLV,
     :                  NE,NUP,NUVP,NFP,NN,NBP,NX3,NY3,NP3,NP,
     :                  IYDIR,MESHNUM,XP1,YP1,XP2,YP2,
     :                  XCMIN,YCMIN,XCMAX,YCMAX,
     :                  IVERBOSE,LABEL,LEN)
C 
C     Subroutine to arrange the strain rate or stress components
C     onto a mesh for contouring (IPLOT=5) or for arrow plots (IPLOT=4)
C     or for profiles (IPLOT=31 for 1D, IPLOT=32 for 2D)
C       The values for IPLOT must match those in plotdefs.h
C     This interpolation to a profile line uses NX3 for both X
C     and Y directions. If this is varied, check the size allocated
C     for the working arrays IHELP2, PROFILTMP, PROFILTMP2
C
      INCLUDE "../basilsrc/indices.parameters"
      INTEGER ARROW
      CHARACTER FNAME*10
      CHARACTER LABEL*30
C     DIMENSION CMESH(NUP),DMESH(NUP),AMESH(NUP),BMESH(NUP)
C     DIMENSION EMESH(NUP)
      DIMENSION SNTRP(NP3),BNTRP(NP3),CNTRP(NP3),DNTRP(NP3)
      DIMENSION EDGEPTS(4,NX3),PROFILTMP(NX3),PROFILTMP2(NX3)
      DIMENSION PTS(NBP)
      DIMENSION INTV(64),RLV(64)
      DIMENSION EX(NUP),EY(NUP),VHB(8,NE),UVP(NUVP)
      DIMENSION SSQ(NUP),IELFIX(NUP)
      DIMENSION LEM(6,NE),NOR(NUP),IBC(NBP),IBNGH(NBP*2)
      DIMENSION IHELP(NP3),IHELP2(NX3)
      ARROW=0
      IF (IPLOT.EQ.4) ARROW=1
C
C    for contours 
C
      IF((IPLOT.EQ.4).OR.(IPLOT.EQ.5))THEN
        CALL STRAIX(IQU,IPLOT,NX3,NY3,NP3,SNTRP,IHELP,BNTRP,
     :              CNTRP,DNTRP,XCMIN,YCMIN,XCMAX,YCMAX,
     :              RLV(IVC),RLV(ISE),
     :              RLV(IARGANP),RLV(IBRGANP),RLV(IHLENSC),
     :              RLV(IBDEPSC),RLV(IBIG),EX,EY,VHB,UVP,SSQ,IELFIX,
     :              LEM,NOR,NE,NUP,NUVP,NFP,INTV(IIVV),INTV(IIVIS),
     :              INTV(INCOMP),INTV(IICR),RLV(ITBXOFF),RLV(ITBYOFF),
     :              0,MESHNUM)
C
C       PSUM = 0
C        CALL PAVG(PSUM,EX,EY,UVP,NN,NUP,NFP,NUVP)
C
C   Output the interpolated values for scalar plot; for arrow plot
C   output is in ARWTWO.  The functionality of the following block
C   of code is moved to routine contour in pref77.c in order that
C   one routine handles different output arrays in a consistent
C   manner.   This block should be deleted if no problems emerge
C
C       IF((IVERBOSE.NE.0).AND.(IPLOT.NE.4))THEN
C         FNAME='sybil.'//LABEL(1:LEN)
C         WRITE(6,10002)LABEL(1:LEN),FNAME,LABEL(1:LEN)
C         OPEN(10,FILE=FNAME)
C         HSX=(XCMAX-XCMIN)/FLOAT(NX3-3)
C         HSY=(YCMAX-YCMIN)/FLOAT(NY3-3)
C         DO 111 J=1,NY3
C           YBP=YCMIN+HSY*FLOAT(J-2)
C           Y=YBP
C           DO 112 I=1,NX3
C             IJ=(J-1)*NX3+I
C             X=XCMIN+HSX*FLOAT(I-2)
C             IF(IHELP(IJ).NE.0)THEN
C               IF(INTV(INCOMP).EQ.-1)THEN
C                 Y=YBP
C                 CALL PROJECTDEG(X,Y,RLV(IXREFM),RLV(IYREFM),1,
C    :                         INTV(INCOMP),IERR)
C               ENDIF
C               WRITE(10,10003)X,Y,SNTRP(IJ)
C             ENDIF
C 112       CONTINUE
C 111     CONTINUE
C         CLOSE(10)
C       ENDIF
C
C        Profile 1_D; only find values on a specific line
C
      ELSE IF(IPLOT.EQ.31)THEN
        CALL STRAIL(IQU,NP,PROFILTMP,IHELP,XP1,YP1,XP2,YP2,
     :              RLV(IVC),RLV(ISE),
     :              RLV(IARGANP),RLV(IBRGANP),
     :              RLV(IHLENSC),RLV(IBDEPSC),RLV(IBIG),
     :              EX,EY,VHB,UVP,SSQ,IELFIX,LEM,NOR,NE,NUP,NUVP,
     :              NFP,INTV(IIVV),INTV(IIVIS),INTV(INCOMP),
     :              INTV(IICR))
C
C   Output the interpolated values if required
C
        IF(IVERBOSE.NE.0)THEN
          FNAME='sybil.'//LABEL(1:LEN)
          WRITE(6,10002)LABEL(1:LEN),FNAME,LABEL(1:LEN)
          OPEN(10,FILE=FNAME)
          WRITE(10,*)FNAME
          WRITE(10,*)'     X       Y       value'
          DXP=(XP2-XP1)/FLOAT(NX3-1)
          DYP=(YP2-YP1)/FLOAT(NX3-1)
          DO 212 I=1,NX3
              X=XP1+DXP*FLOAT(I-1)
              Y=YP1+DYP*FLOAT(I-1)
              IF(INTV(INCOMP).EQ.-1)THEN
                CALL PROJECTDEG(X,Y,RLV(IXREFM),RLV(IYREFM),1,
     :                         INTV(INCOMP),IERR)
              ENDIF
              WRITE(10,10003)X,Y,SNTRP(I)
  212     CONTINUE
          CLOSE(10)
        ENDIF
C
C     removing average pressure may be needed in some cases ?
C
C       PSUM = 0
C        CALL PAVG(PSUM,EX,EY,UVP,NN,NUP,NFP,NUVP)
C
C        Profile 2_D
C
      ELSE IF(IPLOT.EQ.32)THEN
C
C     X direction averages plotted against Y
C
        IF (IYDIR.EQ.1) THEN
          CALL EDGES(EDGEPTS,YP1,YP2,EX,EY,PTS,IBNGH,IBC,
     1                                  NOR,NUP,NE,NBP,NP)
          IX1=1
          IX2=3
          IY1=2
          IY2=4
          NPROF=NP
C         WRITE(6,10001)"X","Y","Y",YP1,YP2
C
C    or Y direction averages plotted against X
C
        ELSE
          CALL EDGES(EDGEPTS,XP1,XP2,EY,EX,PTS,IBNGH,IBC,
     1                                  NOR,NUP,NE,NBP,NP)
          IX1=2
          IX2=4
          IY1=1
          IY2=3
          NPROF=NP
C         WRITE(6,10001)"Y","X","X",XP1,XP2
        END IF
10001   FORMAT('2D profile: the average in the ',A1,'-direct',
     :         'ion is plotted vs ',A1,' for ',A1,' = ',F10.4,
     :         ' to ',F10.4)
C
C     Do the set of integrations required for each level
C
        DO 140 JJ = 1,NPROF
          X = EDGEPTS(IX2,JJ)-EDGEPTS(IX1,JJ)
          Y = EDGEPTS(IY2,JJ)-EDGEPTS(IY1,JJ)
          HLEN=SQRT(X*X+Y*Y)
C
C     the required quantity is interpolated into the regular 
C
          CALL STRAIL(IQU,NP,PROFILTMP2,IHELP2,
     :                EDGEPTS(IX1,JJ),EDGEPTS(IY1,JJ),
     :                EDGEPTS(IX2,JJ),EDGEPTS(IY2,JJ),
     :                RLV(IVC),RLV(ISE),
     :                RLV(IARGANP),RLV(IBRGANP),
     :                RLV(IHLENSC),RLV(IBDEPSC),RLV(IBIG),
     :                EX,EY,VHB,UVP,SSQ,IELFIX,LEM,NOR,NE,NUP,NUVP,
     :                NFP,INTV(IIVV),INTV(IIVIS),INTV(INCOMP),
     :                INTV(IICR))
          CALL INTEGRATTRPZ(PROFILTMP2,IHELP2,NP,HLEN,AREA,LABEL,0)
C
C       following lines were used to confirm integration done in PROFIL
C         IF(JJ.EQ.1)THEN
C           ASUM=AREA*0.5
C         ELSE IF(JJ.EQ.NPROF)THEN
C           ASUM=ASUM+AREA*0.5
C         ELSE
C           ASUM=ASUM+AREA
C         ENDIF
          PROFILTMP(JJ)=AREA
          IHELP(JJ)=1
  140   CONTINUE
C       ASUM=ASUM/FLOAT(NX3-1)
C       WRITE(*,*)'Average value of integral is ',ASUM
C
C    IHELP has the role of identifying points that are outside the solution
C    domain, but is here reset to be used to identify lines that
C    somewhere intersect the domain.
C
      END IF         ! if 2D plot
10002 FORMAT('interpolated values of ',A4,' written to ',A10,
     1 '; 3 columns: X, Y, ',A4)
10003 FORMAT(2F12.5,G15.6)
      RETURN
      END
