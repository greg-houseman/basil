C*--------------------------------------------------------------------
C*    Basil / Sybil:   c3code.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE C3CODE(AMAT,IWORK,N2,M2,NV,NH,IC0,IC2,
     :                  DELX,DELY,TMIN,TMAX,
     :                  XCMIN,XCMAX,YCMIN,YCMAX)
C
C     Draws colour coded 2-dimensional plot of the values in AMAT
C     This routine supersedes C3CODEX which follows later in this file
C     rewritten while tracking down a bug related to misplacement of plot
C     on page. GAH 12/7/21
C
      DIMENSION AMAT(N2,M2),IWORK(N2,M2)
      DIMENSION POLYX(4),POLYY(4)
      LOGICAL FILL
C
      ICD =IC2-IC0
      TSCAL=TMAX-TMIN
      DCDT=FLOAT(ICD)/TSCAL
      IOUTLINE = 0
C
C   NV amd NH correspond to NY2 and NX2 respectively
C   X(Y)CMIN is in position 2, X(Y)CMAX in position NH (NV)
C   For each X-coordinate
C
      DO 18 I=2,NH
        XV=XCMIN+DELX*FLOAT(I-2)
        XVM=XV-DELX*0.5
        XVP=XV+DELX*0.5
        IF(I.EQ.2) XVM=XCMIN
        IF(I.EQ.NH)XVP=XCMAX
        POLYX(1)=XVM
        POLYX(2)=XVP
        POLYX(3)=POLYX(2)
        POLYX(4)=POLYX(1)
C
C   check if first cell occupied
C
        FILL=.FALSE.
        IF(IWORK(2,I).NE.0)THEN         ! start the column fill
          TVAL1=AMAT(2,I)
          ICJ1=IC0+INT(0.5+(TVAL1-TMIN)*DCDT)
          IF(ICJ1.LT.IC0)ICJ1=IC0
          IF(ICJ1.GT.IC2)ICJ1=IC2
          POLYY(1)=YCMIN
          POLYY(2)=POLYY(1)
          FILL=.TRUE.
        ENDIF
C
C     work up each column
C
        DO 17 J=2,NV
          YV=YCMIN+DELY*FLOAT(J-2)
          YVM=YV-DELY*0.5
          YVP=YV+DELY*0.5
          IF(J.EQ.2)YVM=YCMIN
          IF(J.EQ.NV)YVP=YCMAX
          POLYY(3)=YVP
          POLYY(4)=YVP
C
C     close polygon if end of column and still open
C
          IF(FILL.EQV..TRUE.)THEN
            IF((J.EQ.NV))THEN
              CALL FILLPOLY(POLYX,POLYY,4,ICJ1,IOUTLINE)
              FILL=.FALSE.
C
C     close polygon if not end of column but colour is absent
C
            ELSE
              IF(IWORK(J+1,I).EQ.0)THEN
                CALL FILLPOLY(POLYX,POLYY,4,ICJ1,IOUTLINE)
                FILL=.FALSE.
C
C     or if colour changes
C
              ELSE
                TVAL2=AMAT(J+1,I)
                ICJ2=IC0+INT(0.5+(TVAL2-TMIN)*DCDT)
                IF(ICJ2.LT.IC0)ICJ2=IC0
                IF(ICJ2.GT.IC2)ICJ2=IC2
                IF(ICJ2.NE.ICJ1)THEN
                  CALL FILLPOLY(POLYX,POLYY,4,ICJ1,IOUTLINE)
                  POLYY(1)=POLYY(3)       ! set values for next segment
                  POLYY(2)=POLYY(1)
                  ICJ1=ICJ2
                ENDIF
              ENDIF
            ENDIF
C
C    check for start of colour if FILL=.false.
C
          ELSE
            IF(IWORK(J,I).NE.0)THEN
              TVAL1=AMAT(J,I)
              ICJ1=IC0+INT(0.5+(TVAL1-TMIN)*DCDT)
              POLYY(1)=YVM
              POLYY(2)=POLYY(1)
              FILL=.TRUE.
            ENDIF
          ENDIF
   17   CONTINUE
   18 CONTINUE
      CALL SYBFLUSH
      RETURN
      END
C*--------------------------------------------------------------------

      SUBROUTINE C3CODEX(AMAT,IWORK,N2,M2,NV,NH,IC0,IC2,
     :                  DELX,DELY,TMIN,TMAX,
     :                  XCMIN,XCMAX,YCMIN,YCMAX)
C
C     Draws colour coded 2-dimensional plot of the values in AMAT
C
      DIMENSION AMAT(N2,M2),IWORK(N2,M2)
      DIMENSION POLYX(4),POLYY(4)
      LOGICAL FILL
C
      ICD =IC2-IC0
      TSCAL=TMAX-TMIN
C
C   The brush height of CHIT pixels is used on all rows 
C
      CHIT=(XCMAX-XCMIN)/FLOAT(NH-2)
      CH=CHIT
      CALL CONVERTUTOYPTS(CHIT,WDTH)
      IOUTLINE = 0
      X = XCMIN
      XMAX = XCMIN+DELX
C
C   For each X-coordinate
C
      DO 18 I=2,NH-1
      Y=YCMIN
      CALL PLOTU(X,Y,3)
      YMAX = Y+DELY
      BASE = 0.0
      TOP = 0.0
      NUMCNRS = 0
      NUMBASECNRS=0
      NUMTOPCNRS=0
      IF (IWORK(2,I).NE.0) THEN
        NUMBASECNRS = NUMBASECNRS+1
        BASE = AMAT(2,I)
      END IF
      IF (IWORK(2,I+1).NE.0) THEN
        NUMBASECNRS = NUMBASECNRS+1
        BASE = BASE + AMAT(2,I+1)
      END IF
      NUMCNRS = NUMBASECNRS
C     IF (NUMBASECNRS.GT.1) THEN
C       POLYX(1)=X
C       POLYY(1)=Y
C       POLYX(2)=XMAX
C       POLYY(2)=Y
C     END IF
      IF (IWORK(3,I).NE.0) THEN
        NUMTOPCNRS = NUMTOPCNRS+1
        TOP = TOP + AMAT(3,I)
      END IF
      IF (IWORK(3,I+1).NE.0) THEN
        NUMTOPCNRS = NUMTOPCNRS+1
        TOP = (AMAT(3,I+1) + AMAT(3,I+1))
      END IF
      NUMCNRS = NUMCNRS+NUMTOPCNRS
C     IF (NUMTOPCNRS.GT.1) THEN
C       POLYX(3)=XMAX
C       POLYY(3)=YMAX
C       POLYX(4)=X
C       POLYY(4)=YMAX
C       NUMCNRS = NUMCNRS+NUMTOPCNRS
C     END IF
C     IF(NUMCNRS.NE.0)THEN
      IF(NUMCNRS.GT.0)THEN
        T=(BASE + TOP)/FLOAT(NUMCNRS)
        IC=IC0 + INT(0.5+(T-TMIN)*FLOAT(ICD)/TSCAL)
        IF(IC.LT.IC0)IC=IC0
        IF(IC.GT.IC2)IC=IC2
      ELSE
        IC= 0
      END IF
      DO 17 J=3,NV-1
        FILL = .FALSE.
        BASE = TOP
        TOP = 0.0
        NUMBASECNRS = NUMTOPCNRS
        LASTCNRS = NUMCNRS
        NUMCNRS = NUMBASECNRS
        NUMTOPCNRS=0
        IF (IWORK(J+1,I).NE.0) THEN
          NUMTOPCNRS=NUMTOPCNRS+1
          TOP = TOP + AMAT(J+1,I)
        END IF
        IF (IWORK(J+1,I+1).NE.0) THEN
          NUMTOPCNRS=NUMTOPCNRS+1
          TOP = TOP + AMAT(J+1,I+1)
        END IF
        NUMCNRS = NUMCNRS + NUMTOPCNRS
        IF(NUMCNRS.GT.2)THEN
          T=(BASE + TOP)/FLOAT(NUMCNRS)
          NXTC=IC0 + INT(0.5+(T-TMIN)*FLOAT(ICD)/TSCAL)
          IF(NXTC.LT.IC0)NXTC=IC0
          IF(NXTC.GT.IC2)NXTC=IC2
          IF (IC.NE.NXTC.AND.IC.NE.0) THEN
            FILL = .TRUE.
          ELSE
            YMAX = YMAX+DELY
          END IF
        ELSE 
          IF (IC.NE.0) THEN
            FILL = .TRUE.
          ELSE
            Y=Y+DELY
            YMAX = Y+DELY
          END IF
          NXTC = 0
        END IF
        IF (FILL.EQV..TRUE.) THEN
          POLYX(1)=X
          POLYY(1)=Y
          POLYX(2)=XMAX
          POLYY(2)=Y
          POLYX(3)=XMAX
          POLYY(3)=YMAX
          POLYX(4)=X
          POLYY(4)=YMAX
          CALL FILLPOLY(POLYX,POLYY,4,IC,IOUTLINE)
          JNDX=1
C         DO 19 INDX=4,3,-1
C           POLYX(JNDX) = POLYX(INDX)
C           POLYY(JNDX) = POLYY(INDX)
C           JNDX=JNDX+1
C  19     CONTINUE
          Y=YMAX
          YMAX = YMAX+DELY
        END IF
        IC= NXTC
   17 CONTINUE
      IF (IC.NE.0) THEN
        POLYX(1)=X
        POLYY(1)=Y
        POLYX(2)=XMAX
        POLYY(2)=Y
        POLYX(3)=XMAX
        POLYY(3)=YMAX
        POLYX(4)=X
        POLYY(4)=YMAX
        CALL FILLPOLY(POLYX,POLYY,4,IC,IOUTLINE)
      END IF
      X=X+DELX
      XMAX = X+DELX
   18 CONTINUE
      CALL SYBFLUSH
      RETURN
      END
   
      SUBROUTINE DRAWCOLOURBAR(IC0,IC2,BARWDTH,BARLEN,XLEV0,YLEV0,
     1                         XMIN,XMAX,YMIN,YMAX,IORIENTVERT)
C
C      draw colour scale bar - aspect 10:1
C
      ICD =IC2-IC0
      FRAC=0.7
      IF (IORIENTVERT.EQ.1) THEN
        BARLEN=FRAC*(YMAX-YMIN)
      ELSE
        BARLEN=FRAC*(XMAX-XMIN)
      END IF
C     BARWDTH = BARLEN*0.08*(XMAX-XMIN)/(YMAX-YMIN)
      BARWDTH = BARLEN*0.1
      IF (IORIENTVERT.EQ.1) THEN
C          move 1.5xwidth as brush stroke centres on xlev0
        XLEV0=XMAX + BARWDTH*1.5
        YLEV0=YMIN+0.5*(1.0-FRAC)*(YMAX-YMIN)
        CALL CONVERTUTOXPTS(BARWDTH,WDTH)
        CALL CONVERTUTOYPTS(BARLEN,HGT)
      ELSE
C       YLEV0=YMIN - BARWDTH*1.5
        YLEV0=YMIN - BARWDTH*2.0
        XLEV0=XMIN+0.5*(1.0-FRAC)*(XMAX-XMIN)
        CALL CONVERTUTOXPTS(BARLEN,HGT)
        CALL CONVERTUTOYPTS(BARWDTH,WDTH)
      END IF
C      the next line is necessary
      CALL PLOTU(XLEV0,YLEV0,3)
      DINTRVL = FLOAT(ICD)/HGT
      IF (DINTRVL.LT.1.0) THEN
        INTRVL = 1
      ELSE
        INTRVL = DINTRVL
      END IF
      DINCR = BARLEN/FLOAT(ICD/INTRVL)
      YLEV = YLEV0
      XLEV = XLEV0
      CALL SETLINEWIDTH(WDTH)
      DO 800 I=IC0,IC2,INTRVL
        CALL SETPENCOLOR(I)
        IF (IORIENTVERT.EQ.1) THEN
          CALL PLOTU(XLEV0,YLEV,3)
          YLEV=YLEV+DINCR
          CALL PLOTU(XLEV0,YLEV,2)
        ELSE
          CALL PLOTU(XLEV,YLEV0,3)
          XLEV=XLEV+DINCR
          CALL PLOTU(XLEV,YLEV0,2)
        END IF
  800 CONTINUE
      IF (IORIENTVERT.EQ.1) THEN
        BARLEN = YLEV-YLEV0
      ELSE
        BARLEN = XLEV-XLEV0
      END IF
      CALL PLOTU(XLEV0,YLEV0,3)
      CALL SYBFLUSH
      RETURN
      END
C
      SUBROUTINE MARKCOLOURBAR(ZCMIN,ZCMAX,NUMCNTRS,STEP,SLVL,
     1                         IORIENTVERT,XLEV0,YLEV0,BARWDTH,BARLEN)
C
C     add contour levels to colour scale bar.  Allow for ZCMAX < ZCMIN
C
      ZKMIN=ZCMIN
      ZKMAX=ZCMAX
      YLV=YLEV0
      XLV=XLEV0
      NUMCNT=NUMCNTRS
      IF(ZCMIN.GT.ZCMAX)THEN
        ZKMIN=ZCMAX
        ZKMAX=ZCMIN
        YLV=YLEV0+BARLEN
        XLV=XLEV0+BARLEN
      ENDIF
      IF(NUMCNTRS.LT.0)NUMCNT=-NUMCNTRS
      TICLEN = BARWDTH*0.5
      SCFAC=BARLEN/(ZCMAX-ZCMIN)
      IF(SLVL.GT.ZKMIN)THEN
        NSTD=INT((SLVL-ZKMIN)/STEP)
        ZMIN=SLVL-FLOAT(NSTD)*STEP
      ELSE
        NSTD=1+INT((ZKMIN-SLVL)/STEP)
        ZMIN=SLVL+FLOAT(NSTD)*STEP
      ENDIF
      DO 200 J=1,NUMCNT
        CLEV=ZMIN+FLOAT(J-1)*STEP
        IF((CLEV.GT.ZKMIN).AND.(CLEV.LT.ZKMAX))THEN
          IF (IORIENTVERT.EQ.1) THEN
            YLEV=YLV+(CLEV-ZKMIN)*SCFAC
            CALL PLOTU(XLEV0,YLEV,3)
            CALL PLOTU(XLEV0+TICLEN,YLEV,2)
          ELSE
            XLEV=XLV+(CLEV-ZKMIN)*SCFAC
            CALL PLOTU(XLEV,YLEV0,3)
            CALL PLOTU(XLEV,YLEV0+TICLEN,2)
          END IF
        END IF
  200 CONTINUE
      RETURN
      END
