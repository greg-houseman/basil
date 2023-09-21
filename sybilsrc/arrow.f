C*--------------------------------------------------------------------
C*    Basil / Sybil:   arrow.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE ARROW(X,Y,RMAX,NXY,M1,M2,MP,N1,N2,NP,NX3,NY2,
     :                    XCMIN,XCMAX,YCMIN,YCMAX,XMID,YMID,
     :                    IVERBOSE,NCOMP)
      DIMENSION X(NXY),Y(NXY)
      PI = 3.14159
      RNX=(XCMAX-XCMIN)/FLOAT(M2 - M1)
      RNY=(YCMAX-YCMIN)/FLOAT(N2 - N1)
C
C   Verbose output is now provide by DataPrintMeshVelo
C     IF(IVERBOSE.NE.0)THEN
C       OPEN(14,FILE='velg.out',STATUS='UNKNOWN',IOSTAT=IERR)
C       WRITE(14,10002)
C       DO I = M1,M2,MP
C         DO J = N1,N2,NP
C           K = (J-1)*NX3 + I
C           XPT=XCMIN+FLOAT(I-M1)*RNX
C           YPT=YCMIN+FLOAT(J-N1)*RNY
C           IF (NCOMP.EQ.-1) 
C    :        CALL PROJECTDEG(XPT,YPT,XMID,YMID,1,NCOMP,IERR)
C           WRITE(14,10003)I,J,XPT,YPT,X(K),Y(K)
C         ENDDO
C       ENDDO
C       CLOSE(14)
C     END IF
      IF (NCOMP.EQ.-1)
     :  CALL VELOVECTORS(X,Y,XCMIN,YCMIN,RNX,RNY,M1,M2,N1,N2,NXY,NX3)
      RMIN = 0.0
      IF(RMAX.EQ.0.0)CALL RSCALE(X,Y,NXY,M1,M2,N1,N2,RMAX,RMIN,NX3)
      WRITE(6,10001)RMAX,RMIN
10001 FORMAT('Maximum velocity magnitude =',G12.5,8X,'Minimum =',G12.5)
      ARLEN=0.9*FLOAT(MP)*MIN(RNX,RNY)
      ARSC=ARLEN/RMAX
      APEX=10.0
      TAPEX=ARSC*TAN(APEX*0.5*3.141592654/180.0)
C
C     draw velocity vectors
C
      IF(MP.LE.0)MP=1
      IF(NP.LE.0)NP=1
      DO I = M1,M2,MP
        DO J = N1,N2,NP
          K = (J-1)*NX3 + I
          RX=X(K)*X(K) + Y(K)*Y(K)
          IF(RX.GE.0.0005*RMAX*RMAX)THEN
            XPT=XCMIN+FLOAT(I-M1)*RNX
            YPT=YCMIN+FLOAT(J-N1)*RNY
            DX=TAPEX*X(K)
            DY=TAPEX*Y(K)
            XT=XPT - DY
            YT=YPT + DX
            CALL PLOTU(XT,YT,3)
            XT=XPT + X(K)*ARSC
            YT=YPT + Y(K)*ARSC
            CALL PLOTU(XT,YT,2)
            XT=XPT + DY
            YT=YPT - DX
            CALL PLOTU(XT,YT,2)
          ENDIF        ! if arrows large enough
        ENDDO          ! on J
      ENDDO            ! in I
C
C    print out velocity data if required   (?where?)
C
10002 FORMAT('Gridded Velocity vectors follow',/,
     1 '  I    J       X          Y          U          V')
10003 FORMAT(2I5,4F11.5)
      RETURN
      END

      SUBROUTINE ARWTWO(C2T,S2T,WX,WY,RMAX,RMIN,NXY,M1,M2,MP,N1,N2,NP,
     :                    NP3,NY2,NX3,INCI,XCMIN,YCMIN,XCMAX,YCMAX,
     :                    XREFM,YREFM,NCOMP,LABEL,IVERBOSE)
C
      CHARACTER FNAME*10,LABEL*30
      DIMENSION C2T(NXY),S2T(NP3)
      DIMENSION WX(NXY),WY(NXY)
      DIMENSION INCI(2)
      DATA PI/3.141592654/
      SAVE PI
C
      RNX=(XCMAX-XCMIN)/FLOAT(M2 - M1)
      RNY=(YCMAX-YCMIN)/FLOAT(N2 - N1)
      ARLEN=0.8*FLOAT(MP)*MIN(RNX,RNY)
      ARSC=ARLEN/RMAX
      DBL=0.002
      WRITE(6,10000)RMAX,RMIN
10000 FORMAT('Maximum stress/strain-rate magnitude =',G12.5,6X,
     :       'Minimum =',G12.5)
C   verbose output is produced in DataPrintArrow in pref77.c
C     IF(IVERBOSE.NE.0)THEN
C       FNAME='sybil.'//LABEL(1:4)
C       WRITE(6,10002)LABEL(1:4),FNAME
C10002   FORMAT('interpolated values of ',A4,' written to ',A10,
C    1 '; 6 columns: X, Y, s1, s2, s3 theta')
C       OPEN(10,FILE=FNAME)
C       WRITE(10,10001)
C10001   FORMAT('   Centre (xc,yc)      axis 1    axis 2    ',
C    :         'axis c  orientation')
C     ENDIF
C
C    first plot all the relatively positive axes with colour 2, then all the
C    relatively less positive with colour 1
C
      DO KC=1,2
        INC=INCI(3-KC)
        CALL SETPENCOLOR(INC)
C
        DO J = N1,N2,NP
          YC=YCMIN+FLOAT(J-N1)*RNY
          DO I = M1,M2,MP
            XC=XCMIN+FLOAT(I-M1)*RNX
C
C    the two stress components are stored in WX and WY, the
C     orientation of the max principal stress are delivered as
C     sine and cosine by inverting the function set up in STCOMP
C
            K=(J-1)*NX3 + I
            S1=WX(K)
            S2=WY(K)
            IF((S1.NE.0.0).OR.(S2.NE.0.0))THEN
              IF (S2T(K).EQ.0.0.AND.C2T(K).EQ.0.0) THEN
                TH = 0.0
              ELSE
                TH=ATAN2(S2T(K),C2T(K))
              END IF
C
C   verbose output is produced in DataPrintArrow in pref77.c
C             IF((IVERBOSE.NE.0).AND.(KC.EQ.1))THEN
C               CC=-S1-S2
C               THDEG=TH*180.0/PI
C               IF (NCOMP.EQ.-1) THEN
C                 XDEG = XC
C                 YDEG = YC
C                 CALL PROJECTDEG(XDEG,YDEG,XREFM,YREFM,1,
C    :                            NCOMP,IERR)
C                 WRITE(10,10003)XDEG,YDEG,S1,S2,CC,THDEG
C               ELSE
C                 WRITE(10,10003)XC,YC,S1,S2,CC,THDEG
C               ENDIF
C10003          FORMAT(5F10.5,F11.5)
C             END IF
C
C      Draw max principal axis
C
              IF((NCOMP.EQ.1).AND.(S1.GE.0.0))GO TO 3
              IF((KC.EQ.1).AND.(S1.LT.0.0))GO TO 3
              IF((KC.EQ.2).AND.(S1.GE.0.0))GO TO 3
              XM=S1*ARSC*C2T(K)
              YM=S1*ARSC*S2T(K)
              XO=XCMIN + FLOAT(I-M1)*RNX - XM*0.5
              YO=YCMIN + FLOAT(J-N1)*RNY - YM*0.5
              CALL PLOTU(XO,YO,3)
              CALL PLOTU(XO+XM,YO+YM,2)
C
C     draw 2nd stress axis
C
    3         IF((NCOMP.EQ.1).AND.(S2.GE.0.0))GO TO 4
              IF((KC.EQ.1).AND.(S2.LT.0.0))GO TO 4
              IF((KC.EQ.2).AND.(S2.GE.0.0))GO TO 4
              XM=S2*ARSC*S2T(K)
              YM=-S2*ARSC*C2T(K)
              XO=XCMIN + FLOAT(I-M1)*RNX - XM*0.5
              YO=YCMIN + FLOAT(J-N1)*RNY - YM*0.5
              CALL PLOTU(XO,YO,3)
              CALL PLOTU(XO+XM,YO+YM,2)
 4            CONTINUE
            ENDIF           ! if arrow non-zero
          ENDDO             ! on I
        ENDDO               ! on J
      ENDDO                 ! on KC
      CALL SETPENCOLOR(1)
      RETURN
      END

      SUBROUTINE CNTOUR(X,WORK,IWORK,IHELP,NXY,M1,M2,N1,N2,NX3,
     1 NNX,NNY,HMESH,VMESH,ZMIN,ZMAX,XCMIN,XCMAX,YCMIN,YCMAX,SUMINT)
      DIMENSION X(NXY),WORK(N2,M2),IHELP(NXY)
      DIMENSION IWORK(N2,M2)
C
C     find the maximum and minimum values in the region to be contoured
C
      CALL RANGXY(X,IHELP,NXY,M1,M2,1,N1,N2,1,ZMAX,IMX,JMX,
     1                             ZMIN,IMN,JMN,NX3,SUMINT)
C
C    transpose the array into WORK to be passed to CNTR and C3CODE
C    The initial buffer row is retained the final buffer row is stripped
C
      NNX = M2
      NNY = N2
      DO IC = 2,M2
        DO JC = 2,N2
          K = (JC - 1)*NX3 + IC
          WORK(JC,IC) = X(K)
          IWORK(JC,IC)=IHELP(K)
        ENDDO
      ENDDO
      HMESH=(XCMAX-XCMIN)/FLOAT(M2-M1)
      VMESH=(YCMAX-YCMIN)/FLOAT(N2-N1)
      SUMINT=SUMINT*HMESH*VMESH
      RETURN
      END

      SUBROUTINE LINECNTOUR(X,WORK,IWORK,ICON,M2,N2,NNX,NNY,HMESH,VMESH,
     1                ZMIN,ZMAX,XCMIN,YCMIN,ISTIP,
     2                NCONTR,CON,JCON)

      DIMENSION WORK(N2,M2),IWORK(N2,M2),CON(3)
      DIMENSION X(M2*N2)
C
C     write X and Y coordinates into first row and column of WORK

      DO 16 I = 1,NNY
      WORK(I,1)=YCMIN+FLOAT(I-2)*VMESH
   16 CONTINUE
      DO 7 J = 1,NNX
      WORK(1,J)=XCMIN+FLOAT(J-2)*HMESH
    7 CONTINUE
C
C      To decide the contour levels: Check that CLEVEL is in the actual
C      range of values.  If not move it near the midpoint of the range
C      (in multiples of CSTEP).  Determine CMIN and CMAX
C      in multiples of CSTEP removed from the adjusted CLEVEL.
C      ZMAX may be < ZMIN for purpose of inverting colour scale
C      for contours we ensure CMAX > CMIN
C
        CLEVEL=CON(1)
        CSTEP=CON(2)
        EPS=CSTEP*0.01
        NABOVE=0
        NBELOW=0
        IF(((ZMAX-CLEVEL)*(ZMIN-CLEVEL)).GT.0)THEN
          CLEVEL=CSTEP*INT(0.5*(ZMAX+ZMIN)/CSTEP)
        END IF
        CMIN=ZMIN
        CMAX=ZMAX
        IF(ZMAX.LT.ZMIN)THEN
          CMIN=ZMAX
          CMAX=ZMIN
        END IF
        IF(CMAX.GT.CLEVEL)NABOVE=INT((CMAX-CLEVEL)/CSTEP+EPS)
        IF(CMIN.LT.CLEVEL)NBELOW=INT((CLEVEL-CMIN)/CSTEP+EPS)
        CMIN=CLEVEL-CSTEP*NBELOW
        CMAX=CLEVEL+CSTEP*NABOVE
C       WRITE(*,*)'NBELOW,NABOVE',NBELOW,NABOVE
        NCONT=NABOVE+NBELOW+1
        IF(NCONT.GT.50)NCONT=50
C     END IF
C
C    save contour parameters (CSTEP used in labelling), and report
C
      CON(1) = CLEVEL
      CON(2) = CSTEP
      CON(3) = CMAX
      WRITE(6,100) CMIN,CMAX,CSTEP,NCONT,CLEVEL
 100  FORMAT('Contours: ',G11.4,' to 1',G11.4,' step: ',
     :       G11.4,' num:',I3,' level: ',G11.4)
C
C    stipple regions below SLEV if required

   26 IF(ISTIP.NE.0)THEN
        SLEV=CLEVEL
        DO J=2,M2,2
          XST=WORK(1,J)
          DO I=2,N2,2
            YST=WORK(I,1)
            IF(ISTIP.EQ.1.AND.IWORK(I,J).GT.0
     :                                .AND.WORK(I,J).LE.SLEV) THEN
              CALL PLOTU(XST,YST,3)
              CALL PLOT(-0.015,0.0,-3)
              CALL PLOT( 0.030,0.0,-2)
C             CALL SQUARE(XST,YST,0.5,1)
            ELSE IF(ISTIP.EQ.2.AND.IWORK(I,J).GT.0
     :                                .AND.WORK(I,J).GE.SLEV) THEN
              CALL PLOTU(XST,YST,3)
              CALL PLOT(-0.015,0.0,-3)
              CALL PLOT( 0.030,0.0,-2)
C             CALL SQUARE(XST,YST,0.5,1)
            ELSE IF(ISTIP.EQ.3.AND.IWORK(I,J).GT.0
     :                           .AND.ABS(WORK(I,J)).LE.SLEV) THEN
              CALL PLOTU(XST,YST,3)
              CALL PLOT(-0.015,0.0,-3)
              CALL PLOT( 0.030,0.0,-2)
C             CALL SQUARE(XST,YST,0.5,1)
            ELSE IF(ISTIP.EQ.4.AND.IWORK(I,J).GT.0
     :                           .AND.ABS(WORK(I,J)).GE.SLEV) THEN
              CALL PLOTU(XST,YST,3)
              CALL PLOT(-0.015,0.0,-3)
              CALL PLOT( 0.030,0.0,-2)
C             CALL SQUARE(XST,YST,0.5,1)
            END IF
          ENDDO
        ENDDO
      END IF
C
C    CNTR wants CMIN < CMAX and NCONT > 0
C    Draw the contours between levels of CMIN and CMAX
C
C               array,rowdim,coldim,maxdim,locntr,hicntr,maxcntr
      CALL CNTR(WORK,NNY,NNX,N2,CMIN,CMAX,NCONT,JCON,0,0)
C
      RETURN
      END

      SUBROUTINE ROTOR(ANGL,FRMIN,BASE,RMAX,
     :                 NXY,M1,M2,MP,N1,N2,NP,NX3,
     :                 XCMIN,XCMAX,YCMIN,YCMAX)
      DIMENSION ANGL(NXY)
C
C    the parameter BASE allows for solution rotation
C
      PIB180=3.14159265/180.0
      ROMAX=0.0
      RNX=(XCMAX-XCMIN)/FLOAT(M2 - M1)
      RNY=(YCMAX-YCMIN)/FLOAT(N2 - N1)
C
      ARLEN=FLOAT(MP)*RNX
      ARSC = ARLEN/RMAX
      DO I = M1,M2,MP
        DO J = N1,N2,NP
          K = (J-1)*NX3 + I
          XO=XCMIN + FLOAT(I-M1)*RNX
          YO=YCMIN + FLOAT(J-N1)*RNY
          IF(ABS(ANGL(K)).GT.ABS(ROMAX))ROMAX=ANGL(K)
          IF(ABS(ANGL(K)).GE.FRMIN)THEN
            XM=ARSC*COS((ANGL(K)+BASE)*PIB180)
            YM=ARSC*SIN((ANGL(K)+BASE)*PIB180)
            CALL PLOTU(XO-XM*0.5,YO-YM*0.5,3)
            CALL PLOTU(XO+XM*0.5,YO+YM*0.5,2)
          ENDIF       ! if ANGL large enough
        ENDDO         ! on J
      ENDDO           ! on I
      RETURN
      END
C
      SUBROUTINE VPRINT(EX,EY,UVP,NOR,NCOMP,NUP,NN,XREFM)
C
C    routine to print out vertex node coordinates and velocities
C
      DIMENSION EX(NN),EY(NN),UVP(2*NUP),NOR(NUP)
C
      RTODEG=57.29577951
      PION2=1.570796327
      OPEN(14,FILE='veln.out',STATUS='UNKNOWN',IOSTAT=IERR)
      DO 10 J=1,NUP
        NODE=NOR(J)
        IF(NODE.LE.NN)THEN
          XX=EX(NODE)
          YY=EY(NODE)
          IF(NCOMP.EQ.-1)THEN
            XX=XREFM+RTODEG*XX/SIN(YY)
            YY=(YY-PION2)*RTODEG
          END IF
          WRITE(14,10001)NODE,XX,YY,UVP(J),UVP(J+NUP)
10001     FORMAT(I5,4F13.6)
        END IF
   10 CONTINUE
      CLOSE(14)
      RETURN
      END
C
      SUBROUTINE LPRINT(EX,EY,SS,NOR,NCOMP,NUP,NN,XREFM)
C
C    routine to print out vertex node coordinates and layer thickness (or topo)
C
      DIMENSION EX(NN),EY(NN),SS(NUP),NOR(NUP)
C
      RTODEG=57.29577951
      PION2=1.570796327
      OPEN(14,FILE='layer.out',STATUS='UNKNOWN',IOSTAT=IERR)
      DO 10 J=1,NUP
        NODE=NOR(J)
        IF(NODE.LE.NN)THEN
          XX=EX(NODE)
          YY=EY(NODE)
          IF(NCOMP.EQ.-1)THEN
            XX=XREFM+RTODEG*XX/SIN(YY)
            YY=(YY-PION2)*RTODEG
          END IF
          WRITE(14,10001)NODE,XX,YY,EXP(SS(J))
10001     FORMAT(I5,3F13.6)
        END IF
   10 CONTINUE
      CLOSE(14)
      RETURN
      END
