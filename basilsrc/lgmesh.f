C*--------------------------------------------------------------------
C*    Basil / Sybil:   lgmesh.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE LGMESH(NX,NY,XMIN,YMIN,XMAX,YMAX,IFCASE,LGEM,
     :                 LGIBC,LGIBCF,EXLG,EYLG,EX,EY,
     :                 IBC,IBNGH,IBCTYP,NOR,
     :                 NEL,NUL,NNL,NML,NBL,NUP,NBP,IFLT,LUW,IERR)
C
C    A regular mesh of triangular elements is defined in a
C   box which is initially square (or rectangular).  EXLG and EYLG
C   are the X and Y coordinates of the nodes and LGEM gives the
C   nodes associated with each element (in anticlockwise order
C   from the bottom left : 1, 2 & 3).  The nodes are numbered
C   from 1 in the bottom left hand corner, in columns.
C
      DIMENSION NOR(NUP)
      DIMENSION LGEM(6,NEL),EXLG(NUL),EYLG(NUL)
      DIMENSION LGIBC(NBL),LGIBCF(NBL)
      DIMENSION IBC(NBP),IBNGH(NBP*2),IBCTYP(NBP*2)
      DIMENSION EX(NUP),EY(NUP)
C
C   LaGrange mesh is fractionally smaller than the bounding box
C   determined by parameter FRAC
      FRAC = 0.999

C   XMAX already doubled for FLT=2
C   XMAX is reset after first lot of elements for FLT=2
      IF (IFLT.EQ.2) XMAX=XMAX/2.0
      XLMIN=(XMIN+XMAX)/2 - (XMAX-XMIN)/2*FRAC
      XLMAX=(XMIN+XMAX)/2 + (XMAX-XMIN)/2*FRAC
      YLMIN=(YMIN+YMAX)/2 - (YMAX-YMIN)/2*FRAC
      YLMAX=(YMIN+YMAX)/2 + (YMAX-YMIN)/2*FRAC
C
      NXL=NX
      NYL=NY
      NXL1=NXL+1
      NYL1=NYL+1
      DX=(XLMAX-XLMIN)/FLOAT(NXL)
      DY=(YLMAX-YLMIN)/FLOAT(NYL)
C
C     JEL is an index for the elements
C
      JEL=0
C
C     First the nodes
C
      DO 10 I=1,NXL1
      DO 20 J=1,NYL1
         NL=(I-1)*NYL1+J
         EXLG(NL)=FLOAT(I-1)*DX+XLMIN
         EYLG(NL)=FLOAT(J-1)*DY+YLMIN
   20 CONTINUE
   10 CONTINUE
C
C    Now the elements
C
      DO 30 I=1,NXL
      DO 40 J=1,NYL
C
         JEL=JEL+1
C
C vertex nodes
C
         NL1=(I-1)*NYL1+J
         NL2=NL1+1
         NL3=NL1+NYL1
         NL4=NL3+1
C
C midpoint nodes
C
         NM1=NNL+(I-1)*(3*NYL+1)+J
         NM2=NM1+NYL+J-1
         NM3=NM2+1
         NM4=NM2+2
         NM5=NM1+3*NYL+1
C
C    the lower side
C
         LGEM(1,JEL)=NL1
         LGEM(2,JEL)=NL3
         LGEM(3,JEL)=NL2
         LGEM(4,JEL)=NM1
         LGEM(5,JEL)=NM2
         LGEM(6,JEL)=NM3
C
C    the upper side
C
         JELN=JEL+NYL
         LGEM(1,JELN)=NL3
         LGEM(2,JELN)=NL4
         LGEM(3,JELN)=NL2
         LGEM(4,JELN)=NM3
         LGEM(5,JELN)=NM5
         LGEM(6,JELN)=NM4
   40 CONTINUE
         JEL=JEL+NYL
   30 CONTINUE
C
C    If there is an internal fault, duplicate other side of fault
C        and double XMAX for use by LGBND
C
      IF(IFLT.EQ.2) THEN
        DO 50 II=1,NNL/2
           IJ=II+NNL/2
C         EXLG(IJ)=2.0*EXLG(NNL/2)-EXLG(II)
           EXLG(IJ)=EXLG(NNL/2)+EXLG(II)
           EYLG(IJ)=EYLG(II)
   50   CONTINUE
        DO 60 II=1,NEL/2
           IJ=II+NEL/2
           LGEM(1,IJ)=LGEM(1,II)+NNL/2
           LGEM(2,IJ)=LGEM(2,II)+NNL/2
           LGEM(3,IJ)=LGEM(3,II)+NNL/2
           LGEM(4,IJ)=LGEM(4,II)+NML/2
           LGEM(5,IJ)=LGEM(5,II)+NML/2
           LGEM(6,IJ)=LGEM(6,II)+NML/2
   60   CONTINUE
        XMAX=2.0*XMAX
        XLMAX=2.0*XLMAX
        JEL = JEL*2
      END IF
C
C  Find coordinates for the midpoint nodes
C
      DO 120 JJ=1,NEL
        NJ1=LGEM(1,JJ)
        NJ2=LGEM(2,JJ)
        NJ3=LGEM(3,JJ)
        NJ4=LGEM(4,JJ)
        NJ5=LGEM(5,JJ)
        NJ6=LGEM(6,JJ)  
        EXLG(NJ4)=(EXLG(NJ3)+EXLG(NJ1))/2
        EXLG(NJ5)=(EXLG(NJ1)+EXLG(NJ2))/2
        EXLG(NJ6)=(EXLG(NJ2)+EXLG(NJ3))/2
        EYLG(NJ4)=(EYLG(NJ3)+EYLG(NJ1))/2
        EYLG(NJ5)=(EYLG(NJ1)+EYLG(NJ2))/2
        EYLG(NJ6)=(EYLG(NJ2)+EYLG(NJ3))/2
  120 CONTINUE
C
C    Check the number of elements and nodes
C
C      NMT=NM5-NN
      IF(JEL.NE.NEL)WRITE(LUW,10050)JEL,NEL
C      IF(NLC.NE.NN)WRITE(LUW,10051)NLC,NN/2
C      IF(NMT.NE.NMPP)WRITE(LUW,10002)NMT,NMPP/2
C       WRITE(LSC,10003)NUP,NN,NMP
C      DO 130 JJ=1,NUPP
C         WRITE(LSC,10001)JJ,EXLG(JJ),EYLG(JJ)
C  130 CONTINUE
C      CALL MITPRT(LGEM,6,NE,6*NE,LSC)
C
C     WRITE(*,*) 'Finished LaGrangian mesh'
C
C  Set up boundary flags - to be used for tracking
C   deformation of the fault and lagrangian mesh
C
C     CALL LGBND(IFCASE,0.0,XMAX,0.0,YMAX,LGEM,EXLG,EYLG,
      CALL LGBND(IFCASE,XLMIN,XLMAX,YLMIN,YLMAX,LGEM,EXLG,EYLG,
     :                 LGIBC,LGIBCF,EX,EY,IBC,IBNGH,IBCTYP,NOR,
     :                 NEL,NUL,NBL,NUP,NBP,IERR)
C
10050 FORMAT(' NO. OF LG ELEMENTS INCORRECT :',2I6)
10051 FORMAT(' NO. OF LG NODES INCORRECT :',2I6)
10001 FORMAT('I=',I6,' EXLG=',F8.4,' EYLG=',F8.4)
10002 FORMAT(' NO. OF MIDPOINT NODES INCORRECT :',2I6)
10003 FORMAT('NUP=',I6,' NN=',I6,' NMP=',I6)
C
      RETURN
      END

      SUBROUTINE LGMARK(LBC,STELPX,STELPY,
     :                  STELPR,
     :                  NSM,NPM,NRM,LUW,LSC,IERR)
      INCLUDE 'input.parameters'
      INCLUDE 'limits.parameters'
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      DIMENSION STELPX(NPM,NSM),STELPY(NPM,NSM)
      DIMENSION CXX(NSM),CYY(NSM)
      INCLUDE 'input.data'

C    NST is the number of strain markers already in the array
      IEND=0
      NST = 0
 110  CONTINUE
      IVALID=0
      CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
      IF (IEND.EQ.1) GO TO 700
      ITERM=0
      INDX2=1
      CALL GETWORD(INSTR,LENGTH,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0.AND.INDX2.LT.LENGTH) THEN
        ITERM = MATCH(SETUPTERMS,NUMSETUP,WORD)
        IF (ITERM.EQ.INDXMARKERS) THEN
          IVALID=1
          CALL GETSTRMARKDATA(INSTR,LENGTH,INDX2,NROWS,NCOLS,
     :                          STELPR,XMIN,XMAX,YMIN,YMAX,IERR)
        END IF
      END IF
      IF (IVALID.NE.1) GO TO 110
C
      IF (IERR.NE.0) GO TO 700
C If Radius <= 0 , then this is a strain marker line rather than a strain elipse
      IF (STELPR.LT.0.0001) GO TO 400
      IF (NST+NROWS*NCOLS.GT.NSM) THEN
        WRITE(LSC,10010)NSM
        WRITE(LUW,10010)NSM
        IERR=1
        GO TO 700
      END IF
      KK=0
      DX=0
      DY=0
      IF (NCOLS.GT.1) DX = (XMAX-XMIN)/(NCOLS-1)
      IF (NROWS.GT.1) DY = (YMAX-YMIN)/(NROWS-1)
      DO 310 II=1,NCOLS
      DO 320 JJ=1,NROWS
         KK=KK+1
         CXX(KK)=XMIN + (II-1)*DX
         CYY(KK)=YMIN + (JJ-1)*DY
  320 CONTINUE
  310 CONTINUE
C
C  Set up strain elipse markers
C
      DO 200 I=1,NROWS*NCOLS
      DO 210 J=1,NPM
         R=STELPR
         THETA=FLOAT(J-1)/FLOAT(NPM-NRM-1)*2.0*PI
         IF(J.GE.NPM-NRM) THEN
            THETA=0.0
            R=FLOAT(NPM-J)/FLOAT(NRM)*STELPR
         END IF
         STELPX(J,I+NST)=R*COS(THETA)+CXX(I)
         STELPY(J,I+NST)=R*SIN(THETA)+CYY(I)
!The arrays which will indicate if a point is inside or outside the domain
C        ISTELP(J,I+NST)=0
  210 CONTINUE
  200 CONTINUE
      NST = NST + KK
      GO TO 110
C
C  Set up strain line markers
C
 400  CONTINUE
      IF (NST+1.GT.NSM) THEN
        WRITE(LSC,10010)NSM
        WRITE(LUW,10010)NSM
        IERR=1
        GO TO 700
      END IF
      DX = (XMAX-XMIN)/(NPM-1)
      DY = (YMAX-YMIN)/(NPM-1)
      DO 410 II=1,NPM
         STELPX(II,NST+1)=XMIN + (II-1)*DX
         STELPY(II,NST+1)=YMIN + (II-1)*DY
  410 CONTINUE
      NST = NST + 1
      GO TO 110
 700  RETURN
10010 FORMAT('Total number of strain markers exceeds ',I4,' (NSMP)')
      END

      SUBROUTINE COUNTSTRMARKERS(LBC,NSM,LSC,IERR)
      INCLUDE 'input.parameters'
      INCLUDE 'limits.parameters'
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      INCLUDE 'input.data'

      IEND=0
      NSM = 0
      IFOUND=0
 110  CONTINUE
      CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
      IF (IEND.EQ.1) GO TO 700
      ITERM=0
      INDX2=1
      CALL GETWORD(INSTR,LENGTH,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0.AND.INDX2.LT.LENGTH) THEN
        ITERM = MATCH(SETUPTERMS,NUMSETUP,WORD)
        IF (ITERM.EQ.INDXMARKERS) THEN
          IFOUND=1
          NROWS=0
          NCOLS=0
          CALL GETSTRMARKDATA(INSTR,LENGTH,INDX2,NROWS,NCOLS,
     :                          STELPR,DUM2,DUM3,DUM4,DUM5,IERR)
          IF(STELPR.LT.0.0001) THEN
            NROWS = 1
            NCOLS = 1
          END IF
          NSM = NSM + NROWS*NCOLS
        END IF
      END IF
      IF (IFOUND.NE.1) GO TO 110
C
      IF (IERR.NE.0) GO TO 700
      GO TO 110
 700  IF (IFOUND.NE.1) THEN
        WRITE(LSC,10020)
        IERR=1
      END IF
      RETURN
10020 FORMAT('No MARKERS statements found')
      END

      SUBROUTINE GETSTRMARKDATA(INSTR,LEN,INDX2,NROWS,NCOLS,
     2                       STELPR,XMIN,XMAX,YMIN,YMAX,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXROWS=1,
     :                     INDXCOLS=2,
     :                     INDXRADIUS=3,
     :                     INDXLGXMIN=4,
     :                     INDXLGXMAX=5,
     :                     INDXLGYMIN=6,
     :                     INDXLGYMAX=7,
     :                     NUMSTRMARKTERMS=7)
      CHARACTER STRMARKTERMS(1:NUMSTRMARKTERMS)*12
      DATA STRMARKTERMS  /'ROWS        ',
     :                    'COLS        ',
     :                    'R           ',
     :                    'XMIN        ',
     :                    'XMAX        ',
     :                    'YMIN        ',
     :                    'YMAX        '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0.AND.INDX2.LT.LEN) THEN
        ITERM = MATCH(STRMARKTERMS,NUMSTRMARKTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          IERR=10
        ELSE IF (ITERM.EQ.INDXROWS) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.NE.0)
     1      READ(INSTR(INDX1:INDX2-1),*)NROWS
        ELSE IF (ITERM.EQ.INDXCOLS) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.NE.0)
     1      READ(INSTR(INDX1:INDX2-1),*)NCOLS
        ELSE IF (ITERM.EQ.INDXRADIUS) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.NE.0)
     1      READ(INSTR(INDX1:INDX2-1),*)STELPR
        ELSE IF (ITERM.EQ.INDXLGXMIN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.NE.0)
     1      READ(INSTR(INDX1:INDX2-1),*)XMIN
        ELSE IF (ITERM.EQ.INDXLGXMAX) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.NE.0)
     1      READ(INSTR(INDX1:INDX2-1),*)XMAX
        ELSE IF (ITERM.EQ.INDXLGYMIN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.NE.0)
     1      READ(INSTR(INDX1:INDX2-1),*)YMIN
        ELSE IF (ITERM.EQ.INDXLGYMAX) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.NE.0)
     1      READ(INSTR(INDX1:INDX2-1),*)YMAX
        END IF
        IF (INDX2.EQ.0) IERR=10
        IF (IERR.EQ.0) GO TO 5
      END IF
      IF (IERR.NE.0) INDX2=INDX1
      RETURN
      END
C
      SUBROUTINE LGBND(IFCASE,X1,X2,Y1,Y2,LGEM,EXLG,EYLG,
     :                 LGIBC,LGIBCF,EX,EY,IBC,IBNGH,IBCTYP,NOR,
     :                 NEL,NUL,NBL,NUP,NBP,IERR)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION LGEM(6,NEL),EXLG(NUL),EYLG(NUL)
      DIMENSION LGIBC(NBL),LGIBCF(NBL)
      DIMENSION IBC(NBP),IBNGH(NBP*2),IBCTYP(NBP*2)
C
      XMID=(X2+X1)/2.0
      YMID=(Y2+Y1)/2.0
      Y14=(Y2-Y1)*0.25+Y1
      Y34=(Y2-Y1)*0.75+Y1
      DO I=1,NBL
        LGIBCF(I)=0
        LGIBC(I)=0
      ENDDO
C
C    Check the 3 sides of each triangle for a midpoint
C        on an external boundary
C
      DO 20 JEL=1,NEL
      DO 30 J1=1,3
      J2=MOD(J1,3)+1
      J3=MOD(J2,3)+1
      J5=J2+3
      LK1=LGEM(J1,JEL)
      LK2=LGEM(J2,JEL)
      LK3=LGEM(J3,JEL)
      LK5=LGEM(J5,JEL)
      IF(LK5.LT.0)LK5=-LK5
C
      XM=EXLG(LK5)
      YM=EYLG(LK5)
      IF((ABS(XM-X1).LT.1.E-4).OR.(ABS(YM-Y1).LT.1.E-4).OR.
     1      (ABS(XM-X2).LT.1.E-4).OR.(ABS(YM-Y2).LT.1.E-4))THEN
C
C   OK, the boundary segment is J1-J5-J2
C
        DO 100 I=1,NBL
          IB=LGIBC(I)
          IF(IB.EQ.0.OR.IB.EQ.LK5)THEN
            LGIBC(I)=LK5
            LGIBCF(I)=0
            GO TO 30
          END IF
  100   CONTINUE
      ELSE
        DO 80 J=1,NBP
         LVAL=11
         IA=IBC(J)
         IB=IBNGH(J)
         ITA=IBCTYP(J)
         ITB=0
         DO 70 JJ=1,NBP
           IF (IBC(JJ).EQ.IB) THEN
             ITB=IBCTYP(JJ)
             GO TO 60
           END IF
   70    CONTINUE
   60    IF(ITA.GE.10.AND.ITA.LE.14) THEN
           XA=EX(NOR(IA))
           YA=EY(NOR(IA))
           XB=EX(NOR(IB))
           YB=EY(NOR(IB))
           XMAX = XA
           YMAX = YA
           XMIN = XB
           YMIN = YB
           IF (XMIN.GT.XMAX) THEN
             TMP = XMIN
             XMIN = XMAX
             XMAX = TMP
           END IF
           IF (YMIN.GT.YMAX) THEN
             TMP = YMIN
             YMIN = YMAX
             YMAX = TMP
           END IF
           IF(XM.LE.XMAX.AND.XM.GE.XMIN.AND.
     1              YM.LE.YMAX.AND.YM.GE.YMIN) THEN
C  if we reach here - the boundary segment is on a fault!
C
             IF(ITA.EQ.13.OR.ITA.EQ.14.OR.
     1          ITB.EQ.13.OR.ITB.EQ.14) LVAL=13
             DO 90 I=1,NBL
                IB=LGIBC(I)
                IF(IB.EQ.0.OR.IB.EQ.LK5)THEN
                   LGIBC(I)=LK5
                   LGIBCF(I)=LVAL
                   GO TO 30
                END IF
   90        CONTINUE
           END IF
         END IF
   80 CONTINUE

      END IF
C     IERR=1
C     WRITE(LSC,10100)KS,XS,YS
C     WRITE(LUW,10100)KS,XM,YM
C     GO TO 250
   30 CONTINUE
   20 CONTINUE
  250 RETURN
10100 FORMAT('  NO LG.B.C. FOR N =',I5,' X =',G12.5,' Y =',
     1G12.5,//)
      END
C
      SUBROUTINE LGDEF(DT,UVP,UXLG,UYLG,UVELPX,UVELPY,
     :                 EX,EY,EXLG,EYLG,STELPX,STELPY,
     :                 ISTELP,LEM,NOR,
     :                 NUP,NE,NROWS,NUL,NNL,NML,NSM,NPM,
     :                 ILAG,IFLT)
C
C    Routine to calculate the position of the lagrangian nodes.
C    The lagrangian nodes are used to determine the strain or
C    total deformation.
C
C     DOUBLE PRECISION UVP
      DOUBLE PRECISION TRIA2,A0(3),DPEX2,DPEX3,DPEY2,DPEY3
      DIMENSION UVP(NROWS)
      DIMENSION UXLG(NUL),UYLG(NUL)
      DIMENSION EXLG(NUL),EYLG(NUL)
      DIMENSION UVELPX(NPM,NSM),UVELPY(NPM,NSM)
      DIMENSION STELPX(NPM,NSM),STELPY(NPM,NSM)
      DIMENSION ISTELP(NPM,NSM)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION BY0(3),CX0(3),COOL(3)

      DIMENSION BBOX(2,4)
      DOUBLE PRECISION DBBOX(8), DPT(2)
C
      EPS=1.E-8
      CALL FINDBBOX(EX,EY,BBOX,NUP)
      XLEN = BBOX(1,2)-BBOX(1,1)
      XMINR = BBOX(1,2)
      IF (BBOX(1,3).LT.XMINR) THEN
        XMAXR = XMINR
        XMINR = BBOX(1,3)
      ELSE
        XMAXR = BBOX(1,3)
      ENDIF
      XMAXL = BBOX(1,1)
      IF (BBOX(1,4).LT.XMAXL) THEN
        XMINL = XMAXL
        XMAXL = BBOX(1,4)
      ELSE
        XMINL = BBOX(1,4)
      ENDIF
      DBBOX(1) = BBOX(1,1)
      DBBOX(2) = BBOX(2,1)
      DBBOX(3) = BBOX(1,2)
      DBBOX(4) = BBOX(2,2)
      DBBOX(5) = BBOX(1,3)
      DBBOX(6) = BBOX(2,3)
      DBBOX(7) = BBOX(1,4)
      DBBOX(8) = BBOX(2,4)
      NVERT = 4

      IF (ILAG.EQ.1.OR.ILAG.EQ.2) THEN
C
C zero the ISTELP array
C
        DO 165 J=1,NSM
          DO 166 I=1,NPM
            ISTELP(I,J)=0
  166     CONTINUE
  165   CONTINUE
      END IF
C
C   update the velocity of each lagrangian point (uxlg and uylg)
C    first, search through every element in the finite element mesh
C
      DO 100 JEL=1,NE
C
C     First calculate the geometrical coefficients for JEL
C
         TRIA2=0.0
         DO 120 KK=1,3
            K2=MOD(KK,3)+1
            K3=MOD(KK+1,3)+1
            LK=NOR(LEM(KK,JEL))
            LK2=NOR(LEM(K2,JEL))
            LK3=NOR(LEM(K3,JEL))
            X1=EX(LK)
            X2=EX(LK2)
            X3=EX(LK3)
            Y1=EY(LK)
            Y2=EY(LK2)
            Y3=EY(LK3)
            DPEX2 = X2
            DPEX3 = X3
            DPEY2 = Y2
            DPEY3 = Y3
            A0(KK)=DPEX2*DPEY3-DPEX3*DPEY2
            BY0(KK)=Y2-Y3
            CX0(KK)=X3-X2
            TRIA2=TRIA2+A0(KK)
  120    CONTINUE
C
C     Now calculate the natural coordinates of Lagrange pt IL
C      CNL must be between 0 and 1 if NKJ is inside the element
C      We assume no points cross the fault - i.e. if JEL <= NEP/2,
C      we only look at the points <= NNLP/2 or NNLP<=points<=NNLP+NMLP/2
C      and if JEL > NEP/2, we look at the NNLP/2+1<=points<NNLP or
C      points > NNLP+NMLP/2
C
         IF((ILAG.EQ.1).OR.(ILAG.EQ.3))THEN
           DO 160 IL=1,NUL
             LFLT=0
C            these lines only apply to centred internal faults
             IF (IFLT.EQ.2) THEN
               IF(IL.GT.NNL/2.AND.IL.LE.NNL) LFLT=1
               IF(IL.GT.NNL+NML/2) LFLT=1
               IF(JEL.GT.NE/2.AND.LFLT.EQ.0) GOTO 160
               IF(JEL.LE.NE/2.AND.LFLT.EQ.1) GOTO 160
             END IF
             XP=EXLG(IL)
             YP=EYLG(IL)
C
C    these checks for elle runs using lg mesh to track finite strain
             IF ((IFLT.EQ.1.OR.Iflt.EQ.4).AND.XP.GT.XMINR) THEN
C    use XMINL to stop endless loop if test failing
               DPT(1) = XP
               DPT(2) = YP
   10          CALL CROSSINGSTEST(DBBOX,NVERT,DPT,IRES)
               IF (IRES.EQ.0.AND.XP.GT.XMINL) THEN
                 XP = XP-XLEN
                 DPT(1) = XP
                 GOTO  10
               ENDIF
             ELSE IF ((IFLT.EQ.1.OR.IFLT.EQ.4).AND.XP.LT.XMAXL) THEN
               DPT(1) = XP
               DPT(2) = YP
   20          CALL CROSSINGSTEST(DBBOX,NVERT,DPT,IRES)
               IF (IRES.EQ.0.AND.XP.LT.XMAXR) THEN
                 XP = XP+XLEN
                 DPT(1) = XP
                 GOTO  20
               ENDIF
             ENDIF

             DO 130 KK=1,3
               CNL=(A0(KK) + XP*BY0(KK) + YP*CX0(KK))/TRIA2
               IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 160
               COOL(KK)=CNL
  130       CONTINUE
C
C     If we reach here the point is within the triangle, so
C     get the interpolated values
C
            UXT=0.0
            UYT=0.0
            DO 135 KK=1,3
               LK=LEM(KK,JEL)
               UXV=UVP(LK)
               UYV=UVP(LK+NUP)
               K2=MOD(KK,3)+1
               K3=MOD(KK+1,3)+1
               CNK=COOL(KK)*(COOL(KK) - COOL(K2) - COOL(K3))
               UXT=UXT + UXV*CNK
               UYT=UYT + UYV*CNK
  135       CONTINUE
            DO 140 KK=4,6
               LK=IABS(LEM(KK,JEL))
               UXV=UVP(LK)
               UYV=UVP(LK+NUP)
               K2=KK-3
               K3=MOD(KK+1,3)+1
               CNK=4.0*COOL(K2)*COOL(K3)
               UXT=UXT + UXV*CNK
               UYT=UYT + UYV*CNK
  140       CONTINUE
            UXLG(IL)=UXT
            UYLG(IL)=UYT
  160      CONTINUE
         END IF
C
C  Now update the lagrangian strain markers
C
         DO 170 JL=1,NSM
           DO 180 IL=1,NPM
              XP=STELPX(IL,JL)
              YP=STELPY(IL,JL)
              DO 185 KK=1,3
                 CNL=(A0(KK) + XP*BY0(KK) + YP*CX0(KK))/TRIA2
                 IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 180
                 COOL(KK)=CNL
  185         CONTINUE
C
C     If we reach here the point is within the triangle, so
C     get the interpolated values
C     ISTELP becomes 1 if point is within a triangle
              ISTELP(IL,JL)=1
              UXT=0.0
              UYT=0.0
              DO 190 KK=1,3
                 LK=LEM(KK,JEL)
                 UXV=UVP(LK)
                 UYV=UVP(LK+NUP)
                 K2=MOD(KK,3)+1
                 K3=MOD(KK+1,3)+1
                 CNK=COOL(KK)*(COOL(KK) - COOL(K2) - COOL(K3))
                 UXT=UXT + UXV*CNK
                 UYT=UYT + UYV*CNK
  190         CONTINUE
              DO 195 KK=4,6
                 LK=IABS(LEM(KK,JEL))
                 UXV=UVP(LK)
                 UYV=UVP(LK+NUP)
                 K2=KK-3
                 K3=MOD(KK+1,3)+1
                 CNK=4.0*COOL(K2)*COOL(K3)
                 UXT=UXT + UXV*CNK
                 UYT=UYT + UYV*CNK
  195         CONTINUE
              UVELPX(IL,JL)=UXT
              UVELPY(IL,JL)=UYT
  180        CONTINUE
  170      CONTINUE
  100 CONTINUE

C Call subr to test if there are any points which are outside the domain
      CALL MARKADJ(STELPX,STELPY,ISTELP,NSM,NPM)
C
C   now update the position of each lagrangian mesh point (exlg and eylg)
C
      IF((ILAG.EQ.1).OR.(ILAG.EQ.3))THEN
      DO 200 I=1,NUL
         EXLG(I)=EXLG(I)+UXLG(I)*DT
         EYLG(I)=EYLG(I)+UYLG(I)*DT
  200 CONTINUE
      END IF
C
C   now update the position of each marker point (stelpx and stelpy)
C
      IF (ILAG.EQ.1.OR.ILAG.EQ.2) THEN
        DO 220 JL=1,NSM
           DO 210 IL=1,NPM
              STELPX(IL,JL)=STELPX(IL,JL)+UVELPX(IL,JL)*DT
              STELPY(IL,JL)=STELPY(IL,JL)+UVELPY(IL,JL)*DT
  210      CONTINUE
  220    CONTINUE
      END IF
C
      RETURN
      END

      SUBROUTINE MARKADJ(STELPX,STELPY,ISTELP,NSM,NPM)
C
C    Routine to test if any points within the strain markers are outside the
C    domain, in which case they will be relocated inside the domain
C
      DIMENSION STELPX(NPM,NSM),STELPY(NPM,NSM)
      DIMENSION ISTELP(NPM,NSM)

      DO 230 JL=1,NSM
      K=1
        DO WHILE ((K.LE.NPM+1).AND.(ISTELP(K,JL).EQ.0))
            K=K+1
        ENDDO
        IF (K.GT.NPM) THEN
          DO 240 IL=1,NPM
             STELPX(IL,JL)=-99999.0
             STELPY(IL,JL)=-99999.0
  240     CONTINUE
        ELSEIF (K.GT.1.AND.K.LT.NPM) THEN
          DO 250 I=1,K-1
                 STELPX(I,JL)=STELPX(K,JL)
                 STELPY(I,JL)=STELPY(K,JL)
  250     CONTINUE
          DO 260 I=K+1,NPM
             IF (ISTELP(I,JL).EQ.0) THEN
                STELPX(I,JL)=STELPX(I-1,JL)
                STELPY(I,JL)=STELPY(I-1,JL)
             ENDIF
  260     CONTINUE
        ELSEIF (K.EQ.1) THEN
          DO 270 I=K+1,NPM
             IF (ISTELP(I,JL).EQ.0) THEN
               STELPX(I,JL)=STELPX(I-1,JL)
                STELPY(I,JL)=STELPY(I-1,JL)
            ENDIF
  270     CONTINUE
        ELSE !K=NPM
          DO 280 I=1,K-1
                 STELPX(I,JL)=STELPX(K,JL)
                 STELPY(I,JL)=STELPY(K,JL)
  280     CONTINUE
        ENDIF

  230 CONTINUE

      RETURN
      END
