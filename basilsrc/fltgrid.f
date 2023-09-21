C*--------------------------------------------------------------------
C*    Basil / Sybil:   fltgrid.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE REGRID(HLENSC,EX,EY,EXO,EYO,SSQ,FROT,
     :                  NOR,LEM,IFBC,IFEQV,IBCTYP,IENDP1,IENDP2,
C    :                  EX1,EY1,SSQ1,FROT1,JNOALT,JNOAL2,
     :                  NUP,NE,NN,NBP,NFP)
C
C  Subroutine to regrid after each iteration
C
C  Fault nodes are updated first, then outer boundary nodes
C  Next all internal nodes are regridded using the the initial
C  fractional distance between fault nodes and external
C  boundary nodes.
C
C  This routine also reinterpolates the velocity, pressure
C   crustal thickness values etc., for any nodes that have been moved
C
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION EXO(NUP),EYO(NUP)
      DIMENSION IENDP1(NN),IENDP2(NN)
      DIMENSION IBCTYP(NBP*2)
      DIMENSION IFBC(NFP),IFEQV(NFP)
      DIMENSION SSQ(NUP),FROT(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
C      local arrays
      DIMENSION EX1(NUP),EY1(NUP)
      DIMENSION SSQ1(NUP),FROT1(NUP)
      DIMENSION JNOALT(NUP),JNOAL2(NUP)
      DIMENSION A0(3),BY0(3),CX0(3),COOL(3)
      DIMENSION NJ(6)
C
      EPS=1.E-4
      DO 5 JJ=1,NUP
          EX1(JJ)=EX(JJ)
          EY1(JJ)=EY(JJ)
          IF (ISSQACTIVE.EQ.1)SSQ1(JJ)=SSQ(JJ)
          IF (IROTACTIVE.EQ.1)FROT1(JJ)=FROT(JJ)
    5 CONTINUE
      JNONUM=0
      DO 6 JJ=1,NUP
          JNOALT(JJ)=0
          JNOAL2(JJ)=0
    6 CONTINUE

C
C    Determine if fault node must be moved
C     The array JNOALT contains all the node number on
C     and off the fault that have been moved. There are
C     JNONUM of these nodes.
C
      DO 10 I=1,NFP
         IP1=ABS(IFBC(I))
         IP2=ABS(IFBC(IFEQV(I)))
         ISIDE=IP1/IFBC(I)
         IF(ISIDE.GT.0) THEN
             X1=EX(NOR(IP1))
             Y1=EY(NOR(IP1))
             X2=EX(NOR(IP2))
             Y2=EY(NOR(IP2))
             DELX=ABS(X1-X2)
             DELY=ABS(Y1-Y2)
             IF(DELX.GT.EPS.OR.DELY.GT.EPS) THEN
                  X3=0.5*(X1+X2)
                  Y3=0.5*(Y1+Y2)
                  EX(NOR(IP1))=X3
                  EY(NOR(IP1))=Y3
                  EX(NOR(IP2))=X3
                  EY(NOR(IP2))=Y3
                  JNOALT(JNONUM+1)=IP1
                  JNOALT(JNONUM+2)=IP2
                  JNONUM=JNONUM+2
             END IF
         END IF
C corners
         ITYPE = IBCTYP(IP1)
         IF (ITYPE.EQ.10.OR.ITYPE.EQ.12.OR.ITYPE.EQ.14)THEN
            EX(NOR(IP1))=EXO(NOR(IP1))
            EY(NOR(IP1))=EYO(NOR(IP1))
         END IF
         ITYPE = IBCTYP(IP2)
         IF (ITYPE.EQ.10.OR.ITYPE.EQ.12.OR.ITYPE.EQ.14)THEN
            EX(NOR(IP2))=EXO(NOR(IP2))
            EY(NOR(IP2))=EYO(NOR(IP2))
         END IF
   10 CONTINUE
C
C  Calculate all the position of all vertex nodes off the fault
C
      DO 15 I=1,NUP
       IF(NOR(I).LE.NN) THEN
         XO=EXO(NOR(I))
         YO=EYO(NOR(I))
         XO1=EXO(NOR(IENDP1(NOR(I))))
         YO1=EYO(NOR(IENDP1(NOR(I))))
         XO2=EXO(NOR(IENDP2(NOR(I))))
         YO2=EYO(NOR(IENDP2(NOR(I))))
         XN1=XO1
         YN1=YO1
         XN2=EX(NOR(IENDP2(NOR(I))))
         YN2=EY(NOR(IENDP2(NOR(I))))
         DXO=XO1-XO2
         DYO=YO1-YO2
         DTO=SQRT(DXO*DXO+DYO*DYO)
         DXO=XO-XO2
         DYO=YO-YO2
         DFO=SQRT(DXO*DXO+DYO*DYO)
         IF(DTO.GT.0.0) THEN
            XN=DFO/DTO*(XN1-XN2)+XN2
         ELSE
            XN=XN2
         END IF
         IF(DTO.GT.0.0) THEN
            YN=DFO/DTO*(YN1-YN2)+YN2
         ELSE
            YN=YN2
         END IF
         DELX=ABS(XN-EX(NOR(I)))
         DELY=ABS(YN-EY(NOR(I)))
         IF(DELX.GT.EPS.OR.DELY.GT.EPS) THEN 
             EX(NOR(I))=XN
             EY(NOR(I))=YN
             DO 20 JJ=1,NUP
                     K=JNOALT(JJ)
                     IF(K.EQ.I) GOTO 17
                     IF(K.EQ.0) THEN
                         JNOALT(JJ)=I
                         JNONUM=JNONUM+1
                         GOTO 17
                     END IF
   20             CONTINUE
   17    END IF
       END IF
   15 CONTINUE
C
C       Calculate position of midpoint nodes
C  
      DO 25 JEL=1,NE
         NJ(1)=NOR(LEM(1,JEL))
         NJ(2)=NOR(LEM(2,JEL))
         NJ(3)=NOR(LEM(3,JEL))
         DO 30 I=4,6
            IA=I-3
            IB=MOD(I-2,3)+1
            IJ=ABS(LEM(I,JEL))
            NJ(I)=NOR(ABS(LEM(I,JEL)))
            X3=(EX(NJ(IA))+EX(NJ(IB)))*0.5
            Y3=(EY(NJ(IA))+EY(NJ(IB)))*0.5
            DELX=ABS(X3-EX1(NJ(I)))
            DELY=ABS(Y3-EY1(NJ(I)))
            IF(DELX.GT.EPS.OR.DELY.GT.EPS) THEN
                  EX(NJ(I))=X3
                  EY(NJ(I))=Y3
                  DO 50 JJ=1,NUP
                     K=JNOALT(JJ)
                     IF(K.EQ.IJ) GOTO 40
                     IF(K.EQ.0) THEN
                         JNOALT(JJ)=IJ
                         JNONUM=JNONUM+1
                         GOTO 40
                     END IF
   50             CONTINUE
   40       END IF
   30    CONTINUE
   25 CONTINUE
C
C    For each element, check if any of the moved nodes fall within it
C
      DO 60 JEL=1,NE
C
C     If one node in this element has been moved, check all
C     moved nodes, to see whether they are now within this element
C
C         DO 70 KK=1,6
C            NKK=IABS(LEM(KK,JEL))
C            DO 80 II=1,JNONUM
C               IF(NKK.EQ.JNOALT(II))GO TO 100
C   80       CONTINUE
C   70    CONTINUE
C         GO TO 60
C  100    CONTINUE
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
            X1=EX1(LK)
            X2=EX1(LK2)
            X3=EX1(LK3)
            Y1=EY1(LK)
            Y2=EY1(LK2)
            Y3=EY1(LK3)
            A0(KK)=X2*Y3-X3*Y2
            BY0(KK)=Y2-Y3
            CX0(KK)=X3-X2
            TRIA2=TRIA2+A0(KK)
  120    CONTINUE
C
C     Now calculate the natural coordinates of (new) NKJ
C      CNL must be between 0 and 1 if NKJ is inside the element
C
         DO 160 KJN=1,JNONUM
            IF(JNOAL2(KJN).EQ.1) GOTO 160
            NKJ=JNOALT(KJN)
            XP=EX(NOR(NKJ))
            YP=EY(NOR(NKJ))
            DO 130 KK=1,3
               CNL=(A0(KK) + XP*BY0(KK) + YP*CX0(KK))/TRIA2
               IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 160
               COOL(KK)=CNL
  130       CONTINUE
C
C     If we reach here the point is within the triangle, so
C     get the interpolated values
C
            SI=0.0
            FI=0.0
            SK=DFLTSSQ
            FK=0.0
            DO 135 KK=1,3
               LK=LEM(KK,JEL)
               IF(ISSQACTIVE.EQ.1)SK=SSQ1(LK)
               IF(IROTACTIVE.EQ.1)FK=FROT1(LK)
               K2=MOD(KK,3)+1
               K3=MOD(KK+1,3)+1
               CNK=COOL(KK)*(COOL(KK) - COOL(K2) - COOL(K3))
               SI=SI + SK*CNK
               FI=FI + FK*CNK
  135       CONTINUE
            DO 140 KK=4,6
               LK=IABS(LEM(KK,JEL))
               IF(ISSQACTIVE.EQ.1)SK=SSQ1(LK)
               IF(IROTACTIVE.EQ.1)FK=FROT1(LK)
               K2=KK-3
               K3=MOD(KK+1,3)+1
               CNK=4.0*COOL(K2)*COOL(K3)
               SI=SI + SK*CNK
               FI=FI + FK*CNK
  140       CONTINUE
            IF(ISSQACTIVE.EQ.1)SSQ(NKJ)=SI
            IF(IROTACTIVE.EQ.1)FROT(NKJ)=FI
            JNOAL2(KJN)=1
C      WRITE(LSC,10020)NKJ,JEL,XP,YP,SSQ1(NKJ),SI,FROT1(NKJ),FI
C      WRITE(LUW,10020)NKJ,JEL,XP,YP,SSQ1(NKJ),SI,FROT1(NKJ),FI
C
C   go back to check other nodes in this element and other
C    elements
C
  160    CONTINUE
   60 CONTINUE
C
      DO 170 II=1,JNONUM
         JJ=JNOALT(II)
         IF(JNOAL2(II).EQ.0)THEN
             IF(ISSQACTIVE.EQ.1) SSQ(JJ)=-ALOG(HLENSC)
             IF(IROTACTIVE.EQ.1) FROT(JJ)=0.0
C            WRITE(LSC,10010)JNOALT(II)
C            WRITE(LUW,10010)JNOALT(II)
C             STOP
         END IF
  170 CONTINUE
10010 FORMAT(' Node no',I5,' has been moved but SSQ and FROT',
     1' have been set to initial values')
10020 FORMAT(' Node/element ',2I5,' (x,y) =',2F9.5,/,
     1' old/new S',2G12.5,'  old/new F',2G12.5)
      RETURN
      END
C
      SUBROUTINE GRINIT(NX,NY,XLEN,YLEN,EX,EY,NOR,IENDP1,IENDP2,
     :                  NUP,NN,LUW,LSC,IERR)
C
C   This routine initializes the necessary arrays for a general
C    regridding
C    Arrays IENDP1 and IENDP2 are set
C
C     COMMON/AI/LUW,LSC,LBC,LLG
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION IENDP1(NN),IENDP2(NN)
C       local arrays
      DIMENSION M1(NY+1),M2(NY+1),M3(NY+1),M4(NY+1)
C
      EPS=1.E-4
      NX1=NX+1
      NY1=NY+1
      DY=YLEN/FLOAT(NY)
      DX=XLEN/FLOAT(NX)
      DO 5 I=1,NY1
         M1(I)=0
         M2(I)=0
         M3(I)=0
         M4(I)=0
    5 CONTINUE
      DO 6 I=1,NN
        IENDP1(I)=0
        IENDP2(I)=0
    6 CONTINUE
C
      DO 10 I=1,NUP
       IF(NOR(I).LE.NN) THEN
        DELX0=ABS(EX(NOR(I))-0.0)
        DELX1=ABS(EX(NOR(I))-XLEN)
        DELX2=ABS(EX(NOR(I))-2.0*XLEN)
        IF(DELX0.LT.EPS) THEN
          ISET=0
          DO 20 II=1,NY1
            IF(ISET.EQ.0.AND.M1(II).EQ.0) THEN
              M1(II)=I
              ISET=1
C              II=NY1
            END IF 
   20     CONTINUE
          IF(ISET.EQ.0) THEN
              WRITE(LSC,10010)I
              WRITE(LUW,10010)I
              IERR=1
              GO TO 200
          END IF
        ELSE IF(DELX1.LT.EPS) THEN
          ISET=0
          IF(NOR(I).LE.NN/2) THEN
          DO 30 II=1,NY1
            IF(ISET.EQ.0.AND.M2(II).EQ.0) THEN
              M2(II)=I
              ISET=1
C              II=NY1
            END IF 
   30     CONTINUE
          ELSE
          DO 35 II=1,NY1
            IF(ISET.EQ.0.AND.M3(II).EQ.0) THEN
              M3(II)=I
              ISET=1
C              II=NY1
            END IF 
   35     CONTINUE
          END IF
          IF(ISET.EQ.0) THEN
              WRITE(LSC,10010)I
              WRITE(LUW,10010)I
              IERR=1
              GO TO 200
          END IF
        ELSE IF(DELX2.LT.EPS) THEN
          ISET=0
          DO 40 II=1,NY1
            IF(ISET.EQ.0.AND.M4(II).EQ.0) THEN
              M4(II)=I
              ISET=1
C              II=NY1
            END IF 
   40     CONTINUE
          IF(ISET.EQ.0) THEN
              WRITE(LSC,10010)I
              WRITE(LUW,10010)I
              IERR=1
              GO TO 200
          END IF
        END IF
       END IF
   10 CONTINUE
      DO 110 I=1,NUP
       IF(NOR(I).LE.NN) THEN
        YY=EY(NOR(I))
        IF(NOR(I).LE.NN/2) THEN
          DO 120 II=1,NY1
             Y1=EY(NOR(M1(II)))
             DELY=ABS(Y1-YY)
             IF(DELY.LT.EPS) THEN
                IENDP1(NOR(I))=M1(II)
C                II=NY1
             END IF
  120     CONTINUE
          DO 130 II=1,NY1
             Y1=EY(NOR(M2(II)))
             DELY=ABS(Y1-YY)
             IF(DELY.LT.EPS) THEN
                IENDP2(NOR(I))=M2(II)
C                II=NY1
             END IF
  130     CONTINUE
        ELSE
          DO 140 II=1,NY1
             Y1=EY(NOR(M4(II)))
             DELY=ABS(Y1-YY)
             IF(DELY.LT.EPS) THEN
                IENDP1(NOR(I))=M4(II)
C                II=NY1
             END IF
  140     CONTINUE
          DO 150 II=1,NY1
             Y1=EY(NOR(M3(II)))
             DELY=ABS(Y1-YY)
             IF(DELY.LT.EPS) THEN
                IENDP2(NOR(I))=M3(II)
C                II=NY1
             END IF
  150     CONTINUE
        END IF
       END IF
  110 CONTINUE
      DO 160 I=1,NN
         IF(IENDP1(I).LE.0) THEN
            WRITE(LSC,10020)I
            WRITE(LUW,10020)I
            IERR=1
            GO TO 200
         END IF
         IF(IENDP2(I).LE.0) THEN
            WRITE(LSC,10030)I
            WRITE(LUW,10030)I
            IERR=1
            GO TO 200
         END IF
  160 CONTINUE
  200 RETURN
10010 FORMAT(' Boundary point ',I4,' not set in GRINIT')
10020 FORMAT(' IENDP1 point ',I4,' not set in GRINIT')
10030 FORMAT(' IENDP2 point ',I4,' not set in GRINIT')
      END
C
      SUBROUTINE FLTGRD(HLENSC,TBXOFF,TBYOFF,EX,EY,EXO,EYO,SSQ,FROT,
     :                  NOR,LEM,IFBC,IFEQV,
C    :                  EX1,EY1,SSQ1,FROT1,JNOALT,JNOAL2,
     :                  NUP,NE,NN,NFP,LUW,LSC)
C
C  Old subroutine - still used for the tilt block model until
C    REGRID is modified if necessary
C
C    to be called following CRUST to regrid any nodes along
C    the fault and to re-interpolate crustal thickness
C    values etc., for any nodes that have been moved
C
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION EXO(NUP),EYO(NUP)
      DIMENSION IFBC(NFP),IFEQV(NFP)
      DIMENSION SSQ(NUP),FROT(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
C       local arrays
      DIMENSION EX1(NUP),EY1(NUP)
      DIMENSION SSQ1(NUP),FROT1(NUP)
      DIMENSION JNOALT(NFP*2),JNOAL2(NFP*2)
      DIMENSION A0(3),BY(3),CX(3),COOL(3)
      DIMENSION NJ(6)
C
      NFP2 = NFP*2
      EPS=1.E-4
      DO 5 JJ=1,NUP
          EX1(JJ)=EX(JJ)
          EY1(JJ)=EY(JJ)
          IF(ISSQACTIVE.EQ.1)SSQ1(JJ)=SSQ(JJ)
          IF(IROTACTIVE.EQ.1)FROT1(JJ)=FROT(JJ)
    5 CONTINUE
      JNONUM=0
      DO 6 JJ=1,NFP2
          JNOALT(JJ)=0
          JNOAL2(JJ)=0
    6 CONTINUE

C
C    Determine if fault node must be moved
C     The array JNOALT contains all the node number on
C     and off the fault that have been moved. There are
C     JNONUM of these nodes.
C
      DO 10 I=1,NFP
         IP1=ABS(IFBC(I))
         IP2=ABS(IFBC(IFEQV(I)))
         ISIDE=IP1/IFBC(I)
         IF(ISIDE.GT.0) THEN
             X1=EX(NOR(IP1))
             Y1=EY(NOR(IP1))
             X2=EX(NOR(IP2))
             Y2=EY(NOR(IP2))
             DELX=ABS(X1-X2)
             DELY=ABS(Y1-Y2)
             IF(DELX.GT.EPS.OR.DELY.GT.EPS) THEN
                  X3=0.5*(X1+X2)
                  Y3=0.5*(Y1+Y2)
                  EX(NOR(IP1))=X3+TBXOFF/2.0
                  EY(NOR(IP1))=Y3+TBYOFF/2.0
                  EX(NOR(IP2))=X3-TBXOFF/2.0
                  EY(NOR(IP2))=Y3-TBYOFF/2.0
                  JNOALT(JNONUM+1)=IP1
                  JNOALT(JNONUM+2)=IP2
                  JNONUM=JNONUM+2
             END IF
         END IF
   10 CONTINUE

C
C       Calculate position of midpoint nodes
C  
      DO 20 JEL=1,NE
         NJ(1)=NOR(LEM(1,JEL))
         NJ(2)=NOR(LEM(2,JEL))
         NJ(3)=NOR(LEM(3,JEL))
         DO 30 I=4,6
            IA=I-3
            IB=MOD(I-2,3)+1
            IJ=ABS(LEM(I,JEL))
            NJ(I)=NOR(ABS(LEM(I,JEL)))
            X3=(EX(NJ(IA))+EX(NJ(IB)))*0.5
            Y3=(EY(NJ(IA))+EY(NJ(IB)))*0.5
            DELX=ABS(X3-EX1(NJ(I)))
            DELY=ABS(Y3-EY1(NJ(I)))
            IF(DELX.GT.EPS.OR.DELY.GT.EPS) THEN
                  EX(NJ(I))=X3
                  EY(NJ(I))=Y3
                  DO 50 JJ=1,NFP2
                     K=JNOALT(JJ)
                     IF(K.EQ.IJ) GOTO 40
                     IF(K.EQ.0) THEN
                         JNOALT(JJ)=IJ
                         JNONUM=JNONUM+1
                         GOTO 40
                     END IF
   50             CONTINUE
   40       END IF
   30    CONTINUE
   20 CONTINUE
C
C    For each element, check whether any of its nodes has been moved
C
      DO 60 JEL=1,NE
         DO 70 KK=1,6
            NKK=IABS(LEM(KK,JEL))
            DO 80 II=1,JNONUM
               IF(NKK.EQ.JNOALT(II))GO TO 100
   80       CONTINUE
   70    CONTINUE
         GO TO 60
C
C     If one node in this element has been moved, check all
C     moved nodes, to see whether they are now within this element
C
  100    CONTINUE
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
            X1=EX1(LK)
            X2=EX1(LK2)
            X3=EX1(LK3)
            Y1=EY1(LK)
            Y2=EY1(LK2)
            Y3=EY1(LK3)
            A0(KK)=X2*Y3-X3*Y2
            BY(KK)=Y2-Y3
            CX(KK)=X3-X2
            TRIA2=TRIA2+A0(KK)
  120    CONTINUE
C
C     Now calculate the natural coordinates of (new) NKJ
C      CNL must be between 0 and 1 if NKJ is inside the element
C
         DO 160 KJN=1,JNONUM
            IF(JNOAL2(KJN).EQ.1) GOTO 160
            NKJ=JNOALT(KJN)
            XP=EX(NOR(NKJ))
            YP=EY(NOR(NKJ))
            DO 130 KK=1,3
               CNL=(A0(KK) + XP*BY(KK) + YP*CX(KK))/TRIA2
               IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 160
               COOL(KK)=CNL
  130       CONTINUE
C
C     If we reach here the point is within the triangle, so
C     get the interpolated values
C
            SI=0.0
            FI=0.0
            SK=DFLTSSQ
            FK=0.0
            DO 135 KK=1,3
               LK=LEM(KK,JEL)
               IF(ISSQACTIVE.EQ.1)SK=SSQ1(LK)
               IF(IROTACTIVE.EQ.1)FK=FROT1(LK)
               K2=MOD(KK,3)+1
               K3=MOD(KK+1,3)+1
               CNK=COOL(KK)*(COOL(KK) - COOL(K2) - COOL(K3))
               SI=SI + SK*CNK
               FI=FI + FK*CNK
  135       CONTINUE
            DO 140 KK=4,6
               LK=IABS(LEM(KK,JEL))
               IF(ISSQACTIVE.EQ.1)SK=SSQ1(LK)
               IF(IROTACTIVE.EQ.1)FK=FROT1(LK)
               K2=KK-3
               K3=MOD(KK+1,3)+1
               CNK=4.0*COOL(K2)*COOL(K3)
               SI=SI + SK*CNK
               FI=FI + FK*CNK
  140       CONTINUE
            IF(ISSQACTIVE.EQ.1)SSQ(NKJ)=SI
            IF(IROTACTIVE.EQ.1)FROT(NKJ)=FI
            JNOAL2(KJN)=1
C      WRITE(LSC,10020)NKJ,JEL,XP,YP,SSQ1(NKJ),SI,FROT1(NKJ),FI
C      WRITE(LUW,10020)NKJ,JEL,XP,YP,SSQ1(NKJ),SI,FROT1(NKJ),FI
C
C   go back to check other nodes in this element and other
C    elements
C
  160    CONTINUE
   60 CONTINUE
C
      DO 170 II=1,JNONUM
         JJ=JNOALT(II)
         IF(JNOAL2(II).EQ.0)THEN
             IF(ISSQACTIVE.EQ.1)SSQ(JJ)=-ALOG(HLENSC)
             IF(IROTACTIVE.EQ.1)FROT(JJ)=0.0
C             WRITE(LSC,10010)JNOALT(II)
C             WRITE(LUW,10010)JNOALT(II)
C             STOP
         END IF
  170 CONTINUE
10010 FORMAT(' Node no',I5,' has been moved but SSQ and FROT',
     1' have been set to initial values')
10020 FORMAT(' Node/element ',2I5,' (x,y) =',2F9.5,/,
     1' old/new S',2G12.5,'  old/new F',2G12.5)
      RETURN
      END
