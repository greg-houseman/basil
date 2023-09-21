C*--------------------------------------------------------------------
C*    Basil / Sybil:   deform.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE MESHDF(IM,MPE,INAME,EX,EY,EXREF,EYREF,LEM,NOR,
     1                WORK,NE,NUP,NY,YCMIN,YCMAX,IVERBOSE)
C
C     Y1 CONTAINS COORDINATES OF UNDEFORMED ELEMENT
C
C     Y2 CONTAINS COORDINATES OF DEFORMED   ELEMENT
C
C     ORDER OF COORDINATES: X(PT1),Y(PT1), ....... ,X(PT3),Y(PT3)
C
C     N ALLOWS CHOICE OF RELATIVE POSITION VECTORS TO
C
C     BE PLACED IN X1 AND X2
C
C     NNN = 1:      (PT2 - PT1, PT3 - PT2) (PT2 is reference pt)
C
C     NNN = 2:      (PT2 - PT1, PT3 - PT1) (PT1 is reference pt)
C
C     NNN = 3:      (PT3 - PT2, PT3 - PT1) (PT3 is reference pt)
C
C
      CHARACTER INAME*4
      DIMENSION EX(NUP),EY(NUP),EXREF(NUP),EYREF(NUP)
      DIMENSION LEM(6,NE),NOR(NUP)
      DIMENSION WORK(NUP)
C     DIMENSION WORK(NWK)
C     COMMON/CONTRO/LUW,LRBAT,LRINT,LDIN,LWINT,LPLOT,IAUTO,ILABEL
C    1 ,ICOLO
      DIMENSION Y1(6),Y2(6),X1(4),X2(4),F(4),N1(4),N2(4)
      DATA PI /3.14159/
C
      YREF=YCMIN+YCMAX
C     R0=FLOAT(MPE)*YREF/(6.0*FLOAT(LY1-1))
C     R0=FLOAT(MPE)*YREF/(6.0*FLOAT(NY))
      R0=FLOAT(MPE)*(YCMAX-YCMIN)/(6.0*FLOAT(NY))
      IF (IVERBOSE.NE.0) THEN
        WRITE(6,*)'Data for ellipses'
        WRITE(6,*)'Element     x0          y0          a1     ',
     1  '     a2          a3         theta'
      END IF
      DO 60 N=1,NE,MPE
C
C    Put the element coordinates into the Y1 and Y2 matrices
C
      DO I=1,3
        MEL=NOR(LEM(I,N))
        JJ=2*(I-1)+1
        Y1(JJ)=EXREF(MEL)
        Y1(JJ+1)=EYREF(MEL)
        Y2(JJ)=EX(MEL)
        Y2(JJ+1)=EY(MEL)
      ENDDO
      NNN=1
      DO I = 1,4
        N1(I) = I + 2
        N2(I) = I
      ENDDO
C     GOTO (2,3,4) NNN
      IF(NNN.EQ.2)THEN
 3      N2(3) = 1
        N2(4) = 2
C       GOTO 2
      ELSEIF(NNN.EQ.3)THEN
 4      N1(1) = 5
        N1(2) = 6
      ELSE
C2    CONTINUE
        DO I = 1,4
          X1(I) = Y1(N1(I)) - Y1(N2(I))
          X2(I) = Y2(N1(I)) - Y2(N2(I))
        ENDDO
      ENDIF
C
C     PLOT DEFORMED GRID
C
      IF((IM.NE.2).AND.(IM.NE.3))GO TO 10
      CALL PLOTU(Y2(1),Y2(2),3)
      CALL PLOTU(Y2(3),Y2(4),2)
      CALL PLOTU(Y2(5),Y2(6),2)
      CALL PLOTU(Y2(1),Y2(2),2)
C
C     CALCULATE ELEMENTS OF STRAIN TENSOR
C
C       Calculated as:
C       F(1) = Exx = du/dx+1
C       F(2) = Exy = du/dy
C       F(3) = Eyx = dv/dx
C       F(4) = Eyy = dv/dy+1
C       If there is no deformation, strain ellipse is a circle with
C       radius of 1 and F(1) = F(4) = 1
   10 DET = X1(1)*X1(4) - X1(2)*X1(3)
      F(1) = (X2(1)*X1(4) - X1(2)*X2(3))/DET
      F(2) = (X1(1)*X2(3) - X2(1)*X1(3))/DET
      F(3) = (X2(2)*X1(4) - X1(2)*X2(4))/DET
      F(4) = (X1(1)*X2(4) - X2(2)*X1(3))/DET
C
C     CALCULATE ORIENTATION OF PRINCIPAL AXES BEFORE DEFORMATION
C
C       AL1 is STRAIN "in plane" MIN
C       AL2 is STRAIN "in plane" MAX
C       AL3 is STRAIN vertical
      BOT = F(1)*F(1) + F(3)*F(3) - F(2)*F(2) - F(4)*F(4)
      TOP = 2.*(F(1)*F(2) + F(3)*F(4))
      THETA1 = PI/2
      IF(ABS(BOT).GT.1.E-37) THETA1 = 0.5*ATAN(TOP/BOT)
      C = COS(THETA1)
      S = SIN(THETA1)
      AL1 = SQRT((F(1)*C + F(2)*S)**2 + (F(3)*C + F(4)*S)**2)
      TOP = F(3)*C + F(4)*S
      BOT = F(1)*C + F(2)*S
      THETA2 = THETA1 + PI/2
      C = COS(THETA2)
      S = SIN(THETA2)
      AL2 = SQRT((F(1)*C + F(2)*S)**2 + (F(3)*C + F(4)*S)**2)
      AL3 = 1.0/(AL1*AL2)
C
C     Calculate orientation of principal axes after deformation
C
      TH1 = PI/2
      IF(ABS(BOT).GT.1.E-37) TH1 = ATAN(TOP/BOT)
      TH2 = TH1 + PI/2
      ROT = TH1 - THETA1
      IF(AL2.GT.AL1)DEF=ALOG10(AL2/AL1)
      IF(AL2.LE.AL1)DEF=ALOG10(AL1/AL2)
      NC=2*(N-1) + 1
      WORK(NC)=DEF
      WORK(NC+1)=ROT*180.0/PI
      IF(IM.GE.3)GO TO 60
      XO = (Y2(1) + Y2(3) + Y2(5))/3
      YO = (Y2(2) + Y2(4) + Y2(6))/3
      IF (IVERBOSE.NE.0) THEN
        TH1DEG=TH1*180.0/PI
        WRITE(6,10555)N,XO,YO,AL1,AL2,AL3,TH1DEG
10555   FORMAT(I6,6G13.5)
      END IF
      CALL ELLPLT(XO,YO,R0,TH1,AL1,AL2,YREF)
   60 CONTINUE
C     IF(ILABEL.EQ.0)RETURN
      RETURN
      END

      SUBROUTINE ELLPLT(XO,YO,R0,THETA,AL1,AL2,YREF)
      DATA PI /3.14159/
      A1 = AL1*AL1
      A2 = AL2*AL2
      R2 = R0*R0
      IUP = 3
      DO 1 I = 1,46
      THE = PI/22.5*(I-1)
      TH1 = THE - THETA
      C = COS(TH1)
      S = SIN(TH1)
      R = R2/(C*C/A1 + S*S/A2)
      C = COS(THE)
      S = SIN(THE)
      R = SQRT(R)
      X = R*C + XO
      Y = R*S + YO
      CALL PLOTU(X,Y,IUP)
      IUP = 2
 1    CONTINUE
      RETURN
      END

      SUBROUTINE FDPUT(ICHO,S,WORK,EX,EY,NOR,LEM,
     1                NE,NUP,NP3,NX,NX2,NX3,NY2,NY3)
C
C     Routine to transfer data for elements (stored in WORK)
C     onto regular mesh (S).  but transfer is done only for points
C     that fall within some element of the mesh.  No interpolation
C     is done, value assumed constant within an element
C
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION LEM(6,NE),NOR(NUP)
      DIMENSION WORK(NE*2),S(NP3)
      DOUBLE PRECISION BY(3),CX(3),CNL,TRI,X1,X2,X3,Y1,Y2,Y3
C     DOUBLE PRECISION A0(3),COOL(3)
C     DOUBLE PRECISION DPX1,DPX2,DPX3,DPY1,DPY2,DPY3
C     DATA EPS/1.E-4/
C
C    Zero the mesh first
C
      HS=1.0/FLOAT(NX)
      EPS=HS*0.01
      NP=NX3*NY3
      DO I=1,NP3
        S(I)=0.0
      ENDDO
C
C    Look at each element in turn
C
      DO 60 N=1,NE
C
C     Calculate the geometrical coefficients
C
        XMIN=999.
        YMIN=999.
        XMAX=-999.
        YMAX=-999.
        DO 20 K=1,3
          K2=MOD(K,3)+1
          K3=MOD(K+1,3)+1
          LK=NOR(LEM(K,N))
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          X1=EX(LK)
          X2=EX(LK2)
          X3=EX(LK3)
          Y1=EY(LK)
          Y2=EY(LK2)
          Y3=EY(LK3)
          IF(X1.LT.XMIN)XMIN=X1
          IF(X1.GT.XMAX)XMAX=X1
          IF(Y1.LT.YMIN)YMIN=Y1
          IF(Y1.GT.YMAX)YMAX=Y1
          BY(K)=Y2    - Y3
          CX(K)=X3    - X2
   20   CONTINUE
C
C    TRI is twice the area of the triangle element
C
        TRI=CX(3)*BY(2)-CX(2)*BY(3)
C
C     Now look for the points in the mesh that are inside the tri.
C
        DO J=2,NY2
          YP=FLOAT(J-2)*HS
          IF((YP.LE.YMAX+EPS).AND.(YP.GE.YMIN-EPS))THEN
            DO I=2,NX2
              XP=FLOAT(I-2)*HS
              IF((XP.LE.XMAX+EPS).AND.(XP.GE.XMIN-EPS))THEN
C
C    Calculate the natural coordinates.
C
                DO K1=1,3
                  K2=MOD(K1,3)+1
                  LK2=NOR(LEM(K2,N))
                  CNL=(CX(K1)*(YP-EY(LK2)) + BY(K1)*(XP-EX(LK2)))/TRI
                  IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 50
                ENDDO    ! on K1
C
C    If local coords in range the point is within the triangle
C
                IJ=(J-1)*NX3 + I
                NW=(N-1)*2 + ICHO
                S(IJ)=WORK(NW)
              ENDIF    ! within rough X bounds
   50       CONTINUE   ! go to next I, if any local coordinate outside
            ENDDO      ! on I
          ENDIF        ! within rough Y bounds
        ENDDO          ! on J
   60 CONTINUE         ! on N
      RETURN
      END

