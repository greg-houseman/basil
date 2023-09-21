C*--------------------------------------------------------------------
C*    Basil / Sybil:   dncom.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
C      contains subroutines DNCOM, RANGXY, RSCALE, SWAPCO, ROTATE
C
      SUBROUTINE LNCORD(PNI,PLI)
C
C    This routine sets up the local coordinates for the
C    integration points in the triangle.  These are the
C    same for every element, and this routine need only
C    be called once at start of program
C    For each I, 7 values of L(i) and N(i) are stored in
C    consecutive memory slots.
C
C     COMMON/AI/LUW,LSC,LBC,LLG
      DIMENSION PNI(42),PLI(21)
      SAVE ALFA
      DIMENSION CL(3)
      DIMENSION ALFA(10)
      DATA ALFA/0.05971587,0.47014206,0.79742699,0.10128651,
     10.33333333,0.0,0.5,1.0,0.0,0.33333333/
C
      INP=0
      INL=0
      I1=1
      I2=2
      I3=3
      DO 30 K=1,7
C
C   set up coordinates for K near vertex
C
      IF(K.LE.3)THEN
      CL(1)=ALFA(4)
      CL(2)=ALFA(4)
      CL(3)=ALFA(4)
      CL(K)=ALFA(3)
      GO TO 20
      END IF
C
C    set up coordinates for K near midpoint
C
      IF(K.LE.6)THEN
      IX=MOD(K,3)+1
      CL(1)=ALFA(2)
      CL(2)=ALFA(2)
      CL(3)=ALFA(2)
      CL(IX)=ALFA(1)
      GO TO 20
      END IF
C
C    K is the centroid
C
      CL(1)=ALFA(5)
      CL(2)=ALFA(5)
      CL(3)=ALFA(5)
C
C    PLI and PNI are arranged in rows of 7, each row consists
C    of interpolation function Li, Ni values at 7 integration point
C    Row 1 is for I=1, through to row 3 for L, row 6 for N.
C
   20 INP=INP+1
      INL=INL+1
      PLI(INL)   =CL(I1)
      PLI(INL+7) =CL(I2)
      PLI(INL+14)=CL(I3)
      PNI(INP)   =CL(I1)*CL(I1) - CL(I1)*(CL(I2)+CL(I3))
      PNI(INP+7) =CL(I2)*CL(I2) - CL(I2)*(CL(I3)+CL(I1))
      PNI(INP+14)=CL(I3)*CL(I3) - CL(I3)*(CL(I1)+CL(I2))
      PNI(INP+21)=4.0*CL(I3)*CL(I1)
      PNI(INP+28)=4.0*CL(I1)*CL(I2)
      PNI(INP+35)=4.0*CL(I2)*CL(I3)
   30 CONTINUE
      RETURN
      END
C
      SUBROUTINE DNCOM(JJ,BY,CX,DNDP)
      DIMENSION BY(3),CX(3),DNDP(84)
      SAVE ALFA
      DIMENSION ALFA(10)
      DATA ALFA/0.05971587,0.47014206,0.79742699,0.10128651,
     10.33333333,0.0,0.5,1.0,0.0,0.33333333/
C
C     This routine calculates the gradient of the interpolation
C     function N(I) at a location in the element determined by
C     the parameter J.  The calculation for X and Y gradient is
C     the same, though different geometrical elements in B are
C     used.  For numbering the points (I & J), vertices are 1-3
C     midpoints are 4-6, and the centroid (J only) is 7.  The
C     interpolation functions are given in Huebner (p345).
C
C     The gradient values DNDP produced by this routine need to be
C     divided by the factor of 2*(triangle area) to produce actual
C     spatial gradients.  This factor is introduced in matrix assembly.
C
C     The interpolation function is for a vertex point
C
      IND=1
      DO 30 I=1,3
        DO J=1,7
C
C     J  is a vertex point
C
          IF(J.LE.3)THEN
            KN2=4
            IF(J.EQ.I)THEN
              KN1=3
              KN3=4
            ELSE
              KN1=4
              KN3=3
            END IF
            GO TO 20
          ENDIF
C
C    J is a midpoint
C
          IF(J.LE.6)THEN
            IX=MOD(I+1,3)+4
            KN3=2
            IF(J.EQ.IX)THEN
              KN1=1
              KN2=2
            ELSE
              KN1=2
              KN2=1
            END IF
            GO TO 20
          ENDIF
C
C     J is the centroid
C
          KN1=5
          KN2=5
          KN3=5
   20     I2=MOD(I,3)+1
          I3=MOD(I+1,3)+1
          CL1=ALFA(KN1+JJ)
          CL2=ALFA(KN2+JJ)
          CL3=ALFA(KN3+JJ)
          DNDP(IND)  =BY(I)*(2.0*CL1-CL2-CL3)-CL1*(BY(I2)+BY(I3))
          DNDP(IND+1)=CX(I)*(2.0*CL1-CL2-CL3)-CL1*(CX(I2)+CX(I3))
          IND=IND+2
        ENDDO
   30 CONTINUE
C
C    The interpolation function is for a midpoint
C
      DO 230 I=4,6
        DO J=1,7
          I1=MOD(I,3)+1
          I2=MOD(I1,3)+1
          I3=MOD(I1+1,3)+1
C
C    J is a vertex point
C
          IF(J.LE.3)THEN
            IF(J.EQ.I2)THEN
              KN2=3
              KN3=4
            ELSE IF(J.EQ.I3)THEN
              KN2=4
              KN3=3
            ELSE
              KN2=4
              KN3=4
            END IF
            GO TO 220
          END IF
C
C    J is a mid-point node
C
          IF(J.LE.6)THEN
            I5 = 4 + MOD(I,3)
            I6 = 4 + MOD(I+1,3)
            IF(J.EQ.I5)THEN
              KN2=1
              KN3=2
            ELSE IF(J.EQ.I6)THEN
              KN2=2
              KN3=1
            ELSE
              KN2=2
              KN3=2
            END IF
            GO TO 220
          END IF
C
C    J is the centroid
C
          KN2=5
          KN3=5
C
  220     CL2=ALFA(KN2+JJ)
          CL3=ALFA(KN3+JJ)
          DNDP(IND)   = 4.0*(CL3*BY(I2) + CL2*BY(I3))
          DNDP(IND+1) = 4.0*(CL3*CX(I2) + CL2*CX(I3))
          IND=IND+2
        ENDDO
  230 CONTINUE
      RETURN
      END
C
      SUBROUTINE RANGXY(X,IHELP,NXY,M1,M2,MP,N1,N2,NP,XMAX,IM,JM,
     1                  XMIN,IN,JN,NX3,SUMINT)
C
C    This routine determines Min and Max of an array to be contoured
C    It does not look at points for which IHELP=0 (outside the mesh)
C
      DIMENSION X(NXY),IHELP(NXY)
      ISET=0
      IVERB=0
      SUMINT=0.0
      DO I = M1,M2,MP
        DO J = N1,N2,NP
          K = (J - 1)*NX3 + I
          IF(IHELP(K).NE.0)THEN
            IF (ISET.EQ.0) THEN
              XMAX = X(K)
              XMIN = X(K)
              IM = I
              JM = J
              ISET = 1
            ENDIF
            IF(X(K).GT.XMAX)THEN
              XMAX = X(K)
              IM = I
              JM = J
            ELSE IF(X(K).LT.XMIN)THEN
              XMIN = X(K)
              IN = I
              JN = J
            END IF
C
C     boundary points use FACTOR for more precise integral
C     (not precise for irregular boundaries)
C
            FACTOR=1.0
             IF((I.EQ.M1).OR.(I.EQ.M2))FACTOR=FACTOR*0.5
             IF((J.EQ.N1).OR.(J.EQ.N2))FACTOR=FACTOR*0.5
            SUMINT=SUMINT+FACTOR*X(K)
          END IF  ! IHELP non-zero
        ENDDO     ! on J
      ENDDO       ! on I
      IF(IVERB.NE.0)WRITE(LUW,10100)XMAX,IM,JM,XMIN,IN,JN
10100 FORMAT(' MAX VALUE OF FUNCTION IS ',G12.5,' AT I =',I6,
     1', J =',I6,'   MIN VALUE IS ',G12.5,' AT I =',I6,', J =',
     2I6,/)
      RETURN
      END

      SUBROUTINE RSCALE(X,Y,NXY,M1,M2,N1,N2,RMAX,RMIN,NX3)
      DIMENSION X(NXY),Y(NXY)
      K=(N1-1)*NX3+M1
      RMAX=SQRT(X(K)*X(K)+Y(K)*Y(K))
      RMIN=RMAX
      DO I = M1,M2
        DO J = N1,N2
          K = (J - 1)*NX3 + I
          R = SQRT(X(K)*X(K) + Y(K)*Y(K))
          RMAX=MAX(R,RMAX)
          RMIN=MIN(R,RMIN)
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE SWAPCO(EX,EY,UVP,NUP,JROT,JVELROT)
C
C  routine to swap X and Y coordinates for purposes of
C   plotting solution
C
      DIMENSION EX(NUP),EY(NUP),UVP(NUP,2)
C
C    Rotate entire mesh and, if JVELROT=1, velocity field anti - clockwise
C    in steps of 90 degrees (JROT steps)
C
      PIB2=3.141592653/2.0
      DO 50 I=1,NUP
        Y=EY(I)
        X=EX(I)
        RD=SQRT(X*X + Y*Y)
        IF (X.EQ.0.0.AND.Y.EQ.0.0) THEN
          TD=0
        ELSE
          TD=ATAN2(Y,X)
        END IF
        TD=TD+PIB2*FLOAT(JROT)
        EX(I)=RD*COS(TD)
        EY(I)=RD*SIN(TD)
        IF (JVELROT.NE.0) THEN
          U=UVP(I,1)
          V=UVP(I,2)
          RD=SQRT(U*U + V*V)
          IF (U.EQ.0.0.AND.V.EQ.0.0) THEN
            TD=0
          ELSE
            TD=ATAN2(V,U)
          END IF
          TD=TD+PIB2*FLOAT(JROT)
          UVP(I,1)=RD*COS(TD)
          UVP(I,2)=RD*SIN(TD)
        END IF
   50 CONTINUE
      RETURN
      END

      SUBROUTINE ROTATE(FROT,NUP,JROT)
C
C    updates the FROT array (rotation in degrees) by 90.0*JROT
C    if the solution is rotated
C
      DIMENSION FROT(NUP)

      ANGL = 90.0*JROT
      DO 50 I=1,NUP
        FROT(I) = FROT(I) + ANGL
   50 CONTINUE
      RETURN
      END

