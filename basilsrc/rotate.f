C*--------------------------------------------------------------------
C*    Basil / Sybil:   rotate.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE DFORM(NDFORM,IFCASE,X1,X2,Y1,Y2,IFLTTIPS,
     :                 EX,EY,NOR,LEM,NUP,NE,NN)
C
C     This subroutine deforms the mesh to have a higher density
C     of elements around the crack tips.  Note, the outer boundary
C     is left untouched.  This routine deforms the the mesh radially
C     around each tip for radius < RAD
C
C     IMESH=1: one tip at (1,1)
C     IMESH=2: two tips at (1,0.5) and (1,1.5) if INFLT=2
C              two tips at (X1,(Y1+Y2)/2),(X2,(Y1+Y2)/2) if INFLT=1
C
C     If NDFORM=1 the mesh is left square and 
C     if NDFORM=2 the mesh is deformed into a rectangle (x-cord stretched by 2)
C     if NDFORM=3 the mesh is deformed into a rectangle (y-cord stretched by 2)
C     if NDFORM=4 upper half stretched 3 times in y-direction and left
C                 half stretched 2 times in -x direction (iaspei)
C     if NDFORM=9 the mesh is deformed into a circle
C     if NDFORM=10 the mesh inside a certain radius is deformed into a circle
C
C
      INCLUDE "limits.parameters"
C     DIMENSION PX(2),PY(2)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION IFLTTIPS(MAXTIPS*2)
C
C     print *,(IFLTTIPS(I),I=1,MAXTIPS*2)
      IMESH=1
      RAD1=0.5
      IF(NDFORM.EQ.10) RAD1=0.85
      PY=0.5*(Y2+Y1)
C     IF (INFLT.EQ.1) THEN
C       PX(1)=X1
C       PX(2)=X2
C       PY(2)=PY(1)
C       IMESH=2
C     ELSE
C       PX(1)=0.5*(X2+X1)
C     ENDIF
      RAD=0.5*(X2+X1)-X1
      IF(PY-Y1.LT.RAD) RAD=PY-Y1
      RAD=RAD*RAD1
C     IF(IFCASE.EQ.1) THEN
C       IMESH=1
C       PX(1)=0.5*(X2+X1)
C       PY(1)=0.5*(Y2+Y1)
C       RAD=PX(1)-X1
C       IF(PY(1)-Y1.LT.RAD) RAD=PY(1)-Y1
C     ELSE IF(IFCASE.EQ.2) THEN
C       IMESH=2
C       PX(1)=0.25*(X2-X1)+X1
C       PY(1)=0.5*(Y2+Y1)
C       PX(2)=0.75*(X2-X1)+X1
C       PY(2)=0.5*(Y2+Y1)
C       RAD=0.5*(PX(2)-PX(1))
C       IF(PY(1)-Y1.LT.RAD) RAD=PY(1)-Y1
C     END IF
C     IF(NDFORM.EQ.2) THEN
C           X1=2.0*(X1-PX(JJ))+PX(JJ)
C           X2=2.0*(X2-PX(JJ))+PX(JJ)
C     END IF
C     IF(NDFORM.EQ.3) THEN
C           Y1=2.0*(Y1-PY(JJ))+PY(JJ)
C           Y2=2.0*(Y2-PY(JJ))+PY(JJ)
C     END IF
C     IF(NDFORM.EQ.4) THEN
C           Y2=3.0*(Y2-PY(JJ))+PY(JJ)
C           X1=2.0*(X1-PX(JJ))+PX(JJ)
C     END IF
C
C    Change the spacing of the nodes
C
      DO 20 JJ=1,10
         IF (IFLTTIPS(JJ).EQ.0) GO TO 20
         PX = EX(IFLTTIPS(JJ))
         PY = EY(IFLTTIPS(JJ))
      DO 100 II=1,NN
         XX=EX(II)-PX
         YY=EY(II)-PY
         IF(XX.EQ.0.AND.YY.EQ.0) THEN
           THA=0.0
         ELSE
           THA=ATAN2(YY,XX)
         END IF
         RR=SQRT(XX*XX+YY*YY)
         IF(NDFORM.EQ.9) THEN
                IF (XX.NE.0.0.OR.YY.NE.0.0) THEN
                   IF (THA.LE.PI/4.0.AND.THA.GE.-PI/4.0) THEN
                         RR=ABS(XX)
                   ELSE IF (THA.GT.PI/4.0.AND.THA.LT.3.0*PI/4.0) THEN
                         RR=ABS(YY)
                   ELSE IF (THA.LT.-PI/4.0.AND.THA.GT.-3.0*PI/4.0) THEN
                         RR=ABS(YY)
                   ELSE
                         RR=ABS(XX)
                   END IF
                END IF
         END IF
         IF (RR.LE.RAD) THEN
            RR=RR*RR/RAD
            EX(II)=RR*COS(THA)+PX
            EY(II)=RR*SIN(THA)+PY
         END IF
C        IF(NDFORM.EQ.2) THEN
C           EX(II)=2.0*(EX(II)-PX(JJ))+PX(JJ)
C        END IF
C        IF(NDFORM.EQ.3) THEN
C           EY(II)=2.0*(EY(II)-PY(JJ))+PY(JJ)
C        END IF
C        IF(NDFORM.EQ.4) THEN
C           IF(EY(II).GT.PY(JJ)) THEN
C              EY(II)=3.0*(EY(II)-PY(JJ))+PY(JJ)
C           END IF
C           IF(EX(II).LT.PX(JJ)) THEN
C              EX(II)=2.0*(EX(II)-PX(JJ))+PX(JJ)
C           END IF
C        END IF
C
C For Pressure solution case
C
         IF(NDFORM.EQ.10) THEN
                IF (XX.NE.0.0.OR.YY.NE.0.0) THEN
                   IF (THA.LE.PI/4.0.AND.THA.GE.-PI/4.0) THEN
                         RR=ABS(XX)
                   ELSE IF (THA.GT.PI/4.0.AND.THA.LT.3.0*PI/4.0) THEN
                         RR=ABS(YY)
                   ELSE IF (THA.LT.-PI/4.0.AND.THA.GT.-3.0*PI/4.0) THEN
                        RR=ABS(YY)
                   ELSE
                         RR=ABS(XX)
                   END IF
                END IF
                IF (RR.LT.RAD) THEN
                   EX(II)=RR*COS(THA)+PX
                   EY(II)=RR*SIN(THA)+PY
                END IF
         END IF
  100 CONTINUE
   20 CONTINUE
C
C   Calculate position of midpoint nodes
C
      DO 200 JEL=1,NE
            NJ1=NOR(LEM(1,JEL))
            NJ2=NOR(LEM(2,JEL))
            NJ3=NOR(LEM(3,JEL))
            NJ4=NOR(ABS(LEM(4,JEL)))
            NJ5=NOR(ABS(LEM(5,JEL)))
            NJ6=NOR(ABS(LEM(6,JEL)))
            EX(NJ4)=(EX(NJ3)+EX(NJ1))*0.5
            EX(NJ5)=(EX(NJ1)+EX(NJ2))*0.5
            EX(NJ6)=(EX(NJ2)+EX(NJ3))*0.5
            EY(NJ4)=(EY(NJ3)+EY(NJ1))*0.5
            EY(NJ5)=(EY(NJ1)+EY(NJ2))*0.5
            EY(NJ6)=(EY(NJ2)+EY(NJ3))*0.5
  200 CONTINUE
C
C  Deform Lagrangian grid if necessary
C
C     IF(NDFORM.EQ.2.OR.NDFORM.EQ.3.OR.NDFORM.EQ.4) THEN
C       DO 300 II=1,NULP
C          IF(NDFORM.EQ.2) THEN
C              EXLG(II)=2.0*(EXLG(II)-PX(1))+PX(1)
C          END IF
C          IF(NDFORM.EQ.3) THEN
C              EYLG(II)=2.0*(EYLG(II)-PY(1))+PY(1)
C          END IF
C          IF(NDFORM.EQ.4) THEN
C              IF(EXLG(II).LT.PX(1)) THEN
C                   EXLG(II)=2.0*(EXLG(II)-PX(1))+PX(1)
C              END IF
C              IF(EYLG(II).GT.PY(1)) THEN
C                   EYLG(II)=3.0*(EYLG(II)-PY(1))+PY(1)
C              END IF
C          END IF
C 300   CONTINUE
C     END IF
      RETURN
      END
C
      SUBROUTINE TILT(FLTDEG,EX,EY,LEM,NOR,NUP,NE,NN)
C
C     This subroutine deforms the mesh from a square into a tilt
C     block of the desired shape - a parallelagram with an angle of
C     fltdeg in the upper left corner.
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
C
      PI=3.14159265
C
      ANGLE=(90.0-FLTDEG)*PI/180.0
C
C    Change the X coordinate by -Y*TAN(ANGLE)
C
      DO 100 II=1,NN
         XX=EX(II)
         YY=EY(II)
         EX(II)=XX-YY*TAN(ANGLE) 
  100 CONTINUE
C
C   Calculate position of midpoint nodes
C
      DO 200 JEL=1,NE
            NJ1=NOR(LEM(1,JEL))
            NJ2=NOR(LEM(2,JEL))
            NJ3=NOR(LEM(3,JEL))
            NJ4=NOR(ABS(LEM(4,JEL)))
            NJ5=NOR(ABS(LEM(5,JEL)))
            NJ6=NOR(ABS(LEM(6,JEL)))
            EX(NJ4)=(EX(NJ3)+EX(NJ1))*0.5
            EX(NJ5)=(EX(NJ1)+EX(NJ2))*0.5
            EX(NJ6)=(EX(NJ2)+EX(NJ3))*0.5
            EY(NJ4)=(EY(NJ3)+EY(NJ1))*0.5
            EY(NJ5)=(EY(NJ1)+EY(NJ2))*0.5
            EY(NJ6)=(EY(NJ2)+EY(NJ3))*0.5
  200 CONTINUE
      RETURN
      END
