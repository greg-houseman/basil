C*--------------------------------------------------------------------
C*    Basil / Sybil:   cform.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE CFORM(IDEFTYP,BANGL,ERA,XLEN,YLEN,DEFV,
     :                 NCOMP,XZERO,YZERO,EX,EY,LEM,NOR,
     :                 NUP,NE,NN,LUW,LSC)
C
C    Routine to modify the initial shape of the mesh, after the
C    boundary conditions have been applied
C
C      IDEFTYP=0 : no change to the rectangular mesh
C      IDEFTYP=1 : circular arc added to southern boundary
C      IDEFTYP=2 : arc added and southern boundary rotated
C      IDEFTYP=3 : southern boundary rotated, no arc
C      IDEFTYP=4 : plate bending algorithm implemented
C                to provide initial geometry of subducted slab
C      IDEFTYP=5 : perturbs y = 0 boundary
C      IDEFTYP=6 : is used to modify the layer thickness (see SSINIT)
C      IDEFTYP=7 : stretch mesh to fit weak zone
C      IDEFTYP=8 : similar to TYPE=5, but allows a different
C                  stretch factor in upper part of layer (above Moho)
C      IDEFTYP=9 : stretch an initially rectangular region by
C                  arbitrarily moving one corner to a new location
C      IDEFTYP=10: cause mesh concentration in a circular region around
C                  one point on the mesh
C      IDEFTYP=11: cause rotation of a circular patch of the mesh
C      IDEFTYP=12: perturbs Moho and lithosphere base
C      IDEFTYP=16: used to modify the layer thickness using a special
C                  function of position(see SSINIT2)
C      IDEFTYP=100: deform rectangular mesh into a parallelogram with an
C                  angle of TILTDEG in the upper left corner(see TILT)
C      IDEFTYP=101: cause mesh concentration around FLTTIPS (see DFORM)
C      IDEFTYP=110: project to spherical coordinates (see PROJECTXY)
C      IDEFTYP=111: rotate to equator in spherical coordinates (see GETPOLE)
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION DEFV(5)
C
      WRITE(LSC,10011)IDEFTYP
      WRITE(LUW,10011)IDEFTYP
10011 FORMAT('Mesh deformation of type ',I4,' is now implemented')
      IF(IDEFTYP.LE.0)RETURN
      BANGLR=BANGL*3.141592653/180.0
C
C    The boundary X = X0 is deformed into the arc of a circle, with
C    radius of curvature ARA, and centre at X = X0 + ARA - ERA,
C    Y = 1. The maximum deformation (ERA) is at (X0,1).
C     Y coordinates are not altered in this section
C
      IF(IDEFTYP.LE.2)THEN
      X0=0.0
      XB=1.0
      WRI=0.25
C  quick and dirty to allow concave curvature
      ARA=(WRI*WRI+ERA*ERA)/(2.0*ABS(ERA))
      DO 60 I=1,NN
        Y=ABS(1.0-EY(I))
        IF(Y.GE.WRI)GO TO 60
        X=EX(I)
        IF(X.GE.XB)GO TO 60
        F=ABS(ERA) - ARA + SQRT(ARA*ARA - Y*Y)
        IF (ERA.LT.0.0) F=-F
        G=(XB-X)/(XB-X0)
        EX(I)=X - F*G
   60 CONTINUE
      WRITE(LUW,10555)ERA,ARA
      WRITE(LSC,10555)ERA,ARA
10555 FORMAT(' Boundary on X=X0 has been deformed into the arc',
     1' of a circle, with ERA, ARA =',2G12.5,/)
      IF(XB.NE.1.0)WRITE(LUW,10556)
      IF(XB.NE.1.0)WRITE(LSC,10556)
10556 FORMAT(' Mesh to the North of the Altyn Tagh remains',
     1' undistorted',/)
      END IF
C
C    Now rotate the southern boundary of the mesh so that its
C    new orientation is BANGL degrees clockwise of its former
C     orientation
C
      IF((IDEFTYP.GE.2).AND.(IDEFTYP.LE.3))THEN
      DO 80 I=1,NN
      X=EX(I)
      Y=EY(I)
      RRY=Y-1.0
      RRX=X-X0
C
C   For points with X <= X0, 0.5 <= Y <= 1.5
C    rigid body rotation around X=X0, Y=1.
C
      IF((RRX.LE.1.E-4).AND.(ABS(RRY).LE.0.5+1.E-4))THEN
      RSS=SQRT(RRX*RRX + RRY*RRY)
      IF (RRX.EQ.0.AND.RRY.EQ.0) THEN
        TH=0
      ELSE
        TH=ATAN2(RRX,RRY)
      END IF
      ANG=BANGLR + TH
      XN=X0 + RSS*SIN(ANG)
      YN=1.0 + RSS*COS(ANG)
      END IF
C
C   For points with X >= X0 and 0.5 <= Y <= 1.5
C     stretch down towards rotated southern margin
C
      IF((RRX.GT.1.E-4).AND.(ABS(RRY).LE.0.5+1.E-4))THEN
      XMOD=(XLEN-X)/(XLEN-X0)
      DX= RRY*SIN(BANGLR)
      DY=-RRY*(1.0-COS(BANGLR))
      XN=X+DX*XMOD
      YN=Y+DY*XMOD
      END IF
C
C   For points with Y <= 0.5 or Y >= 1.5
C
      IF(ABS(RRY).GT.0.5+1.E-4)THEN
      XMOD=(XLEN-X)/(XLEN-X0)
      YMOD=(1.0-ABS(RRY))*2.0
      DX1= RRY*SIN(BANGLR)
      DX2= RRY*TAN(BANGLR)
      DY1=-RRY*(1.0-COS(BANGLR))
      DY2=0.0
      DX=YMOD*DX1 + (1.0-YMOD)*DX2
      DY=YMOD*DY1 + (1.0-YMOD)*DY2
      XN=X+DX*XMOD
      YN=Y+DY*XMOD
      END IF
C
      EX(I)=XN
      EY(I)=YN
   80 CONTINUE
      WRITE(LUW,10557)BANGL
      WRITE(LSC,10557)BANGL
10557 FORMAT(' Southern boundary has been rotated',F6.1,' degrees',/)
      END IF
C
C    For X  < X1 the plate is not moved
C    For X1 < X < X2 the plate is distorted into a circular arc
C      with radius of curvature = R1
C    For X  > X2 the plate is rotated and translated down to a constant dip
C
      IF((IDEFTYP.EQ.4).AND.(BANGL.NE.0.0))THEN
      X1=DEFV(1)
      X2=DEFV(2)
      RC=(X2-X1)/BANGLR
      YC=0.5*YLEN-RC
      XC=X1
C
C   Center of curvature for curved segment at (XC,YC)
C
      DO 160 I=1,NN
      Y=EY(I)
      X=EX(I)
      IF(X.GE.X2)THEN
        RR=SQRT((X-X2)*(X-X2)+(Y-YC)*(Y-YC))
        IF ((X-X2).EQ.0.AND.(Y-YC).EQ.0) THEN
          TH=0
        ELSE
          TH=ATAN2((Y-YC),(X-X2))
        END IF
        TH=TH-BANGLR
        EX(I)=XC+RR*COS(TH)
        EY(I)=YC+RR*SIN(TH)
      ELSE IF(X.LE.X1)THEN
        GO TO 160
      ELSE
        RR=RC+Y-0.5*YLEN
        TH=(X-X1)/RC
        EX(I)=XC+RR*SIN(TH)
        EY(I)=YC+RR*COS(TH)
      END IF
  160 CONTINUE
      WRITE(LUW,10558)X1,X2,BANGL
      WRITE(LSC,10558)X1,X2,BANGL
10558 FORMAT('Mesh has been distorted into circular arc between X=',
     1F8.4,' and X =',F8.4,/,'Rotation angle = ',F8.4,' degrees',/)
      END IF
C
C    stretch the mesh below Y = 0
C
      IF(IDEFTYP.EQ.5)THEN
C       IF(DEFV(3).NE.0.0)THEN
C         DO 238 I=1,NN
C           EY(I)=EY(I)-DEFV(3)
C 238     CONTINUE
          YMIN=0.0
          DO 239 I=1,NN
            Y=EY(I)
C           IF(Y.LT.0.0)EY(I)=-DEFV(4)*Y*Y
            IF(Y.LT.0.0)THEN
              DY =-DEFV(4)*Y
              EY(I)=EY(I)-DY
              IF(EY(I).LT.YMIN)YMIN=EY(I)
            END IF
  239     CONTINUE
          WRITE(LUW,10559)1.0+DEFV(4)
          WRITE(LSC,10559)1.0+DEFV(4)
10559     FORMAT('Stretch factor ',F6.3,' applied below Y = 0')
C         WRITE(LUW,10559)DEFV(3),DEFV(4)
C         WRITE(LSC,10559)DEFV(3),DEFV(4)
C10559     FORMAT('Mesh has been translated down ',F6.3,' unit,',
C     1   ' and stretch factor ',F6.3,' applied, below Y=0')
C       END IF
C
C    perturb layer centered on y = DEFV(3) with displacement
C    in y-direction. Displacement drops to zero at y = DEFV(4)
C    and y = DEFV(5) above and below.  Amplitude is DEFV(1), relative
C    wavenumber is DEFV(2)
C
        PI=3.141592653
        IF(YMIN.EQ.0)YMIN=1.0
        IF(DEFV(1).NE.0.0)THEN
        IF(DEFV(4).EQ.DEFV(3))DEFV(4)=DEFV(3)-0.1
        IF(DEFV(5).EQ.DEFV(3))DEFV(5)=DEFV(3)+0.1
          DO 240 I=1,NN
            X=EX(I)
            Y=EY(I)
C           DY=DEFV(1)*(1.0+COS(PI*DEFV(2)*(X-DEFV(5))))*0.5
C           IF(X.LT.DEFV(5))DY=DEFV(1)
C           IF(X.GT.(DEFV(5)+1.0/DEFV(2)))DY=0.0
            DY=DEFV(1)*COS(PI*DEFV(2)*X/XLEN)
            SCF=0.0
            IF((Y.LE.DEFV(3)).AND.(Y.GT.DEFV(4)))
     .             SCF = ABS((Y-DEFV(4))/(DEFV(3)-DEFV(4)))
            IF((Y.GE.DEFV(3)).AND.(Y.LT.DEFV(5)))
     .             SCF = ABS((DEFV(5)-Y)/(DEFV(5)-DEFV(3)))
C           IF(Y.LT.DEFV(3))SCF=ABS((DEFV(3)-ABS(Y))/DEFV(3))
C           IF(Y.LT.0.0)SCF = 1.0-ABS(Y/YMIN)
C           SCF=(YLEN-DEFV(3)-ABS(Y))
C           IF(ABS(Y).GT.(YLEN-DEFV(3)))SCF=0.0
            EY(I)=Y-DY*SCF
  240     CONTINUE
          WRITE(LUW,10560)DEFV(1),DEFV(2)
          WRITE(LSC,10560)DEFV(1),DEFV(2)
10560     FORMAT('Boundary (Y=0) has been perturbed: amplitude =',
     1    F6.4,'  relative wavenumber =',F6.3)
        END IF
      END IF
C
C    perturb y = 0 boundary and Moho
C
      IF(IDEFTYP.EQ.8)THEN
        PI=3.141592653
        IF(DEFV(1).NE.0.0)THEN
          YMOHO=1.0-DEFV(4)
          IF(NCOMP.LE.1)THEN
            DO 250 I=1,NN
              X=EX(I)
              Y=EY(I)
              DY=DEFV(1)*COS(PI*DEFV(2)*X/XLEN)
              SCF=0.0
              STRETCH=1.0
                IF((Y.LT.(YLEN+YZERO)).AND.(Y.GT.YMOHO))THEN !apr20,2004
                 SCF=((YLEN+YZERO-Y)/(YLEN+YZERO-YMOHO))*DEFV(3) !apr20,2004
               ELSE IF((Y.LE.YMOHO).AND.(Y.GE.0.0))THEN  !apr20,2004
                SCF=1.0+(Y/YMOHO)*(DEFV(3)-1.0)
               ELSE IF(Y.LT.0.0)THEN  !apr20,2004
                SCF=1.0-(Y/YZERO)
                STRETCH=1.0-DEFV(5)*Y
              END IF
              EY(I)=Y*STRETCH-DY*SCF
  250       CONTINUE   !I=1,NN
            WRITE(LUW,10561)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
            WRITE(LSC,10561)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
10561       FORMAT('Boundary (Y=0) has been perturbed: amplitude =',
     1      F6.4,'  relative wavenumber =',F6.3,/,
     2      '   Moho deflection =',F7.4,' at y = ',F7.4)
            WRITE(LUW,10566)DEFV(5)
            WRITE(LSC,10566)DEFV(5)
10566       FORMAT('Box below Y = 0, stretched, DEFV(5) =',F6.3)
C
C   Bessel function boundary perturbation for cylindrical geometry
C  (NCOMP=2). J1 zeros occure at x = 3.8317, 7.0156, 10.173
C   J1 zero is selected by putting value in DEFV(2)
C
          ELSE IF(NCOMP.EQ.2)THEN  !cylindrical geometry
C           YMOHO=DEFV(3)
            DO 251 I=1,NN
              X=EX(I)
              Y=EY(I)
              DY=DEFV(1)*BESSJ0(DEFV(2)*X/XLEN)
              SCF=0.0
              STRETCH=1.0
              IF((Y.LT.(YLEN+YZERO)).AND.(Y.GT.YMOHO))THEN  !apr20,2004
                SCF=((YLEN+YZERO-Y)/(YLEN+YZERO-YMOHO))*DEFV(3) !apr20,2004
C               SCF=((YLEN+YZERO-Y)/(YLEN+YZERO-YMOHO))
              ELSE IF((Y.LE.YMOHO).AND.(Y.GE.0.0))THEN  !apr20,2004
                SCF=1.0+(Y/YMOHO)*(DEFV(3)-1.0)
C               SCF=(Y/YMOHO)
              ELSE IF(Y.LT.0.0)THEN  !apr20,2004
                SCF=1.0-(Y/YZERO)
                STRETCH=1.0-DEFV(5)*Y
              END IF
              EY(I)=Y*STRETCH-DY*SCF
C             EY(I)=Y*STRETCH+DY*SCF
  251       CONTINUE   !I=1,NN
            WRITE(LUW,10565)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
            WRITE(LSC,10565)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
10565       FORMAT('Boundary has been perturbed (BESSJ0): amplitude =',
     1      F6.4,'  relative wavenumber =',F6.3,/,
     2      '   Moho deflection =',F7.4,' at y = ',F7.4)
            WRITE(LUW,10566)DEFV(5)
            WRITE(LSC,10566)DEFV(5)
          END IF  !NCOMP.EQ.2
        END IF  !DEFV(1).NE.0.0
      END IF  !IDEFTYP.EQ.8
C
C    stretch mesh on one side of the central Aust. weak zone
C
      IF((IDEFTYP.EQ.7) .AND. (DEFV(1).GT.0.0)) THEN
        DO 320 I=1,NN
          X=EX(I)
          Y=EY(I)
          IF (X.LE.DEFV(1))THEN 
            IF ((Y.LT.DEFV(2)) .AND. (DEFV(2).GT.0.0)) THEN
              DY=DEFV(4)*(Y/DEFV(2))
            ELSE IF ((Y.GT.DEFV(3)) .AND. (DEFV(3).LT.YLEN))THEN
              DY=DEFV(4)*(YLEN-Y)/(YLEN-DEFV(3))
            ELSE 
              DY=DEFV(4)
            END IF
            DYX=DY*(DEFV(1)-X)/(DEFV(1))
C
            EY(I)=Y+DYX
          END IF
  320   CONTINUE
        WRITE(LUW,10569)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
        WRITE(LSC,10569)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
10569   FORMAT('Mesh has been stretched for X < ',F6.4,'and',F6.4,
     1  ' < Y <',F6.4,'DY = ',F6.4)
      END IF  !IDEFTYP.EQ.7
C
C    Stretch an initially rectangular region by moving the point
C    at (DEFV(1), DEFV(2)) to (DEFV(3), DEFV(4)), ensuring that
C    the other corner points are not moved.
C
      IF(IDEFTYP.EQ.9)THEN
        XC1=DEFV(1)
        YC1=DEFV(2)
        IF(ABS(XC1).LT.1.E-4)THEN
          XOPP=XLEN
          XPTM=0.0
        ELSE IF(ABS(XC1-XLEN).LT.1.E-4)THEN
          XOPP=0.0
          XPTM=XLEN
        ELSE
          WRITE(LSC,10572)XC1,YC1
          WRITE(LUW,10572)XC1,YC1
10572     FORMAT('Coordinates of point ',2F6.3,
     1           ' differ from corner points')
          STOP
        END IF
        IF(ABS(YC1).LT.1.E-4)THEN
          YOPP=YLEN
          YPTM=0.0
        ELSE IF(ABS(YC1-YLEN).LT.1.E-4)THEN
          YOPP=0.0
          YPTM=YLEN
        ELSE
          WRITE(LSC,10572)XC1,YC1
          WRITE(LUW,10572)XC1,YC1
          STOP
        END IF
        XDISP=DEFV(3)-XC1
        YDISP=DEFV(4)-YC1
C
C   each point in the mesh is moved a distance proportional to XDISP
C   in the X-direction and YDISP in the Y-direction
C
      DO 70 I=1,NN
        X=EX(I)
        Y=EY(I)
        RRY=(Y-YOPP)/(YPTM-YOPP)
        RRX=(X-XOPP)/(XPTM-XOPP)
C
C   the stretch factor is scaled back in proportion to distance
C   from the opposite boundaries.  This is a bilinear interpolation.
C
        DX=XDISP*RRY*RRX
        DY=YDISP*RRY*RRX
        EX(I)=X+DX
        EY(I)=Y+DY
   70 CONTINUE
C
      WRITE(LUW,10562)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
      WRITE(LUW,10562)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
10562 FORMAT(' Node at (',F6.3,',',F6.3,') moved to (',
     1 F6.3,',',F6.3,')',/)
      END IF  !IDEFTYP.EQ.9
C
C    mesh localisation option v.1 : recommend use v.2 (IDEFTYP=20)
C
      IF(IDEFTYP.EQ.10)THEN
        XC1=DEFV(1)
        YC1=DEFV(2)
        RAD=DEFV(3)
        SMIN=DEFV(4)
C    
C   each vertex point within radius RAD of (XC1, YC1) is moved closer
C   to (XC1,YC1) by a scaling factor that decreases from 1
C   at radius RAD to SMIN at radius zero using a linear (radius)
C   dependence.
C
        DO 71 I=1,NN
          X=EX(I)
          Y=EY(I)
          RR=SQRT((Y-YC1)*(Y-YC1)+(X-XC1)*(X-XC1))
          IF(RR.LT.RAD)THEN
            SFAC=(RR+SMIN*(RAD-RR))/RAD
            XN=XC1+(X-XC1)*SFAC
            YN=YC1+(Y-YC1)*SFAC
            EX(I)=XN
            EY(I)=YN
          END IF
   71   CONTINUE
C
        WRITE(LSC,10563)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
        WRITE(LUW,10563)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
10563   FORMAT(' Mesh has been distorted to increase node ',
     1  'concentration around (',F6.3,',',F6.3,')',/,' within radius',
     2  F6.3,' node separations scaled by factor of up to:',F6.3,/)
      END IF  !IDEFTYP.EQ.10
C
C    mesh localisation option v.2: focus mesh near r=RADA at the
C    expense of mesh near radius r=RADB, relative to (XC1,YC1)
C    ALPHA between 0.5 and 1 is recommended.
C
      IF(IDEFTYP.EQ.20)THEN
        RADA=DEFV(1)
        RADB=DEFV(2)
        IF(ABS(RADB-RADA).LT.1.E-10)THEN
          WRITE(LSC,*)'TYPE 20 deformation requires DEFV(1) and ',
     :                'DEFV(2) are not the same; no action taken'
        ELSE
        ALPHA1=1.0+DEFV(3)
        XC1=DEFV(4)
        YC1=DEFV(5)
        BTERM=RADB**ALPHA1
        ATERM=RADA**ALPHA1
        DNORM=1.0/(BTERM-ATERM)
C    
C   each vertex point within RADA < r < RAD (XC1, YC1) is moved 
C   relative to (XC1,YC1) by a scaling factor that decreases mesh
C   discretization closer to r = RADB.  Midpoints are moved at
C   end of this routine.
C
        DO 73 I=1,NN
          X=EX(I)
          Y=EY(I)
          RR=SQRT((Y-YC1)*(Y-YC1)+(X-XC1)*(X-XC1))
          IF((RR.GT.RADA).AND.(RR.LT.RADB))THEN
            IF(RR.GT.0)THEN
              SINT=Y/RR
              COST=X/RR
              RTERM=RR**ALPHA1
              RPR=(RADA*(BTERM-RTERM)+RADB*(RTERM-ATERM))*DNORM
              XN=XC1+RPR*COST
              YN=YC1+RPR*SINT
              EX(I)=XN
              EY(I)=YN
            ENDIF
          END IF
   73   CONTINUE
C
        WRITE(LSC,10567)DEFV(1),DEFV(2),DEFV(4),DEFV(5),DEFV(3)
        WRITE(LUW,10567)DEFV(1),DEFV(2),DEFV(4),DEFV(5),DEFV(3)
10567   FORMAT('Mesh resolution has been increased near r = ',
     :  F6.3,' at the expense of mesh',/,'near r = ',F6.3,
     :  ' relative to centre (',F6.3,',',F6.3,') using alpha = ',
     :  F6.3,/)
        ENDIF
      END IF  !IDEFTYP.EQ.20
C
C      set up rotation parameters for mesh rotation
C
      IF(IDEFTYP.EQ.11)THEN
        PIB180=3.141592653/180.0
        RANGL=PIB180*DEFV(1)
        XC1=DEFV(2)
        YC1=DEFV(3)
        RAD1=DEFV(4)
        RAD2=DEFV(5)
C
C   each vertex point within radius RAD1 is rotated about the  point
C  (XC1, YC1)C   by angle RANGL.  The rotation is tapered linearly to zero
C  out at radius RAD2.
C
        DO 72 I=1,NN
          X=EX(I)
          Y=EY(I)
          RR=SQRT((Y-YC1)*(Y-YC1)+(X-XC1)*(X-XC1))
          IF(RR.LT.RAD1)THEN
            RANGLE=RANGL
          ELSE IF((RR.GE.RAD1).AND.(RR.LT.RAD2))THEN
            RANGLE=RANGL*(RR-RAD2)/(RAD1-RAD2)
          ELSE
            RANGLE=0.0
          END IF
            COSQ=COS(RANGLE)
            SINQ=SIN(RANGLE)
            XN=XC1+(X-XC1)*COSQ-(Y-YC1)*SINQ
            YN=YC1+(X-XC1)*SINQ+(Y-YC1)*COSQ
            EX(I)=XN
            EY(I)=YN
   72   CONTINUE
C
        WRITE(LSC,10564)DEFV(1),DEFV(2),DEFV(3),DEFV(4),DEFV(5)
        WRITE(LUW,10564)DEFV(1),DEFV(2),DEFV(3),DEFV(4),DEFV(5)
10564   FORMAT(' Mesh has been distorted to rotate by angle ',F6.3,
     1  ' degrees, the region around (',F6.3,',',F6.3,')',/,
     2  ' within radius',F6.3,
     2  ' and rotation tapered linearly to zero at radius ',F6.3,/)
      END IF   !IDEFTYP=11
      IF(IDEFTYP.EQ.12)THEN
         write(LSC,*)'begin IDEFTYP=12'
         write(LSC,*)'DEFV:'
         write(LSC,*)DEFV(1),DEFV(2),DEFV(3),DEFV(4),DEFV(5)
         IF(DEFV(1).NE.1.0)THEN
            YMOHO=1.0-DEFV(4)
            write(LSC,*)'YMOHO.YLEN',YMOHO,YLEN
            DO I=1,NN   !loop over all nodes
               X=EX(I)  !x-coordinate
               Y=EY(I)  !y-coordinate (height over lithosphere base)
C              width of transition zone
               TRANSWIDTH=DEFV(2)*0.1
               IF(X.LE.(DEFV(2)-TRANSWIDTH))THEN  !thickened region
                  DY=DEFV(1) !downward displacement of base lithosphere
               ELSE  !half Gaussian function
                  X1=DEFV(2)-TRANSWIDTH
                  DY=DEFV(1)*exp(-(X-X1)*(X-X1)/(TRANSWIDTH*TRANSWIDTH))
               ENDIF

               SCF=0.0
               STRETCH=1.0
               IF((Y.LT.YLEN+YZERO).AND.(Y.GT.YMOHO))THEN  !crust
                  SCF=(YLEN+YZERO-Y)/(YLEN+YZERO-YMOHO)*DEFV(3)
               ELSEIF((Y.LE.YMOHO).AND.(Y.GE.0.0))THEN !mantle lithosphere
                  SCF=1.0+Y/YMOHO*(DEFV(3)-1.0)
               ELSEIF(Y.LT.0.0)THEN
                  SCF=1.0-(Y/YZERO)
                  STRETCH=1.0-DEFV(5)*Y
               ENDIF
               EY(I)=Y*STRETCH-DY*SCF
            ENDDO  !I
         ENDIF
         WRITE(LUW,10570)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
         WRITE(LSC,10570)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
10570    FORMAT('Boundary has been perturbed: amplitude =',
     1   F6.4,'  width =',F6.3,/,
     2   '   Lithosphere deflection =',F7.4,' at y = ',F7.4)
         WRITE(LUW,10566)DEFV(5)
         WRITE(LSC,10566)DEFV(5)
      ENDIF  !IDEFTYP=12

      IF(IDEFTYP.EQ.13)THEN
         write(LSC,*)'begin IDEFTYP=13'
         write(LSC,*)'DEFV:'
         write(LSC,*)DEFV(1),DEFV(2),DEFV(3),DEFV(4),DEFV(5)
         IF(DEFV(1).NE.1.0)THEN
            YMOHO=1.0-DEFV(4)
            write(LSC,*)'YMOHO.YLEN',YMOHO,YLEN
            DO I=1,NN   !loop over all nodes
               X=EX(I)  !x-coordinate
               Y=EY(I)  !y-coordinate (height over lithosphere base)
C              width of transition zone
               TRANSWIDTH=DEFV(2)*0.1
               IF(X.LE.(DEFV(2)-TRANSWIDTH))THEN  !thickened region
                  DYL=DEFV(1)*DEFV(3) !downward displacement of base lithosphere
                  DYM=DEFV(1)
               ELSE  !half Gaussian function
                  X1=DEFV(2)-TRANSWIDTH
                  DYL=DEFV(1)*DEFV(3)*exp(-(X-X1)*(X-X1)/
     &                       (TRANSWIDTH*TRANSWIDTH))
                  DYM=DEFV(1)*exp(-(X-X1)*(X-X1)/
     &                       (TRANSWIDTH*TRANSWIDTH))
               ENDIF

               SCF=0.0
               STRETCH=1.0
               IF((Y.LT.YLEN+YZERO).AND.(Y.GT.YMOHO))THEN  !crust
                  SCF=(YLEN+YZERO-Y)/(YLEN+YZERO-YMOHO)
               ELSEIF((Y.LE.YMOHO).AND.(Y.GE.0.0))THEN !mantle lithosphere
                  SCF=1.0+(1.0-Y/YMOHO)*(DEFV(3)-1.0)
               ELSEIF(Y.LT.0.0)THEN
                  SCF=(1.0-(Y/YZERO))*DEFV(3)
                  STRETCH=1.0-DEFV(5)*Y
               ENDIF
               EY(I)=Y*STRETCH-DYM*SCF
            ENDDO  !I
         ENDIF
         WRITE(LUW,10571)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
         WRITE(LSC,10571)DEFV(1),DEFV(2),DEFV(3),DEFV(4)
10571    FORMAT('Boundary has been perturbed: amplitude =',
     1   F6.4,'  width =',F6.3,/,
     2   '   Moho deflection =',F7.4,' at y = ',F7.4)
         WRITE(LUW,10566)DEFV(5)
         WRITE(LSC,10566)DEFV(5)

      ENDIF  !IDEFTYP=13

C     just for testing oregano buoyancy routines. Remove later
C     perturb interface at a given level. Harmonic perturbation
      IF(IDEFTYP.EQ.14)THEN
         PI=3.141592653
         DO I=1,NN   !loop over all nodes
            X=EX(I)  !x-coordinate
            Y=EY(I)  !y-coordinate (height over lithosphere base)
            IF (Y.GT.(DEFV(2)-0.001).AND.Y.LT.(DEFV(2)+0.001))THEN
               EY(I)=DEFV(2)+DEFV(1)*COS(X*PI/DEFV(3))
            ENDIF
         ENDDO
      ENDIF  !IDEFTYP=14
C     stretch the entire mesh in the X or Y direction by a factor
C     DEFV(1) - X direction and DEFV(2) - Y direction
C
      IF(IDEFTYP.EQ.15)THEN
         IF((DEFV(1).LT.0.001).OR.(DEFV(2).LT.0.001))THEN
           WRITE(LUW,*)'DEFORM command DEFV(1) or DEFV(2) is too small'
           WRITE(LSC,*)'DEFORM command DEFV(1) or DEFV(2) is too small'
           STOP
         END IF
         DO I=1,NN   !loop over all nodes
           X=EX(I)  !x-coordinate
           EX(I)=X*DEFV(1)
           Y=EY(I)  !y-coordinate
           EY(I)=Y*DEFV(2)
         ENDDO
         WRITE(LSC,10015)DEFV(1),DEFV(2)
         WRITE(LUW,10015)DEFV(1),DEFV(2)
10015    FORMAT('Mesh has been stretched in X direction by factor',
     :   F9.4,/,'Mesh has been stretched in Y direction by factor',
     :   F9.4)
      ENDIF  !IDEFTYP=15
C
C       Reposition midpoint nodes
C
      DO JEL=1,NE
        NJ1=NOR(LEM(1,JEL))
        NJ2=NOR(LEM(2,JEL))
        NJ3=NOR(LEM(3,JEL))
        NJ4=NOR(IABS(LEM(4,JEL)))
        NJ5=NOR(IABS(LEM(5,JEL)))
        NJ6=NOR(IABS(LEM(6,JEL)))    
        EX(NJ4)=(EX(NJ3)+EX(NJ1))*0.5
        EX(NJ5)=(EX(NJ1)+EX(NJ2))*0.5
        EX(NJ6)=(EX(NJ2)+EX(NJ3))*0.5
        EY(NJ4)=(EY(NJ3)+EY(NJ1))*0.5
        EY(NJ5)=(EY(NJ1)+EY(NJ2))*0.5
        EY(NJ6)=(EY(NJ2)+EY(NJ3))*0.5
      ENDDO
      RETURN
      END
C**********************************************************
      SUBROUTINE SSINIT(LFLAT,HLENSC,AMX1,AMY1,AMH,PHIX,PHIY,
     :                  XLEN,YLEN,EX,EY,SSQ,NOR,NUP,LUW,LSC,IERR)
C
C     subroutine to initialise the SSQ array with a non-uniform
C     distribution of crustal thickness
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION SSQ(NUP)
      PI=3.141592653
C
C    for initial perturbation to crustal thickness
C
      WVNBR=PI
C     WVNBRX=PI/XLEN
C     WVNBRY=PI/YLEN
      WVNBRX=PHIX/XLEN
      WVNBRY=PHIY/YLEN
      PHIX=0.0
      PHIY=0.0
      SQ3=SQRT(3.0)
      DO I=1,NUP   
        XX=EX(NOR(I))
        YY=EY(NOR(I))
C
C    for rolls, squares, rectangles
C
        DIFF=AMX1*COS(WVNBRX*XX+PHIX)+AMY1*COS(WVNBRY*YY+PHIY)
C
C    for hexagons
C
        DIFF=DIFF + AMH*( COS(WVNBR*XX) +
     1      2.0*COS(0.5*WVNBR*XX)*COS(0.5*SQ3*WVNBR*YY) )
        IF(DIFF.LE.-1.0)THEN
          WRITE(LSC,*)'DIFF wrong in SSINIT: modify AMX1,AMY1, AMH'
          IERR = 1
          GO TO 60
        END IF
        SSQ(I)=SSQ(I) + ALOG(REAL(1.0+DIFF))
      ENDDO
      WRITE(LUW,10109)AMX1, AMY1, AMH, PHIX, PHIY
10109 FORMAT('Layer thickness perturbed using the ',
     1'parameters: AMX1, AMY1, AMH, PHIX, PHIY =',5F10.4)
C
   60 RETURN
      END
C***********************************************************
      SUBROUTINE SSINIT2(LFLAT,HLENSC,AMX1,AMY1,AMH,PHIX,PHIY,
     :                  XLEN,YLEN,EX,EY,SSQ,NOR,NUP,LUW,LSC,IERR)
C
C     subroutine to initialise the SSQ array with a non-uniform
C     distribution of crustal thickness, using a special funtion of position
C     in the interior of an ellipse of axis a=AMX1; b=AMY1;
C     f=AMH=declared thickening ratio;
C     C0=PHIX,C1=PHIY- inner/outer parameters i nside the ellipse
C     THG-thickening ratio

      REAL EP,THG,X0,Y0,C0,C1,L0
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION SSQ(NUP)
      PI=3.141592653
      X0=0.85   !The centre of the ellipse: in the middle of the SW boundary: (0 .85,0.0)
      Y0=-0.1
C      L0=0.35
C
C      PHIX=0.0
C      PHIY=0.0
C
      DO I=1,NUP
        XX=EX(NOR(I))
        YY=EY(NOR(I))

        EP=(XX-X0)**2/AMX1**2+(YY-Y0)**2/AMY1**2
C     Define thickening ratio- depending on the position relative to the ellipse

        IF(EP.LE.PHIX) THEN
          THG=AMH
        ELSEIF ((PHIX.LE.EP).AND.(EP.LE.PHIY)) THEN
          THG=1.0+(AMH-1.0)*((EP-PHIY)/(PHIX-PHIY))
        ELSE
          THG=1.0
        ENDIF
C       WRITE(LUW,*)I,XX,YY,THG
C
C    for rolls, squares, rectangles
C       DIFF=AMX1*COS(WVNBRX*XX+PHIX)+AMY1*COS(WVNBRY*YY+PHIY)
C
C    Define crustal thickness: SSQ=ln(THG)+ln(L0); L0-reference layer thickness
C       SSQ(I)=SSQ(I) + ALOG(REAL(1.0+DIFF))
C
        SSQ(I)=SSQ(I) + ALOG(THG)
      ENDDO

      WRITE(LUW,10109)AMX1, AMY1, AMH, PHIX, PHIY
10109 FORMAT('Layer thickness perturbed using the ',
     1'parameters: AMX1, AMY1, AMH, PHIX, PHIY =',5F10.4)
C
   60 RETURN
      END
C**********************************************************************
C
C    Subroutine to impose bessel function perturbation for (cylindrical
C    geometry ie.NCOMP=2
C
      FUNCTION BESSJ0(x)
      REAL BESSJ0,x
      REAL ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     *.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     *651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,
     *s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     *59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        BESSJ0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2 
        xx=ax-.785398164
        BESSJ0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.

