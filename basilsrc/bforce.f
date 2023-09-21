C*--------------------------------------------------------------------
C*    Basil / Sybil:   bforce.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE BFORCE(IGRAV,RHOG,NCOMP,EX,EY,LEM,NOR,PNI,
     :                  PLI,DENS,QLOAD,NUP,NE,NROWS)
C
C     Routine to calculate the load vector due to body force
C     Force acts in x (IGRAV=+/-1) or y (IGRAV=+/-2) direction
C     Not implemented yet: IGRAV=3, gravity acts in arbitrary direction
C      theta, relative to the x-axis
C     Values put into QLOAD.
C     the DENS array is here assumed to contain the density array
C     and RHOG contains the normalisation constant rho*grav
C
      DOUBLE PRECISION PNI,PLI,BY,CX,DNDP,XX,XRAD,TRIA
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION PNI(42)
      DIMENSION PLI(7,3)
      DIMENSION LEM(6,NE)
      DIMENSION DENS(7,NE)
      DIMENSION QLOAD(NROWS)
      DIMENSION BY(3),CX(3),DNDP(84)
      SAVE W
      DIMENSION W(7)
      DIMENSION XX(3)
      DIMENSION XRAD(7)
C     DATA W/0.05,0.05,0.05,0.13333333,0.13333333,0.13333333,
C    10.45/
      DATA W/0.12593918,0.12593918,0.12593918,0.13239415,
     10.13239415,0.13239415,0.225/
C
      IDIR=IABS(IGRAV)
      ISGN=ISIGN(1,IGRAV)
      IF((IDIR.LT.1).OR.(IDIR.GT.2))RETURN
      IOFF=NUP*(IDIR-1)
      FNORM=ISGN*RHOG
C
C     WRITE(*,*)'check DENS array initialised in BFORCE'
C     CALL MATPRT(DENS,7,NE,7*NE,6)
C
C    Clear the load array
C
      DO I=1,NROWS
        QLOAD(I)=0.0
      ENDDO
C
C    Cycle over the elements
C
      DO 60 N=1,NE
C
C    Set up the geometrical coefficients for the natural
C    coordinates used in the triangle elements.  TRIA is the
C    area of the element
C
      TRIA=0.0
      DO 10 K1=1,3
      K2=MOD(K1,3)+1
      K3=MOD(K1+1,3)+1
      LK1=NOR(LEM(K1,N))
      LK2=NOR(LEM(K2,N))
      LK3=NOR(LEM(K3,N))
      XX(K1)=EX(LK1)
      BY(K1)=EY(LK2)-EY(LK3)
      CX(K1)=EX(LK3)-EX(LK2)
      TRIA=TRIA + (EX(LK2)*EY(LK3)-EX(LK3)*EY(LK2))
   10 CONTINUE
      TRIA=TRIA*0.5
C
      DO 25 K7=1,7
        XRAD(K7)=XX(1)*PLI(K7,1)+XX(2)*PLI(K7,2)+XX(3)*PLI(K7,3)
   25 CONTINUE
C
C    Calculate the values and gradients of the interpolation
C    function at the seven integration points
C
      CALL DNCOM(0,BY,CX,DNDP)
C
C     for node number I
C
      DO I=1,6
        LI=IABS(LEM(I,N))
        SUM1=0.0
C
C     Integrate over the element
C
        DO K=1,7
C
C    perform integration using density at integraton points
C
          KIN=(I-1)*7 + K
          PSUM=W(K)*DENS(K,N)*PNI(KIN)
C
C    For axisymmetric case (around x = 0)
C
          IF(NCOMP.EQ.2)PSUM=PSUM*XRAD(K)
          SUM1=SUM1+PSUM
C
        ENDDO
          LID=IOFF+LI
          QLOAD(LID)=QLOAD(LID)+TRIA*SUM1
      ENDDO
   60 CONTINUE
C
C     normalise using sign and magnitude of gravity
C
      DO 70 I=IOFF+1,IOFF+NUP
        QLOAD(I)=FNORM*QLOAD(I)
   70 CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE DENSIT(DENS,EX,EY,NOR,LEM,NE,NUP,
     :                  RHOBG,RHOP,YLEV,LUW,LSC,IDBUG)
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION DENS(7,NE)
C
C    The DENS array here contains the dimensionless density distribution
C     as used in BFORCE.  the buoyancy terms are scaled by RHOG in BFORCE
C
      DO 100 J=1,NE
C
C    set background level
C
      DO K=1,7
        DENS(K,J)=RHOBG
      ENDDO
C
C    get central coordinate of element and reset according to depth
C
      YC=0.0
      DO K1=1,3
        YC=YC+EY(NOR(LEM(K1,J)))
      ENDDO
      YC=YC/3.0
      YVRT=0.30385953*YC
      YMID=0.17914761*YC
C
C    for vertex node,   (X,Y)=0.69614048(XV,YV)+0.30385953(XC,YC)
C    for midpoint node, (X,Y)=0.82085238(XM,YM)+0.17914761(XC,YC)
C    (obtained from manipulations of data on p421 of Huebner)
C
      DO K=1,3
        NP=ABS(LEM(K,J))
        YDEP=0.69614048*EY(NOR(NP)) + YVRT
        IF(YDEP.LT.YLEV)DENS(K,J)=RHOP
      ENDDO
      DO K=4,6
        NP=ABS(LEM(K,J))
        YDEP=0.82085238*EY(NOR(NP)) + YMID
        IF(YDEP.LT.YLEV)DENS(K,J)=RHOP
      ENDDO
      IF(YC.LT.YLEV)DENS(7,J)=RHOP
  100 CONTINUE
      IF(IDBUG.NE.0)THEN
        WRITE(LSC,*)'If y <',YLEV,',  rho =',RHOP,' else rho =',RHOBG
        CALL MATPRT(DENS,7,NE,7*NE,LUW)
      END IF
      RETURN
      END
      SUBROUTINE MENDBC(ITYPE,BCV,TIMEL,XLEN,YLEN,BIG,EX,EY,NOR,
     :                  IBC,IBCTYP,QBND,NUP,NBP,LSC,LUW)
C
C    this routine adjusts the boundary condition, depending on
C    spatial coordinates of the mesh
C    vertical boundary condition is changed from fixed velocity to
C    fixed stress for those nodes with x > XPASS, Y > YPASS
C
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION IBC(NBP),IBCTYP(NBP*2)
      DIMENSION QBND(NBP*2)
      DIMENSION BCV(5)
C
C    TYPE=4,5 are used in the subducted slab problem
C    TYPE=8 is used to set the velocity bc for the Tranverse Ranges
C    TYPE=14 is used to stop an extending blob at a given depth by
C      changing the boundary condition from TY=0 to UY=0
C
      EPS=0.001
      CONST=3.1415926
      IF(ITYPE.EQ.4)XPASS=BCV(2)
      IF(ITYPE.EQ.5)XPASS=BCV(2)+TIMEL*BCV(4)
      IF(ITYPE.EQ.8)THEN
        XPASS=BCV(2)
        XNOPAS=XLEN-XPASS
        IF(XNOPAS.NE.XPASS)CONST=3.1415926/(XNOPAS-XPASS)
        IF(XNOPAS.LE.XPASS*0.99)THEN
          XNOPAS=XLEN
          IF(XNOPAS.NE.XPASS)CONST=3.1415926*0.5/(XNOPAS-XPASS)
        END IF
      END IF
      YPASS=-2.0
      ISET=0
C
C    test all boundary nodes
C
      DO 10 I=1,NBP
        IXP=I
        IYP=I+NBP
        NDB=IBC(I)
        LDB=NOR(NDB)
        XCM=EX(LDB)
        YCM=EY(LDB)
C
C   Reset the velocity function on y=YLEN for a convergent zone of fixed
C    width.  Assume symmetry about the midpoint if XPASS is to the left of
C    the midpoint, or else the right boundary provides the symmetry plane.
C
        IF(ITYPE.EQ.8)THEN
          IF(ABS(YCM-YLEN).LT.EPS)THEN
            IF((XCM.GE.XPASS).AND.(XCM.LE.XNOPAS))THEN
              QBND(I)=BCV(1)*COS((XCM-XPASS)*CONST)
            ELSE IF(XCM.LE.XPASS)THEN
              QBND(I)=BCV(1)
            ELSE
              QBND(I)=-BCV(1)
            END IF
          END IF
C
C  the following block resets x-component velocity at every timestep
C  as used in the Western US (CVSN) project
C
        ELSE IF(ITYPE.EQ.21)THEN
          IF((YCM.LT.YLEN*0.9).AND.(XCM.GT.BCV(2)).AND.
     1       (XCM.LT.BCV(5)))THEN
              HVEL=0.0
              IF((XCM.LE.BCV(4)).AND.(XCM.GE.BCV(3)))THEN
                HVEL=BCV(1)
              ELSE IF((XCM.LT.BCV(3)).AND.(XCM.GT.BCV(2)))THEN
                FRAC=(XCM-BCV(2))/(BCV(3)-BCV(2))
                SFR=SIN(CONST*0.5*FRAC)
                HVEL=BCV(1)*SFR*SFR
              ELSE IF((XCM.LT.BCV(5)).AND.(XCM.GT.BCV(4)))THEN
                FRAC=(XCM-BCV(5))/(BCV(4)-BCV(5))
                SFR=SIN(CONST*0.5*FRAC)
                HVEL=BCV(1)*SFR*SFR
              END IF
              QBND(I)=HVEL
          END IF
C
C   If a velocity condition has previously been set on a
C    part of the boundary that moves past x=XPASS, modify it
C
        ELSE IF((ITYPE.EQ.4).OR.(ITYPE.EQ.5))THEN
          IF((XCM.GT.XPASS).AND.(YCM.GT.YPASS))THEN
C
C   Change velocity to stress condition
C
            IF(BCV(1).GT.0.0)THEN
              IF(IBCTYP(IYP).EQ.0)THEN
                IBCTYP(IYP)=1
                QBND(IYP)=BCV(3)
                WRITE(LSC,10001)NDB,XCM,YCM,BCV(3)
10001     FORMAT('Node ',I6,', at X =',F8.3,'   Y =',F8.3,'  velocity',
     1    ' -> traction =',F8.3)
              END IF
C
C   Change velocity to friction condition
C
            ELSE
              IF(IBCTYP(IYP).NE.2)THEN
                IBCTYP(IYP)=2
                IBCTYP(I)=2
                QBND(IYP)=BCV(3)
                QBND(I)=BCV(3)
                WRITE(LSC,10002)NDB,XCM,YCM,BCV(3)
10002     FORMAT('Node ',I6,', at X =',F8.3,'   Y =',F8.3,'  velocity',
     1    ' -> friction =',F8.3)
              END IF
            END IF
C
C      eveywhere below YPASS has a frictional condition set, but
C      friction increases linearly from 0 at a depth of Y = YPASS
C      to BCV(3) at depth YLEV
C
C     IF(YCM.LT.YPASS)THEN
C     note that YLEV is not yet defined here
C       DVAL=BCV(3)*(YPASS-YCM)/(YPASS-YLEV)
C       IBCTYP(IYP)=2
C       IBCTYP(I)=2
C       QBND(IYP)=DVAL
C       QBND(I)=DVAL
C     END IF

C
C      eveywhere below YPASS (BCV(1)) a vertical stress condition is reset to
C      a vertical velocity condition (BCV(2))
C      This boundary modification is specific to the RT calculation with a
C      reflecting boundary at x=0.  It assumes that TY=0 on x=0 and that 
C      the boundary condition on the lower boundary is initially TY=0 and 
C      can change to UY=0 if points fall below a specified level and may
C      be reset to TY=0 if a point rebounds
C
          END IF
        ELSE IF(ITYPE.EQ.14)THEN
           YTEST = BCV(1)*0.8
           IF(YCM.LT.YTEST) THEN
             IF(YCM.LT.BCV(1))THEN
               IF (IBCTYP(IYP).EQ.1) THEN
                 IBCTYP(IYP)=0
                 QBND(IYP)=BCV(2)
                 ISET=ISET+1
               END IF
             ELSE IF (IBCTYP(IYP).EQ.0) THEN
               IBCTYP(IYP)=1
               QBND(IYP)=0
               ISET=ISET+1
             END IF
          END IF
C
C      re-set the friction coefficient to increase from BCV(3)
C      to BCV(4) on a friction boundary condition for Y-coordinate
C      in the range     BCV(1)-0.5*BCV(2) < Y < BCV(1)+0.5*BCV(2)
C      as used when a sinking drip encounters the 660 km transition
C      the change is implemented as a weighted geometric mean
C
        ELSE IF(ITYPE.EQ.15)THEN
          YTEST = BCV(1)+0.5*BCV(2)
          IF(YCM.LT.YTEST) THEN
            WGHT=(YTEST-YCM)/BCV(2)
            IF(WGHT.GT.1.0)WGHT=1.0
            CWGHT=1.0-WGHT
            FNEW=EXP(WGHT*LOG(BCV(4))+CWGHT*LOG(BCV(3)))
            IF (IBCTYP(IYP).EQ.2)QBND(IYP)=FNEW
            IF (IBCTYP(IXP).EQ.2)QBND(IXP)=FNEW
            ISET=ISET+1
          END IF
        END IF
   10 CONTINUE
      IF(ISET.GT.0)THEN
        WRITE(LUW,10003)ITYPE,ISET
        WRITE(LSC,10003)ITYPE,ISET
      ENDIF
10003 FORMAT('Nodal Boundary conditions modified, ITYPE =',I3,
     :       '.  nodes altered: ',I5)
C       WRITE(LUW,*)'QBND array after MENDBC'
C       CALL MATPRT(QBND,NBP,2,NBP*2,LUW)
      RETURN
      END
      SUBROUTINE EXYLIM(EX,EY,NOR,NUP,NN,UVP,NROWS,
     1 XMIN,XMAX,YMIN,YMAX)
C
C     Routine to locate position of extreme mesh nodes
C
C     DOUBLE PRECISION UVP
      DIMENSION EX(NUP),EY(NUP),NOR(NUP),UVP(NROWS)
      DIMENSION XMIN(2),XMAX(2),YMIN(2),YMAX(2)
C
      XMINT=EX(1)
      YMINT=EY(1)
      XMAXT=EX(1)
      YMAXT=EY(1)
      NPOSXA=1
      NPOSYA=1
      NPOSXI=1
      NPOSYI=1
      DO 100 N=2,NUP
      K=NOR(N)
      IF(K.LE.NN)THEN
        X=EX(K)
        Y=EY(K)
        IF(X.GT.XMAXT)THEN
          NPOSXA=N
          XMAXT=X
        END IF
        IF(X.LT.XMINT)THEN
          NPOSXI=N
          XMINT=X
        END IF
        IF(Y.GT.YMAXT)THEN
          NPOSYA=N
          YMAXT=Y
        END IF
        IF(Y.LT.YMINT)THEN
          NPOSYI=N
          YMINT=Y
        END IF
      END IF
  100 CONTINUE
      KP=NOR(NPOSXI)
      XMIN(1)=EX(KP)
      XMIN(2)=EY(KP)
      KP=NOR(NPOSXA)
      XMAX(1)=EX(KP)
      XMAX(2)=EY(KP)
      KP=NOR(NPOSYI)
      YMIN(1)=EX(KP)
      YMIN(2)=EY(KP)
      VLEFT=UVP(NUP+NPOSYI)
      KP=NOR(NPOSYA)
      YMAX(1)=EX(KP)
      YMAX(2)=EY(KP)
      RETURN
      END
      SUBROUTINE SSQLIM(SSQ,EX,EY,NOR,NUP,SMIN,SMAX)
C
C     Routine to locate positions of extreme thicknesses
C
      DIMENSION SSQ(NUP),EX(NUP),EY(NUP),NOR(NUP)
      DIMENSION SMIN(3),SMAX(3)
C
      SMINT=SSQ(1)
      SMAXT=SSQ(1)
      NPOSSA=1
      NPOSSI=1
      DO 100 N=2,NUP
        SS=SSQ(N)
        IF(SS.GT.SMAXT)THEN
          NPOSSA=N
          SMAXT=SS
        END IF
        IF(SS.LT.SMINT)THEN
          NPOSSI=N
          SMINT=SS
        END IF
  100 CONTINUE
      KP=NOR(NPOSSI)
      SMIN(2)=EX(KP)
      SMIN(3)=EY(KP)
      KP=NOR(NPOSSA)
      SMAX(2)=EX(KP)
      SMAX(3)=EY(KP)
      SMIN(1)=EXP(SMINT)
      SMAX(1)=EXP(SMAXT)
      RETURN
      END
      SUBROUTINE ROTLIM(FROT,EX,EY,NOR,NUP,ROMIN,ROMAX)
C
C     Routine to locate positions of extreme rotation
C
      DIMENSION FROT(NUP),EX(NUP),EY(NUP),NOR(NUP)
      DIMENSION ROMIN(3),ROMAX(3)
C
      RMINT=FROT(1)
      RMAXT=FROT(2)
      NPOSSA=1
      NPOSSI=1
      DO 100 N=2,NUP
        SS=FROT(N)
        IF(SS.GT.RMAXT)THEN
          NPOSSA=N
          RMAXT=SS
        END IF
        IF(SS.LT.RMINT)THEN
          NPOSSI=N
          RMINT=SS
        END IF
  100 CONTINUE
      KP=NOR(NPOSSI)
      ROMIN(2)=EX(KP)
      ROMIN(3)=EY(KP)
      KP=NOR(NPOSSA)
      ROMAX(2)=EX(KP)
      ROMAX(3)=EY(KP)
      ROMIN(1)=RMINT
      ROMAX(1)=RMAXT
      RETURN
      END
      SUBROUTINE XYUVP(PVALUE,EX,EY,UVP,NOR,THDINT,
     :                 NUP,NN,NCOMP,NROWS,
     :                 FROT,SSQ,SMIN,SMAX,
     :                 UMAX,VMAX,UMIN,VMIN,
     :                 XMIN,XMAX,YMIN,YMAX,
     :                 ED2IMIN,ED2IMAX,VISCMIN,VISCMAX,
     :                 THDIMIN,THDIMAX,WKINTMIN,WKINTMAX,
     :                 TIME,PTLOC,VHB,RLV,
     :                 LEM,IELFIX,NE,NFP,IVV,ICR,
     :                 MSINDX,MEASUR,MSNODE,LUW,LSC,IPR)
C
C     Routine to record position, velocity and thickness at nodes
C
C     DOUBLE PRECISION UVP
      INCLUDE 'indices.parameters'
      INCLUDE 'input.parameters'
      INCLUDE 'limits.parameters'
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      CHARACTER FORM*9
      DIMENSION EX(NUP),EY(NUP),SSQ(NUP),FROT(NUP)
      DIMENSION VHB(8,NE),THDINT(7,NE),RLV(64)
      DIMENSION LEM(6,NE),NOR(NUP),UVP(NROWS)
      DIMENSION MEASUR(MAXMEAS),MSNODE(MAXMEAS)
      DIMENSION PVALUE(MAXMEAS),PTLOC(3,MAXMEAS)
      DIMENSION XMIN(2),XMAX(2),YMIN(2),YMAX(2),SMIN(3),SMAX(3)
C    local arrays
      DIMENSION PMIN(3),PMAX(3),ROMIN(3),ROMAX(3)
      DIMENSION COOL(3)
C
      CALL EXYLIM(EX,EY,NOR,NUP,NN,UVP,NROWS,XMIN,XMAX,YMIN,YMAX)
      IF(ISSQACTIVE.EQ.1) CALL SSQLIM(SSQ,EX,EY,NOR,NUP,SMIN,SMAX)
      IF(IROTACTIVE.EQ.1) CALL ROTLIM(FROT,EX,EY,NOR,NUP,ROMIN,ROMAX)
      IF(IPR.NE.0)THEN
        WRITE(LSC,10002)XMIN(1),XMIN(2),XMAX(1),XMAX(2),
     1                YMIN(1),YMIN(2),YMAX(1),YMAX(2)
        WRITE(LUW,10002)XMIN(1),XMIN(2),XMAX(1),XMAX(2),
     1                YMIN(1),YMIN(2),YMAX(1),YMAX(2)
10002   FORMAT('    XMINx    XMINy    XMAXx    XMAXy',
     1'    YMINx    YMINy    YMAXx    YMAXy',/,8F9.4)
C        WRITE(LSC,10003)SMIN(2),SMIN(3),SMIN(1),SMAX(2),
C     :                  SMAX(3),SMAX(1)
C        WRITE(LUW,10003)SMIN(2),SMIN(3),SMIN(1),SMAX(2),
C     :                  SMAX(3),SMAX(1)
C10003   FORMAT('    SMINx    SMINy     SMIN    SMAXx'
C     :         '    SMAXy     SMAX',/,6F9.4)
      END IF
      IF (NN.GT.0) THEN
        CALL PRESLIM(UVP(NROWS-NN+1),EX,EY,NN,NOR,NUP,PMIN,PMAX,
     :               LUW,LSC)
        IF(IPR.NE.0)THEN
          WRITE(LSC,10004)PMIN(2),PMIN(3),PMIN(1),PMAX(2),
     :                  PMAX(3),PMAX(1)
          WRITE(LUW,10004)PMIN(2),PMIN(3),PMIN(1),PMAX(2),
     :                  PMAX(3),PMAX(1)
10004     FORMAT('    PMINx      PMINy       PMIN      PMAXx',
     :         '      PMAXy       PMAX',/,6F11.5)
        END IF
      END IF
      DO 110 J=1,MSINDX
        MEAS=MEASUR(J)
C       WRITE(LSC,*)'XYUVP: J =',J,'  MEAS = ',MEAS
        PVALUE(J)=0
      IF(MEAS.EQ.INDXTIME)THEN
        PVALUE(J)=TIME
      ELSE IF(MEAS.EQ.INDXXMIN)THEN
        PVALUE(J)=XMIN(1)
      ELSE IF(MEAS.EQ.INDXXMAX)THEN
        PVALUE(J)=XMAX(1)
      ELSE IF(MEAS.EQ.INDXYMIN)THEN
        PVALUE(J)=YMIN(2)
      ELSE IF(MEAS.EQ.INDXYMAX)THEN
        PVALUE(J)=YMAX(2)
      ELSE IF(MEAS.EQ.INDXUMIN)THEN
        PVALUE(J)=UMIN
      ELSE IF(MEAS.EQ.INDXUMAX)THEN
        PVALUE(J)=UMAX
      ELSE IF(MEAS.EQ.INDXVMIN)THEN
        PVALUE(J)=VMIN
      ELSE IF(MEAS.EQ.INDXVMAX)THEN
        PVALUE(J)=VMAX
      ELSE IF(MEAS.EQ.INDXLTMIN)THEN
        PVALUE(J)=SMIN(1)
      ELSE IF(MEAS.EQ.INDXLTMAX)THEN
        PVALUE(J)=SMAX(1)
      ELSE IF(MEAS.EQ.INDXPRESMIN)THEN
        PVALUE(J)=PMIN(1)
      ELSE IF(MEAS.EQ.INDXPRESMAX)THEN
        PVALUE(J)=PMAX(1)
      ELSE IF(MEAS.EQ.INDXROMIN)THEN
        PVALUE(J)=ROMIN(1)
      ELSE IF(MEAS.EQ.INDXROMAX)THEN
        PVALUE(J)=ROMAX(1)
      ELSE IF(MEAS.EQ.INDXED2IMIN)THEN
        PVALUE(J)=ED2IMIN
      ELSE IF(MEAS.EQ.INDXED2IMAX)THEN
        PVALUE(J)=ED2IMAX
      ELSE IF(MEAS.EQ.INDXTHDIMIN)THEN
        PVALUE(J)=THDIMIN
      ELSE IF(MEAS.EQ.INDXTHDIMAX)THEN
        PVALUE(J)=THDIMAX
      ELSE IF(MEAS.EQ.INDXWKINTMIN)THEN
        PVALUE(J)=WKINTMIN
      ELSE IF(MEAS.EQ.INDXWKINTMAX)THEN
        PVALUE(J)=WKINTMAX
      ELSE IF(MEAS.EQ.INDXVISCMIN)THEN
        PVALUE(J)=VISCMIN
      ELSE IF(MEAS.EQ.INDXVISCMAX)THEN
        PVALUE(J)=VISCMAX
      ELSE IF((MEAS.GE.INDXXX).AND.(MEAS.LE.INDXRO))THEN
        IF(MSNODE(J).GT.0)THEN
          K=NOR(MSNODE(J))
          PTLOC(1,J) = EX(K)
          PTLOC(2,J) = EY(K)
          IF(MEAS.EQ.INDXXX)PVALUE(J)=EX(K)
          IF(MEAS.EQ.INDXYY)PVALUE(J)=EY(K)
          IF(MEAS.EQ.INDXUX)PVALUE(J)=UVP(MSNODE(J))
          IF(MEAS.EQ.INDXUY)PVALUE(J)=UVP(NUP+MSNODE(J))
          IF(MEAS.EQ.INDXLT)PVALUE(J)=EXP(SSQ(MSNODE(J)))
          IF(MEAS.EQ.INDXPRES)PVALUE(J)=UVP(NROWS-NN+K)
          IF(MEAS.EQ.INDXRO)PVALUE(J)=FROT(MSNODE(J))
        ENDIF
      ELSE IF(MEAS.EQ.INDXWKINT)THEN
C
C       need something new to deal with this one
C
C     ELSE IF(MEAS.GE.INDXEDXX.AND.MEAS.LE.INDXFOLT.AND.
C    :        MSNODE(J).GT.0)THEN
      ELSE
        IF(MSNODE(J).GT.0)THEN
          DO N=1,3
            COOL(N)=PTLOC(N,J)
          END DO
          CALL STRAIE(MEAS,MSNODE(J),PVALUE(J),COOL,
     :             RLV(IVC),RLV(ISE),RLV(IARGANP),RLV(IBRGANP),
     :             RLV(IHLENSC),RLV(IBDEPSC),RLV(IBIG),RLV(IGAMMA),
     :             RLV(ITREF),EX,EY,VHB,THDINT,UVP,SSQ,IELFIX,
     :             LEM,NOR,NE,NUP,NROWS,NFP,IVV,NCOMP,ICR)
        END IF
      END IF
  110 CONTINUE
C      WRITE(FORM,10002)MSINDX
C10002 FORMAT('(',I2,'G15.7)')
C      WRITE(LWR,FORM)(PVALUE(J),J=1,MSINDX)
      RETURN
      END
      SUBROUTINE UVPLIM(UVMAX,UVP,EX,EY,NOR,
     :                  UMAX,VMAX,UMIN,VMIN,NROWS,NUP,
     :                  LUW,LSC,IPR)
C
C     Routine to locate position of extreme velocity
C
C     DOUBLE PRECISION UVP(NROWS)
      DIMENSION UVP(NROWS)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
C
      UMAX=UVP(1)
      UMIN=UMAX
      VMAX=UVP(NUP+1)
      VMIN=VMAX
      NUMAX=1
      NVMAX=1
      NUVMAX=1
      NUMIN=1
      NVMIN=1
      UVMAX=SQRT(UMAX*UMAX+VMAX*VMAX)
      DO 100 N=1,NUP
      U=UVP(N)
      V=UVP(N+NUP)
      UV=SQRT(U*U+V*V)
      IF(U.GT.UMAX)THEN
        NUMAX=N
        UMAX=U
      ELSE IF(U.LT.UMIN) THEN
        NUMIN=N
        UMIN=U
      END IF
      IF(V.GT.VMAX)THEN
        NVMAX=N
        VMAX=V
      ELSE IF(V.LT.VMIN) THEN
        NVMIN=N
        VMIN=V
      END IF
      IF(UV.GT.UVMAX)THEN
        NUVMAX=N
        UVMAX=UV
      END IF
  100 CONTINUE
      U = UMAX
      IF (ABS(UMIN).GT.ABS(UMAX)) THEN
        U=UMIN
        NUMAX = NUMIN
      END IF
      V = VMAX
      IF (ABS(VMIN).GT.ABS(VMAX)) THEN
        V=VMIN
        NVMAX = NVMIN
      END IF
      KP=NOR(NUMAX)
      XUMAX=EX(KP)
      YUMAX=EY(KP)
      KP=NOR(NVMAX)
      XVMAX=EX(KP)
      YVMAX=EY(KP)
      KP=NOR(NUVMAX)
      XUVMAX=EX(KP)
      YUVMAX=EY(KP)

      IF (IPR.GT.0) THEN
        WRITE(LSC,10002)XUMAX,YUMAX,U,XVMAX,YVMAX,V,XUVMAX,
     1                 YUVMAX,UVMAX
        WRITE(LUW,10002)XUMAX,YUMAX,U,XVMAX,YVMAX,V,XUVMAX,
     1                 YUVMAX,UVMAX
      END IF
10002 FORMAT('    XUMAX    YUMAX     UMAX    XVMAX    YVMAX     VMAX',
     1'   XUVMAX   YUVMAX    UVMAX',/,9F9.4)
      RETURN
      END

      SUBROUTINE PRESLIM(PRES,EX,EY,NN,NOR,NUP,PMIN,PMAX,LUW,LSC)
C
C     Routine to locate positions of extreme pressure
C
      DIMENSION PRES(NN),EX(NUP),EY(NUP),NOR(NUP)
      DIMENSION PMIN(3),PMAX(3)
C
      PMINT=PRES(1)
      PMAXT=PRES(1)
      NPOSSA=1
      NPOSSI=1
      DO 100 N=2,NN
        P=PRES(N)
        IF(P.GT.PMAXT)THEN
          NPOSSA=N
          PMAXT=P
        END IF
        IF(P.LT.PMINT)THEN
          NPOSSI=N
          PMINT=P
        END IF
  100 CONTINUE
      KP=NOR(NPOSSI)
      PMIN(2)=EX(KP)
      PMIN(3)=EY(KP)
      KP=NOR(NPOSSA)
      PMAX(2)=EX(KP)
      PMAX(3)=EY(KP)
      PMIN(1)=PMINT
      PMAX(1)=PMAXT
      RETURN
      END

      SUBROUTINE GETNODENO(MEASUR,PTLOC,MSNODE,MSINDX,
     :                     EX,EY,NOR,NUP,NN,LUW,LSC,IERR)
C
C    This subroutine identifies a node number from its
C     present cartesian coordinates in PTLOC
C
      INCLUDE 'input.parameters'
      INCLUDE 'limits.parameters'
      DIMENSION EX(NUP),EY(NUP),NOR(NUP)
      DIMENSION MEASUR(MSINDX),PTLOC(3,MSINDX),MSNODE(MSINDX)
      DIMENSION PTDIFF(2,MAXMEAS),BBOX(2,4),DIST(MAXMEAS)
C
      EPS=1.E-4
      IERR=0
C     CALL FINDBBOX(EX,EY,BBOX,NUP)
      CALL FINDMLIM(EX,EY,NUP,XMIN,XMAX,YMIN,YMAX)
      PTDIFF(1,1) = ABS(XMAX - XMIN) * 0.1
      PTDIFF(2,1) = ABS(YMAX - YMIN) * 0.1
      DIST(1)=PTDIFF(1,1)*PTDIFF(1,1) + PTDIFF(2,1)*PTDIFF(2,1)
      DO N=2,MSINDX
        DIST(N)=DIST(1)
      ENDDO
      DO 100 N=1,NUP
        K=NOR(N)
        IF(K.LE.NN)THEN
          DO 80 M=1,MSINDX
C           IF((MEASUR(M).GE.INDXXX.AND.MEASUR(M).LE.INDXUY).OR.
C    1       (MEASUR(M).EQ.INDXLT).OR.(MEASUR(M).EQ.INDXPRES).OR.
C    :       (MEASUR(M).EQ.INDXRO))THEN
            IF((MEASUR(M).GE.INDXXX).AND.(MEASUR(M).LE.INDXRO))THEN
              XDIF=ABS(EX(K)-PTLOC(1,M))
              IF(XDIF.LE.PTDIFF(1,1))THEN
                YDIF=ABS(EY(K)-PTLOC(2,M))
                IF(YDIF.LE.PTDIFF(2,1).AND.
     :             ((YDIF*YDIF+XDIF*XDIF).LT.DIST(M)))THEN
                  DIST(M) = YDIF*YDIF+XDIF*XDIF
C                 IF(MSNODE(M).NE.0)THEN
C                   WRITE(*,*)'Warning: multiple node id in GETNODENO'
C                   IERR=1
C                 END IF
                  MSNODE(M)=N
                END IF
              END IF
            END IF
   80     CONTINUE
        END IF
  100 CONTINUE
      DO 120 I=1,MSINDX
C       IF((MEASUR(I).GE.INDXXX.AND.MEASUR(I).LE.INDXUY).OR.
C    :     (MEASUR(I).EQ.INDXLT).OR.(MEASUR(I).EQ.INDXPRES).OR.
C    :     (MEASUR(I).EQ.INDXRO))THEN
        IF((MEASUR(I).GE.INDXXX).AND.(MEASUR(I).LE.INDXRO))THEN
          IF(MSNODE(I).EQ.0)THEN
            WRITE(*,*)'Warning: failed to locate ',I,'th variable ',
     :                'in SERIES'
            IERR=1
          ELSE IF (MEASUR(I).EQ.INDXPRES) THEN
            WRITE(LSC,1001)MSNODE(I),EX(NOR(MSNODE(I))),
     :                    EY(NOR(MSNODE(I))),PTLOC(1,I),PTLOC(2,I)
            WRITE(LUW,1001)MSNODE(I),EX(NOR(MSNODE(I))),
     :                    EY(NOR(MSNODE(I))),PTLOC(1,I),PTLOC(2,I)
          ELSE IF (MEASUR(I).EQ.INDXRO) THEN
            WRITE(LSC,1002)MSNODE(I),EX(NOR(MSNODE(I))),
     :                    EY(NOR(MSNODE(I))),PTLOC(1,I),PTLOC(2,I)
            WRITE(LUW,1002)MSNODE(I),EX(NOR(MSNODE(I))),
     :                    EY(NOR(MSNODE(I))),PTLOC(1,I),PTLOC(2,I)
          END IF
        END IF
  120 CONTINUE
1001  FORMAT('SERIES: PRES using node ',I6,',',F9.4,',',F9.4,
     :       ' at',F9.4,',',F9.4)
1002  FORMAT('SERIES: RO using node ',I6,',',F9.4,',',F9.4,
     :       ' at',F9.4,',',F9.4)
      RETURN
      END

      SUBROUTINE GETELEMENTNO(MEASUR,PTLOC,MSNODE,MSINDX,
     :                     EX,EY,NOR,NUP,NN,
     :                     LEM,NE,LUW,LSC,IERR)
C
C    This subroutine identifies the relevant element and calculates the
C    local coordinates of the points that are identified in the SERIES
C    command. Initial cartesian coordinates are in PTLOC.  These are 
C    replaced by natural coordinates in this subroutine.
C    The element number is stored in MSNODE
C
      INCLUDE 'input.parameters'
      INCLUDE 'limits.parameters'
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE),NOR(NUP)
      DIMENSION MEASUR(MSINDX),PTLOC(3,MSINDX),MSNODE(MSINDX)
      DIMENSION COOL(3)
C
C  EPS reduced by a factor of 10 relative to previous value
C
      EPS=5.E-4
      IERR=0
      DO 18 M=1,MSINDX
        IF (MEASUR(M).GE.INDXEDXX.AND.MEASUR(M).LE.INDXFOLT)THEN
          MSNODE(M)=0
        END IF
  18  CONTINUE
C
C     for each element in the mesh
C
      DO 80 N=1,NE
C
C     Calculate the geometrical coefficients
C
        LK1=NOR(LEM(1,N))
        LK2=NOR(LEM(2,N))
        LK3=NOR(LEM(3,N))
        X1=EX(LK1)
        X2=EX(LK2)
        X3=EX(LK3)
        Y1=EY(LK1)
        Y2=EY(LK2)
        Y3=EY(LK3)
C
C  in this method we effectively use (X1,Y1) as a local origin, for
C     computing local coordinates and determinants
C
        Y3M1=Y3-Y1
        X3M1=X3-X1
        Y1M2=Y1-Y2
        X1M2=X1-X2
        TRI=Y1M2*X3M1-X1M2*Y3M1
C
C    TRI is twice the area of the triangle element
C    For each point in the SERIES list, which requires local coordinates:
C
        DO 180 M=1,MSINDX
          IF(MSNODE(M).EQ.0.AND.
     :      MEASUR(M).GE.INDXEDXX.AND.MEASUR(M).LE.INDXFOLT)THEN
C
C   amendment to next 2 lines shifts local origin to (X1,Y1)
C
            XP=PTLOC(1,M)-X1
            YP=PTLOC(2,M)-Y1
C
C     Calculate the natural coordinates, allow points to be just outside
C      element with error EPS (drag them back in)
C
            COOL(2)=XP*Y3M1-YP*X3M1
            COOL(3)=XP*Y1M2-YP*X1M2
            COOL(1)=TRI-COOL(2)-COOL(3)
            DO K=1,3
              CNL=COOL(K)/TRI
              IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 180
              IF (CNL.GT.1.0) CNL=1.0
              IF (CNL.LT.0.0) CNL=0.0
              COOL(K)=CNL
            ENDDO
C
C     If we reach here the point is within the triangle
C
            PTLOC(1,M)=COOL(1)
            PTLOC(2,M)=COOL(2)
            PTLOC(3,M)=COOL(3)
            MSNODE(M)=N
          END IF
 180    CONTINUE
 80   CONTINUE
C
C   check that all points could be located in the mesh
C
      DO 120 I=1,MSINDX
        IF( MEASUR(I).GE.INDXEDXX.AND.MEASUR(I).LE.INDXFOLT)THEN
          IF(MSNODE(I).EQ.0) THEN
            WRITE(*,*)'Warning: failed to locate position ',
     :        PTLOC(1,I),' ',PTLOC(2,I),' ',I,' in SERIES'
C           IERR=1
          ELSE
            LK1=NOR(LEM(1,MSNODE(I)))
            LK2=NOR(LEM(2,MSNODE(I)))
            LK3=NOR(LEM(3,MSNODE(I)))
            X1=EX(LK1)
            X2=EX(LK2)
            X3=EX(LK3)
            Y1=EY(LK1)
            Y2=EY(LK2)
            Y3=EY(LK3)
            XP=PTLOC(1,I)*X1 + PTLOC(2,I)*X2 + PTLOC(3,I)*X3
            YP=PTLOC(1,I)*Y1 + PTLOC(2,I)*Y2 + PTLOC(3,I)*Y3
            WRITE(LSC,1001)I,MEASUR(I),MSNODE(I),XP,YP
            WRITE(LUW,1001)I,MEASUR(I),MSNODE(I),XP,YP
          END IF
        END IF
  120 CONTINUE
 1001 FORMAT('Series entry ',I3,', measure ',I3,', element ',I6,
     :       ', mesh coords ',F9.4,',',F9.4)
      RETURN
      END
C
      SUBROUTINE RESETSERIES(LEM,NOR,NE,NUP,MEASUR,MSNODE,MSINDX,
     :                       PTLOC,EX,EY)
      INCLUDE 'input.parameters'
      INCLUDE 'limits.parameters'
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE),NOR(NUP)
      DIMENSION MEASUR(MSINDX),PTLOC(3,MSINDX),MSNODE(MSINDX)
      DIMENSION COOL(3)

      DO I=1,MSINDX
        IF( MEASUR(I).GE.INDXEDXX.AND.MEASUR(I).LE.INDXFOLT.AND.
     :      MSNODE(I).GT.0)THEN
          N=MSNODE(I)
          DO K=1,3
            COOL(K) = PTLOC(K,I)
          END DO
          LK=NOR(LEM(1,N))
          LK2=NOR(LEM(2,N))
          LK3=NOR(LEM(3,N))
          X1=EX(LK)
          X2=EX(LK2)
          X3=EX(LK3)
          Y1=EY(LK)
          Y2=EY(LK2)
          Y3=EY(LK3)
          XP=COOL(1)*X1 + COOL(2)*X2 + COOL(3)*X3
          YP=COOL(1)*Y1 + COOL(2)*Y2 + COOL(3)*Y3
          PTLOC(1,I) = XP
          PTLOC(2,I) = YP
          PTLOC(3,I) = 0.0
        END IF
      END DO
      RETURN
      END
C---------------------------------------------------------------
C
      SUBROUTINE PMINMAX(QLOAD,NDFREE,IRANK)
C
C    this is a general purpose routine to assist with debugging;
C    advises how many non-zero elements, and what are min and max of
C    array entriesr. IRANK is processor number in parallel environment
C    (ignored in serial environment)
C
      DIMENSION QLOAD(NDFREE)
      NZSUM=0
      QMIN=QLOAD(1)
      QMAX=QLOAD(1)
      DO KK=1,NDFREE
        IF(QLOAD(KK).NE.0.0)NZSUM=NZSUM+1
        IF(QLOAD(KK).GT.QMAX)THEN
          QMAX=QLOAD(KK)
        END IF
        IF(QLOAD(KK).LT.QMIN)THEN
          QMIN=QLOAD(KK)
        END IF
      ENDDO
      WRITE(*,*)'IRANK,NZSUM,QMIN,QMAX',IRANK,NZSUM,QMIN,QMAX
      RETURN
      END
