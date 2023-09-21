C*--------------------------------------------------------------------
C*    Basil / Sybil:   tvspres.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE TVSPRES(EX,EY,NOR,LEM,UVP,
     :                   NUP,NE,NROWS,NN,NFP,SE)
C
C    Routine to calculate different components of the strain or
C    stress tensor (or derived quantities) from the velocity field
C
      DOUBLE PRECISION UVP,XMESH,BMESH
      DOUBLE PRECISION BY,CX,DNDP
      DOUBLE PRECISION VF,TRI,X2,X3,Y2,Y3,ANGL,DB2,DA2,SANGL,CSANGL,PI
      DOUBLE PRECISION XP,YP,DUDX,DUDY,DVDX,DVDY,DWDZ,UI,VI,EDXY,ED2I
      DOUBLE PRECISION SI,SELOC
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION UVP(NROWS)
      DIMENSION NOR(NUP)
C      local arrays
      DIMENSION XMESH(NUP)
      DIMENSION BMESH(NUP)
      DIMENSION BY(3),CX(3),DNDP(84)
      DATA PI/3.14159265358979/
      SAVE PI
C
C    Zero the mesh first
C
      DO 10 I=1,NUP
        BMESH(I)=0.0
        XMESH(I)=0.0
   10 CONTINUE
      NUP2=NUP*2
      SELOC=SE
C
C    Look at each element in turn
C
      DO 62 N=1,NE
C
C     Check if the viscosity of the element is enhanced
C     *** to be fixed: Must check IVIS here and then adjust VF using VHB ***
C
        VF=1.D0
C
C     Calculate the geometrical coefficients
C
        DO 20 K=1,3
          K2=MOD(K,3)+1
          K3=MOD(K+1,3)+1
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          X2=EX(LK2)
          X3=EX(LK3)
          Y2=EY(LK2)
          Y3=EY(LK3)
          BY(K)=Y2-Y3
          CX(K)=X3-X2
   20   CONTINUE
        TRI=CX(3)*BY(2)-CX(2)*BY(3)
C
C    TRI is twice the area of the triangle element
C    Get the derivatives of the interpolation function at the nodes
C
        CALL DNCOM(5,BY,CX,DNDP)
C
C     Find the internal angle made by the boundaries of the
C     triangle at this node
C
        DO K=1,7
          IF(K.EQ.7)GO TO 29
          LK=LEM(K,N)
          ANGL=0.5
          IF(K.GE.4)GO TO 29
          K2=MOD(K,3)+1
          K3=MOD(K+1,3)+1
          DB2= BY(K2)*BY(K2) + CX(K2)*CX(K2)
          DA2= BY(K3)*BY(K3) + CX(K3)*CX(K3)
          SANGL=TRI/SQRT(DB2*DA2)
          CSANGL=1.0-SANGL*SANGL
          IF(CSANGL.LT.0.)CSANGL=0.0
          CSANGL=SQRT(CSANGL)
          IF (SANGL.EQ.0.AND.CSANGL.EQ.0) THEN
            ANGL=0
          ELSE
            ANGL=0.5*ATAN2(SANGL,CSANGL)/PI
          END IF
C
C    Location of node K
C
   29     IF(K.LE.3)THEN
            ILK=NOR(LEM(K,N))
            XP=EX(ILK)
            YP=EY(ILK)
          ELSE IF((K.GE.4).AND.(K.LE.6))THEN
            KK1=MOD(K+1,3)+1
            KK2=MOD(K+2,3)+1
            ILK1=NOR(LEM(KK1,N))
            ILK2=NOR(LEM(KK2,N))
            XP=0.5*(EX(ILK1)+EX(ILK2))
            YP=0.5*(EY(ILK1)+EY(ILK2))
          ELSE
            ILK1=NOR(LEM(1,N))
            ILK2=NOR(LEM(2,N))
            ILK3=NOR(LEM(3,N))
            XP=(EX(ILK1)+EX(ILK2)+EX(ILK3))/3.0
            YP=(EY(ILK1)+EY(ILK2)+EY(ILK3))/3.0
          END IF
C
C   Calculate the velocity derivatives at node K
C
   25     DUDX=0.0
          DUDY=0.0
          DVDX=0.0
          DVDY=0.0
C
C    Sum the interpolation functions
C
          DO 55 I=1,6
            NI=LEM(I,N)
            KIN=(I-1)*14 + (K-1)*2 + 1
            UI=UVP(NI)
            VI=UVP(NI+NUP)
C
C    du/dx, du/dy, dv/dx, dv/dy at node K
C
            DUDX=DUDX + UI*DNDP(KIN)
            DUDY=DUDY + UI*DNDP(KIN+1)
            DVDX=DVDX + VI*DNDP(KIN)
            DVDY=DVDY + VI*DNDP(KIN+1)
   55     CONTINUE
          DUDX=DUDX/TRI
          DVDX=DVDX/TRI
          DUDY=DUDY/TRI
          DVDY=DVDY/TRI
          DWDZ=-DUDX-DVDY
          EDXY=0.5D0*(DUDY+DVDX)
C
C    Find and calculate the desired quantity and put it in XMESH
C
          ED2I=SQRT(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ + 2.D0*EDXY*EDXY)
C
C      effective viscosity
C
          IF(ED2I.LE.0.D0)THEN
            SI=1.D5
          ELSE
            SI = VF*(ED2I**(1.D0/SELOC - 1.D0))
            IF(SI.GT.1.0D5)SI=1.0D5
          END IF
C
C    To obtain stress multiply by viscosity
C
          SI=2.D0*SI*DWDZ
          XMESH(LK)=XMESH(LK) + SI*ANGL
          BMESH(LK)=BMESH(LK) + ANGL
        ENDDO        ! on K
   62 CONTINUE       ! on N
C
C    Normalise by the total angle around the node
C
      DO 60 I=1,NUP
        XMESH(I)=XMESH(I)/BMESH(I)
   60 CONTINUE
      DO 70 I=1,NUP
         IN=NOR(I)
         IF(IN.LE.NN.AND.IN.GE.1) THEN
           UVP(NUP2+NFP+IN)=-XMESH(I)
         END IF
   70 CONTINUE
  900 CONTINUE
      RETURN
      END
