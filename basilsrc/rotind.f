C*--------------------------------------------------------------------
C*    Basil / Sybil:   rotind.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE LITHBC(HLENSC,EX,EY,NOR,IBC,IBNGH,IBCTYP,
     :                  QBND,ARGAN,NUP,NBP,NN,ITYPE,
     :                  NCOMP,CENTLNG,BCV,LUW)
C
C    This routine sets the traction components on the boundary
C    to simulate a traction which is normal to the boundary.  It
C    may be called at each timestep in order to adjust the traction
C    vector as the boundary is deformed.
C
C    In thin sheet mode (ITYPE.LE.11) it is assumed that the normal
C    stress outside the boundary is due to a lithostatic column with
C    original crustal thickness (designed for the lithostatic boundary
C    condition in the collision models)   These traction values can be
C    modified in an adhoc way, by scaling the x and y cmponents by BCV(1,2)
C
C    In vertical section mode (ITYPE.EQ.12) it is assumed that the normal
C    stress is proportional to the Y-coordinate: p = -BCV(1)*(BCV(2)-Y)
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION IBC(NBP),IBNGH(NBP*2),IBCTYP(NBP*2)
      DIMENSION BCV(5)
C
      SAVE PI,PION2,DTOR
      DATA PI,PION2,DTOR/3.141592653,1.570796327,1.7453292519943295E-2/
C
C   CT0 is the hydrostatic part of the stress field due
C   due to the weight of the medium outside the stress-free
C    boundary.  It is set to zero here because a reference load
C    based on ARGAN and the reference crustal thickness has been
C    subtracted from the load matrix calculated by PLOAD
C
      CT0=ARGAN/(HLENSC*HLENSC)
      PVERT=1.0/6.0
      PMID=2.0/3.0
C
C     For each node on the boundary for which:
C             both components of traction are fixed
C
      DO 60 JX=1,NBP
        JY=JX+NBP
        IF((IBCTYP(JX).EQ.1).AND.(IBCTYP(JY).EQ.1))THEN
          NODE=IBC(JX)
          NNODE=NOR(NODE)
          XN=EX(NNODE)
          YN=EY(NNODE)
C
C      NODEA is next node clockwise, NODEB is next node anticlockwise
C
          NODEA=IBNGH(JX)
          NODEB=IBNGH(JY)
          NNODEA=NOR(NODEA)
          NNODEB=NOR(NODEB)
          XNA=EX(NNODEA)
          YNA=EY(NNODEA)
          XNB=EX(NNODEB)
          YNB=EY(NNODEB)
          IF(NNODE.GT.NN)THEN
            PC=PMID
            XDIF=XNB-XNA
            YDIF=YNB-YNA
          ELSE
            PC=PVERT
            XDIFA=XN-XNA
            YDIFA=YN-YNA
            XDIFB=XNB-XN
            YDIFB=YNB-YN 
            XDIF=XDIFA+XDIFB
            YDIF=YDIFA+YDIFB
          END IF
C
C          integral{Pnx.ds} = integral{P(-dy)}
C          integral{Pny.ds} = integral{P dx}
C   for vertices, int{P dx}+ int{P dx'} = int{P (dx+dx')}
C
C   these conditions apply only to boundary nodes which are fixed
C   traction in both directions (first for lithostatic boundary in
C    thin sheet problem, ITYPE.le.3)
C
          IF(ITYPE.LE.3)THEN
            QBND(JX)=0.0
            QBND(JY)=0.0
C
C    ITYPE = 11 applies tractions of separate magnitude to X and Y
C    directions (for application to arc problem in thin sheet case)
C    These traction components are directed outward relative to the
C    boundary orientation, multiplied by cosine of angle between outward
C    normal and principal coordinate direction (X or Y), scaled by
C    the lithostatic normal stress
C
          ELSE IF(ITYPE.EQ.11) THEN
            RMAG=PC*CT0
            QBND(JX)=BCV(1)*YDIF*RMAG
            QBND(JY)=-BCV(2)*XDIF*RMAG
C
C   the following block stops the 'bulge' spreading outside of X-range
C    0 to XLIM, boundary condition changed into reflecting boundary
C
            IF((XN.LE.BCV(3)).OR.(XN.GE.BCV(4)))THEN
              IBCTYP(JX)=0
              QBND(JX)=0.0
              QBND(JY)=0.0
            END IF
C
C    ITYPE = 12 assumes hydrostatic stress acting on any fixed traction
C     boundary node. BCV(1)=rho.g for the external fluid, BCV(2) is surface level.
C
          ELSE IF(ITYPE.EQ.12)THEN
            IF(YN.LE.BCV(2))THEN
              RMAG=-PC*BCV(1)*(BCV(2)-YN)
              QBND(JX)=YDIF*RMAG
              QBND(JY)=-XDIF*RMAG
            END IF
C
C   ITYPE = 13: same as 11, but allows two different external tractions
C   to be specified - dependent on x-ccordinate
C
          ELSE IF(ITYPE.EQ.13)THEN
            RMAG=PC*CT0
            IF(NCOMP.EQ.-1)THEN
              XND=CENTLNG + XN/(DTOR*SIN(YN))
            ELSE
              XND=XN
            END IF
            IF(XND.LE.BCV(3))THEN
              QBND(JX)=BCV(1)*YDIF*RMAG
              QBND(JY)=-BCV(1)*XDIF*RMAG
            ELSE
              QBND(JX)=BCV(2)*YDIF*RMAG
              QBND(JY)=-BCV(2)*XDIF*RMAG
            END IF
C
C   ITYPE = 14: same as 11,12 but allows 3 different external tractions
C   to be specified - dependent on x-ccordinate 
C    (enter 4,5 in degrees if spherical)
C
          ELSE IF(ITYPE.EQ.14)THEN
            RMAG=PC*CT0
            IF(NCOMP.EQ.-1)THEN
              XND=CENTLNG + XN/(DTOR*SIN(YN))
            ELSE
              XND=XN
            END IF
            IF(XND.LE.BCV(4))THEN
              QBND(JX)=BCV(1)*YDIF*RMAG
              QBND(JY)=-BCV(1)*XDIF*RMAG
            ELSE IF(XND.GE.BCV(5))THEN
              QBND(JX)=BCV(3)*YDIF*RMAG
              QBND(JY)=-BCV(3)*XDIF*RMAG
            ELSE
              QBND(JX)=BCV(2)*YDIF*RMAG
              QBND(JY)=-BCV(2)*XDIF*RMAG
            END IF
          END IF
        END IF
   60 CONTINUE

C
C     Output the boundary condition matrices
C
      IDBUG=0
      IF(IDBUG.NE.0)THEN
        LSC=6
        WRITE(LSC,10101)
10101 FORMAT(' ',/,'Array data from LITHBC:',/,
     1'       JX     IBC     X    Y    IBNGH1  IBNGH2 IBCTYPX   QBND X',
     2'    IBCTYPY   QBND Y')
        DO 600 JX=1,NBP
          JY=JX+NBP
          NNODE=NOR(IBC(JX))
          EXX=EX(NNODE)
          EYY=EY(NNODE)
          WRITE(LSC,10102)JX,IBC(JX),EXX,EYY,IBNGH(JX),IBNGH(JY),
     1                    IBCTYP(JX),QBND(JX),IBCTYP(JY),QBND(JY)
10102     FORMAT(2I8,2F8.4,3I8,G12.5,I8,G12.5)
  600 CONTINUE
      END IF
      RETURN
      END
C
      SUBROUTINE GENSHBC(NUP,NBP,BIG,QBND,EX,EY,IBC,IBCTYP,NOR,
     :                   VELXO,VELYO,TBXOFF,DUDX,DUDY,DVDX,DVDY)
C
C     This routine updates the boundary conditions to apply a constant 
C     general shear strain rate defined by DUDX, DUDY, DVDX, DVDY
C 
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION IBC(NBP),IBCTYP(NBP*2)
C
C    Update variables used in the external fault case 
C     This only considers external fault on y-boundaries
C
      VELXO=DUDX*TBXOFF
      VELYO=DVDX*TBXOFF
C
C    Loop through all the boundary points
C
C     WRITE(9,10059)
      DO 10 I=1,NBP
        NODE=IBC(I)
        NODUT=IBCTYP(I)
        NODVT=IBCTYP(I+NBP)
        NORN=NOR(IABS(NODE))
        XX=EX(NORN)
        YY=EY(NORN)
        U=DUDX*XX+DUDY*YY
        V=DVDX*XX+DVDY*YY
        IF(NODUT.EQ.0.OR.NODUT.EQ.10)QBND(I)=U
        IF(NODVT.EQ.0.OR.NODVT.EQ.10)QBND(I+NBP)=V
C
C     WRITE(9,10001)I,NODE,NODUT,NODVT,XX,YY,U,V
C
   10 CONTINUE
10059 FORMAT(' GENSHBC ')
10001 FORMAT('I=',I5,' NODE=',I5,' UXT=',I5,' UYT=',I5,
     1' X=',F5.3,' y=',F5.3,' UX=',G12.5,' UY=',G12.5)
      END
