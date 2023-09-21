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
      SUBROUTINE ROTIND(BIG,OMDOT,OMEGA,XLEN,YLEN,
     :                  EX,EY,QBND,IBC,IBCTYP,NOR,
     :                  NUP,NBP,NY1,LUW,IERR)
C
C    Adjusts the QBND matrix during execution to simulate the
C    rotation of the indenter.  Rotation of the arc is about the 
C    midpoint of the two cusps.  The midpoint moves at the constant
C    indentation velocity, so it is located at y = 1, x = x0+u0*t
C    It also rotates rigidly.  
C    The two wings have a cosine**2 taper on both velocity components.  
C    The rotation rate is OMDOT, in units of radians/dtu.
C    (XA,YA), (UA,VA) are the coordinates and velocity of the
C    eastern node
C    (XB,YB), (UB,VB) are the coordinates and velocity of the
C    western node
C    This routine is written specifically for the '33' mesh, i.e.
C    it refers to the node numbers of the two cusps, and relies
C    on the nodes appearing in the IBC matrix in sequential order
C    of increasing Y ordinate.
C    Probably the best way to set this up is to identify the two
C    cusp nodes and the two 'wing' nodes on the first call and
C    keep referring to those nodes on subsequent calls.  This remains to
C    be implemented.
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION IBC(NBP),IBCTYP(NBP*2)
      SAVE EPS,PI
      DATA EPS,PI/1.0E-4,3.14159265/
C
C   Note that 2 cusp nos below are for modified '33' mesh
C
      ICSP1=1022
      ICSP2=1702
      IWNGR=681
      IWNGL=2041
C
      XA=EX(NOR(ICSP1))
      YA=EY(NOR(ICSP1))
      XB=EX(NOR(ICSP2))
      YB=EY(NOR(ICSP2))
      XMID=0.5*(XA+XB)
      YMID=0.5*(YA+YB)
      HARC=0.5*SQRT((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB))
      IF ((XA-XB).EQ.0.AND.(YB-YA).EQ.0) THEN
        OMEGA=0
      ELSE
        OMEGA=ATAN2(XA-XB,YB-YA)
      END IF
      CS=COS(OMEGA)
      SN=SIN(OMEGA)
      UA=1.0+HARC*OMDOT*CS
      VA=HARC*OMDOT*SN
      UB=1.0-HARC*OMDOT*CS
      VB=-HARC*OMDOT*SN
C
C  Velocities are only set on a velocity-type boundary node
C
      INDER=0
      INDEL=(NY1-1)/4
      YINC=YLEN*0.5/FLOAT(NY1-1)
C
C    OMEGA is the total rotation in radians, so far
C    The algorithm searches through the boundary nodes, but
C    relies on the fact that they are found sequentially in
C    order of increasing Y, for setting the velocity on the
C    wings
C
      DO 10 I=1,NBP
      NODU=IBC(I)
      NODUT=IBCTYP(I)
      NODVT=IBCTYP(I+NBP)
      NORN=NOR(IABS(NODU))
      X=EX(NORN)
      Y=EY(NORN)
      IF((NODU.LE.IWNGR).OR.(NODU.GE.IWNGL).OR.
     1 (X.GT.0.9*XLEN))GO TO 10
C
C   The following is for the right hand wing
C
      IF(NODU.GE.ICSP1)GO TO 20
      INDER=INDER+1
      IF(NODUT.GT.0)GO TO 10
      YPR=2.0*PI*YINC*INDER
      SFAC=SIN(YPR)
      SFAC=SFAC*SFAC
      IF(NODUT.EQ.0)QBND(I)=UA*SFAC
      IF(NODVT.EQ.0)QBND(I+NBP)=VA*SFAC
      IF(INDER.LT.INDEL)GO TO 10
      WRITE(LUW,10001)INDER,NODU,X,Y
      IERR=1
      GO TO 500
C
C    The following is for the central arc
C
   20 IF(NODU.GT.ICSP2)GO TO 30
      RAD=SQRT((X-XMID)*(X-XMID) + (Y-YMID)*(Y-YMID))
      IF ((XMID-X).EQ.0.AND.(Y-YMID).EQ.0) THEN
        PHI=0
      ELSE
        PHI=ATAN2(XMID-X,Y-YMID)
      END IF
      U=1.0-RAD*OMDOT*COS(PHI)
      V=-RAD*OMDOT*SIN(PHI)
      QBND(I)=U
      QBND(I+NBP)=V
      GO TO 10
C
C   The following is for the left hand wing
C
   30 INDEL=INDEL-1
      YPR=2.0*PI*YINC*INDEL
      SFAC=SIN(YPR)
      SFAC=SFAC*SFAC
      IF(NODUT.EQ.0)QBND(I)=UB*SFAC
      IF(NODVT.EQ.0)QBND(I+NBP)=VB*SFAC
      IF(INDEL.GT.0)GO TO 10
      WRITE(LUW,10001)INDEL,NODU,X,Y
      IERR=1
      GO TO 500
C
   10 CONTINUE
C     WRITE(LUW,10059)
C10059 FORMAT(' QBND matrix at end of ROTIND')
C     CALL MATPRT(QBND,NBPP,2,NBP2P,LUW)
 500  RETURN
10001 FORMAT(' ROTIND GOING WRONG, IND =',I5,' NOD =',I5,
     1' X =',G12.5,' Y =',G12.5)
10002 FORMAT(' INADMISSABLE NODE NUMBER (',I5,') IN',
     1' ROUTINE ROTIND')
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
