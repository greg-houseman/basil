C*--------------------------------------------------------------------
C*    Basil / Sybil:   strain.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE STRAIX(IQU,IPLOT,NX3,NY3,NP3,SNTRP,IHELP,
     :               BNTRP,CNTRP,DNTRP,XCMIN,YCMIN,XCMAX,YCMAX,
     :               VC,SE,ARGAN,BRGAN,HLENSC,BDEPSC,BIG,EX,EY,
     :               VHB,UVP,SSQ,IELFIX,LEM,NOR,NE,NUP,NUVP,
     :               NFP,IVV,IVIS,NCOMP,ICR,TBXOFF,TBYOFF,VOFFSET,
     :               MESHNUM)
C
C    Routine to obtain interpolated values of the strain or
C    stress tensor (or derived quantities) from the velocity field
C    Note that option numbers are tied to names of physical quantities
C    by definitions in strain.h  Values are returned in array 
C    SNTRP(NX3,NY3), with IHELP(NX3,NY3) indicating location in or out
C    of solution domain.
C    STRAIX is intended to replace STRAIN, using a different more
C    accurate means of interpolating quantities based on the spatial
C    derivatives. STRAIX finds components on a rectangular array of 
C    interpolation points.
C
      DIMENSION SNTRP(NP3),BNTRP(NP3),CNTRP(NP3),DNTRP(NP3),IHELP(NP3)
      DIMENSION SSQ(NUP),IELFIX(NUP)
      DIMENSION EX(NUP),EY(NUP),VHB(8,NE),UVP(NUVP)
      DIMENSION LEM(6,NE),NOR(NUP)
C     DIMENSION A0(3),BY(3),CX(3),DLDX(3),DLDY(3)
      DIMENSION PNI(7,6),PLI(7,3),XX(3),YY(3),COOL(3)
      DOUBLE PRECISION TRI,X1,X2,X3,Y1,Y2,Y3,CNL
      DOUBLE PRECISION A0(3),BY(3),CX(3),DLDX(3),DLDY(3)
C
      NUP2=NUP*2
      EPS=1.0E-4
      CALL LNCORD(PNI,PLI)
      GPEZERO=-0.5*ARGAN/(HLENSC*HLENSC)
      NX2=NX3-1
      NY2=NY3-1
      HSX=(XCMAX-XCMIN)/FLOAT(NX3-3)
      HSY=(YCMAX-YCMIN)/FLOAT(NY3-3)
      EPSX=HSX*0.5
      EPSY=HSY*0.5
      UUI=0.0
      VVI=0.0
      SXP=SE
C
C   Zero the arrays.  IHELP is an index array which records whether a
C   grid point is actually within the finite element mesh.  It contains
C   the relevant element number - or else zero.
C
      DO I=1,NP3
        IHELP(I)=0
        SNTRP(I)=0.0
        BNTRP(I)=0.0
        CNTRP(I)=0.0
        DNTRP(I)=0.0
      ENDDO
C
C    BIGV is defined in STCOMP consistent with usage in VISK (cginit.f)
C
      BIGV=BIG/10000.0
C
C    Look at each element in turn
C
      DO 80 N=1,NE
      IF(IVV.GE.3)SXP=VHB(8,N)
      VF=VC
C
C     Calculate the geometrical coefficients
C
   12 XMIN=999.
      YMIN=999.
      XMAX=-999.
      YMAX=-999.
      DO 20 K1=1,3
        K2=MOD(K1,3)+1
        K3=MOD(K1+1,3)+1
        LK1=NOR(LEM(K1,N))
        LK2=NOR(LEM(K2,N))
        LK3=NOR(LEM(K3,N))
        X1=EX(LK1)
        X2=EX(LK2)
        X3=EX(LK3)
        Y1=EY(LK1)
        Y2=EY(LK2)
        Y3=EY(LK3)
        IF(X1.LT.XMIN)XMIN=X1
        IF(X1.GT.XMAX)XMAX=X1
        IF(Y1.LT.YMIN)YMIN=Y1
        IF(Y1.GT.YMAX)YMAX=Y1
        A0(K1)=X2*Y3 - X3*Y2
        BY(K1)=Y2    - Y3
        CX(K1)=X3    - X2
   20 CONTINUE
C
C    TRI is twice the area of the triangle element
C    spatial gradients of the Li interpolation functions are constant
C
      TRI=X2*Y3-X3*Y2+X3*Y1-X1*Y3+X1*Y2-X2*Y1
      DO KP=1,3
        DLDX(KP)=BY(KP)/TRI
        DLDY(KP)=CX(KP)/TRI
      ENDDO
      XOFF=0.0
      YOFF=0.0
      VOFF=0.0
      DO 75 IDUP=1,MESHNUM
        IF (IDUP.GT.1) THEN
          XOFF=XOFF+TBXOFF
          YOFF=YOFF+TBYOFF
          XMIN = XMIN + XOFF
          XMAX = XMAX + XOFF
          YMIN = YMIN + YOFF
          YMAX = YMAX + YOFF
          VOFF = VOFF + VOFFSET
        END IF
C
C     Now look for the points in the mesh that are inside the tri.
C     define a rectangular patch in SNTRP that contains this element
C      IMIN to IMAX, JMIN to JMAX.  If element is outside defined
C      region covered by grid, go to the next element (80)
C      J=2 is for YCMIN; J=NY2 is for YCMAX, etc
C
      JMAX=INT((YMAX+EPSY-YCMIN)/HSY) +2
      IF (JMAX.LT.2) GO TO 75
      IF (JMAX.GT.NY2) JMAX=NY2
      JMIN=INT((YMIN-EPSY-YCMIN)/HSY) +2
      IF (JMIN.LT.2) JMIN=2
      IF (JMIN.GT.NY2) GO TO 75
      IMAX=INT((XMAX+EPSX-XCMIN)/HSX) +2
      IF (IMAX.LT.2) GO TO 75
      IF (IMAX.GT.NX2) IMAX=NX2
      IMIN=INT((XMIN-EPSX-XCMIN)/HSX) +2
      IF (IMIN.LT.2) IMIN=2
      IF (IMIN.GT.NX2) GO TO 75
C
C    get the coordinates (XP,YP) of point to be interpolated
C
      YP=YCMIN + (JMIN-3)*HSY - YOFF
      DO 70 J=JMIN,JMAX
        YP=YP+HSY
        XP=XCMIN + (IMIN-3)*HSX - XOFF
        DO 50 I=IMIN,IMAX
          XP=XP+HSX
C
C     Calculate the natural coordinates of (XP,YP)
C
          DO 25 K=1,3
            CNL=(A0(K) + XP*BY(K) + YP*CX(K))/TRI
            IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 50
            COOL(K)=CNL
   25     CONTINUE
C
C     If we reach here the point is within the triangle, so
C         Interpolate all strain-rate components
C
          DUDX=0.0
          DUDY=0.0
          DVDX=0.0
          DVDY=0.0
          DO 35 KI=1,6
            LK=LEM(KI,N)
            UK=UVP(LK)
            VK=UVP(LK+NUP)
            IF(KI.LE.3)THEN
              DNDX=(4.0*COOL(KI)-1.0)*DLDX(KI)
              DNDY=(4.0*COOL(KI)-1.0)*DLDY(KI)
            ELSE
              KIF=KI-3
              KIB=MOD(KI+1,3)+1
              DNDX=4.0*(COOL(KIF)*DLDX(KIB)+COOL(KIB)*DLDX(KIF))
              DNDY=4.0*(COOL(KIF)*DLDY(KIB)+COOL(KIB)*DLDY(KIF))
            END IF
            DUDX=DUDX+DNDX*UK
            DUDY=DUDY+DNDY*UK
            DVDX=DVDX+DNDX*VK
            DVDY=DVDY+DNDY*VK
   35     CONTINUE
C
C   if required interpolate the velocity components
C
          IF((NCOMP.GE.2).OR.(NCOMP.LE.-1))THEN
            UUI=0.0
            VVI=0.0
            DO KI=1,6
              LK=LEM(KI,N)
              IF(KI.LE.3)THEN
                QLI=COOL(KI)*(2.0*COOL(KI)-1.0)
              ELSE
                KIF=KI-3
                KIB=MOD(KI+1,3)+1
                QLI=4.0*COOL(KIF)*COOL(KIB)
              END IF
              UUI=UUI+UVP(LK)*QLI
              VVI=VVI+UVP(LK+NUP)*QLI
            END DO
          END IF
C
C   if required, get the viscosity coefficient VF by interpolation
C
          IF(IVIS.GE.1)CALL BINTIP(VHB(1,N),COOL,VF)
C
C   if required, interpolate the pressure field (NCOMP.GE.1)
C
          IF(IQU.GE.25)THEN
            SIZZ=0.0
            IF(NCOMP.GE.1)THEN
              PRES=0.0
              DO K=1,3
                LK=NUP2+NFP+NOR(LEM(K,N))
                PRES=PRES+UVP(LK)*COOL(K)
              END DO
            ELSE
C
C   obtain pressure from SIZZ for thin sheet calculations (NCOMP.LE.0)
C   (from interpolated crustal thickness distribution and Argand number)
C
              IF(ICR.NE.0)THEN
                SFAC=-0.5*ARGAN
                SSQIP=0.0
                ELF=0.0
                DO KI=1,6
                  LK=LEM(KI,N)
                  IF(KI.LE.3)THEN
                    QLI=COOL(KI)*(2.0*COOL(KI)-1.0)
                  ELSE
                    KIF=KI-3
                    KIB=MOD(KI+1,3)+1
                    QLI=4.0*COOL(KIF)*COOL(KIB)
                  END IF
                  SSQIP=SSQIP+SSQ(LK)*QLI
                  ELF=ELF+FLOAT(IELFIX(LK))*QLI
                END DO
                IF(ELF.GE.0.5)SFAC=-0.5*(ARGAN+BRGAN)
                SIZZ=SFAC*EXP(2.0*SSQIP) - GPEZERO
              END IF
            END IF
          END IF  
C
C   move required quantity(s) to interpolation grid for scalar colour plot
C
          IJ=(J-1)*NX3 + I
          IHELP(IJ)=N
          CALL STCOMP(IPLOT,IQU,SVAL,BVAL,CVAL,DVAL,
     :                DUDX,DUDY,DVDX,DVDY,PRES,SIZZ,
     :                NCOMP,UUI,VVI,XP,YP,IVV,SXP,VF,BIGV)
          SNTRP(IJ)=SVAL
C
C   move two scalars and direction to interpolation grids for arrow plots
C
          IF(IPLOT.EQ.4)THEN
            BNTRP(IJ)=BVAL
            CNTRP(IJ)=CVAL
            DNTRP(IJ)=DVAL
          ENDIF
   50   CONTINUE
   70 CONTINUE
   75 CONTINUE
   80 CONTINUE
      IF(IPLOT.NE.5)RETURN
C
C   store ccordinates (x, y) in edges of SNTRP for use in contours
C   This is not required for arrow plots indeed causes a problem there
C   possible PROBLEMS with ad hoc change here: check contours
C   lines commented out were overwriting wrong part of array
C
C     DO I=1,NX3-2
C       SNTRP(NX3*(NY3-2)+I)=XCMIN+FLOAT(I-1)*HSX
      DO I=2,NX3-1
        SNTRP(I)=XCMIN+FLOAT(I-2)*HSX
      ENDDO
C     DO J=1,NY3-2
C       SNTRP(J*NX3-2)=YCMIN+FLOAT(J-1)*HSY
      DO J=2,NY2-1
        SNTRP((J-1)*NX3+1)=YCMIN+FLOAT(J-2)*HSY
      ENDDO
      RETURN
      END
      SUBROUTINE STRAIL(IQU,NL,RLINE,ILINE,XP1,YP1,XP2,YP2,
     1              VC,SE,ARGAN,BRGAN,HLENSC,BDEPSC,BIG,EX,EY,
     2              VHB,UVP,SSQ,IELFIX,LEM,NOR,NE,NUP,NUVP,
     3              NFP,IVV,IVIS,NCOMP,ICR)
C
C    Routine to obtain interpolated values of the strain or
C    stress tensor (or derived quantities) from the velocity field
C    Note that option numbers are tied to names of physical quantities
C    by definitions in strain.h
C    Similar to STRAIX, STRAIL finds components on a trajectory
C    across the solution region
C
      DIMENSION RLINE(NL),ILINE(NL)
      DIMENSION SSQ(NUP),IELFIX(NUP)
      DIMENSION EX(NUP),EY(NUP),VHB(8,NE),UVP(NUVP)
      DIMENSION LEM(6,NE),NOR(NUP)
      DIMENSION A0(3),BY(3),CX(3),DNDP(84),COOL(3)
      DIMENSION PNI(7,6),PLI(7,3),XX(3),YY(3)
      DIMENSION DLDX(3),DLDY(3)
C     DOUBLE PRECISION TRI,X1,X2,X3,Y1,Y2,Y3
C     DOUBLE PRECISION A0(3),BY(3),CX(3),COOL(3)
C
      NUP2=NUP*2
      EPS=1.0E-4
      CALL LNCORD(PNI,PLI)
      GPEZERO=-0.5*ARGAN/(HLENSC*HLENSC)
      UUI=0.0
      VVI=0.0
C
C   Zero the arrays
C   RLINE(NL) contains the intrpolated values on the line
C   ILINE(NL) contains a switch indicating a value in RLINE
C   Line runs between (XP1,YP1) and (XP2,YP2)
C
      DO I=1,NL
        ILINE(I)=0
        RLINE(I)=0.0
      ENDDO
      DXP=(XP2-XP1)/FLOAT(NL-1)
      DYP=(YP2-YP1)/FLOAT(NL-1)
C     EPSX=DXP*0.1
C     EPSY=DYP*0.1
      EPS=0.01*SQRT(DXP*DXP+DYP*DYP)
      SXP=SE
C
C    BIGV is defined in STCOMP consistent with usage in VISK (cginit.f)
C
      BIGV=BIG/10000.0
C
C    Look at each element in turn
C
      DO 80 N=1,NE
      IF(IVV.GE.3)SXP=VHB(8,N)
      VF=VC
C
C     Calculate the geometrical coefficients
C
   12 XMIN=999.
      YMIN=999.
      XMAX=-999.
      YMAX=-999.
      DO 20 K1=1,3
        K2=MOD(K1,3)+1
        K3=MOD(K1+1,3)+1
        LK1=NOR(LEM(K1,N))
        LK2=NOR(LEM(K2,N))
        LK3=NOR(LEM(K3,N))
        X1=EX(LK1)
        X2=EX(LK2)
        X3=EX(LK3)
        Y1=EY(LK1)
        Y2=EY(LK2)
        Y3=EY(LK3)
        XMIN=MIN(XMIN,X1)
        XMAX=MAX(XMAX,X1)
        YMIN=MIN(YMIN,Y1)
        YMAX=MAX(YMAX,Y1)
        A0(K1)=X2*Y3 - X3*Y2
        BY(K1)=Y2    - Y3
        CX(K1)=X3    - X2
   20 CONTINUE
C
C    TRI is twice the area of the triangle element
C    spatial gradients of the Li interpolation functions are constant
C
      TRI=X2*Y3-X3*Y2+X3*Y1-X1*Y3+X1*Y2-X2*Y1
      DO KP=1,3
        DLDX(KP)=BY(KP)/TRI
        DLDY(KP)=CX(KP)/TRI
      ENDDO
C
C     Now look for the points on the line that are inside the tri.
C
      DO 70 J=1,NL
        XP=XP1+DXP*FLOAT(J-1)
        YP=YP1+DYP*FLOAT(J-1)
        IF((YP.GT.YMAX+EPS).OR.(YP.LT.YMIN-EPS).OR.
     1     (XP.GT.XMAX+EPS).OR.(XP.LT.XMIN-EPS))GO TO 70
C
C     Calculate the natural coordinates of (XP,YP)
C
          DO 25 K=1,3
            CNL=(A0(K) + XP*BY(K) + YP*CX(K))/TRI
            IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 70
            COOL(K)=CNL
   25     CONTINUE
C
C     If we reach here the point is within the triangle, so
C         Interpolate all strain-rate components
C
          DUDX=0.0
          DUDY=0.0
          DVDX=0.0
          DVDY=0.0
          DO 35 KI=1,6
            LK=LEM(KI,N)
            UK=UVP(LK)
            VK=UVP(LK+NUP)
            IF(KI.LE.3)THEN
              DNDX=(4.0*COOL(KI)-1.0)*DLDX(KI)
              DNDY=(4.0*COOL(KI)-1.0)*DLDY(KI)
            ELSE
              KIF=KI-3
              KIB=MOD(KI+1,3)+1
              DNDX=4.0*(COOL(KIF)*DLDX(KIB)+COOL(KIB)*DLDX(KIF))
              DNDY=4.0*(COOL(KIF)*DLDY(KIB)+COOL(KIB)*DLDY(KIF))
            END IF
            DUDX=DUDX+DNDX*UK
            DUDY=DUDY+DNDY*UK
            DVDX=DVDX+DNDX*VK
            DVDY=DVDY+DNDY*VK
   35     CONTINUE
C
C   if required interpolate the velocity components
C
          IF((NCOMP.GE.2).OR.(NCOMP.LE.-1))THEN
            UUI=0.0
            VVI=0.0
            DO KI=1,6
              LK=LEM(KI,N)
              IF(KI.LE.3)THEN
                QLI=COOL(KI)*(2.0*COOL(KI)-1.0)
              ELSE
                KIF=KI-3
                KIB=MOD(KI+1,3)+1
                QLI=4.0*COOL(KIF)*COOL(KIB)
              END IF
              UUI=UUI+UVP(LK)*QLI
              VVI=VVI+UVP(LK+NUP)*QLI
            END DO
          END IF
C
C   if required, get the viscosity coefficient by interpolation
C
          IF(IVIS.GE.1)CALL BINTIP(VHB(1,N),COOL,VF)
C
C   if required, interpolate the pressure field (NCOMP.GE.1)
C
          IF(IQU.GE.25)THEN
            SIZZ=0.0
            IF(NCOMP.GE.1)THEN
              PRES=0.0
              DO K=1,3
                LK=NUP2+NFP+NOR(LEM(K,N))
                PRES=PRES+UVP(LK)*COOL(K)
              END DO
            ELSE
C
C   obtain pressure from SIZZ for thin sheet calculations (NCOMP.LE.0)
C   (from interpolated crustal thickness distribution and Argand number)
C
              IF(ICR.NE.0)THEN
                SFAC=-0.5*ARGAN
                SSQIP=0.0
                ELF=0.0
                DO KI=1,6
                  LK=LEM(KI,N)
                  IF(KI.LE.3)THEN
                    QLI=COOL(KI)*(2.0*COOL(KI)-1.0)
                  ELSE
                    KIF=KI-3
                    KIB=MOD(KI+1,3)+1
                    QLI=4.0*COOL(KIF)*COOL(KIB)
                  END IF
                  SSQIP=SSQIP+SSQ(LK)*QLI
                  ELF=ELF+FLOAT(IELFIX(LK))*QLI
                END DO
                IF(ELF.GE.0.5)SFAC=-0.5*(ARGAN+BRGAN)
                SIZZ=SFAC*EXP(2.0*SSQIP) - GPEZERO
              END IF
            END IF
          END IF  
C
C   place the required quantity in the interpolation grid
C
          ILINE(J)=N
          IPLOT=0
          CALL STCOMP(IPLOT,IQU,SVAL,BVAL,TVAL,UVAL,
     :                DUDX,DUDY,DVDX,DVDY,PRES,SIZZ,
     :                NCOMP,UUI,VVI,XP,YP,IVV,SXP,VF,BIGV)
          RLINE(J)=SVAL
   70   CONTINUE
   80 CONTINUE
      RETURN
      END
      SUBROUTINE BINTIP(VHB,COOL,VINT)
C
C    Interpolate a quantity defined on the 7 integration points
C    using a best fit quadratic interpolation function
C
      DIMENSION VHB(7),VNOD(6),ARR(7,6),COOL(3),QOOL(6)
C
C   the following matrix is the 7x6 matrix that projects from the
C   7 interpolation points to a best fit set of 6 nodal values
C
      SAVE ARR
      DATA ARR/1.9743924,0.1434449,0.1434449,-0.4126757,-0.4126757,
     : 0.2563768,-0.6923076,0.1434449,1.9743924,0.1434449,0.2563768,
     :-0.4126757,-0.4126757,-0.6923076,0.1434449,0.1434449,1.9743924,
     :-0.4126757,0.2563768,-0.4126757,-0.6923076,0.09812581,0.1562660,
     :0.09812581,1.0345624,-0.2822974,-0.2822974,0.1775148,0.09812581,
     :0.09812581,0.1562660,-0.2822974,1.0345624,-0.2822974,0.1775148,
     :0.1562660,0.09812581,0.09812581,-0.2822974,-0.2822974,1.0345624,
     :0.1775148/
C
C    Feb 2015: amended by GAH, so that interpolation is done on logarithm
C    of VHB, inverted by EXP operator at end of intepolation - to ensure
C    no negative values produced by interpolation of variable field
C
C    obtain the best-fit equivalent nodal values
C
      DO KE=1,6
        VNOD(KE)=0.0
        DO KPP=1,7
          VNOD(KE)=VNOD(KE)+ARR(KPP,KE)*ALOG(VHB(KPP))
        ENDDO
      ENDDO
C
C     compute the quadratic interpolation function values at
C     the given local coordinates.
C
      DO KE=1,3
        QOOL(KE)=COOL(KE)*(2.0*COOL(KE)-1.0)
      ENDDO
      DO KE=4,6
        KPF=KE-3
        KPB=MOD(KE+1,3)+1
        QOOL(KE)=4.0*COOL(KPF)*COOL(KPB)
      ENDDO
C
C   the interpolated fuction value then is inner product of VNOD and QOOL
C
      SUM=0.0
      DO KE=1,6
        SUM=SUM+QOOL(KE)*VNOD(KE)
      ENDDO
      VINT=EXP(SUM)
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE STCOMP(IPLOT,IQU,SVAL,BVAL,CVAL,DVAL,
     :                  DUDX,DUDY,DVDX,DVDY,PRES,SIZZ,
     :                  NCOMP,UUI,VVI,XLA,YLA,IVV,SXP,VF,BIGV)
      SAVE PI
      DATA PI/3.141592653/
C
C    A maximum viscosity of BIGV is permitted, smoothly tapered in the
C    range BIGV/100 to BIGV
C
C     BIGV=1.0E+20
C     BIGV10=BIGV*0.01
C     BIGINV=1.0/BIGV
C
C    this routine assigns the relevant component of the stress or
C    strain-rate tensor for use by STRAIX
C
      EDXX=DUDX
      EDYY=DVDY
      EDXY=0.5*(DUDY+DVDX)
      VORT=DVDX-DUDY
C
C    Corrections to the strain-rate expressions for the
C     thin spherical sheet
C
      IF(NCOMP.EQ.-1)THEN
        YCOTK=1.0/TAN(YLA)
        EDXX=EDXX+VVI*YCOTK
        EDYY=EDYY+XLA*YCOTK*DVDX
        EDXY=EDXY+0.5*(XLA*DUDX-UUI)*YCOTK
        VORT=VORT-(XLA*DUDX-UUI)*YCOTK
      END IF
      EDZZ=-(EDXX+EDYY)
C
C   strain-rates perpendicular to the plane in cylindrical axisymmetry
C
      IF(NCOMP.EQ.2)THEN
        IF(XLA.GT.0.)EDZZ=UUI/XLA
      END IF
C
C    Second invariant of the strain-rate tensor, thermal dissipation
C
      ED2I=SQRT(EDXX*EDXX + EDYY*EDYY + EDZZ*EDZZ + 2.0*EDXY*EDXY)
      THDI=0.0
      IF(ED2I.GT.0.0)THDI=2.0*VF*ED2I**(1.0/SXP + 1.0)
C
C     viscosity; interpolated VF from BINTIP should always be > 0
C
      VISC=VF
      IF(IVV.GE.2)THEN
        VISC=BIGV
        IF(ED2I.GT.0.0)VISC=VF*ED2I**(1.0/SXP - 1.0)
C
C     a maximum viscosity is permitted but is approached gradually
C
C       IF(VISC.GT.BIGV10)VISC=1.D0/(BIGINV+(1.D0/VISC))
        IF(VISC.GT.BIGV)VISC=BIGV
      END IF
C
C    Change pressure for thin sheet
C
      IF(NCOMP.LE.0)PRES=SIZZ-2.0*EDZZ*VISC
C
C    Diagonalise the strain-rate tensor if required
C    THETA and THETA1 are the angles of the principal axes
C
      IF((EDXY.EQ.0.0).AND.(EDXX-EDYY).EQ.0.0)THEN
        THETA = 0.0
      ELSE
        THETA=0.5*ATAN2(2.0*EDXY,(EDXX-EDYY))
      END IF
      THETA1=THETA+PI*0.5
      TCOS=COS(THETA)
      TCS1=COS(THETA1)
      TSIN=SIN(THETA)
      TSN1=SIN(THETA1)
      S01=EDXX*TCOS*TCOS + EDYY*TSIN*TSIN + 2.0*EDXY*TSIN*TCOS
      S02=EDXX*TCS1*TCS1 + EDYY*TSN1*TSN1 + 2.0*EDXY*TSN1*TCS1
      IF(S01.GE.S02)THEN
        PSR1=S01
        PSR2=S02
        CANG=THETA
        TANG=THETA1
      ELSE
        PSR1=S02
        PSR2=S01
        CANG=THETA1
        TANG=THETA
      END IF
      PSRM=0.5*(PSR1-PSR2)
      SANG=0.5*(CANG+TANG)
C
C    Type of faulting : sum of two double couples :
C       0 to 1/4 : thrust + thrust ; 1/4 to 1/2 : thrust + strike slip ;
C       1/2 to 3/4 : normal + strike slip ; 3/4 to 1 : normal + normal
C       Note:(1/4 to 3/8 thrust dominates strike slip faulting
C             3/8 to 1/2 strike slip dominates thrust faulting
C             1/2 to 5/8 strike slip dominates normal faulting)
C             5/8 to 3/4 normal dominates strike slip faulting
C
        DBLC=0.5
        IF((PSR1.NE.0.0).OR.(PSR2.NE.0.0))
     :              DBLC=(0.75*PI+ATAN2(PSR2,PSR1))/PI
C
C    default values if not set
C
      SVAL=0.0
      BVAL=0.0
      CVAL=0.0
      DVAL=0.0
C
C    following block assigns required values for scalar plot
C    Values for arrow plots follow this block
C
      IF(IPLOT.NE.4)THEN
      IF(IQU.EQ.1)THEN
        SVAL=EDXX
        RETURN
      ELSE IF(IQU.EQ.2)THEN
        SVAL=EDYY
        RETURN
      ELSE IF(IQU.EQ.3)THEN
        SVAL=EDZZ
        RETURN
      ELSE IF(IQU.EQ.4)THEN
        SVAL=EDXY
        RETURN
      ELSE IF(IQU.EQ.5)THEN
        SVAL=PSR1
        RETURN
      ELSE IF(IQU.EQ.6)THEN
        SVAL=PSR2
        RETURN
      ELSE IF(IQU.EQ.7)THEN
        SVAL=PSRM
        RETURN
      ELSE IF(IQU.EQ.8)THEN
        SVAL=CANG
        RETURN
      ELSE IF(IQU.EQ.9)THEN
        SVAL=TANG
        RETURN
      ELSE IF(IQU.EQ.10)THEN
        SVAL=SANG
        RETURN
      ELSE IF(IQU.EQ.11)THEN
        SVAL=DBLC
        RETURN
C
C    vorticity VORT
C
      ELSE IF(IQU.EQ.12)THEN
        SVAL=VORT
        RETURN
      ELSE IF(IQU.EQ.13)THEN
        SVAL=ED2I
        RETURN
C
C     thermal dissipation
C
      ELSE IF(IQU.EQ.14)THEN
        SVAL=THDI
        RETURN
      ELSE IF(IQU.EQ.15)THEN
        SVAL=ALOG10(VISC)
        RETURN
C
C    TAXX, TAYY, TAZZ, TAXY
C
      ELSE IF(IQU.EQ.16)THEN
        SVAL=2.0*EDXX*VISC
        RETURN
      ELSE IF(IQU.EQ.17)THEN
        SVAL=2.0*EDYY*VISC
        RETURN
      ELSE IF(IQU.EQ.18)THEN
        SVAL=2.0*EDZZ*VISC
        RETURN
      ELSE IF(IQU.EQ.19)THEN
        SVAL=2.0*EDXY*VISC
        RETURN
      ELSE IF(IQU.EQ.20)THEN
        SVAL=2.0*VISC*PSR1
        RETURN
      ELSE IF(IQU.EQ.21)THEN
        SVAL=2.0*VISC*PSR2
        RETURN
      ELSE IF(IQU.EQ.22)THEN
        SVAL=2.0*VISC*PSRM
        RETURN
C
C   vorticity to shear stress ratio: VOTA
C
      ELSE IF(IQU.EQ.24)THEN
        SVAL=0.0
        IF(ABS(PSR1-PSR2).GT.1.E-6)RVAL=VORT/(PSR1-PSR2)
        RETURN
C
C    PRES, SIXX, SIYY, SIZZ
C
      ELSE IF(IQU.EQ.25)THEN
        SVAL=2.0*EDXX*VISC+PRES
        RETURN
      ELSE IF(IQU.EQ.26)THEN
        SVAL=2.0*EDYY*VISC+PRES
        RETURN
      ELSE IF(IQU.EQ.27)THEN
        SVAL=2.0*EDZZ*VISC+PRES
        RETURN
      ELSE IF(IQU.EQ.28)THEN
        SVAL=2.0*VISC*PSR1+PRES
        RETURN
      ELSE IF(IQU.EQ.29)THEN
        SVAL=2.0*VISC*PSR2+PRES
        RETURN
      ELSE IF(IQU.EQ.30)THEN
        SVAL=PRES
        RETURN
      END IF
C
C    Not clear if options 31 to 33 are still valid and useful
C
C    Type of brittle failure: AMU is the coefficient of friction
C    APSS for strike-slip, APN for normal, APT for thrust, latter
C    two only relevant for thin sheet
C
      AMU=0.85
      SAVG=0.5*(PSR1+PSR2)
      APSS= -PSRM*SQRT(1.0+1.0/AMU/AMU)-SAVG
      TMAX=0.5*(PSR1-EDZZ)
      SAVG=0.5*(PSR1+EDZZ)
      APN= -TMAX*SQRT(1.0+1.0/AMU/AMU)-SAVG
      TMAX=0.5*(EDZZ-PSR2)
      SAVG=0.5*(EDZZ+PSR2)
      APT= -TMAX*SQRT(1.0+1.0/AMU/AMU)-SAVG
C
C     Take the minimum of APSS, APN and APT
C
      IF(IQU.EQ.31)THEN
        IF(APSS.LE.APN.AND.APSS.LE.APT) TI=APSS
        IF(APN.LE.APSS.AND.APN.LE.APT) TI=APN
        IF(APT.LE.APSS.AND.APT.LE.APN) TI=APT
        SVAL=2.0*VISC*TI    ! BRIT
        RETURN
C
C Type of faulting (SS=0.5, N=-0.5, T=1.5)
C    (Contour with level=0.0 and interval=1.0)
C
      ELSE IF(IQU.EQ.32)THEN
        IF(APN.LE.APSS.AND.APN.LE.APT) TI= -0.5
        IF(APT.LE.APSS.AND.APT.LE.APN) TI=1.5
        IF(APSS.LE.APN.AND.APSS.LE.APT) TI=0.5
        SVAL=TI             ! BRI2
        RETURN
C
C  Orientation of the intermediate deviatoric stress
C   Vertical=0.5, PSR1=-0.5, PSR2=1.5
C    (Contour with level=0.0 and interval=1.0)
C    
      ELSE IF(IQU.EQ.33)THEN  !may not be active any more
C       FOLT=0.5
C       IF((EDZZ.GT.PSR1).AND.(EDZZ.GT.PSR2)) FOLT=-0.5
C       IF((EDZZ.LT.PSR1).AND.(EDZZ.LT.PSR2)) FOLT=1.5
        SVAL=0.0
        RETURN
      END IF
      RETURN
      END IF       ! end of IPLOT.ne.4 condition
C
C    values required for arrow plots (IPLOT=4) follow
C    IQU values have distinct interpretation for IPLOT = 4
C
C    Amplitude and direction of principal strain axes
C
      IF(IQU.EQ.1) THEN
        SVAL=PSR1
        BVAL=PSR2
        CVAL=COS(CANG)
        DVAL=SIN(CANG)
C
C     Direction only of principal strain axes
C
      ELSE IF(IQU.EQ.2) THEN
        SVAL=1.0
        BVAL=-1.0
        CVAL=COS(CANG)
        DVAL=SIN(CANG)
C
C     Amplitude and direction of principal deviatoric stresses
C
      ELSE IF(IQU.EQ.3) THEN
        SVAL=2.0*VISC*PSR1
        BVAL=2.0*VISC*PSR2
        CVAL=COS(CANG)
        DVAL=SIN(CANG)
C
C     amplitude and direction of maximum shear strain rate
C
      ELSE IF(IQU.EQ.4) THEN
        SVAL=PSRM
        BVAL=-PSRM
        CVAL=COS(CANG + PI*0.25)
        DVAL=SIN(CANG + PI*0.25)
C
C     Direction only of maximum shear stresses
C
      ELSE IF(IQU.EQ.5) THEN
        SVAL=1.0
        BVAL=-1.0
        CVAL=COS(CANG + PI*0.25)
        DVAL=SIN(CANG + PI*0.25)
C
C     Direction and magnitude of maximum shear stress
C
      ELSE IF(IQU.EQ.22) THEN
        SVAL=2.0*VISC*PSRM
        BVAL=-2.0*VISC*PSRM
        CVAL=COS(CANG + PI*0.25)
        DVAL=SIN(CANG + PI*0.25)
C
C     Direction of likely strike-slip faulting 
C
      ELSE IF(IQU.EQ.7) THEN
        SVAL=1.0
        BVAL=0.0
        ARRANG=CANG+PI/3.0
        CAT=COS(ARRANG)
        SAT=SIN(ARRANG)
        ROTRA1=(DVDY-DUDX)*SAT*CAT + DVDX*CAT*CAT - DUDY*SAT*SAT
        ARRANG=CANG+PI*2.0/3.0
        CAT=COS(ARRANG)
        SAT=SIN(ARRANG)
        ROTRA2=(DVDY-DUDX)*SAT*CAT + DVDX*CAT*CAT - DUDY*SAT*SAT
        IF(ABS(ROTRA1).LE.ABS(ROTRA2))THEN
          CVAL=COS(CANG+PI/3.0)
          SVAL=SIN(CANG+PI/3.0)
        ELSE
          CVAL=COS(CANG+PI*2.0/3.0)
          SVAL=SIN(CANG+PI*2.0/3.0)
        END IF
      END IF
      RETURN
      END
C
      SUBROUTINE PAVG(SUM,EX,EY,UVP,NN,NUP,NFP,NUVP)
C
C  Routine to calculate the mean of the pressure if the pressure
C  is indeterminate by a constant amount in absence of normal
C  stress condition on some part of the boundary
C
      DIMENSION EX(NUP),EY(NUP),UVP(NUVP)
C     INCLUDE 'cg.parameters'
C     COMMON/B1/EX1(NUPP)
C     COMMON/C1/EY1(NUPP)
C     COMMON/AD/UVP(NROWSP)
C     COMMON/CONTRO/LUW,LRBAT,LRINT,LDIN,LWINT,LPLOT,IAUTO,ILABEL,ICOLO
      SUM=0.0
      ITOT=0
      DO 10 I=1,NN
        AX=EX(I)-1.0
        AY=EY(I)-1.0
        RR=SQRT(AX*AX+AY*AY)
        NI=I+NUP+NUP+NFP
        IF(RR.GT.0.05.AND.RR.LT.0.5) THEN
           SUM=SUM+UVP(NI)
           ITOT=ITOT+1
        END IF
   10 CONTINUE
      IF (ITOT.NE.0) SUM=SUM/ITOT
C      WRITE(LSC,10001)SUM
C10001 FORMAT("The mean pressure is ",G12.5)
      RETURN
      END
C
      SUBROUTINE  MATPRT(Z,NX,NY,NK,LUW)
C 
C    MATPRT and MITPRT are designed for outputting matrix data in
C    tabulated form, either for purpose of debugging, or for verbose
C    mode output. They are respectively for real and integer variables
C
      DIMENSION Z(NK)
      DIMENSION L(10)
C
C            FIND NUMBER OF BLOCKS
C
      JB=10
      NB=NX/JB+1
      MM=MOD(NX,JB)
      IF (MM.EQ.0) NB=NB-1
C
C            CYCLE OVER BLOCKS
C
      DO I=1,NB
        JA=(I-1)*JB + 1
        JC=I*JB
        IF (I.EQ.NB) JC=NX
        JJ=JC-JA+1
C
        DO J=1,JJ
          L(J)=JA+J-1
        ENDDO
        WRITE (LUW,101) (L(J),J=1,JJ)
  101   FORMAT (7X,10('J=',I5,3X))
C
C            CYCLE OVER ROWS
C
        DO K1=1,NY
          WRITE (LUW,102) K1, (Z((K1-1)*NX+J),J=JA,JC)
        ENDDO
  102   FORMAT (1X,'K=',I4,1X,10G11.5)
C
        WRITE(LUW,103)
      ENDDO
  103 FORMAT(/)
C
      RETURN
      END
C
      SUBROUTINE  MITPRT(JZ,NX,NY,NK,LUW)
C
C    MATPRT and MITPRT are designed for outputting matrix data in
C    tabulated form, either for purpose of debugging, or for verbose
C    mode output. They are respectively for real and integer variables
C
      DIMENSION JZ(NK)
      DIMENSION L(10)
C
C            FIND NUMBER OF BLOCKS
C
      JB=10
      NB=NX/JB+1
      MM=MOD(NX,JB)
      IF (MM.EQ.0) NB=NB-1
C
C            CYCLE OVER BLOCKS
C
      DO I=1,NB
        JA=(I-1)*JB + 1
        JC=I*JB
        IF (I.EQ.NB) JC=NX
        JJ=JC-JA+1
C
        DO J=1,JJ
          L(J)=JA+J-1
        ENDDO
        WRITE (LUW,101) (L(J),J=1,JJ)
  101   FORMAT (10X,10('J=',I5,3X))
C
C            CYCLE OVER ROWS
C
        DO K1=1,NY
          WRITE (LUW,102) K1, (JZ((K1-1)*NX+J),J=JA,JC)
        ENDDO
  102   FORMAT (1X,'K=',I4,1X,10I11)
C
        WRITE(LUW,103)
      ENDDO
  103 FORMAT(/)
C
      RETURN
      END
