C*--------------------------------------------------------------------
C*    Basil / Sybil:   strain.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE STRAIE(IQU,IDELE,VALUE,COOL,VC,SE,ARGAN,
     :                  BRGAN,HLENSC,BDEPSC,BIG,GAMM,TREF,
     :                  EX,EY,VHB,THDINT,UVP,SSQ,IELFIX,
     :                  LEM,NOR,NE,NUP,NUVP,NFP,IVV,NCOMP,ICR)
C
C    Routine to obtain interpolated values of the strain or
C    stress tensor (or derived quantities) from the velocity field
C    Note that option numbers are tied to names of physical quantities
C    by definitions in strain.h
C    Similar to STRAIX and STRAIL, STRAIE finds components for the natural
C    coordinates (COOL) of a given element (IDELE)
C    The order in COOL assumes the vertex coords in order 1-3
C
      DIMENSION SSQ(NUP),IELFIX(NUP)
      DIMENSION EX(NUP),EY(NUP),VHB(8,NE),UVP(NUVP),THDINT(7,NE)
      DIMENSION LEM(6,NE),NOR(NUP)
      DIMENSION COOL(3)
C    local arrays
      DIMENSION A0(3),BY(3),CX(3),DNDP(84)
      DIMENSION XX(3),YY(3)
      DIMENSION DLDX(3),DLDY(3)
C     DOUBLE PRECISION TRI,X1,X2,X3,Y1,Y2,Y3
C     DOUBLE PRECISION A0(3),BY(3),CX(3),COOL(3)
C
      NUP2=NUP*2
      EPS=1.0E-4
      GPEZERO=-0.5*ARGAN/(HLENSC*HLENSC)
      UUI=0.0
      VVI=0.0
      SXP=SE
C
C    BIGV should be the same as that calculated in VISK (cginit.f)
C
      BIGV=BIG/10000.0
      IF(IVV.GE.3)SXP=VHB(8,IDELE)
      VF=VC
      TEMDEP=1.0
C
C     Calculate the geometrical coefficients
C
      DO 20 K1=1,3
        K2=MOD(K1,3)+1
        K3=MOD(K1+1,3)+1
        LK1=NOR(LEM(K1,IDELE))
        LK2=NOR(LEM(K2,IDELE))
        LK3=NOR(LEM(K3,IDELE))
        X1=EX(LK1)
        X2=EX(LK2)
        X3=EX(LK3)
        Y1=EY(LK1)
        Y2=EY(LK2)
        Y3=EY(LK3)
C    The order in COOL assumes the vertex coords in order 1-3
        IF (K1.EQ.1) THEN
          XP=COOL(1)*X1 + COOL(2)*X2 + COOL(3)*X3
          YP=COOL(1)*Y1 + COOL(2)*Y2 + COOL(3)*Y3
        END IF
        A0(K1)=X2*Y3 - X3*Y2
        BY(K1)=Y2    - Y3
        CX(K1)=X3    - X2
   20 CONTINUE
C
C    TRI is twice the area of the triangle element
C    spatial gradients of the Li interpolation functions are constant
C
      TRI=A0(1)+A0(2)+A0(3)
      DO KP=1,3
        DLDX(KP)=BY(KP)/TRI
        DLDY(KP)=CX(KP)/TRI
      ENDDO
C
C         Interpolate all strain-rate components
C
      DUDX=0.0
      DUDY=0.0
      DVDX=0.0
      DVDY=0.0
      DO 35 KI=1,6
        LK=LEM(KI,IDELE)
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
   35 CONTINUE
C
C   if required interpolate the velocity components
C
      IF((NCOMP.GE.2).OR.(NCOMP.LE.-1))THEN
        UUI=0.0
        VVI=0.0
        DO KI=1,6
          LK=LEM(KI,IDELE)
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
      IF(IVV.GE.1)CALL BINTIP(VHB(1,IDELE),COOL,VF,1)
      IF(GAMM.NE.0.0)THEN
        CALL BINTIP(THDINT(1,IDELE),COOL,DAMG,0)
        WRKTRM=GAMM*DAMG
        IF(TREF.EQ.0.D0)THEN       ! if GAMMA is constant
          TEMDEP=EXP(-WRKTRM)
          VF=VF*TEMDEP
        ELSE                        ! if GAMMA is temperature dependent
          BLOG=LOG(VF)
          DNUMER=BLOG+(TREF*BLOG-1.0)*WRKTRM
          DENOMR=1.0+TREF*WRKTRM
          VF=EXP(DNUMER/DENOMR)
        ENDIF
      ENDIF
C
C   if required, interpolate the pressure field (NCOMP.GE.1)
C
      IF(IQU.GE.INDXSIXX)THEN
        SIZZ=0.0
        IF(NCOMP.GE.1)THEN
          PRES=0.0
          DO K=1,3
            LK=NUP2+NFP+NOR(LEM(K,IDELE))
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
              LK=LEM(KI,IDELE)
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
      CALL STCOMP(IQU,VALUE,DUDX,DUDY,DVDX,DVDY,PRES,SIZZ,
     :                NCOMP,UUI,VVI,XP,YP,IVV,SXP,VF,BIGV)
      RETURN
      END
      SUBROUTINE BINTIP(VHB,COOL,VINT,ILOG)
C
C    Interpolate a quantity defined on the 7 integration points
C    using a best fit quadratic interpolation function.  ILOG is 
C    a switch that determines if the interpolation is done on the
C    function itself (ILOG=0) or its logarithm (ILOG.ne.0).  The
C    latter ensures that the interpolation will not produce negative
C    values on quantities like viscosity.
C
      DIMENSION VHB(7),COOL(3)
      DIMENSION VNOD(6),ARR(7,6),QOOL(6)
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
C    obtain the best-fit equivalent nodal values
C
      IF(ILOG.NE.0)THEN
        DO KE=1,6
          VNOD(KE)=0.0
          DO KPP=1,7
            VNOD(KE)=VNOD(KE)+ARR(KPP,KE)*ALOG(VHB(KPP))
          ENDDO
        ENDDO
      ELSE
        DO KE=1,6
          VNOD(KE)=0.0
          DO KPP=1,7
            VNOD(KE)=VNOD(KE)+ARR(KPP,KE)*VHB(KPP)
          ENDDO
        ENDDO
      ENDIF
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
      IF(ILOG.NE.0)THEN
        VINT=EXP(SUM)
      ELSE
        VINT=SUM
      ENDIF
C
      RETURN
      END
C
C*******************************************************************************
C
      SUBROUTINE STCOMP(IQU,VALUE,DUDX,DUDY,DVDX,DVDY,PRES,SIZZ,
     :                  NCOMP,UUI,VVI,XLA,YLA,IVV,SXP,VF,BIGV)
      INCLUDE "input.parameters"
      SAVE PI
      DATA PI/3.141592653/
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
        EDZZ=UUI/XLA
      END IF
C
      IF(IQU.EQ.INDXEDXX)THEN
        VALUE=EDXX
        RETURN
      ELSE IF(IQU.EQ.INDXEDYY)THEN
        VALUE=EDYY
        RETURN
      ELSE IF(IQU.EQ.INDXEDZZ)THEN
        VALUE=EDZZ
        RETURN
      ELSE IF(IQU.EQ.INDXEDXY)THEN
        VALUE=EDXY
        RETURN
C
C    vorticity VORT
C
      ELSE IF(IQU.EQ.INDXVORT)THEN
        VALUE=VORT
        RETURN
      END IF   ! end of IF BLOCK A
C
C    Second invariant of the strain-rate tensor
C
      IF(IQU.GE.INDXED2I)THEN
        ED2I=SQRT(EDXX*EDXX + EDYY*EDYY + EDZZ*EDZZ + 2.0*EDXY*EDXY)
        IF(IQU.EQ.INDXED2I)THEN
          VALUE=ED2I
          RETURN
C
C     thermal dissipation
C
        ELSE IF(IQU.EQ.INDXTHDI)THEN
          THDI=0.0
          IF(ED2I.GT.0.0)THDI=VF*ED2I**(1.0/SXP + 1.0)
          VALUE=THDI
          RETURN
        END IF
C
C     viscosity
C
        VISC=VF
        IF(IVV.GE.2)THEN
          VISC=BIGV
          IF(ED2I.GT.0.0)VISC=0.5*VF*ED2I**(1.0/SXP - 1.0)
          IF(VISC.GT.BIGV)VISC=BIGV
        END IF
      END IF    !end of IF BLOCK B
      IF(IQU.EQ.INDXVISC)THEN
C       VALUE=ALOG10(VISC)
        VALUE=VISC
        RETURN
C
C    TAXX, TAYY, TAZZ, TAXY
C
      ELSE IF(IQU.EQ.INDXTAXX)THEN
        VALUE=2.0*EDXX*VISC
        RETURN
      ELSE IF(IQU.EQ.INDXTAYY)THEN
        VALUE=2.0*EDYY*VISC
        RETURN
      ELSE IF(IQU.EQ.INDXTAZZ)THEN
        VALUE=2.0*EDZZ*VISC
        RETURN
      ELSE IF(IQU.EQ.INDXTAXY)THEN
        VALUE=2.0*EDXY*VISC
        RETURN
      END IF    ! end of IF BLOCK C
C
C    Change pressure for thin sheet
C
      IF(NCOMP.LE.0)PRES=SIZZ-2.0*EDZZ*VISC
C
C    PRES, SIXX, SIYY, SIZZ
C
      IF(IQU.EQ.INDXPRES)THEN
        VALUE=PRES
        RETURN
      ELSE IF(IQU.EQ.INDXSIXX)THEN
        VALUE=2.0*EDXX*VISC+PRES
        RETURN
      ELSE IF(IQU.EQ.INDXSIYY)THEN
        VALUE=2.0*EDYY*VISC+PRES
        RETURN
      ELSE IF(IQU.EQ.INDXSIZZ)THEN
        VALUE=2.0*EDZZ*VISC+PRES
        RETURN
      END IF     ! end of IF BLOCK D
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
      IF(IQU.EQ.INDXPSR1)THEN
        VALUE=PSR1
        RETURN
      ELSE IF(IQU.EQ.INDXPSR2)THEN
        VALUE=PSR2
        RETURN
      ELSE IF(IQU.EQ.INDXMSST)THEN
        VALUE=PSRM
        RETURN
      ELSE IF(IQU.EQ.INDXCANG)THEN
        VALUE=CANG
        RETURN
      ELSE IF(IQU.EQ.INDXTANG)THEN
        VALUE=TANG
        RETURN
      ELSE IF(IQU.EQ.INDXSANG)THEN
        VALUE=SANG
        RETURN
C
C    Type of faulting : sum of two double couples :
C       0 to 1/4 : thrust + thrust ; 1/4 to 1/2 : thrust + strike slip ;
C       1/2 to 3/4 : normal + strike slip ; 3/4 to 1 : normal + normal
C       Note:(1/4 to 3/8 thrust dominates strike slip faulting
C             3/8 to 1/2 strike slip dominates thrust faulting
C             1/2 to 5/8 strike slip dominates normal faulting)
C             5/8 to 3/4 normal dominates strike slip faulting
C
      ELSE IF(IQU.EQ.INDXDBLC)THEN
        DBLC=0.5
        IF((PSR1.NE.0.0).OR.(PSR2.NE.0.0))
     :              DBLC=(0.75*PI+ATAN2(PSR2,PSR1))/PI
        VALUE=DBLC
        RETURN
C
C   vorticity to shear stress ratio: VOTA
C
      ELSE IF(IQU.EQ.INDXVOTA)THEN
        VALUE=0.0
        IF(ABS(PSR1-PSR2).GT.1.E-6)VALUE=VORT/(PSR1-PSR2)
        RETURN
      ELSE IF(IQU.EQ.INDXTAU1)THEN
        VALUE=2.0*VISC*PSR1
        RETURN
      ELSE IF(IQU.EQ.INDXTAU2)THEN
        VALUE=2.0*VISC*PSR2
        RETURN
      ELSE IF(IQU.EQ.INDXTAUM)THEN
        VALUE=2.0*VISC*PSRM
        RETURN
      ELSE IF(IQU.EQ.INDXSIG1)THEN
        VALUE=2.0*VISC*PSR1+PRES
        RETURN
      ELSE IF(IQU.EQ.INDXSIG2)THEN
        VALUE=2.0*VISC*PSR2+PRES
        RETURN
      END IF     ! end of IF BLOCK E
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
      IF(IQU.EQ.INDXBRIT)THEN
        IF(APSS.LE.APN.AND.APSS.LE.APT) TI=APSS
        IF(APN.LE.APSS.AND.APN.LE.APT) TI=APN
        IF(APT.LE.APSS.AND.APT.LE.APN) TI=APT
        VALUE=2.0*VISC*TI    ! BRIT
        RETURN
C
C Type of faulting (SS=0.5, N=-0.5, T=1.5)
C    (Contour with level=0.0 and interval=1.0)
C
      ELSE IF(IQU.EQ.INDXBRI2)THEN
        IF(APN.LE.APSS.AND.APN.LE.APT) TI= -0.5
        IF(APT.LE.APSS.AND.APT.LE.APN) TI=1.5
        IF(APSS.LE.APN.AND.APSS.LE.APT) TI=0.5
        VALUE=TI             ! BRI2
        RETURN
C
C  Orientation of the intermediate deviatoric stress
C   Vertical=0.5, PSR1=-0.5, PSR2=1.5
C    (Contour with level=0.0 and interval=1.0)
C    
      ELSE IF(IQU.EQ.INDXFOLT)THEN  !may not be active any more
C       FOLT=0.5
C       IF((EDZZ.GT.PSR1).AND.(EDZZ.GT.PSR2)) FOLT=-0.5
C       IF((EDZZ.LT.PSR1).AND.(EDZZ.LT.PSR2)) FOLT=1.5
        VALUE=0.0
        RETURN
      END IF
C
C   default value for any unrecognised value of IQU
C
      VALUE=0.0
      RETURN
      END
C
