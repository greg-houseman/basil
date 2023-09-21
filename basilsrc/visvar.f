C*--------------------------------------------------------------------
C*    Basil / Sybil:   visvar.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE VISVAR(LTMP,SE,HLENSC,BDEPSC,REFLEV,RISOST,EX,EY,
     :                  LEM,NOR,VHB,SSQ,IMAT,DENS,VC,CENTLNG,RPT2,
     :                  IDEFTYP,IVIS,IDEN,IVV,IMREG,NUP,NE,LUW,LSC,
     :                  IDBUG,IERR)
C
C     This routine specifies the relative lithospheric viscosity
C     coefficient VC in the array entries VHB(1-7,NEP), keyed to the 
C     integration points of the triangular elements.
C
C     Setting and use of viscosity parameters should be consistent with
C     definition of IVV as follows
C
C       * ivv    visc (vhb[i,1-7])   se (vhb[i,8])
C       *  0       constant           constant(1)    No vhb array
C       *  1       variable           constant(1)
C       *  2       constant           constant(>1)   
C       *  3       variable           constant(>1)
C       *  4       variable           variable
C
      INCLUDE "limits.parameters"
      COMMON/FILES/NAMDUM,LGFILE,PFILE
      COMMON/V01S/NAMEW,NDATE,NAMER,COMMEN
      PARAMETER (MAXPTS=80)
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER(16) NAMDUM,NAMER,NAMEW,NDATE
      CHARACTER(80) LGFILE,PFILE,COMMEN
      CHARACTER(1) IPE
      DOUBLE PRECISION ROTIX(3,3)
      DOUBLE PRECISION PLON8,PROT8
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION SSQ(NUP)
      DIMENSION VHB(8,NE)
      DIMENSION IMAT(NE)
      DIMENSION DENS(7,NE)
      DIMENSION IVX(100),POLYG(2,MAXPTS)
      DIMENSION RXVC(10),RXRO(10),QIFN(6,7)
      DIMENSION RPT2(4)
C
C  the array QIFN is used to convert from data on nodes to data on
C  interpolation points
C
      DATA QIFN/
     1   0.474352628,-0.0807685927,-0.0807685927, 0.323074371,
     ;   0.323074371, 0.0410358272,-0.0807685927, 0.474352628,
     :  -0.0807685927, 0.0410358272, 0.323074371, 0.323074371,
     :  -0.0807685927,-0.0807685927, 0.474352628, 0.323074371,
     :   0.0410358272,  0.323074371, -0.0280749425,-0.0525838993,
     :  -0.0280749425, 0.884134233, 0.11229977, 0.11229977,
     :  -0.0280749425,-0.0280749425,-0.0525838993, 0.11229977,
     ;   0.884134233, 0.11229977, -0.0525838993,-0.0280749425,
     :  -0.0280749425, 0.11229977, 0.11229977, 0.884134233,
     :  -0.111111119,-0.111111119,-0.111111119, 0.444444478,
     :   0.444444478, 0.444444478/
C
C    if spherical shell coordinate rotation is invoked work
C    out rotation operators to apply to coordinates
C
      IF(IDEFTYP.EQ.111)THEN
        CALL GETPOLE(RPT2,PLON8,PROT8,ROTIX) ! find equatorial pole for rotation
        PTALO=RPT2(1)
        PTALA=RPT2(2)
        CALL ROTATE(ROTIX,PTALO,PTALA,1,0.0) ! forward rotation to get DPHI
        DPHI=PTALO
        CALL MATSET(PLON8,-PROT8,ROTIX)    ! modify ROTIX for inverse rotation
      ENDIF
C
C     In method E, the anomalous regions are defined by a list of element numbers
C
C     In method P, we define a polygon(s) by a list of vertices in anti-clockwise
C      order.  POLYG contains the actual (x,y) coordinates of the vertices
C
C     Read special region definitions, allowing for continuation lines
C
      EPS=1.E-4
  460 CONTINUE
      CALL READLINE(INSTR,LENGTH,LTMP,LSC,IEND)
      IF (IEND.EQ.1) GO TO 700
      CALL SKIP(INSTR,LENGTH,1,'n',K1)
      IF(INSTR(K1:K1+2).NE.'REG')GO TO 460
      CALL PARSEV(INSTR,LENGTH,IPE,NPTS,POLYG,MAXPTS,LTMP,
     :            IVX,IVRULE,IDRULE,IMRULE,ISERULE,RXSE,RXVC,RXRO,
     :            IVCCNT,IROCNT,IP,IERR)
C     WRITE(*,*)'PARSEV: IPE,IP,IVRULE,IDRULE,IMRULE,NPTS,RXVC =',
C    :IPE,IP,IVRULE,IDRULE,IMRULE,NPTS,RXVC
      IF(IP.EQ.0)GO TO 460
C
C   check for unphysical parameters that would cause computational problems
C   and check IVV set to appropriate value
C
      IF(RXSE.LE.EPS)THEN
        RXSE=1.0
        WRITE(LUW,*)'Warning RXSE <= 0 is unphysical: reset to 1.0'
        WRITE(LSC,*)'Warning RXSE <= 0 is unphysical: reset to 1.0'
      END IF
      IF(((IPE.EQ.'E').OR.(IVRULE.EQ.0)).AND.(RXVC(1).LE.0.0))THEN
        WRITE(LUW,*)'Warning RXVC <= 0 is unphysical: reset to 1.0'
        WRITE(LSC,*)'Warning RXVC <= 0 is unphysical: reset to 1.0'
        IERR=1
        RETURN
      END IF
      IF(IVV.LT.1)IVV=1
      IF((ABS(RXSE-1.0).GT.EPS).AND.(IVV.LT.3)) IVV=3
      IF((IPE.EQ.'P').AND.(NPTS.LE.2))THEN
        WRITE(LUW,*)'Warning REG P statement requires > 2 points'
        WRITE(LSC,*)'Warning REG P statement requires > 2 points'
        IERR=1
        RETURN
      ENDIF
C
C    Method 1. - use list of element numbers from REG E {...} statement
C
      IF((IPE.EQ.'E').AND.(NPTS.GT.0))THEN
        DO 110 I=1,NPTS
          LNO=IVX(I)
          IF(LNO.GT.NE)THEN
            WRITE(LSC,*)'Warning (VISVAR): Element no.',LNO,' > NE'
            WRITE(LUW,*)'Warning (VISVAR): Element no.',LNO,' > NE'
          ELSE IF((LNO.GT.0).AND.(LNO.LE.NE))THEN
            IF(RXVC(1).GT.0.0)THEN
              IF(IVV.LT.1)IVV=1
              DO K=1,7
                VHB(K,LNO)=RXVC(1)
              ENDDO
            END IF
            IF((ISERULE.GE.0).AND.(IVIS.GT.0))THEN
              VHB(8,LNO)=RXSE
            END IF
          END IF
  110   CONTINUE
      END IF
C
C    VR = 10 maps the viscosity from a poly file triangulation
C
      IF(IPE.EQ.'A'.AND.IVRULE.EQ.10)THEN
        CALL VCFROMPOLY(VHB,EX,EY,LEM,NOR,IVV,NUP,NE,
     :                  PFILE,COMMEN,IERR)
C
C    Method 2. - use corners of polygon from REG P {...} statement
C                note corners in anticlockwise order
C
      ELSE IF(IPE.EQ.'P'.OR.(IPE.EQ.'A'.AND.IVRULE.NE.10))THEN
        IF (IPE.EQ.'P') THEN
C
C    get min and max of X and Y in this polygon
C
          XMIN=POLYG(1,1)
          YMIN=POLYG(2,1)
          XMAX=POLYG(1,1)
          YMAX=POLYG(2,1)
          DO 25 J=2,NPTS
            IF(POLYG(1,J).LT.XMIN)XMIN=POLYG(1,J)
            IF(POLYG(1,J).GT.XMAX)XMAX=POLYG(1,J)
            IF(POLYG(2,J).LT.YMIN)YMIN=POLYG(2,J)
            IF(POLYG(2,J).GT.YMAX)YMAX=POLYG(2,J)
   25     CONTINUE
        END IF
C
C  values required for IVRULE = 11
C
        IF(IVRULE.EQ.11)THEN
          RFL=REFLEV*0.001
          VR10B=-(RXVC(1)-RXVC(2))/(RXVC(3)-RXVC(4))
          VR10A=RXVC(1) + VR10B*RXVC(3)
        END IF
C
C    for each element
C
        DO 40 IJ=1,NE
C
C     find centroid coordinates of element
C
          NCROSSINGS=0
          XP=0.0
          YP=0.0
          DO JJ=1,3
            XP=XP+EX(NOR(LEM(JJ,IJ)))
            YP=YP+EY(NOR(LEM(JJ,IJ)))
          ENDDO
          XP=XP/3.0
          YP=YP/3.0
C
C    if spherical projection / coordinate rotation invoked
C    XP, YP rescaled to long, lat in degrees
C
          IF(IDEFTYP.GE.110)THEN
            PHI0=CENTLNG
            IF(IDEFTYP.EQ.111)PHI0=DPHI
            XP=PHI0+XP/(DTOR*COS(YP-PION2))
            YP=(YP-PION2)/DTOR
            IF(IDEFTYP.EQ.111)THEN
              CALL ROTATE(ROTIX,XP,YP,1,0.0)
            ENDIF
          ENDIF
C
C   if centroid of element is not within polygon, go to next element
C
          IF(IPE.EQ.'P')THEN
            IF((XP.GT.XMAX).OR.(XP.LT.XMIN).OR.(YP.GT.YMAX).OR.
     1         (YP.LT.YMIN))GO TO 40
C
C    check the no. of polygon edge crossings made by a ray along
C    the +X axis. See Graphics Gems IV, 1.4, Ed. Paul S. Heckbert, 1994
            XP1=POLYG(1,NPTS)
            YP1=POLYG(2,NPTS)
            XP2=POLYG(1,1)
            YP2=POLYG(2,1)
            IYFLG0=0
            IF (YP1.GE.YP) IYFLG0=1
            DO 35 J=1,NPTS
C    check if segment endpts are on opposite sides of (XP,YP)
              IF (YP2.GE.YP) THEN
                IYFLG1=1
              ELSE
                IYFLG1=0
              END IF
              IF (IYFLG0.NE.IYFLG1) THEN
                IF (XP1.GE.XP) THEN
                  IXFLG0=1
                ELSE
                  IXFLG0=0
                END IF
C    check if segment endpts are on the same side of y=YP
                IF ((IXFLG0.EQ.0.AND.XP2.LT.XP).OR.
     1              (IXFLG0.EQ.1.AND.XP2.GE.XP)) THEN
                  IF (IXFLG0.EQ.1) THEN
                    IF (IYFLG0.EQ.1) NCROSSINGS = NCROSSINGS-1
                    IF (IYFLG0.EQ.0) NCROSSINGS = NCROSSINGS+1
                  END IF
                ELSE
C    compute intersection of segment and +X ray
C    if >= XP then ray hits
                  IF ((XP2-(YP2-YP)*(XP1-XP2)/(YP1-YP2)).GE.XP)
     1                                                    THEN
                    IF (IYFLG0.EQ.1) NCROSSINGS = NCROSSINGS-1
                    IF (IYFLG0.EQ.0) NCROSSINGS = NCROSSINGS+1
                  END IF
                END IF
              END IF
              IYFLG0 = IYFLG1
              XP1 = XP2
              YP1 = YP2
              IF (J.LE.NPTS) THEN
                XP2=POLYG(1,J+1)
                YP2=POLYG(2,J+1)
              END IF
   35       CONTINUE
          END IF
C
C    element lies within polygon - set viscosity parameters
C    VR = 0 gives a constant viscosity distribution
C
          IF((IPE.EQ.'A').OR.(NCROSSINGS.NE.0))THEN
            IF((ISERULE.GE.0).AND.(IVIS.GT.0))THEN
              VHB(8,IJ)=RXSE
            END IF

            IF(IVRULE.EQ.0)THEN
              IF(RXVC(1).GT.0.0)THEN
                IF(IVV.LT.1)IVV=1
                DO KK=1,7
                  VHB(KK,IJ)=RXVC(1)
                ENDDO
              END IF
C
C    Set parameters which depend on depth if necessary
C
C    for vertex node,   (X,Y)=0.69614048(XV,YV)+0.30385953(XC,YC)
C    for midpoint node, (X,Y)=0.82085238(XM,YM)+0.17914761(XC,YC)
C    (obtained from manipulations of data on p421 of Huebner)
C
C   set viscosity values if applicable
C   VR = 1 gives linear variation with depth
C
            ELSE IF(IVRULE.EQ.1)THEN
              XVRT=0.30385953*XP
              XMID=0.17914761*XP
              YVRT=0.30385953*YP
              YMID=0.17914761*YP
              DO KK=1,3
                NZ=NOR(LEM(KK,IJ))
                XX=0.69614048*EX(NZ)+XVRT
                YY=0.69614048*EY(NZ)+YVRT
                VHB(KK,IJ)=RXVC(1) + RXVC(2)*XX + RXVC(3)*YY
              ENDDO
              DO KK=4,6
                NZ=NOR(LEM(KK,IJ))
                XX=0.82085238*EX(NZ)+XMID
                YY=0.82085238*EY(NZ)+YMID
                VHB(KK,IJ)=RXVC(1) + RXVC(2)*XX + RXVC(3)*YY
              ENDDO
              VHB(7,IJ)=RXVC(1) + RXVC(2)*XP + RXVC(3)*YP
C
C   VR = 2 gives exponential variation with depth
C
            ELSE IF(IVRULE.EQ.2)THEN
              XVRT=0.30385953*XP
              XMID=0.17914761*XP
              YVRT=0.30385953*YP
              YMID=0.17914761*YP
              DO KK=1,3
                NZ=NOR(LEM(KK,IJ))
                XX=0.69614048*EX(NZ)+XVRT
                YY=0.69614048*EY(NZ)+YVRT
                VHB(KK,IJ)=RXVC(1)*EXP(RXVC(2)*XX+RXVC(3)*YY)
              ENDDO
              DO KK=4,6
                NZ=NOR(LEM(KK,IJ))
                XX=0.82085238*EX(NZ)+XMID
                YY=0.82085238*EY(NZ)+YMID
                VHB(KK,IJ)=RXVC(1)*EXP(RXVC(2)*XX+RXVC(3)*YY)
              ENDDO
              VHB(7,IJ)=RXVC(1)*EXP(RXVC(2)*XP+RXVC(3)*YP)
C
C    VR = 3 scales the current viscosity by a factor RXVC(1)
C
            ELSE IF(IVRULE.EQ.3)THEN
              IF(RXVC(1).GT.0.0)THEN
                DO KK=1,7
                  VHB(KK,IJ)=VHB(KK,IJ)*RXVC(1)
                ENDDO
              END IF
C
C    VR = 5 produces a smooth variation from V1 at X1/Y1 to V2 at X2/Y2
C           VC contains: V1, X1, V2, X2, SW; SW = 0 for X, /= 0 for Y
C
            ELSE IF(IVRULE.EQ.5)THEN
              XVRT=0.30385953*XP
              XMID=0.17914761*XP
              YVRT=0.30385953*YP
              YMID=0.17914761*YP
              DO KK=1,3
                NZ=NOR(LEM(KK,IJ))
                XX=0.69614048*EX(NZ)+XVRT
                YY=0.69614048*EY(NZ)+YVRT
                IF(RXVC(5).EQ.0)THEN
                  WW=0.5*(1.0+COS(PI*(XX-RXVC(2))/(RXVC(4)-RXVC(2))))
                ELSE
                  WW=0.5*(1.0+COS(PI*(YY-RXVC(2))/(RXVC(4)-RXVC(2))))
                ENDIF
                VHB(KK,IJ)=RXVC(1)*WW+RXVC(3)*(1.0-WW)
              ENDDO
              DO KK=4,6
                NZ=NOR(LEM(KK,IJ))
                XX=0.82085238*EX(NZ)+XMID
                YY=0.82085238*EY(NZ)+YMID
                IF(RXVC(5).EQ.0)THEN
                  WW=0.5*(1.0+COS(PI*(XX-RXVC(2))/(RXVC(4)-RXVC(2))))
                ELSE
                  WW=0.5*(1.0+COS(PI*(YY-RXVC(2))/(RXVC(4)-RXVC(2))))
                ENDIF
                VHB(KK,IJ)=RXVC(1)*WW+RXVC(3)*(1.0-WW)
              ENDDO
              IF(RXVC(5).EQ.0)THEN
                WW=0.5*(1.0+COS(PI*(XP-RXVC(2))/(RXVC(4)-RXVC(2))))
              ELSE
                WW=0.5*(1.0+COS(PI*(YP-RXVC(2))/(RXVC(4)-RXVC(2))))
              ENDIF
              VHB(7,IJ)=RXVC(1)*WW+RXVC(3)*(1.0-WW)
C
C    VR = 6 is similar to VR = 5, but it enables an oblique band of
C    varying VC value, defined by 4 points that define a quadrilateral,
C    producing a smooth variation from V1 at (C1,D1) to V2 at (C2,D2)
C    the same gradation is imposed from (C3,D3) to (C4,D4).
C    Interpolation is linear in C and cos^2 in D; points 1 and 2
C    are presumed on the same C coordinate, as are points 3 and 4
C    and VC contains: V1, V2, C1, D1, D2, C3, D3, D4, SW
C    if SW!=0 or ==0  (C,D) is interpreted as (X,Y) or (Y,X)
C
            ELSE IF(IVRULE.EQ.6)THEN
              XVRT=0.30385953*XP
              XMID=0.17914761*XP
              YVRT=0.30385953*YP
              YMID=0.17914761*YP
              SW=RXVC(9)
              IF(SW.NE.0.0)THEN
                CC1=RXVC(3)
                DD1=RXVC(4)
                CC2=RXVC(3)
                DD2=RXVC(5)
                CC3=RXVC(6)
                DD3=RXVC(7)
               CC4=RXVC(6)
                DD4=RXVC(8)
              ELSE
                DD1=RXVC(3)
                CC1=RXVC(4)
                DD2=RXVC(3)
                CC2=RXVC(5)
                DD3=RXVC(6)
                CC3=RXVC(7)
                DD4=RXVC(6)
                CC4=RXVC(8)
              ENDIF
              DO KK=1,7
                NZ=NOR(LEM(KK,IJ))
                IF(KK.LE.3)THEN
                  XX=0.69614048*EX(NZ)+XVRT    !vertices
                  YY=0.69614048*EY(NZ)+YVRT
                ELSE IF(KK.EQ.7)THEN
                  XX=XP                        !centroid
                  YY=YP
                ELSE
                  XX=0.82085238*EX(NZ)+XMID    !midpoints
                  YY=0.82085238*EY(NZ)+YMID
                ENDIF
                IF(SW.EQ.0)THEN
                  DTP1=DD1+(YY-CC1)*(DD3-DD1)/(CC3-CC1)
                  DTP2=DD2+(YY-CC2)*(DD4-DD2)/(CC4-CC2)
                  WW=0.5*(1.0+COS(PI*(XX-DTP1)/(DTP2-DTP1)))
                ELSE
                  DTP1=DD1+(XX-CC1)*(DD3-DD1)/(CC3-CC1)
                  DTP2=DD2+(XX-CC2)*(DD4-DD2)/(CC4-CC2)
                  WW=0.5*(1.0+COS(PI*(YY-DTP1)/(DTP2-DTP1)))
                ENDIF
                VHB(KK,IJ)=RXVC(1)*WW+RXVC(2)*(1.0-WW)
              ENDDO
C
C    VR=11 set the strength coefficient using the topographic value in SSQ array
C     never fully tested this option
C
            ELSE IF(IVRULE.EQ.11)THEN
              DO KK=1,7
                EDPSUM=0.0
                DO KX=1,6
                  EDPSUM=EDPSUM+QIFN(KX,KK)*SSQ(LEM(KX,IJ))
                ENDDO
                EDP=RFL+BDEPSC*(HLENSC*EXP(EDPSUM)-1.0)/RISOST
                VHBTMP=VR10A+VR10B*EDP
                IF(EDP.LT.-RXVC(3))VHB(KK,IJ)=VHBTMP
              ENDDO
            END IF
C
C    DR = 0 gives a constant density distribution
C
            IF(IDRULE.EQ.0)THEN
              DO KK=1,7
                DENS(KK,IJ)=RXRO(1)
              ENDDO
C
C    Set parameters which depend on depth if necessary
C    DR = 1 gives linear dependence on coordinate (X or Y)
C
            ELSE IF(IDRULE.EQ.1)THEN
              XVRT=0.30385953*XP
              XMID=0.17914761*XP
              YVRT=0.30385953*YP
              YMID=0.17914761*YP
              DO KK=1,3
                NZ=NOR(LEM(KK,IJ))
                XX=0.69614048*EX(NZ)+XVRT
                YY=0.69614048*EY(NZ)+YVRT
                DENS(KK,IJ)=RXRO(1) + RXRO(2)*XX + RXRO(3)*YY
              ENDDO
              DO KK=4,6
                NZ=NOR(LEM(KK,IJ))
                XX=0.82085238*EX(NZ)+XMID
                YY=0.82085238*EY(NZ)+YMID
                DENS(KK,IJ)=RXRO(1) + RXRO(2)*XX + RXRO(3)*YY
              ENDDO
              DENS(7,IJ)=RXRO(1) + RXRO(2)*XP + RXRO(3)*YP
C
C   DR = 2 gives exponential dependence on coordinate (X or Y)
C
            ELSE IF(IDRULE.EQ.2)THEN
              XVRT=0.30385953*XP
              XMID=0.17914761*XP
              YVRT=0.30385953*YP
              YMID=0.17914761*YP
              DO KK=1,3
                NZ=NOR(LEM(KK,IJ))
                XX=0.69614048*EX(NZ)+XVRT
                YY=0.69614048*EY(NZ)+YVRT
                DENS(KK,IJ)=RXRO(1)*EXP(RXRO(2)*XX+RXRO(3)*YY)
              ENDDO
              DO KK=4,6
                NZ=NOR(LEM(KK,IJ))
                XX=0.82085238*EX(NZ)+XMID
                YY=0.82085238*EY(NZ)+YMID
                DENS(KK,IJ)=RXRO(1)*EXP(RXRO(2)*XX+RXRO(3)*YY)
              ENDDO
              DENS(7,IJ)=RXRO(1)*EXP(RXRO(2)*XP+RXRO(3)*YP)
C
C   DR = 3 gives sinusoidal variation with depth
C
            ELSE IF(IDRULE.EQ.3)THEN
              XVRT=0.30385953*XP
              XMID=0.17914761*XP
              YVRT=0.30385953*YP
              YMID=0.17914761*YP
              DO KK=1,3
                NZ=NOR(LEM(KK,IJ))
                XX=0.69614048*EX(NZ)+XVRT
                YY=0.69614048*EY(NZ)+YVRT
                DENS(KK,IJ)=RXRO(1)*COS(RXRO(2)*XX+RXRO(3))*
     :                      COS(RXRO(4)*YY+RXRO(5))
              ENDDO
              DO KK=4,6
                NZ=NOR(LEM(KK,IJ))
                XX=0.82085238*EX(NZ)+XMID
                YY=0.82085238*EY(NZ)+YMID
                DENS(KK,IJ)=RXRO(1)*COS(RXRO(2)*XX+RXRO(3))*
     :                      COS(RXRO(4)*YY+RXRO(5))
              ENDDO
              DENS(7,IJ)=RXRO(1)*COS(RXRO(2)*XP+RXRO(3))*
     1           COS(RXRO(4)*YP+RXRO(5))
            END IF
C
C  MR defines material property type (default is 1)
C
            IF(IMRULE.NE.-1)IMAT(IJ)=IMRULE
          ENDIF
   40   CONTINUE
      END IF
C
C   Set viscosity and density flags
C
      IF (IVRULE.NE.-1) IVIS = 1
      IF (ISERULE.NE.-1) IVIS = 1
      IF (IDRULE.NE.-1) IDEN = 1
      IF (IMRULE.NE.-1) IMREG = 1
      GO TO 460
C
C    Finished reading instructions
C
  700 CONTINUE
C
      IF(IDBUG.NE.0)THEN
        WRITE(LUW,10123)
        CALL MATPRT(VHB,8,NE,8*NE,LUW)
10123   FORMAT('Dimensionless viscosity ratios are determined by',
     1  ' the VHB array (8,NE)',/)
      END IF
 500  CONTINUE
      END
      SUBROUTINE PARSEV(INSTR,LENGTH,IPE,NPTS,POLYG,MAXPTS,LTMP,
     :            IVX,IVRULE,IDRULE,IMRULE,ISERULE,RXSE,RXVC,RXRO,
     :            IVCCNT,IROCNT,IP,IERR)
C
C     This subroutine parses the special region definition statements
C     Used by VISVAR.  IP = 1 is returned if parsing is successful
C     else IP = 0 is returned
C
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER(1) IPE
      DIMENSION POLYG(2,MAXPTS),IVX(100)
      DIMENSION RXVC(10),RXRO(10)
      COMMON/AI/LUW,LSC,LBC,LLG
      IP=0
C
C     set default parameter values
C
      EPS=1.E-4
      IPE=' '
      NPTS=0
      IVRULE=-1
      IDRULE=-1
      IMRULE=-1
      ISERULE=-1
      RXSE=1.0
      DO I=1,10
        RXVC(I)=0.0
        RXRO(I)=0.0
      ENDDO
      IVCCNT=0
      IROCNT=0
C
C     check for keyword 'REG'
C
      CALL SKIP(INSTR,LENGTH,1,'n',J1)
      IF(J1.EQ.0)RETURN
      IF(J1.GE.80)RETURN
      IF(INSTR(J1:J1+2).NE.'REG')RETURN
C
C     A region definition command is now recognised
C     Problems in parsing the string go to 200 below
C
C     find the next non-zero character and set the IPE value
C
      K1=J1+3
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IPE=INSTR(J1:J1)
      IF((IPE.NE.'P').AND.(IPE.NE.'E').AND.(IPE.NE.'A'))GO TO 200
C
C     'P' is the polynomial definition using pairs of reals
C     'E' is the element number definition using list of integers
C
      IF (IPE.NE.'A') THEN
C     read '{' to signal data value input
C
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'{')GO TO 200
      ELSE
        NPTS=2
      END IF
C
C     read ordered pairs of coordinates, determine no.of corners
C      of polygon
C
      IF(IPE.EQ.'P')THEN
        K1=J1+1
        NPTS=0
   90   CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF (J1.EQ.0) THEN
          CALL READLINE(INSTR,LENGTH,LTMP,LSC,IEND)
          IF (IEND.EQ.1) GO TO 200
          K1=1
          CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        END IF
        IF(INSTR(J1:J1).EQ.'}')GO TO 20
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J2,'n',J3)
        IF(J3.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J3,'b',J4)
        IF(J4.EQ.0)GO TO 200
        NPTS=NPTS+1
        IF(NPTS.GT.MAXPTS)THEN
          WRITE(LSC,*)'Internal limit on NPTS passed in VISVAR: ',IPE
          GO TO 500
        END IF
        READ(INSTR(J1:J4-1),*,ERR=198)POLYG(1,NPTS),POLYG(2,NPTS)
        K1=J4
        GO TO 90
  198   WRITE(LSC,199)
  199   FORMAT('NB: do not split ordered pairs across lines of input',/,
     1  'leave spaces after commas and before and after :, }, {')
        GO TO 200
C
C     read list of element numbers and determine number of elements
C
      ELSE IF(IPE.EQ.'E')THEN
        K1=J1+1
        NPTS=0
   92   CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).EQ.'}')GO TO 20
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        NPTS=NPTS+1
        IF(NPTS.GT.100)THEN
          WRITE(LSC,*)'Internal limit on NPTS passed in VISVAR: ',IPE
          GO TO 500
        END IF
        READ(INSTR(J1:J2-1),*)IVX(NPTS)
        K1=J2
        GO TO 92
      END IF
C
C     read ':' to identify start of variables, keywords may be in 
C     any order but appropriate values must follow keyword
C
   20 K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.':')GO TO 200
      KB=J1+1
      LENGS=LENGTH-J1
C
C     read ':' and then Look for keyword VR
C
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'VR',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*,ERR=200)IVRULE
      END IF
C
C    Look for keyword DR and read parameter value
C
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'DR',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*,ERR=200)IDRULE
      END IF
C
C    Look for keyword MR and read parameter value
C
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'MR',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*,ERR=200)IMRULE
      END IF
C
C    Look for keyword SE and read parameter value
C
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'SE',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*,ERR=200)RXSE
        ISERULE=0
      END IF
C
C    Look for keyword RO and read parameter value
C
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'RO',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
  96    CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).EQ.';')GO TO 30
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        IROCNT=IROCNT+1
        IF(IROCNT.GT.10)THEN
          WRITE(LSC,*)'Too many RO values in input for VISVAR'
          GO TO 500
        END IF
        READ(INSTR(J1:J2-1),*,ERR=200)RXRO(IROCNT)
        K1=J2
        GO TO 96
      END IF
C
C    Look for keyword VC and read parameter value
C
   30 K1=J1+1
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'VC',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
  97    CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).EQ.';')GO TO 150
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        IVCCNT=IVCCNT+1
        IF(IVCCNT.GT.10)THEN
          WRITE(LSC,*)'Too many VC values in input for VISVAR'
          GO TO 500
        END IF
        READ(INSTR(J1:J2-1),*,ERR=200)RXVC(IVCCNT)
        K1=J2
        GO TO 97
      END IF
C
C     reading successfully completed
C
  150 IP=1
      IF(IPE.EQ.'A')THEN
        WRITE(LUW,10000)
        WRITE(LSC,10000)
10000   FORMAT('All elements:')
      ELSE IF(IPE.EQ.'P')THEN
        WRITE(LUW,10001)
        WRITE(LSC,10001)
10001   FORMAT('Within the polygonal region whose vertices are:')
        WRITE(LUW,10002)(POLYG(1,J),POLYG(2,J),J=1,NPTS)
        WRITE(LSC,10002)(POLYG(1,J),POLYG(2,J),J=1,NPTS)
10002   FORMAT('(',F10.4,',',F10.4,'),  ','(',F10.4,',',F10.4,')')
      ELSE IF(IPE.EQ.'E')THEN
        WRITE(LUW,10003)
        WRITE(LSC,10003)
10003   FORMAT('Within the region defined by the following list',
     1  ' of elements:')
        WRITE(LUW,10004)(IVX(J),J=1,NPTS)
        WRITE(LSC,10004)(IVX(J),J=1,NPTS)
10004   FORMAT(10I8)
      END IF
      IF((IVRULE.EQ.0).AND.(RXVC(1).NE.0.0))THEN
        WRITE(LUW,10005)RXVC(1)
        WRITE(LSC,10005)RXVC(1)
10005 FORMAT('Effective viscosity coefficient is set to VC =',
     1 G12.5)
      END IF
      IF((ABS(RXSE-1.).GE.EPS).AND.ISERULE.GE.0)THEN
        WRITE(LUW,10006)RXSE
        WRITE(LSC,10006)RXSE
      END IF
      IF(IDRULE.EQ.0)THEN
        WRITE(LUW,10007)RXRO(1)
        WRITE(LSC,10007)RXRO(1)
      END IF
10006 FORMAT('Stress vs strain-rate exponent is set to  SE =',
     1 G12.5)
10007 FORMAT('Constant density = ',G12.5)
      IF(IVRULE.EQ.1)THEN
        WRITE(LUW,10008)'Viscosity',RXVC(1),RXVC(2),RXVC(3)
        WRITE(LSC,10008)'Viscosity',RXVC(1),RXVC(2),RXVC(3)
      END IF
      IF(IDRULE.EQ.1)THEN
        WRITE(LUW,10008)'Density',RXRO(1),RXRO(2),RXRO(3)
        WRITE(LSC,10008)'Density',RXRO(1),RXRO(2),RXRO(3)
      END IF
10008 FORMAT(A10,' is a linear function of coordinates',/,
     1F8.3,' + ',F8.3,'x + ',F8.3,'y)')
      IF(IVRULE.EQ.2)THEN
        WRITE(LUW,10009)'Viscosity',RXVC(1),RXVC(2),RXVC(3)
        WRITE(LSC,10009)'Viscosity',RXVC(1),RXVC(2),RXVC(3)
      END IF
      IF(IVRULE.EQ.11)THEN
        WRITE(LUW,10123)RXVC(1),RXVC(3),RXVC(2),RXVC(4)
        WRITE(LSC,10123)RXVC(1),RXVC(3),RXVC(2),RXVC(4)
10123   FORMAT('Strength coefficient has been set so that VHB =',
     :        F8.3,' when depth = ', F8.3,'  linearly increasing',
     :        ' to ',F8.3,' when depth =',F8.3,/)
      END IF
      IF(IDRULE.EQ.2)THEN
        WRITE(LUW,10009)'Density',RXRO(1),RXRO(2),RXRO(3)
        WRITE(LSC,10009)'Density',RXRO(1),RXRO(2),RXRO(3)
      END IF
10009 FORMAT(A10,' is an exponential function of coordinates',/,
     1F9.4,' * exp(',F9.4,'x + ',F9.4,'y)')
      IF(IVRULE.EQ.3)THEN
        WRITE(LUW,10011)'Viscosity',RXVC(1)
        WRITE(LSC,10011)'Viscosity',RXVC(1)
      END IF
      IF(IDRULE.EQ.3)THEN
        WRITE(LUW,10013)'Density',RXRO(1),RXRO(2),RXRO(3),RXRO(4),
     1                   RXRO(5)
        WRITE(LSC,10013)'Density',RXRO(1),RXRO(2),RXRO(3),RXRO(4),
     1                   RXRO(5)
      END IF
10013 FORMAT(A10,' is a cosine function of coordinates',/,
     1F9.4,' * cos(',F9.4,'x +',F9.4,') * cos(',F9.4,'y + ',F9.4,')')
10011 FORMAT(A10,' is rescaled by the factor ',G12.4)
      IF(IMRULE.NE.-1)THEN
         WRITE(LUW,10012)IMRULE
         WRITE(LSC,10012)IMRULE
      END IF
10012 FORMAT('Material property parameters set to',I2)
      WRITE(LUW,*)
      WRITE(LSC,*)
      RETURN
C
C     problem in parsing the instruction
C
  200 CONTINUE
      WRITE(LUW,*)'PARSEV - problem interpreting command'
      WRITE(LSC,*)'PARSEV - problem interpreting command'
C     LPROB=((K1-1)/80)+1
C     LPINT=K1-(LPROB-1)*80
C     DO 201 J=1,80
C 201 PNTR(J:J)=' '
C     PNTR(LPINT:LPINT)='^'
C     WRITE(LUW,10010)
C     WRITE(LSC,10010)
C0010 FORMAT('The following instruction could not be ',
C    1'successfully parsed:')
C     DO 202 J=1,LPROB
C     J1=80*(J-1)+1
C     J2=80*J
C     WRITE(LUW,'(A80)')INSTR(J1:J2)
C     WRITE(LSC,'(A80)')INSTR(J1:J2)
C 202 CONTINUE
C     WRITE(LUW,'(A80)')PNTR
C     WRITE(LSC,'(A80)')PNTR
      CALL WRITEINSTR(INSTR,LENGTH,K1)
  500 IERR = 1
      RETURN
      END

      SUBROUTINE RHEOL(EX,EY,NOR,LEM,VHB,VOLD,IMAT,RHEO,IV,IRTYP,NE,NUP,
     :NN,IRK,IRHEOXC)
C
C     MIB2002
C     VHB = VOLD*B(x,y), VOLD is starting background viscosity.
C     RHEOTYP = 1:  Gaussian X dependence
C                   WHERE B(X) = 1/(1 + (A-1)*EXP(-((Xc-X)/b)**2))
C                   max(B(X)) = 1; min(B(X)) = 1/A
C     Allows for different A and b value in Crust and lithosphere
C     Uses IMAT to determine crust and lithosphere.
C     ===================================================================
C     RHEOTYP = 2:  Exponential decay in x direction
C     RHEO(1-3): 1-Xc, 2-b, 3-A
C     WHERE B(X) = exp((xc-x)/b)/(A+exp((xc-x)/b))
C     min(B(x))  = 1/A
C     max(B(x))  = ->1
C     ===================================================================
C     RHEO(1-5) are RealVars (RLV,RV) IVISP1-IVISP5
C     RHEO(1): Xc              [0 - default]
C     RHEO(2): b (lithosphere) [0]
C     RHEO(3): A (lithosphere) [1]
C     RHEO(4): b (crust)       [0]
C     RHEO(5): A (crust)       [1]
C     IRHEOXC: Xc' (needed for Xc not at reflecting boundary-- see below).
C     Note, to have no weak zone in the crust or lithosphere,
C     set A=1,b=0 for that layer.
C     ===================================================================
C     IRHEOXC: Need to keep track of center of convergence for asymmetric
C     models so that XC gets advected with overall motion of the lithosphere.
C     First time: get node for XC (XCNO) to use in each succesive time.
C     Save this in IRHEOXC. Then XC at each time step is:
C     XC = RHEO(1) + (EX(K)-RHEO(1)), K = NOR(IRHEOXC).
C     Replaced all RHEO(1) in equations by XC
C     Portion of code to get XCNO copied from GETNODENO.
C     Added IRK to RHEOL: IRK=1 for first time step.
C     ===================================================================

      DIMENSION IV(64)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION VHB(8,NE)
      DIMENSION VOLD(8,NE)
      DIMENSION RHEO(5)
      DIMENSION IMAT(NE)
C   local
      DIMENSION PTLOC(2,1),PTDIFF(2,1),BBOX(2,4)

C     Get node number for initial location of node at center of convergence.
      IF (IRK.EQ.1) THEN
         CALL FINDBBOX(EX,EY,BBOX,NUP)
         PTDIFF(1,1) = ABS(BBOX(1,2) - BBOX(1,1)) * 0.1
         PTDIFF(2,1) = ABS(BBOX(2,4) - BBOX(2,1)) * 0.1
         DIST = PTDIFF(1,1)*PTDIFF(1,1) + PTDIFF(2,1)*PTDIFF(2,1)

C     X location center of convergence
         PTLOC(1,1) = RHEO(1)
C     Y location top of box
         PTLOC(2,1) = BBOX(2,2)

         DO 30 N=1,NUP
            K=NOR(N)
            IF(K.LE.NN)THEN
               XDIF=ABS(EX(K)-PTLOC(1,1))
               IF(XDIF.LE.PTDIFF(1,1))THEN
                  YDIF=ABS(EY(K)-PTLOC(2,1))
                  IF(YDIF.LE.PTDIFF(2,1).AND.
     :                 ((YDIF*YDIF+XDIF*XDIF).LT.DIST))THEN
                     DIST = YDIF*YDIF+XDIF*XDIF
                     IRHEOXC = N
                  END IF
               END IF
            END IF
 30      CONTINUE
      END IF
C     Get current x-position of center node IRHEOXC
      K = NOR(IRHEOXC)
      DXC = (EX(K) - RHEO(1))
      XC = RHEO(1) + DXC

C     RHEOTYP=1, gaussian dependence in X.
      IF (IV(IRTYP).EQ.1) THEN
        DO 40 IJ=1,NE
          XP=0.0
          YP=0.0
C     crust
          IF (IMAT(IJ).EQ.2) THEN
            DO JJ=1,3
              XP=XP+EX(NOR(LEM(JJ,IJ)))
              YP=YP+EY(NOR(LEM(JJ,IJ)))
            ENDDO
            XP=XP/3.0
            YP=YP/3.0
            XVRT=0.30385953*XP
            XMID=0.17914761*XP
            YVRT=0.30385953*YP
            YMID=0.17914761*YP
            DO KK=1,3
              NZ=NOR(LEM(KK,IJ))
              XX=0.69614048*EX(NZ)+XVRT
              YY=0.69614048*EY(NZ)+YVRT
              VHB(KK,IJ)= VOLD(KK,IJ)/((RHEO(5)-1)*
     :                    EXP(-((XC-XX)/RHEO(4))**2) + 1)
            ENDDO
            DO KK=4,6
              NZ=NOR(LEM(KK,IJ))
              XX=0.82085238*EX(NZ)+XMID
              YY=0.82085238*EY(NZ)+YMID
              VHB(KK,IJ)= VOLD(KK,IJ)/((RHEO(5)-1)*
     :                    EXP(-((XC-XX)/RHEO(4))**2) + 1)
            ENDDO
            VHB(7,IJ)= VOLD(7,IJ)/((RHEO(5)-1)*
     :                 EXP(-((XC-XP)/RHEO(4))**2) + 1)
C     lithosphere
          ELSE IF (IMAT(IJ).EQ.1) THEN
            DO JJ=1,3
              XP=XP+EX(NOR(LEM(JJ,IJ)))
              YP=YP+EY(NOR(LEM(JJ,IJ)))
            ENDDO
            XP=XP/3.0
            YP=YP/3.0
            XVRT=0.30385953*XP
            XMID=0.17914761*XP
            YVRT=0.30385953*YP
            YMID=0.17914761*YP
            DO KK=1,3
              NZ=NOR(LEM(KK,IJ))
              XX=0.69614048*EX(NZ)+XVRT
              YY=0.69614048*EY(NZ)+YVRT
              VHB(KK,IJ)= VOLD(KK,IJ)/((RHEO(3)-1)*
     :                    EXP(-((XC-XX)/RHEO(2))**2) + 1)
            ENDDO
            DO KK=4,6
              NZ=NOR(LEM(KK,IJ))
              XX=0.82085238*EX(NZ)+XMID
              YY=0.82085238*EY(NZ)+YMID
              VHB(KK,IJ)= VOLD(KK,IJ)/((RHEO(3)-1)*
     :                    EXP(-((XC-XX)/RHEO(2))**2) + 1)
            ENDDO
            VHB(7,IJ)= VOLD(7,IJ)/((RHEO(3)-1)*
     :                 EXP(-((XC-XP)/RHEO(2))**2) + 1)
          ENDIF
 40     CONTINUE

C  RHEOTYP=2 is exponential variation in x.
      ELSE IF (IV(IRTYP).EQ.2) THEN
        DO 50 IJ=1,NE
          XP=0.0
          YP=0.0
C     crust
          IF (IMAT(IJ).EQ.2) THEN
            DO JJ=1,3
              XP=XP+EX(NOR(LEM(JJ,IJ)))
              YP=YP+EY(NOR(LEM(JJ,IJ)))
            ENDDO
            XP=XP/3.0
            YP=YP/3.0
            XVRT=0.30385953*XP
            XMID=0.17914761*XP
            YVRT=0.30385953*YP
            YMID=0.17914761*YP
            DO KK=1,3
              NZ=NOR(LEM(KK,IJ))
              XX=0.69614048*EX(NZ)+XVRT
              YY=0.69614048*EY(NZ)+YVRT
              YY=0.82085238*EY(NZ)+YMID
              EPT=EXP((XC-XX)/RHEO(4))
              VHB(KK,IJ)= VOLD(KK,IJ)*EPT/(RHEO(5)+EPT)
            ENDDO
            DO KK=4,6
              NZ=NOR(LEM(KK,IJ))
              XX=0.82085238*EX(NZ)+XMID
              YY=0.82085238*EY(NZ)+YMID
              EPT=EXP((XC-XX)/RHEO(4))
              VHB(KK,IJ)= VOLD(KK,IJ)*EPT/(RHEO(5)+EPT)
            ENDDO
            EPT=EXP((XC-XX)/RHEO(4))
            VHB(7,IJ)= VOLD(7,IJ)*EPT/(RHEO(5)+EPT)
C     lithosphere
          ELSE IF (IMAT(IJ).EQ.1) THEN
            DO JJ=1,3
              XP=XP+EX(NOR(LEM(JJ,IJ)))
              YP=YP+EY(NOR(LEM(JJ,IJ)))
            ENDDO
            XP=XP/3.0
            YP=YP/3.0
            XVRT=0.30385953*XP
            XMID=0.17914761*XP
            YVRT=0.30385953*YP
            YMID=0.17914761*YP
            DO KK=1,3
              NZ=NOR(LEM(KK,IJ))
              XX=0.69614048*EX(NZ)+XVRT
              YY=0.69614048*EY(NZ)+YVRT
              YY=0.82085238*EY(NZ)+YMID
              EPT=EXP((XC-XX)/RHEO(2))
              VHB(KK,IJ)= VOLD(KK,IJ)*EPT/(RHEO(3)+EPT)
            ENDDO
            DO KK=4,6
              NZ=NOR(LEM(KK,IJ))
              XX=0.82085238*EX(NZ)+XMID
              YY=0.82085238*EY(NZ)+YMID
              EPT=EXP((XC-XX)/RHEO(2))
              VHB(KK,IJ)= VOLD(KK,IJ)*EPT/(RHEO(3)+EPT)
            ENDDO
            EPT=EXP((XC-XX)/RHEO(2))
            VHB(7,IJ)= VOLD(7,IJ)*EPT/(RHEO(3)+EPT)
          ENDIF         ! on value of IMAT(IJ)
 50     CONTINUE
      END IF
      RETURN
      END

      SUBROUTINE VISVARFLAGS(LTMP,RXSE,IVIS,IDEN,IMREG,
     :                       LUW,LSC,IDBUG,IERR)
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*MAXLINES)
C
      EPS=1.E-4
  460 CONTINUE
      CALL READLINE(INSTR,LENGTH,LTMP,LSC,IEND)
      IF (IEND.EQ.1) GO TO 700
      CALL SKIP(INSTR,LENGTH,1,'n',K1)
      IF(INSTR(K1:K1+2).NE.'REG')GO TO 460
      PREVSE=RXSE
      CALL PARSEVFLAGS(INSTR,LENGTH,LTMP,RXSE,ISERULE,
     :                 IVRULE,IDRULE,IMRULE,IP,IERR)
      IF(IP.EQ.0)GO TO 460
C
C   Set viscosity and density flags
C
C     IF (IVRULE.NE.-1.AND.ABS(RXSE-1.0).LT.EPS) IVV = 1
C     IF (IVRULE.LT.1.AND.ABS(RXSE-1.0).GT.EPS) IVV = 2
C     IF (IVRULE.NE.-1.AND.ABS(RXSE-1.0).GT.EPS) IVV = 3
C   Check for IVV=4 after vhb array set up
      IF (ABS(RXSE-1.0).GT.EPS) IVIS = 1
      IF (IVRULE.NE.-1) IVIS = 1
      IF (ISERULE.NE.-1) IVIS = 1
      IF (IDRULE.NE.-1) IDEN = 1
      IF (IMRULE.NE.-1) IMREG = 1
      GO TO 460
C
C    Finished reading instructions
C
  700 CONTINUE
      RETURN
      END

      SUBROUTINE PARSEVFLAGS(INSTR,LENGTH,LTMP,RXSE,ISERULE,
     :                       IVRULE,IDRULE,IMRULE,IP,IERR)
C
C     This subroutine parses the special region definition statements
C     Used by VISVAR. 
C
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER IPE*1
      COMMON/AI/LUW,LSC,LBC,LLG
C
C     set default parameter values
C
      IPE=' '
      IP=0
      IVRULE=-1
      IDRULE=-1
      IMRULE=-1
      ISERULE=-1
C
C     check for keyword 'REG'
C
      CALL SKIP(INSTR,LENGTH,1,'n',J1)
      IF(J1.EQ.0)RETURN
      IF(J1.GE.80)RETURN
      IF(INSTR(J1:J1+2).NE.'REG')RETURN
C
C     A region definition command is now recognised
C     Problems in parsing the string go to 200 below
C
C     find the next non-zero character and set the IPE value
C
      K1=J1+3
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IPE=INSTR(J1:J1)
      IF((IPE.NE.'P').AND.(IPE.NE.'E').AND.(IPE.NE.'A'))GO TO 200
C
C     'P' is the polynomial definition using pairs of reals
C     'E' is the element number definition using list of integers
C
      IF (IPE.NE.'A') THEN
C     read '{' to signal data value input
C
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'{')GO TO 200
C
C     read ordered pairs of coordinates, determine no.of corners
C      of polygon
C
        K1=J1+1
   90   CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF (J1.EQ.0) THEN
          CALL READLINE(INSTR,LENGTH,LTMP,LSC,IEND)
          IF (IEND.EQ.1) GO TO 200
          K1=1
          CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        END IF
        IF(INSTR(J1:J1).EQ.'}')GO TO 20
        K1=J1+1
        GO TO 90
      END IF
C
C     read ':'
C
   20 K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.':')GO TO 200
C
C     look for keyword VR and read value
C
      KB=J1+1
      LENGS=LENGTH-J1
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'VR',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*,ERR=200)IVRULE
      END IF
C
C    Look for keyword DR and read parameter value
C
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'DR',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*,ERR=200)IDRULE
      END IF
C
C    Look for keyword MR and read parameter value
C
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'MR',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*,ERR=200)IMRULE
      END IF
C
C    Look for keyword SE and read parameter value
C
      CALL LOOK(INSTR(KB:LENGTH),LENGS,'SE',JF)
      IF(JF.NE.0)THEN
        K1=KB+JF+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*,ERR=200)RXSE
      END IF
C
C     reading successfully completed
C
      IP=1
      RETURN
C
C     problem in parsing the instruction
C
  200 CONTINUE
      CALL WRITEINSTR(INSTR,LENGTH,K1)
  500 IERR = 1
      RETURN
      END

      SUBROUTINE SETIVV(IVV,IVIS,SE,VHB,NE)
      DIMENSION VHB(8,NE)
C
C   set flag IVV:
C      0 for SE = 1 and constant viscosity
C      1 for SE = 1 and spatially variable viscosity
C      2 for SE > 1 (but constant) and constant viscosity
C      3 for SE > 1 (but constant) and spatially variable viscosity
C      4 for SE >= 1 and variable and spatially variable viscosity
C
      EPS=1.E-4
      IF (IVIS.GT.0) THEN
        VHAVE=0.0
        SEAVE=0.0
        DO I=1,NE
          DO K=1,7
            VHAVE=VHAVE+VHB(K,I)
          ENDDO
          SEAVE=SEAVE+VHB(8,I)
        ENDDO
        VHAVE=VHSUM/(FLOAT(7*NE))
        SEAVE=SEAVE/FLOAT(NE)
        ISEVAR=0
        IVCVAR=0
        DO I=1,NE
          DO K=1,7
            IF(ABS(VHB(K,I)-VHAVE).GT.EPS*VHAVE) IVCVAR=1
          ENDDO
          IF  (ABS(VHB(8,I)-SEAVE).GT.EPS*SEAVE) ISEVAR=1
        ENDDO
C       WRITE(*,*)'SETIVV: IVCVAR,ISEVAR =',IVCVAR,ISEVAR
        IF(ISEVAR.EQ.0)THEN               ! constant SE
          IF(ABS(SEAVE-1.0).LT.EPS)THEN   ! SE = 1
            IF(IVCVAR.EQ.0)THEN           ! VHB variable
              IVV=0
            ELSE
              IVV=1
            ENDIF
          ELSE                            ! SE != 1
            IF(IVCVAR.EQ.0)THEN           ! constant VHB
              IVV=2
            ELSE
              IVV=3
            ENDIF
          ENDIF
        ELSE                               ! SE variable
          IVV=4
        ENDIF
      ELSE     ! IVIS <= 0
        IF(ABS(SE-1.0).LT.EPS)THEN
          IVV=0
        ELSE
          IVV=2
        ENDIF
      ENDIF
C     WRITE(*,*)'SETIVV: IVV =',IVV,' IVIS =',IVIS
      RETURN
      END
