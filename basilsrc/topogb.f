C*--------------------------------------------------------------------
C*    Basil / Sybil:   topog.f  1.2  1 October 2002
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C  
C    Routines for input of data in longitude, latitude, elevation(m)
C      TOPODATA,TOPOINTRP,PROJECTXY,TRANSLATEMESH,GETPOLE,ROTATE
C      Function BILINEAR
C
      SUBROUTINE TOPODATA(FILER,ITYPE,GXMIN,GYMIN,GXMAX,GYMAX,
     :                    EX,EY,SSQ,NOR,NUP,
     :                    HLENSC,BDEPSC,RISOST,REFLEV,IERR)
C
C    Routine to read data from a digital elevation file of type
C    ITYPE: remains to be defined for different types.  At present
C    this routine assumes an ascii format xyz file with triplets 
C    of values (long, lat, elev) on a regular grid, with constant
C    and equal spacing in both lat and long directions. Extreme
C    values of long, lat, and elev are determined from the file.
C
C    The total number of data points should be less than MAXTOPO
C    as this is the length of the X,Y,ELEV arrays
C    MAXTOPO is set in limits.parameters and can be changed as 
C    these arrays are local to topog.f
C
      INCLUDE "limits.parameters"
      COMMON/AI/LUW,LSC,LBC,LLG
      CHARACTER FILER*60,INSTR*80,WORD*12
      DIMENSION X(MAXTOPO),Y(MAXTOPO),ELEV(MAXTOPO)
      DIMENSION EX(NUP),EY(NUP),SSQ(NUP),NOR(NUP)
      DIMENSION GRDDAT(3)
C
      LDIN=2
      EPS=1e-4
C
C    find the min and max values of latitude (GY) and
C    longitude (GX) in the dataset
C
      GXMIN=999.999
      GYMIN=999.999
      GXMAX=-999.999
      GYMAX=-999.999
      OPEN(LDIN,FILE=FILER,STATUS='old',IOSTAT=IOS)
      IF(IOS.NE.0) THEN
        WRITE(LSC,10125)FILER
        WRITE(LUW,10125)FILER
        IERR=1
        RETURN
      END IF
C
C   read data, line by line
C
      NREC=0
  102 CONTINUE
      READ(LDIN,*,END=10010,IOSTAT=IOS)(GRDDAT(K),K=1,3)
      IF (IOS.EQ.0.OR.IOS.EQ.112) THEN
        IF(GRDDAT(1).GT.GXMAX)GXMAX=GRDDAT(1)
        IF(GRDDAT(1).LT.GXMIN)GXMIN=GRDDAT(1)
        IF(GRDDAT(2).GT.GYMAX)GYMAX=GRDDAT(2)
        IF(GRDDAT(2).LT.GYMIN)GYMIN=GRDDAT(2)
        NREC=NREC+1
      END IF
      GO TO 102
10010 AREA=(GXMAX-GXMIN)*(GYMAX-GYMIN)
      SUM=(GXMAX+GYMAX)-(GXMIN+GYMIN)
      DINC=0.5*(SUM+SQRT(SUM*SUM+4*AREA*(NREC-1)))/(NREC-1)
      NCOL=INT((GXMAX-GXMIN+EPS)/DINC)+1
      NROW=INT((GYMAX-GYMIN+EPS)/DINC)+1
C
C    print parameters for data file, assuming regularly spaced data,
C    with equal grid spacing in each of X and Y directions,  NCOL and NROW
C    are inferred number of rows and colums in the dataset
C
        WRITE(LSC,1005)NREC,FILER
        WRITE(LUW,1005)NREC,FILER
        WRITE(LSC,1003)GXMIN,GXMAX,DINC
        WRITE(LUW,1003)GXMIN,GXMAX,DINC
        WRITE(LSC,1004)GYMIN,GYMAX,DINC
        WRITE(LUW,1004)GYMIN,GYMAX,DINC
 1003   FORMAT('Dataset longitude range:',F10.2,' to ',F10.2,' step :',
     1          F12.4)
 1004   FORMAT('Dataset latitude range:',F10.2,' to ',F10.2,' step :',
     1          F12.4)
 1005   FORMAT('Reading topography data: ',I7,' records from file ',A60)
      REWIND(LDIN)
C
C    check enough space allocated to store topodata
C
      IF (NREC.GT.MAXTOPO) THEN
        WRITE(LSC,1001)NREC,MAXTOPO
        WRITE(LUW,1001)NREC,MAXTOPO
 1001   FORMAT('Data points ',I6,' exceeds array limit (MAXTOPO) ',I6) 
        IERR=1
        RETURN
      END IF
C
C    Read the xyz data, after setting null values on grid
C
      EMIN=999.999
      EMAX=-999.999
      DO J=1,NREC
        X(J)=EMAX
        Y(J)=EMAX
      ENDDO
  101 CONTINUE
      READ(LDIN,*,END=10011,IOSTAT=IOS)(GRDDAT(K),K=1,3)
C
C   check for undefined values in topo dataset
C
      IF (IOS.EQ.112) THEN
        BACKSPACE(LDIN)
        READ(LDIN,'(A80)',IOSTAT=IOS)INSTR
        WRITE(LSC,*)'Problem point: ',INSTR
        N=INDEX(INSTR,'NaN')
        IF (N.EQ.0) IERR=1
      ELSE IF (IOS.EQ.0) THEN
C
C   get point indices, careful to avoid rounding toward zero
C   put elevation and coordinates into rasterised grid, stored
C   in rows of constant latitude, increasing E and N from SW corner
C
        GRDDIF=GRDDAT(1)-GXMIN
        IF(GRDDIF.GT.0)THEN
          IPTX = INT((GRDDIF+EPS)/DINC)+1
        ELSE
          IPTX = INT((GRDDIF-EPS)/DINC)+1
        END IF
        GRDDIF=GRDDAT(2)-GYMIN
        IF(GRDDIF.GT.0)THEN
          IPTY = INT((GRDDIF+EPS)/DINC)+1
        ELSE
          IPTY = INT((GRDDIF-EPS)/DINC)+1
        END IF
        IPT=NCOL*(IPTY-1)+IPTX
        X(IPT)=GRDDAT(1)
        Y(IPT)=GRDDAT(2)
        ELEV(IPT)=GRDDAT(3)
        IF (GRDDAT(3).GT.EMAX) EMAX=GRDDAT(3)
        IF (GRDDAT(3).LT.EMIN) EMIN=GRDDAT(3)
      END IF
      GO TO 101
10011 CONTINUE
C
C     check if any points have not been set
C
        DO J=1,NROW
          DO I=1,NCOL
            K=(J-1)*NCOL+I
            IF((ABS(X(K)).GT. 999.).OR.(ABS(Y(K)).GT.999))THEN
            WRITE(*,*)'Check TOPODATA: point (I,J) = (',I,',',J,')',
     1                ' not set'
            END IF
          ENDDO
        ENDDO
        WRITE(LSC,1002) EMIN,EMAX
        WRITE(LUW,1002) EMIN,EMAX
C
C    proceed to the interpolation
C
        CALL TOPOINTRP(X,Y,ELEV,NCOL,NROW,NPTS,
     :                 EX,EY,SSQ,NOR,NUP,GXMIN,GYMIN,DINC,
     :           HLENSC,BDEPSC,RISOST,REFLEV,IERR)
      CLOSE(LDIN)
      RETURN
 1002 FORMAT('Dataset elevation range           ',F10.3,' -> ',F10.3)
10125 format(' The read file could not be opened ',a60,/)
      END

      SUBROUTINE TOPOINTRP(X,Y,ELEV,NCOL,NROW,NPTS,
     :                     EX,EY,SSQ,NOR,NUP,GXMIN,GYMIN,DINC,
     :         HLENSC,BDEPSC,RISOST,REFLEV,IERROR)
C
C   this interpolation routine puts the finite element mesh on top
C   of the topo data grid read in TOPOREAD and makes a bilinear
C   interpolation for any node that is in the domain of the dataset
C   nodes that are outside that domain are set to the value REFLEV
C   At this time coordinates are still in raw degrees of long and lat
C   with elevations assumed to be in metres.
C
      INCLUDE "limits.parameters"
      DIMENSION X(MAXTOPO)
      DIMENSION Y(MAXTOPO)
      DIMENSION ELEV(MAXTOPO)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION SSQ(NUP)
      DIMENSION NOR(NUP)
      COMMON/AI/LUW,LSC,LBC,LLG
      DIMENSION EVAL(4)
C
      EPS=1.E-5
      EMIN=REFLEV
      EMAX=REFLEV
      ISTOP=0
      NNOT=0
      DO 10 I=1,NUP
        SSQQ=REFLEV
        XPT=EX(NOR(I))
        YPT=EY(NOR(I))
C
C   find offset into arrays
C
        IXOFF=INT((XPT-GXMIN)/DINC)+1
        IYOFF=INT((YPT-GYMIN)/DINC)+1
        IF(XPT.LT.GXMIN)IXOFF=INT((XPT-GXMIN)/DINC)
        IF(YPT.LT.GYMIN)IYOFF=INT((YPT-GYMIN)/DINC)
        IOFF=NCOL*(IYOFF-1)+IXOFF
C
C   check that point to be interpolated is within available data grid
C
        IF((IXOFF.GE.1).AND.(IXOFF.LT.NCOL).AND.
     :     (IYOFF.GE.1).AND.(IYOFF.LT.NROW)) THEN
           XB0=GXMIN+DINC*FLOAT(IXOFF-1)
           XB1=GXMIN+DINC*FLOAT(IXOFF)
           YB0=GYMIN+DINC*FLOAT(IYOFF-1)
           YB1=GYMIN+DINC*FLOAT(IYOFF)
C
C    EVAL holds the 4 elevation values in order BL,BR,TL,TR
C    on the grid surrounding the node to be interpolated
C
           EVAL(1)=ELEV(IOFF)
           EVAL(2)=ELEV(IOFF+1)
           EVAL(3)=ELEV(IOFF+NCOL)
           EVAL(4)=ELEV(IOFF+NCOL+1)
        IF((XPT.LT.XB0).OR.(XPT.GT.XB1).OR.(YPT.LT.YB0).OR.(YPT.GT.YB1))
     :  WRITE(*,*)XB0,XPT,XB1,YB0,YPT,YB1
C
C   do the bilinear interpolation
C
          SSQQ=BILINEAR(XB0,XB1,YB0,YB1,XPT,YPT,EVAL,REFLEV)
          IF(SSQQ.GT.EMAX) EMAX=SSQQ
          IF(SSQQ.LT.EMIN) EMIN=SSQQ
C
C   if points are not in the data grid they are set to REFLEV
C   following message appears as a warning that at least 1 is external
C
        ELSE
          NNOT=NNOT+1
          IF(NNOT.EQ.1)THEN
            WRITE(*,*)'X =',XPT,' Y =',YPT,' appears to be external'
          ENDIF
        END IF
C
C     convert topography into dimensionless layer thickness
C     of which the natural logarithm is stored.
C
        SSQR=(1.0+0.001*(SSQQ-REFLEV)*RISOST/BDEPSC)/HLENSC
        SSQ(I)=ALOG(ABS(SSQR))
C     ABS is probably an unnecessary precaution, but it stops a
C     crash in case of rubbish input.
  10  CONTINUE
C
      WRITE(LSC,1002) EMIN,EMAX
      WRITE(LUW,1002) EMIN,EMAX
      WRITE(LSC,1003)NNOT
      WRITE(LUW,1003)NNOT
      WRITE(LSC,10111)REFLEV,RISOST,BDEPSC,HLENSC
      WRITE(LUW,10111)REFLEV,RISOST,BDEPSC,HLENSC
      RETURN
 1002 FORMAT('Interpolated elevation data range ',F10.3,' -> ',F10.3)
 1003 FORMAT('Interpolation: ',I5,' nodes lie outside domain of data')
10111 FORMAT('Scaling Parameters for Isostatic calculation: ',
     :'REFLEV = ',F12.3,/,'RISOST = ',F12.3,'  BDEPSC = ',F12.3,
     :'  HLENSC = ',F12.3,/)
      END

      SUBROUTINE PROJECTXY(EX,EY,CENTREL,YMID,IVERB,NUP,IERROR)
C
C    Routine to project data in longitude, latitude into x, y space
C    Radius is assumed to be 1
C    Co-Latitude values use the range 0->180 (0 at S.Pole)
C    Note that projection depends on the central longitude: CENTREL
C    Degrees are converted to radians (DTOR=PI/180)
C
      INCLUDE "limits.parameters"
      COMMON/AI/LUW,LSC,LBC,LLG
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION BBOX(2,4)
C
      EXMIN=999.999
      EYMIN=999.999
      EXMAX=-999.999
      EYMAX=-999.999
      DO JE=1,NUP
        IF(EX(JE).LT.EXMIN)EXMIN=EX(JE)
        IF(EY(JE).LT.EYMIN)EYMIN=EY(JE)
        IF(EX(JE).GT.EXMAX)EXMAX=EX(JE)
        IF(EY(JE).GT.EYMAX)EYMAX=EY(JE)
      ENDDO
      YMID=0.5*(EYMIN+EYMAX)
      IF(IVERB.EQ.1)THEN
        WRITE(LSC,10123)CENTREL
        WRITE(LUW,10123)CENTREL
10123   FORMAT('Sinusoidal equal area projection is applied, using ',
     :         'center longitude',F9.2,/)
 
        WRITE(LSC,10131)NUP,EXMIN,EXMAX,EYMIN,EYMAX
        WRITE(LUW,10131)NUP,EXMIN,EXMAX,EYMIN,EYMAX
10131   FORMAT('Before spherical projection:'
     : /,      'For ',I7,' nodes, longitude (X) :',F9.3,' to ',F9.3,
     : /,      '                   latitude  (Y) :',F9.3,' to ',F9.3,/)
        IF((CENTREL.LT.EXMIN).OR.(CENTREL.GT.EXMAX))
     :  WRITE(LSC,10132)CENTREL,EXMIN,EXMAX
10132   FORMAT('Warning: Central longitude: ',F8.3,' lies outside ',
     :         'range defined by min and max longitudes:',2F8.3)
      END IF
C
C     IF (NCOMP.LT.0) THEN
        DO 10 I=1,NUP
          EY(I)=EY(I)*DTOR+PION2
          EX(I)=(EX(I)-CENTREL)*DTOR*SIN(EY(I))
  10    CONTINUE
C
C    locally Mercator projection, no longer used
C
C     ELSE
C       YVAL=SIN(YMID*DTOR+PION2)
C       DO 20 I=1,NUP
C         EY(I)=EY(I)*DTOR+PION2
C         EX(I)=(EX(I))*DTOR*YVAL
C 20    CONTINUE
C     END IF
      IF(IVERB.EQ.1)THEN
        EXMIN=999.999
        EYMIN=999.999
        EXMAX=-999.999
        EYMAX=-999.999
        DO JE=1,NUP
          IF(EX(JE).LT.EXMIN)EXMIN=EX(JE)
          IF(EY(JE).LT.EYMIN)EYMIN=EY(JE)
          IF(EX(JE).GT.EXMAX)EXMAX=EX(JE)
          IF(EY(JE).GT.EYMAX)EYMAX=EY(JE)
        ENDDO
        WRITE(LSC,10133)NUP,EXMIN,EXMAX,EYMIN,EYMAX
        WRITE(LUW,10133)NUP,EXMIN,EXMAX,EYMIN,EYMAX
10133   FORMAT('After spherical projection:'
     : /,      'For ',I7,' nodes, X coordinates:',F9.5,' to ',F9.5,
     : /,      '                   Y coordinates:',F9.5,' to ',F9.5,/)
        END IF
      RETURN
      END
      SUBROUTINE PROJECTDEG(EX,EY,XMID,YMID,NUP,NCOMP,IERROR)
C
C    Routine to project data in x, y  into longitude, latitude space
C    i.e. inverse operation to that defined by PROJECTXY
C    Latitude values adjusted to the range 0->180 (0 at S.Pole)
C    Longitude values centred on the longitude at centre of projection
C    Convert to degrees (RTOD=180/PI)
C
      INCLUDE 'limits.parameters'
      DIMENSION EX(NUP),EY(NUP)

      IERROR=0
C     WRITE(6,1000)EX(1),EY(1)
C
C    this block for sinusoidal equal area projection
C
      IF (NCOMP.LT.0) THEN
        DO 10 I=1,NUP
          EX(I)=EX(I)/SIN(EY(I))/DTOR+XMID
          EY(I)=(EY(I)-PION2)/DTOR
  10    CONTINUE
C
C    this block uses invariant X-scale; actually used ?
C
      ELSE
        YVAL=ASIN(YMID-PION2)/DTOR
        DO 20 I=1,NUP
          EY(I)=(EY(I)-PION2)/DTOR
          EX(I)=EX(I)/SIN(YVAL)/DTOR+XMID
  20    CONTINUE
      END IF
C     WRITE(6,1001)EX(1),EY(1)
      RETURN
 1000 FORMAT('Projecting x,y (',F6.2,', ',F6.2,')')
 1001 FORMAT('to lat, long (',F6.3,', ',F6.3,')')
      END

      REAL FUNCTION BILINEAR(XB0,XB1,YB0,YB1,XPT,YPT,EVAL,DFLT)
      DIMENSION EVAL(4)
      BILINEAR=DFLT
      IF((EVAL(1).NE.DFLT).AND.(EVAL(2).NE.DFLT).AND.
     :   (EVAL(3).NE.DFLT).AND.(EVAL(4).NE.DFLT))THEN
        A=(XPT-XB0)/(XB1-XB0)
        C=(YPT-YB0)/(YB1-YB0)
        BILINEAR=EVAL(3)*(1.0-A)*C +       EVAL(4)*A*C
     :          +EVAL(1)*(1.0-A)*(1.0-C) + EVAL(2)*A*(1.0-C)
      END IF
      END
 
      SUBROUTINE TRANSLATEMESH(XOFF,XSCALE,YOFF,YSCALE,
     :                         EX,EY,ISIZE)
      INCLUDE "limits.parameters"
      DIMENSION EX(ISIZE),EY(ISIZE)
      EPS=1e-5
      IF (ABS(XOFF).GT.EPS) THEN
        DO I=1,ISIZE
          EX(I) = AMOD(EX(I)+XOFF,RLONGMAX)
        ENDDO
      END IF
      IF (ABS(YOFF).GT.EPS) THEN
        DO I=1,ISIZE
          EY(I) = EY(I)+YOFF
        ENDDO
      END IF
      RETURN
      END
C
      SUBROUTINE ROTATE(ROTIX,RLON,RLAT,NP,DPHI)
C
C    routine ROTATE: computes new coordinates of a set of points
C    (RLON,RLAT) on the surface of the sphere, given a rotation
C    defined by the matrix ROTIX computed in GETPOLE
C
C    Input Coordinates and rotation in degrees
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*4 RLON(NP),RLAT(NP),DPHI
      DIMENSION ROTIX(3,3),RHS(3),RES(3)
C
C   constant factors
C
      PI=3.14159265358979323846D0
      DEGRAD=PI/180.D0
      PION2=0.5d0*PI
      DO I=1,NP
C
C   evaluate the right hand side using the long/lat of the point
C   and effect the rotation using a matrix multiply of ROTIX
C
        GLA=RLAT(I)*DEGRAD
        GLO=RLON(I)*DEGRAD
        CLATG=DCOS(GLA)
        CLONG=DCOS(GLO)
        SLONG=DSIN(GLO)
        SLATG=DSIN(GLA)
        RHS(1)=CLATG*CLONG
        RHS(2)=CLATG*SLONG
        RHS(3)=SLATG
        DO K=1,3
          RES(K)=0.d0
          DO J=1,3
            RES(K)=RES(K)+RHS(J)*ROTIX(K,J)
          ENDDO
        ENDDO
C
C    get new latitude and longitude from RES(1-3), convert to degrees
C
        SLON=DATAN2(RES(2),RES(1))
        IF(SLON.LT.-PI)SLON=SLON+2.0*PI
        IF(SLON.GT.PI)SLON=SLON-2.0*PI
        RLON(I)=SLON/DEGRAD-DPHI
        RLAT(I)=(DASIN(RES(3)))/DEGRAD
C
C    flattening of sphere using sinusoidal equal area projection
C    is done with subsequent call to PROJEXTXY if needed
C
      ENDDO                          ! I = 1 to NP
C
      RETURN
      END
C
      SUBROUTINE GETPOLE(PT2,PLON,PROT,ROTIX)
C
C    routine to (1) determine longitude of equatorial pole of rotation 
C    and angle of rotation given two points (PT2) that define a great 
C    circle which is to be rotated into coincidence with the equator.  
C    (2) evaluate the entries of the matrix to be applied (ROTIX)
C    in routine ROTATE.  PT2 input in degrees, PLON,PROT returned
C    in radians (double precision).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*4 PT2
      DIMENSION PT2(4),ROTIX(3,3)
C
C   constant factors
C
      PI=3.14159265358979323846D0
      DEGRAD=PI/180.D0
C
C     convert degrees to radians and get sines and cosines
C
      ALA=PT2(2)*DEGRAD
      ALO=PT2(1)*DEGRAD
      BLA=PT2(4)*DEGRAD
      BLO=PT2(3)*DEGRAD
      CLONA=DCOS(ALO)
      SLONA=DSIN(ALO)
      CLATA=DCOS(ALA)
      SLATA=DSIN(ALA)
      CLONB=DCOS(BLO)
      SLONB=DSIN(BLO)
      CLATB=DCOS(BLA)
      SLATB=DSIN(BLA)
C
C      assemble and solve the 2x2 matrix
C
      A11=CLATA*CLONA
      A12=CLATA*SLONA
      A21=CLATB*CLONB
      A22=CLATB*SLONB
      DETM=A11*A22-A21*A12
C
C     there could be a problem if both points have the same longitude
C     if so the required rotation would be 90 degrees about that longitude
C
      IF(DABS(DETM).LE.1.D-60)THEN
        WRITE(*,*)'Singular matrix in GETPOLE, check parameters'
        PLONR=ALO
        PROTR=PI*0.5D0
      ELSE
        PX=(-A22*SLATA+A12*SLATB)/DETM
        PY=( A21*SLATA-A11*SLATB)/DETM
C
C     normalize the vector and compute pole coordinates and rotation
C
        PMAG=DSQRT(1.D0+PX*PX+PY*PY)
        PX=PX/PMAG
        PY=PY/PMAG
        SINOM=DSQRT(PX*PX+PY*PY)
        PROTR=DASIN(SINOM)
        PLONR=DATAN2(-PX,PY)
      ENDIF
      IF(PLONR.LT.-PI)PLONR=PLONR+2.0*PI
      IF(PLONR.GT.PI)PLONR=PLONR-2.0*PI
C
C     return PLON, PROT in degrees
C
      PLON=PLONR/DEGRAD
      PROT=PROTR/DEGRAD
C
C     get components of vector q
C
      QX=PY/SINOM
      QY=-PX/SINOM
C
C     evaluate the entries for the Rotation matrix
C
      COSOM=DSQRT(1.d0-SINOM*SINOM)
      ROTIX(1,1)=QX*QX*(1.D0-COSOM)+COSOM
      ROTIX(2,2)=QY*QY*(1.D0-COSOM)+COSOM
      ROTIX(3,3)=COSOM
      ROTIX(2,1)=QX*QY*(1.D0-COSOM)
      ROTIX(1,3)=QY*SINOM
      ROTIX(3,2)=QX*SINOM
      ROTIX(1,2)=ROTIX(2,1)
      ROTIX(3,1)=-ROTIX(1,3)
      ROTIX(2,3)=-ROTIX(3,2)
C
      RETURN
      END
      SUBROUTINE MATSET(PLON,PROT,ROTIX)
C
C    routine to set up rotation matrix, given pole of rotation
C    specified by longitude PLON (lat assumed zero) and rotation
C    angle PROT.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ROTIX(3,3)
C
C   constant factors
C
      PI=3.14159265358979323846D0
      DEGRAD=PI/180.D0
C
C     convert degrees to radians and get components of vector q
C
      PLONR=PLON*DEGRAD
      QX=DCOS(PLONR)
      QY=DSIN(PLONR)
C
C     evaluate the entries for the Rotation matrix
C
      PROTR=PROT*DEGRAD
      COSOM=DCOS(PROTR)
      SINOM=DSIN(PROTR)
      ROTIX(1,1)=QX*QX*(1.D0-COSOM)+COSOM
      ROTIX(2,2)=QY*QY*(1.D0-COSOM)+COSOM
      ROTIX(3,3)=COSOM
      ROTIX(2,1)=QX*QY*(1.D0-COSOM)
      ROTIX(1,3)=QY*SINOM
      ROTIX(3,2)=QX*SINOM
      ROTIX(1,2)=ROTIX(2,1)
      ROTIX(3,1)=-ROTIX(1,3)
      ROTIX(2,3)=-ROTIX(3,2)
C
      RETURN
      END
