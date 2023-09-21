C
C          Program: circles
C
C    routine to calculate coordinates of points defined on either
C    a great circle segment, or on a small circle segment on the sphere.
C    or on the Cartesian plane
C
C    multiple modes of operation are defined by a single command line
C    that specifies two Keycodes, followed by 4 or 5 real values and 
C    an integer, whose purpose is defined as in the following examples.
C    The ingeger NL in each case specifies how many points are determined
C    on each line segment.  NL includes both end-points, so set to 2 if
C    only endpoints are required. The arc/line is subdivided into (NL-1)
C    equal increments.  
C
C    Modes of operation
C
C     GRTC PTS lon1 lat1 lon2 lat2 NL
C       points are computed on a great circle defined by two endpoints
C
C     GRTC DAZ lon1 lat1 dist azim NL
C       points are computed on a great circle defined by one endpoint, a
C       distance increment between points and an azimuth.  Distance in 
C       degrees, or use KAZ for distance in km.
C
C     SMLC PAZ plon plat dist azim1 azim2 NL
C       points are computed on a small circle at a specified distance
C       from a pole between two azimuth limits.
C
C     SMLC PLA plon plat dist lat1 lat2 NL
C       points are computed on a small circle at a specified distance
C       from a pole, between two latitude limits
C
C     longitudes may use 0 to 360 degrees, or -180 to +180 degrees
C     for latitudes, positive is north, negative is south.
C
C     CIRC x_o y_o radius theta_1 theta_2 NL
C       points on a Cartesian plane relative to centre (x_o, y_o)
C       radius rad between angles (degrees) theta_1 and theta_2
C       angles in degrees, relative to zero on the x-axis.
C       A closed circle is obtained if theta_2-theta_1 = 360.0
C
C     LINE x_0 y_0 x_1 y_1 NL
C       Points on a straight line between (x_0, y_0) and (x_1, Y_1)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER LINE*120
      CHARACTER ICRC*4,ITYP*3
      PARAMETER (NLINEP=200)
      DIMENSION DISTDEG(NLINEP)
      DIMENSION XPAR(5),AZIMS(2)
C
      PI=3.141592653589793D0
      D2R=PI/180.D0
      ONED=10000.d0/90.d0
      LUR=3
      LUW=6
C
C    open file to read command line
C
      OPEN(LUR,FILE='circles.in')
C
C   read successive command lines to end of file
C
   20 CONTINUE
      DO K=1,12
        LINE((K-1)*10+1:K*10)='          '
      ENDDO
      READ(LUR,'(A120)',END=30)LINE
      IF(LINE(1:3).EQ.'GRT')THEN
        READ(LINE,*,ERR=30)ICRC,ITYP,(XPAR(J),J=1,4),NLINE
      ELSEIF(LINE(1:3).EQ.'SML')THEN
        READ(LINE,*,ERR=30)ICRC,ITYP,(XPAR(J),J=1,5),NLINE
      ELSEIF(LINE(1:3).EQ.'LIN')THEN
        READ(LINE,*,ERR=30)ICRC,(XPAR(J),J=1,4),NLINE
      ELSEIF(LINE(1:3).EQ.'CIR')THEN
        READ(LINE,*,ERR=30)ICRC,(XPAR(J),J=1,5),NLINE
      ELSE
        WRITE(*,*)'First word of circles.in should be either ',
     :            'LINE, CIRC, GRTC or SMLC'
        STOP
      ENDIF
      IF(NLINE.LT.2)NLINE=2
      IF(NLINE.GT.NLINEP)NLINE=NLINEP
C
C      Entering Lat and Long limits for profile
C
      IF(ICRC.EQ.'GRTC')THEN
        XLON=XPAR(1)
        XLAT=XPAR(2)
        IF(ITYP.EQ.'PTS')THEN
          ELON=XPAR(3)
          ELAT=XPAR(4)
          WRITE(LUW,10001)XLON,XLAT,ELON,ELAT,NLINE
10001     FORMAT('Point A (Lon,Lat) : ',F12.5,',',F12.5,/,
     :       'Point B (Lon,Lat) : ',F12.5,',',F12.5,/,
     :       'Number of points along profile :',I6)
C
C      calculate distance and bearing
C
          DLONG=ELON-XLON
          IF(DLONG.LT.-180.d0)DLONG=DLONG+360.D0
          IF(DLONG.GT.180.d0)DLONG=DLONG-360.D0
          CDE=DSIN(ELAT*D2R)*DSIN(XLAT*D2R) + 
     :        DCOS(ELAT*D2R)*DCOS(XLAT*D2R)*DCOS(DLONG*D2R)
          DELT=DACOS(CDE)
          CBE=(DSIN(ELAT*D2R)*DCOS(XLAT*D2R) - DSIN(XLAT*D2R)*
     :        DCOS(ELAT*D2R)*DCOS(DLONG*D2R))/DSIN(DELT)
          IF (CBE.GT.1.d0)CBE=1.0d0
          IF (CBE.LT.-1.d0)CBE=-1.0d0
C
          AZIM=DACOS(CBE)
          BEAR=AZIM/D2R
          IF(DLONG.LT.0.0)BEAR=-BEAR
          DEDE=DELT/D2R
          DEKM=DEDE*ONED
          WRITE(LUW,10100)XLON,XLAT,ELON,ELAT
10100     FORMAT('Great circle from ', F10.5,' E, ',F10.5,' N to ',
     :          F10.5,' E, ',F10.5,' N: ')
          WRITE(LUW,10101)DEDE,DEKM,BEAR
10101     FORMAT('Distance (deg) =',F10.5,'  Distance (km)  =',F10.1,
     :           '  Bearing (deg) =',F10.5)
          DINC=DELT/DFLOAT(NLINE-1)
          ITV=1
          IF(ELON.LT.XLON)ITV=-1
          JTV=1
C
C     circle defined by distance and bearing
C
        ELSEIF(ITYP.EQ.'DAZ'.OR.ITYP.EQ.'KAZ')THEN
          IF(ITYP.EQ.'KAZ')THEN
            DEKM=XPAR(3)
            DEDE=DEKM/ONED
          ELSE
            DEDE=XPAR(3)
            DEKM=DEDE*ONED
          ENDIF
          BEAR=XPAR(4)
          WRITE(LUW,10102)XLAT,XLON,DEKM,DEDE,BEAR
10102     FORMAT('Great circle from ', F10.5,' N, ',F10.5,' E, ',
     :          ' Distance (km)',F10.1,', (deg)',F10.5,'  Bearing: ',
     :          F10.5,' E ')
          DINC=DEDE*D2R
          DELT=DINC*FLOAT(NLINE-1)
          AZIM=BEAR*D2R
          CBE=DCOS(AZIM)
          CDE=DCOS(DELT)
          ITV=1
          IF((BEAR.GT.180.0).OR.(BEAR.LT.0.0))ITV=-1
        END IF     ! ITYP = PTS or DAZ
C
C    for a Great Circle: increment distance 
C
        CTX=DCOS(XLAT*D2R)
        STX=DSIN(XLAT*D2R)
        CFX=DCOS(XLON*D2R)
        SFX=DSIN(XLON*D2R)
        DO 200 J=1,NLINE
          DIST=DFLOAT(J-1)*DINC
          DISTDEG(J)=DIST/D2R
          TVA=DSIN(DIST)*CTX*CBE
          TVB=DCOS(DIST)*STX
          TVLATR=DASIN(TVA + TVB)
          TVC=DCOS(DIST) - STX*DSIN(TVLATR)
          TVD=DCOS(TVLATR)*CTX
          TVE=TVC/TVD
          IF (TVE.GT.1.0)TVE=1.0
          IF (TVE.LT.-1.0)TVE=-1.0
          TVLONR=(XLON*D2R) + FLOAT(ITV)*DACOS(TVE)
          ELON=TVLONR/D2R
          ELAT=TVLATR/D2R
          WRITE(LUW,10050)ELON,ELAT
10050     FORMAT(2F11.5)
  200   CONTINUE
        CTE=DCOS(ELAT*D2R)
        STE=DSIN(ELAT*D2R)
        CFE=DCOS(ELON*D2R)
        SFE=DSIN(ELON*D2R)
C
C    calculate and print pole parameters (location and back-azimuths)
C    using cross product of cartesian vectors, azimuths from cosine law
C
        XCP=STE*CTX*SFX - STX*CTE*SFE
        YCP=STX*CTE*CFE - STE*CTX*CFX
        ZCP=CTX*CFX*CTE*SFE - CTX*SFX*CTE*CFE
        HMAG=DSQRT(XCP*XCP+YCP*YCP)
        POLON=(DATAN2(YCP,XCP))/D2R
        POLATR=DATAN2(ZCP,HMAG)
        COSTP=DCOS(POLATR)
        POLAT=POLATR/D2R
C
C     DACOS is ambiguous in sign, determined by relative longitude
C
        AZIMX=(DACOS(STX/COSTP))/D2R
        DPX=POLON-XLON
        IF((DPX.GT.0.d0).OR.(DPX.LT.-180.d0))AZIMX=-AZIMX
        AZIME=(DACOS(STE/COSTP))/D2R
        DPE=POLON-ELON
        IF((DPE.GT.0.d0).OR.(DPE.LT.-180.d0))AZIME=-AZIME
        WRITE(*,10104)POLON,POLAT,AZIMX,AZIME
10104   FORMAT('For pole at ', F10.5,' E  ',F10.5,' N, ',
     :         'Azimuths of end points are: ',
     :          F10.5,' and ',F10.5,' E of N')
C
C    for a Small Circle: define path relative to pole
C
      ELSEIF(ICRC.EQ.'SMLC')THEN
        POLONR=XPAR(1)*D2R
        POLATR=XPAR(2)*D2R
        DIST=XPAR(3)*D2R
        CPOLATR=DCOS(POLATR)
        SPOLATR=DSIN(POLATR)
        CDIST=DCOS(DIST)
        SDIST=DSIN(DIST)
        IF(ITYP.EQ.'PAZ')THEN     ! azimuth range specified
          AZIMS(1)=XPAR(4)*D2R
          AZIMS(2)=XPAR(5)*D2R
        ELSEIF(ITYP.EQ.'PLA')THEN ! latitude range specified
          DO K=1,2
            XVC=(DSIN(XPAR(3+K)*D2R)-SPOLATR*CDIST)/(CPOLATR*SDIST)
            IF(DABS(XVC).LE.1.D0)THEN
              AZIMS(K)=DACOS(XVC)
            ELSE
              WRITE(*,*)'range of latitudes as specified appears ',
     :                'inconsistent with pole location'
              STOP
            ENDIF
          ENDDO
        ENDIF
        AZIMAD=AZIMS(1)/D2R
        AZIMBD=AZIMS(2)/D2R
        WRITE(LUW,10103)XPAR(3),XPAR(1),XPAR(2),AZIMAD,AZIMBD
10103   FORMAT('Small circle distant ', F10.5,' deg from pole ',
     :          F10.5,' E, ',F10.5,', N; Azimuth ',
     :          F10.5,' to ',F10.5,' E')
        DAZIM=(AZIMS(2)-AZIMS(1))/FLOAT(NLINE-1)
        DO 300 J=1,NLINE
          AZIM=AZIMS(1)+FLOAT(J-1)*DAZIM
          ITV=1
          IF((AZIM.GT.PI).OR.(AZIM.LT.-PI))ITV=-1
          XVA=SDIST*CPOLATR*DCOS(AZIM)+CDIST*SPOLATR
          XVLATR=DASIN(XVA)
          XVB=(CDIST-SPOLATR*XVA)/(CPOLATR*DCOS(XVLATR))
          XVLONR=ITV*DACOS(XVB)+POLONR
          XVLON=XVLONR/D2R
          XVLAT=XVLATR/D2R
          WRITE(LUW,10050)XVLON,XVLAT
  300 CONTINUE
C
C   for a circle on the Cartesian plane
C
      ELSEIF(ICRC.EQ.'CIRC')THEN
        DIST=XPAR(3)
        AZIMS(1)=XPAR(4)*D2R
        AZIMS(2)=XPAR(5)*D2R
        DAZIM=(AZIMS(2)-AZIMS(1))/FLOAT(NLINE-1)
        WRITE(LUW,10105)XPAR(3),XPAR(1),XPAR(2),XPAR(4),XPAR(5)
10105   FORMAT('Cartesian circle radius ', F10.5,' from centre ',
     :          F10.5,' , ',F10.5,', between angles ',
     :          F10.5,' to ',F10.5,' from x = 0')
        DO 400 J=1,NLINE
          DAZ=DAZIM*FLOAT(J-1)
          XV=XPAR(1)+XPAR(3)*COS(DAZ)
          YV=XPAR(2)+XPAR(3)*SIN(DAZ)
          WRITE(LUW,10050)XV,YV
  400   CONTINUE
C
C   for a line on the Cartesian plane
C
      ELSEIF(ICRC.EQ.'LINE')THEN
        FRL=1./FLOAT(NLINE-1)
        DX=(XPAR(3)-XPAR(1))*FRL
        DY=(XPAR(4)-XPAR(2))*FRL
        WRITE(LUW,10106)XPAR(1),XPAR(2),XPAR(3),XPAR(4)
10106   FORMAT('Cartesian line between (',F10.5,' , ',
     :          F10.5,' ) and (',F10.5,', ',F10.5,')')
        DO 500 J=1,NLINE
          XV=XPAR(1)+DX*FLOAT(J-1)
          YV=XPAR(2)+DY*FLOAT(J-1)
          WRITE(LUW,10050)XV,YV
  500   CONTINUE
      ENDIF
      GO TO 20                  ! return for next input line
C
   30 CONTINUE
      STOP
      END
