C
C    program polyfix to modify an existing basic polyfile by adding
C    in a set of points where gps data are available
C
C    written by: G. Houseman, August 2015. amended March 2022
C
      PARAMETER (NBND=1000,NDAT=10000)
      CHARACTER STRING*120,STRINGN*120
      CHARACTER FNAMED*64,FNAMEB*64
      CHARACTER SITE*5
      CHARACTER(20) THRK
      DOUBLE PRECISION XLON,YLAT,XC,YC,XSEG,YSEG,DEL1,AZIM,DMIN
      DOUBLE PRECISION XEXT,YEXT,XLO,YLA
      DIMENSION XLON(NDAT),YLAT(NDAT),UEAST(NDAT),UNORTH(NDAT)
      DIMENSION NIX(NDAT),LINEX(NDAT)
      DIMENSION XC(NBND),YC(NBND),XEXT(NBND),YEXT(NBND)
      DIMENSION IND(NBND),JND(NBND),JC1(NBND),JC2(NBND)
      DIMENSION XSEG(2),YSEG(2)
C
C    100 m error on location is considered the same, consistent
C    with locations written with 3 decimal places in degrees
C
      THRES2=0.2              ! (km) threshold for location match
      THRES1=10.0             ! (km) threshold for mesh proximity
      DEG2KM=111.11111
C
C    usage: polyfix filename thres1 boundary
C    where filename is the file containing GPS data (default: selvect.out)
C    thres1 is the threshold (km) for removing points too close together (default 10)
C    boundary is the name of the boundary file used to exclude data (default: boundary)
C
      DO K=1,64
        FNAMED(K:K) =' '
      ENDDO
      FNAMEB=FNAMED
      DO K=1,NDAT
        LINEX(K)=0
      ENDDO
C
C    read data file name if provided
C
      CALL GETARG(1,FNAMED)
      JFND=0
        DO J=1,64
          IF(FNAMED(J:J).NE.' ')JFND=J
        ENDDO
      IF(JFND.EQ.0)THEN                ! default file name
        JFND=11
        FNAMED(1:JFND)='selvect.out'
      ENDIF
      OPEN(8,FILE=FNAMED(1:JFND),STATUS='OLD',ERR=1500)
C
C    read threshold (km) if provided
C
      CALL GETARG(2,THRK)
      IF(THRK.NE."")READ(THRK,*)THRES1
      THRESH1=THRES1/DEG2KM
      THRESH2=THRES2/DEG2KM
C
C    read boundary filename if provided
C
      CALL GETARG(3,FNAMEB)
      JFNB=0
        DO J=1,64
          IF(FNAMEB(J:J).NE.' ')JFNB=J
        ENDDO
      IF(JFNB.EQ.0)THEN                ! default file name
        JFNB=8
        FNAMEB(1:JFNB)='boundary'
      ENDIF
C
C    first read the boundary file and use it to exclude points outside the
C      boundary when data are read.
C
      OPEN(4,FILE=FNAMEB(1:JFNB),STATUS='OLD')
      NEXT=0
      DMIN=THRES2/DEG2KM
  110 CONTINUE
      STRING=STRINGN
      READ(4,'(A120)',END=120)STRING
      IF(STRING(1:1).NE.'#')THEN
        NEXT=NEXT+1
        IF(NEXT.GT.NBND)THEN
          WRITE(6,*)'Insufficient memory, increase NBND and recompile'
          STOP
        ENDIF
        READ(STRING(1:120),*,END=120)XEXT(NEXT),YEXT(NEXT)
        IF(NEXT.GT.1)THEN
          DEL1=ABS(XEXT(NEXT)-XEXT(NEXT-1))
          DEL2=ABS(YEXT(NEXT)-YEXT(NEXT-1))
          IF(DEL2.GT.DEL1)DEL1=DEL2
          IF(DEL1.LT.0.00099)NEXT=NEXT-1   ! assumed to be a repeat point
        ENDIF
      ENDIF
      GO TO 110
  120 CONTINUE
      CLOSE(4)
C PCE      write(6,'(a20)')FNAMEB(1:JFNB)
C
C    expect last point to duplicate first point
C
      DEL1=ABS(XEXT(NEXT)-XEXT(1))
      DEL2=ABS(YEXT(NEXT)-YEXT(1))
      IF(DEL2.GT.DEL1)DEL1=DEL2
      IF(DEL1.LT.0.00099)NEXT=NEXT-1
      WRITE(*,*)'External boundary used to exclude data has ',NEXT,
     :          'points'
C
C    next read the set of points to be included
C
      NADD=0
      NOTIN=0
      NLIN=0
   10 CONTINUE
      STRING=STRINGN
      READ(8,'(A120)',END=20)STRING
      IF(STRING(1:1).NE.'#')THEN 
        NLIN=NLIN+1
        READ(STRING(1:120),*,END=20)XLO,YLA,UE,UN
C PCE        write(6,'(12f12.3)') XLO,YLA,UE,UN
C
C   before saving point, check it lies within the external boundary
C
        CALL INOROUT(XEXT,YEXT,NEXT,XLO,YLA,INXOUT)
C       write(6,*)NEXT,XLO,YLA,INXOUT
        IF(INXOUT.EQ.0)THEN
          NOTIN=NOTIN+1
          GO TO 10
        ENDIF
C
C   if duplicate points present in file, only keep the first
C
        IF(NADD.GT.0)THEN
          DO NDU=1,NADD
            DEL1=ABS(XLON(NDU)-XLO)
            DEL2=ABS(YLAT(NDU)-YLA)
            IF(DEL2.GT.DEL1)DEL1=DEL2
            IF(DEL1.LT.0.00099)THEN   ! < 0.001 degree difference in lat or long
              NOTIN=NOTIN+1
              GO TO 10
            ENDIF
          ENDDO
        ENDIF
C
C    ok, we keep the point
C
        NADD=NADD+1
        IF(NADD.GT.NDAT)THEN
          WRITE(6,*)'Insufficient memory, increase NDAT and recompile'
          STOP
        ENDIF
        XLON(NADD)=XLO
        YLAT(NADD)=YLA
        UEAST(NADD)=UE
        UNORTH(NADD)=UN
        LINEX(NLIN)=NADD
      ENDIF
      GO TO 10
   20 CONTINUE
      WRITE(6,2010)NADD,FNAMED(1:JFND)
 2010 FORMAT(I5,' data points from ',A16,' lie within the boundary')
      WRITE(6,2011)NOTIN
 2011 FORMAT(I5,' points are excluded as external to the region')
      DO K=1,NADD
        NIX(K)=0
      ENDDO
C
C    now open existing poly file and modify header
C
      OPEN(4,FILE='x.poly',STATUS='OLD')
      OPEN(10,FILE='new.poly')
      STRING=STRINGN
      READ(4,'(A120)',END=30)STRING
      JBEG=120
      DO J=2,120
        IF((STRING(J-1:J-1).NE.' ').AND.(STRING(J:J).EQ.' '))THEN
          JEND=J-1
          GO TO 30
        ENDIF
      ENDDO
   30 READ(STRING(1:JEND),*)NPTS
      WRITE(6,2030)NPTS
 2030 FORMAT(I5,' points already present in poly file')
C
C     read and store coordinates of existing points
C
      DO K=1,NPTS
        READ(4,*)IND(K),XC(K),YC(K)
      ENDDO
C
C     now read connection segments
C
      READ(4,*)NCON
      WRITE(6,2031)NCON
 2031 FORMAT(I5,' connecting segments are defined')
      DO K=1,NCON
        READ(4,*,IOSTAT=IOS)JND(K),JC1(K),JC2(K)
        IF(IOS.NE.0)THEN
          WRITE(*,*)'There is an inconsistency in the poly file ',
     :    'between number and list of segments'
          STOP
        ENDIF
      ENDDO
C
C     first identify points that are already present in x.poly
C     arising from 'boundary' or 'internal' files
C
      OPEN(8,FILE='polyfix_ok')
      NUMNIX=0
      DO K=1,NADD
        DMIN=180.D0
        DO L=1,NPTS
          CALL DISTAZ(XC(L),YC(L),XLON(K),YLAT(K),DEL1,AZIM)
          IF(DEL1.LT.DMIN)THEN
            DMIN=DEL1
            LMIN=L
          ENDIF
        ENDDO
        IF(DMIN.LT.THRESH2)THEN
          NIX(K)=1
          NUMNIX=NUMNIX+1
          DO J=1,NDAT
            IF(LINEX(J).EQ.K)LINEX(J)=0
          ENDDO
          WRITE(6,2032)K,XLON(K),YLAT(K)
 2032     FORMAT('Point ',I5,' at ',F12.3,',',F12.3,
     :                 ' already present, not added')
          WRITE(8,2034)XLON(K),YLAT(K),UEAST(K),UNORTH(K),K
        ENDIF
      ENDDO
      CLOSE(8)
      IF(NUMNIX.EQ.0)THEN
        WRITE(6,2035)
 2035   FORMAT('None of the new points are present in the existing '
     :         'poly file')
      ELSE
        WRITE(6,2042)NUMNIX
 2042   FORMAT(I5,' points already present in the polyfile, not added')
      ENDIF
C
C     now check proximity to existing segments
C
      WRITE(6,2044)THRES1
 2044 FORMAT('Points that are closer to existing points than ',F5.1,
     : ' km will not be used')
      OPEN(7,FILE='polyfix_probs')
      WRITE(6,10050)
10050 FORMAT('The following points (see also polyfix_probs) are too ',
     :'close to existing line segments',/'Modified point definition ',
     :'in boundary or internal files may allow inclusion.')
      DO K=1,NADD
        IF(NIX(K).EQ.0)THEN
          DMIN=180.D0
          DO M=1,NCON
            XSEG(1)=XC(JC1(M))
            YSEG(1)=YC(JC1(M))
            XSEG(2)=XC(JC2(M))
            YSEG(2)=YC(JC2(M))
            CALL BDIST(XLON(K),YLAT(K),XSEG,YSEG,2,0,DEL1)
            IF(DEL1.LT.DMIN)THEN
              DMIN=DEL1
              MMIN=M
            ENDIF
          ENDDO
          IF(DMIN.LT.THRESH1)THEN
            NIX(K)=1
            NUMNIX=NUMNIX+1
            DO J=1,NDAT
              IF(LINEX(J).EQ.K)LINEX(J)=0
            ENDDO
            WRITE(6,2033)K,XLON(K),YLAT(K),MMIN,DMIN*DEG2KM
C           WRITE(6,2033)K,XLON(K),YLAT(K),MMIN,XSEG(1),YSEG(1),
C    :                   XSEG(2),YSEG(2)
 2033       FORMAT('Point ',I5,' at ',F12.3,',',F12.3,'  Distance ',
     :           ' to segment ',I5,' is ',F12.3)
 2034       FORMAT(2F12.5,2F12.3,I5)
            WRITE(7,2034)XLON(K),YLAT(K),UEAST(K),UNORTH(K),K
          ENDIF
        ENDIF
      ENDDO
      WRITE(*,*)NUMNIX,' points excluded on proximity'
C
C    now construct modified poly file
C
      NADDED=NADD-NUMNIX
      NNPTS=NPTS+NADDED
      REWIND(4)
      READ(4,'(A120)',END=30)STRING
      KEND=120
      DO K=1,120
        IF(STRING(K:K).NE.' ')KEND=K
      ENDDO
      WRITE(10,*)NNPTS,STRING(JEND+1:KEND)
C
C    read and write the block of nodes already present
C
      DO J=1,NPTS
        STRING=STRINGN
        READ(4,'(A120)')STRING
        JEND=120
        DO K=1,120
          IF(STRING(K:K).NE.' ')KEND=K
        ENDDO
        WRITE(10,*)STRING(1:KEND)
      ENDDO
C
C     now write the new set of nodes
C
      NUMB=NPTS
      DO J=1,NADD
        IF(NIX(J).EQ.0)THEN
          NUMB=NUMB+1
          WRITE(10,2020)NUMB,XLON(J),YLAT(J)
        ENDIF
      ENDDO
 2020 FORMAT(I4,2F10.5)
C
C    finally copy everything else in the file
C
   40 CONTINUE
      STRING=STRINGN
      READ(4,'(A120)',END=50)STRING
      JEND=120
      DO J=1,120
        IF(STRING(J:J).NE.' ')JEND=J
      ENDDO
      WRITE(10,*)STRING(1:JEND)
      GO TO 40
   50 CONTINUE
      WRITE(6,2043)NPTS,NADDED
 2043 FORMAT('File: new.poly now includes ',I5,' original plus',I5,
     :       ' additional points ')
C
C    copy an edited set of data that have satisfied all the criteria
C
      OPEN(8,FILE=FNAMED(1:JFND),STATUS='OLD')     ! original data file
      OPEN(9,FILE='data_used')
      K=0
      NUSED=0
      NTUSED=0
      DO J=1,NDAT
        STRING=STRINGN
        READ(8,'(A120)',END=320)STRING(1:120)
        IF(STRING(1:1).NE.'#')THEN
          K=K+1
          IF(LINEX(K).NE.0)THEN
            NUSED=NUSED+1
            WRITE(9,'(A120)')STRING(1:120)
          ELSE
            NTUSED=NTUSED+1
          ENDIF
        ENDIF
      ENDDO
  320 CONTINUE
      WRITE(*,*)NUSED,' data values have been written to data_used'
      WRITE(*,*)NTUSED,' data values are not used as external or ',
     :         'too close to existing points'
      CLOSE(9)
      CLOSE(8)
      STOP
 1500 CONTINUE
        WRITE(*,*)'Data file: ',FNAMED(1:JFND),' is not present'
        WRITE(*,*)'Usage: polyfix datafname thresh boundaryfname ',
     :          '(all optional)'
      STOP
      END
C
C    routine to examine proximity of a point (XLOC,YLOC) to a boundary
C    (XLON(NPTS), YLAT(NPTS)), used to reject points that are too close.
C
      SUBROUTINE BDIST(XLOC,YLOC,XLON,YLAT,NPTS,ICLOSE,DMINX)
      DOUBLE PRECISION XLOC,YLOC,XLON,YLAT,AZIM,ARE2,XL,YL
      DOUBLE PRECISION DEL1,DEL2,DEL3,CS1,CS2,SNLAT,DMINX
      DIMENSION XLON(NPTS),YLAT(NPTS)
      DMINX=180.D0
      DEG2R=2.D0*DASIN(1.D0)/DMINX
      IF(NPTS.GT.1)THEN
C
C    ICLOSE=1 means the polygon is closed by connecting first and last points
C
        DO J=1,NPTS-1+ICLOSE
          IF((ICLOSE.EQ.1).AND.(J.EQ.NPTS))THEN
            XL=XLON(1)
            YL=YLAT(1)
          ELSE
            XL=XLON(J+1)
            YL=YLAT(J+1)
          ENDIF
C
C     DEL1,2,3 are the sides of the triangle
C  
          CALL DISTAZ(XLOC,YLOC,XLON(J),YLAT(J),DEL1,AZIM)
          CALL DISTAZ(XLOC,YLOC,XL,YL,DEL2,AZIM)
          CALL DISTAZ(XLON(J),YLAT(J),XL,YL,DEL3,AZIM)
          CS1=DEL2*DEL2 - (DEL1*DEL1 + DEL3*DEL3)
          CS2=DEL1*DEL1 - (DEL2*DEL2 + DEL3*DEL3)
C
C     Internal angles are both acute if (XLOC,YLOC) within shadow of segment
C
          IF((CS1*CS2).LT.0.)THEN
            DMIN=MIN(DEL1,DEL2)
C
C     perpendicular distance is 2*AREA/baseline length, units in degrees
C     but assuming locally flat as distances are small
C
          ELSE
            SNLAT=DABS(DSIN(YLOC*DEG2R))
            ARE2=DABS(((XLON(J)-XLOC)*(YL-YLOC)-
     :                 (XL-XLOC)*(YLAT(J)-YLOC)))
            DMIN=SNLAT*ARE2/DEL3
          ENDIF
          IF(DMIN.LT.DMINX)DMINX=DMIN
        ENDDO
      ENDIF
      RETURN
      END
C
C    routine to make distance and azimuth calculations on a sphere
C
      SUBROUTINE DISTAZ(RLONG,RLATI,PLONG,PLATI,DELTA,AZIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C   constant factors
C
      PI=3.14159265358979323846D0
      DEGRAD=PI/180.D0
C
C     convert degrees to radians
C
      PLAT=PLATI*DEGRAD
      PLON=PLONG*DEGRAD
      RLAT=RLATI*DEGRAD
      RLON=RLONG*DEGRAD
C
C    distance between two points (in radians)
C
      CD=DCOS(PLAT)*DCOS(RLAT)*DCOS(PLON-RLON) + DSIN(PLAT)*DSIN(RLAT)
      DELTA=DACOS(CD)
C
C    azimuth of P at R
C
      CA=1.0D0
      DENOM=DCOS(RLAT)*DSIN(DELTA)
      IF(DENOM.NE.0.D0)CA=(DSIN(PLAT) - DSIN(RLAT)*CD)/DENOM
      DLON=PLON-RLON
      IF(DLON.LT.-PI)DLON=DLON+2.D0*PI
      IF(DLON.GT. PI)DLON=DLON-2.D0*PI
      IF(CA.GT.1.D0)CA=1.D0
      IF(CA.LT.-1.D0)CA=-1.D0
      IF((0.D0.LE.DLON).AND.(DLON.LE.PI))THEN
        AZIM=DACOS(CA)
      ELSE
        AZIM=-DACOS(CA)
      END IF
C
C     convert radians to degrees
C
      DELTA=DELTA/DEGRAD
      AZIM=AZIM/DEGRAD
C
      RETURN
      END
C
      SUBROUTINE INOROUT_GREG(XEXT,YEXT,NEXT,XC,YC,INXOUT)
C
C    this routine checks if a point (XC,YC) is within the external
C    perimeter defined by the set of points (XEXT(J),YEXT(J)), (J=1,NEXT)
C  
      DOUBLE PRECISION XEXT,YEXT,XC,YC,XINC,FR,YINT
      DIMENSION XEXT(NEXT),YEXT(NEXT)
C
      NCROSS=0
      DO K=1,NEXT-1
        XINC=XEXT(K+1)-XEXT(K)
        IF(ABS(XINC).GT.0)THEN
          FR=(XC-XEXT(K))/XINC
          IF((FR.GT.0).AND.(FR.LE.1.D0))THEN
            YINT=YEXT(K)+FR*(YEXT(K+1)-YEXT(K))
            IF(YINT.GT.YC)NCROSS=NCROSS+1
          ENDIF
        ENDIF
      ENDDO
      INXOUT=MOD(NCROSS,2)
      RETURN
      END

      SUBROUTINE INOROUT ( x, y, n, x0, y0, m)
!-----------------------------------------------------------------------
! GIVEN A POLYGONAL LINE CONNECTING THE VERTICES (X(I),Y(I)) (I = 1,...,N)
! TAKEN IN THIS ORDER.  IT IS ASSUMED THAT THE POLYGONAL PATH IS A LOOP,
! WHERE (X(N),Y(N)) = (X(1),Y(1)) OR THERE IS AN ARC FROM (X(N),Y(N)) TO
! (X(1),Y(1)).  N.B. The polygon may cross itself any number of times.

! (X0,Y0) IS AN ARBITRARY POINT AND L AND M ARE VARIABLES.
! On output, L AND M ARE ASSIGNED THE FOLLOWING VALUES ...

!    L = -1   IF (X0,Y0) IS OUTSIDE THE POLYGONAL PATH
!    L =  0   IF (X0,Y0) LIES ON THE POLYGONAL PATH
!    L =  1   IF (X0,Y0) IS INSIDE THE POLYGONAL PATH

! M = 0 IF (X0,Y0) IS ON OR OUTSIDE THE PATH.  IF (X0,Y0) IS INSIDE THE
! PATH THEN M IS THE WINDING NUMBER OF THE PATH AROUND THE POINT (X0,Y0).

! Fortran 66 version by A.H. Morris
! Converted to ELF90 compatibility by Alan Miller, 15 February 1997
!Converted back into F77 by Philip England 16 Oct 2007
!-----------------------

      IMPLICIT NONE
      REAL*8 x0, y0, x(n), y(n)
      INTEGER n, l, m

!     Local variables
      INTEGER i, n0
      REAL*8 angle, eps, pi, pi2, sum, theta, theta1, thetai, tol, u, v

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0

      eps = 0.433681d-19
c      eps = EPSILON(1.0)

!-----------------------------------------------------------------------
      n0 = n
      IF (x(1).eq.x(n) .AND. y(1).eq.y(n)) n0 = n - 1
      pi = ATAN2(0.0, -1.0)
      pi2 = 2.0*pi
      tol = 4.0*eps*pi
      l = -1
      m = 0

      u = x(1) - x0
      v = y(1) - y0
      IF (u.eq.0.0 .AND. v.eq.0.0) GO TO 20
      IF (n0.lt.2) RETURN
      theta1 = ATAN2(v, u)

      sum = 0.0
      theta = theta1
      DO i = 2, n0
        u = x(i) - x0
        v = y(i) - y0
        IF (u.eq.0.0 .AND. v.eq.0.0) GO TO 20
        thetai = ATAN2(v, u)
  
        angle = ABS(thetai - theta)
        IF (ABS(angle - pi).lt.tol) GO TO 20
        IF (angle.gt.pi) angle = angle - pi2
        IF (theta.gt.thetai) angle = -angle
        sum = sum + angle
        theta = thetai
c	write(6,'(i8,12g16.4)') i,thetai,angle,theta,theta1
      END DO

      angle = ABS(theta1 - theta)
      IF (ABS(angle - pi).lt.tol) GO TO 20
      IF (angle.gt.pi) angle = angle - pi2
      IF (theta.gt.theta1) angle = -angle
      sum = sum + angle

!     SUM = 2*PI*M WHERE M IS THE WINDING NUMBER

      m = ABS(sum)/pi2 + 0.2
      IF (m.eq.0) RETURN
      l = 1
      IF (sum.lt.0.0) m = -m
      RETURN

!     (X0, Y0) IS ON THE BOUNDARY OF THE PATH

 20   l = 0
      RETURN
      END
