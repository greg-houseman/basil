C
C   program COROTATE applies a rotation defined by pole coordinates
C   (lon,lat) and angle (degrees) to a set of locations (lon,lat)
C   The rotation here is specified by providing two locations
C   for which the great circle is rotated to the equator
C   also enabled: rotation of surface vectors and covariance matrices
C
C   program written by Greg Houseman, Dec 2015, revised August 2020
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER FNAME*64
      CHARACTER COMM*4
      CHARACTER LINE*140,TLINE*140,RLINE*140
      CHARACTER SITE*8,FORMT*8,FORMIT*6
      DIMENSION XPAR(10),RPAR(10)
      DIMENSION NWORD0(30),NWORD1(30)
      DIMENSION ROTIX(3,3)
      LLN=140
C
      WRITE(6,10111)
10111 FORMAT(' Program COROTATE is controlled by the parameters ',
     :'provided in the file ''corotate_command''')
      PI=3.14159265358979323846D0
      DEGRAD=PI/180.0D0
      EARTH=6371.0   ! earth radius in km
      NPAR=2
      NSITE=0
      DO K=1,8
        K1=(K-1)*8
        FNAME(K1+1:K1+8)='        '
      ENDDO
C
C    open and read command line
C
      OPEN(3,FILE='corotate_command')
C
C    this program allows two modes of operation:
C    COMM is a 4-byte word = 2PTS or UNDO that determines mode
C    PAR1 to PAR4 are parameters used in either mode, see below
C    FNAME is character variable giving name of data file to convert
C    NPAR is the number of numbers on each line of the data file
C    if there are more entries on a line they are treated as text
C    if there are fewer entries on a line, run will fail.
C
      READ(3,*)COMM,PAR1,PAR2,PAR3,PAR4,FNAME,NSITE,NPAR
      IF(NSITE.GT.8)NSITE=8       ! bytes allowed for site name
      CLOSE(3)
      NLEN=0
      DO K=1,64
        IF(FNAME(K:K).NE.' ')NLEN=K
      ENDDO
      WRITE(6,*)'Input file is ',FNAME(1:NLEN),'; Output file is ',
     :           FNAME(1:NLEN)//'.rotated'
C
C    In '2PTS mode:
C    you provide two points which will define the new equator
C    P1 = A longitude, P2 = A latitude
C    P3 = B longitude, P4 = B latitude
C    FNAME is the name of the input data file, NPAR is number of parameters
C    per line and the transformation involves:
C    1. determination of the pole to the great circle of the 2 points
C    2. rotation to bring that pole into coincidence with the N pole (PROT)
C    3. a rotation about the N pole by angle DPHI.
C
C    program can be tested by:
C    1. apply 2PTS mode to input FNAME to determine PLON,PLAT,PROT,DPHI
C    2. apply UNDO mode to output of first step, using parameters output
C       in preceding 2PTS application.
C    note that corotate.out will need to be renamed between steps 1 and 2
C
      IF(COMM.EQ.'2PTS')THEN
        WRITE(6,120)PAR1,PAR2,PAR3,PAR4
  120   FORMAT(' New equator defined by p1 =',2F9.3,'  p2 =',
     :         2F9.3)
        CALL GETPOLE(PAR1,PAR2,PAR3,PAR4,PLON,PLAT,PROT)
C
C    In 'UNDO' mode:
C    you provide rotation parameters in the command line
C    that have been previously obtained in a '2PTS' application
C    P1 = pole longitude, P2 = pole latitude (=0)
C    P3 = rotation angle, P4 = longitude shift
C    FNAME is the name of the input data file, NPAR is number of parameters
C    per line and the transformation involves:
C    1. a rotation about the N pole by angle DPHI, followed by
C    2. a rotation about (PLON,PLAT) by angle PROT
C
      ELSE IF(COMM.EQ.'UNDO')THEN
        PLON=PAR1
        PLAT=PAR2
        PROT=PAR3
        DPHI=PAR4
        WRITE(6,122)PLON,PLAT,PROT,DPHI
  122   FORMAT('UNDO transformation applied: PLON, PLAT, PROT, DPHI = ',
     :         4F10.4)
      ELSE     ! value of COMM
        WRITE(6,10115)
10115   FORMAT('The rotate_command must begin with UNDO or 2PTS')
        STOP
      ENDIF
C
C     Set up the matrix for the rotation, verify action on two points
C     and determine DPHI for use in inverse operation
C
      CALL MATSET(PLON,PROT,ROTIX)
      IF(COMM.EQ.'2PTS')THEN
C       CALL ROTATE(ROTIX,PAR1,PAR2,DPHI,TLA,0.d0,0.d0,VENEW,VNNEW,0)
C       CALL ROTATE(ROTIX,PAR3,PAR4,TOB,TLB,0.d0,0.d0,VENEW,VNNEW,0)
        XPAR(1)=PAR1
        XPAR(2)=PAR2
        CALL ROTATE(ROTIX,XPAR,RPAR,2,SITE)
        DPHI=RPAR(1)
        TLA=RPAR(2)
        XPAR(1)=PAR3
        XPAR(2)=PAR4
        CALL ROTATE(ROTIX,XPAR,RPAR,2,SITE)
        TLB=RPAR(2)
C
C     latitude of rotated points should be zero
C
        IF((DABS(TLA).GT.1.d-3).OR.(DABS(TLB).GT.1.d-3))THEN
          WRITE(6,10113)TLA,TLB
10113     FORMAT('TLA = ',F10.5,' TLB = ',F10.5, 
     :           ' but should be near zero; check!')
          STOP
        ENDIF
        WRITE(6,121)PLON,PLAT,-PROT,DPHI
  121   FORMAT(' For REVERSE transformation: PLON, PLAT, PROT, DPHI = ',
     :         4F10.4)
      ENDIF
C
C     Open the data file and read the data one line at a time
C     data line should begin with long, lat and follow with
C     other parameters as needed. Program expects NPAR numbers
C     followed by any relevant text.
C     
      OPEN(8,FILE=FNAME(1:NLEN))
      OPEN(9,FILE=FNAME(1:NLEN)//'.rotated')
      MNP=0
      MNQ=0
C
C      read the next line, check length of line (bytes) and no. of words
C
   50 CONTINUE
        DO K=1,LLN
          LINE(K:K)=' '
        ENDDO
        TLINE=LINE
        RLINE=LINE
        READ(8,'(A140)',END=55,ERR=52)LINE
C       WRITE(*,*)MNP,LINE
C
C      identify separate words on a line
C
        LENG=0
        NWORD=0
        IF(LINE(1:1).NE.' ')THEN
          NWORD=1
          NWORD0(NWORD)=1
        ENDIF
        DO K=2,LLN
C
C     end of a word
C
          IF((LINE(K:K).EQ.' ').AND.(LINE(K-1:K-1).NE.' '))THEN
            NWORD1(NWORD)=K-1
            LENG=K-1
C
C     start of a new word
C
          ELSE IF((LINE(K:K).NE.' ').AND.(LINE(K-1:K-1).EQ.' '))THEN
            NWORD=NWORD+1
            IF(NWORD.GT.30)THEN
              WRITE(*,*)'Insufficient space defined for NWORD0, NWORD1'
              STOP
            ENDIF
            NWORD0(NWORD)=K
          ENDIF
        ENDDO
C
C      if no content or too few words ignore line, assumed header
C
        SITE='        '
        IF((LENG.EQ.0).OR.(NWORD.LT.NPAR))GO TO 50
C
C      read sitename if NSITE is non zero
C
        IF(NSITE.GT.0)THEN
          KW1=NWORD0(NSITE)
          KW2=NWORD1(NSITE)
          IF(KW1.EQ.KW2)GO TO 51
          READ(LINE(KW1:KW2),*,ERR=51)SITE(1:8)
        ENDIF
C
C      read the numbers
C
        N1=1
        IF(NSITE.EQ.1)N1=2
        KW1=NWORD0(N1)
        KW2=NWORD1(N1+NPAR-1)
C       WRITE(*,*)'KW1,KW2,LINE(KW1:KW2)=',KW1,KW2,LINE(KW1:KW2)
        IF(KW1.EQ.KW2)GO TO 51
        READ(LINE(KW1:KW2),*,ERR=51)(XPAR(J),J=1,NPAR)
C
C     read whatever else is on line as text
C
        NPEND=NWORD1(NPAR)
        IF(NSITE.NE.0)NPEND=NWORD1(NPAR+1)
        LENT=LENG-NPEND
        IF(LENT.GT.0)THEN
          READ(LINE(NPEND+1:LENG),*,ERR=51)TLINE(1:LENT)
        ENDIF
C
C     if reading successful increment the data count
C 
        MNP=MNP+1
C
C   Apply the rotation transformation to one point at a time
C   if NPAR.gt.2 also rotate surface vectors in (XPAR(3),XPAR(4))
C   the sequence of 2 rotations is reversed for 2PTS and UNDO
C
        IF(COMM.EQ.'UNDO')XPAR(1)=XPAR(1)+DPHI
        CALL ROTATE(ROTIX,XPAR,RPAR,NPAR,SITE)
        IF(COMM.EQ.'2PTS')RPAR(1)=RPAR(1)-DPHI
C
C   reconstitute the output line (coordinates in degrees)
C
        LEND=0
        IF(MNP.EQ.1)WRITE(FORMT(1:8),124)NPAR
  124   FORMAT('(',I1,'F12.5)')
        IF(NSITE.EQ.1)THEN
          RLINE(1:8)=SITE
          LEND=8
        ENDIF
        MEND=12*NPAR+LEND
        WRITE(RLINE(LEND+1:MEND),FORMT(1:8))(RPAR(J),J=1,NPAR)
        IF(NSITE.EQ.NPAR+1)THEN
          RLINE(MEND+2:MEND+9)=SITE(1:8)
          MEND=MEND+9
        ENDIF
        IF(LENT.GT.0)THEN
          RLINE(MEND+2:MEND+LENT+1)=TLINE(1:LENT)
        ENDIF
        WRITE(FORMIT(1:6),125)MEND+LENT+1
  125   FORMAT('(A',I3.3,')')
        WRITE(9,FORMIT(1:6))RLINE(1:MEND+LENT)
        GO TO 50      ! go for next line
C
C    line appears to not have numerical format; treat as header
C
   51   CONTINUE
          WRITE(9,*)LINE(1:LENG)
          MNQ=MNQ+1
          IF(MNQ.GT.10)THEN
            WRITE(*,*)'More than 5 line errors, stopping read'
            STOP
          ENDIF
          GO TO 50
C
C   end of file reached
C
   55 CONTINUE
      WRITE(6,100)MNP,NPAR
  100 FORMAT(' Rotation has been applied to ',I8,' data points, ',
     :       'with ',I2,' parameters per line')
      STOP
C
   52 CONTINUE
      WRITE(*,*)'Problem encountered reading data file, NP =',NP
      STOP
      END
      SUBROUTINE ROTATE(ROTIX,XPAR,RPAR,NPAR,SITE)
C
C    routine ROTATE: computes new coordinates in RPAR of a point 
C    on the surface of the sphere provided in XPAR, using a rotation
C    transformation defined by matrix ROTIX, defined in MATSET.
C    this routine replaces XROTATE.
C
C    Input Coordinates and rotation in degrees
C    Assumed order of parameters: lon,lat,ve,vn,se,sn,cen
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER SITE*8
      DIMENSION ROTIX(3,3),RHS(3),RES(3)
      DIMENSION XPAR(NPAR),RPAR(NPAR)
C
C   constant factors
C
      PI=3.14159265358979323846D0
      DEGRAD=PI/180.D0
      DO I=1,NPAR
        RPAR(I)=0.d0
      ENDDO
C
C   evaluate the right hand side using the long/lat of the point
C   and effect the rotation
C
      GLA=XPAR(2)*DEGRAD
      GLO=XPAR(1)*DEGRAD
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
C    get new latitude and longitude from RES(3)
C
      SLAT=DASIN(RES(3))/DEGRAD
      SLON=DATAN2(RES(2),RES(1))/DEGRAD
      IF(SLON.LT.-180.D0)SLON=SLON+360.D0
      IF(SLON.GT.180.D0)SLON=SLON-360.D0
      RPAR(1)=SLON
      RPAR(2)=SLAT
      IF(NPAR.LE.2)RETURN
C
C   rotate vector whose E and N components are XPAR(3), XPAR(4)
C
      CALL VELROT(ROTIX,RES,GLA,GLO,XPAR(3),XPAR(4),RPAR(3),RPAR(4))
      IF(NPAR.LE.4)RETURN
C
C    get the eigenvectors of the covariance matrix
C
      AEE=XPAR(5)*XPAR(5)
      CNN=XPAR(6)*XPAR(6)
      BEN=DSIGN(1.d0,XPAR(7))*XPAR(7)*XPAR(7)
C
      IF((BEN*BEN).LT.0.5d0*(AEE*CNN))THEN
        CALL DSYEV2(AEE,BEN,CNN,EV1,EV2,CSV,SNV)
C
C    negative eigenvalues should not be present; they indicate a
C    problem with data source, but warning message and ad hoc change
C    to covariance enable continuation of calculation without NaNs.
C
      ELSE
        BEN=0.d0
        CALL DSYEV2(AEE,BEN,CNN,EV1,EV2,CSV,SNV)
        WRITE(*,101)SITE,XPAR(7)
101     FORMAT('For site: ',A8,'  covariance of ',G14.5,'  set to zero')
      ENDIF
      IF(EV1*EV2.LT.0.d0)THEN
        WRITE(*,102)SITE,AEE,BEN,CNN
102     FORMAT('Negative eigenvalue for ',A8,'  Cd = ',3G12.5)
        STOP
      ENDIF
C
C    rotate orthogonal eigenvectors (CSV,SNV)^T and (-SNV,CSV)^T
C
      CALL VELROT(ROTIX,RES,GLA,GLO,CSV,SNV,CSVR,SNVR)
C
C     [AEE  BEN]    [CSVR  -SNVR][EV1  0 ][ CSVR  SNVR]
C     [BEN  CNN] =  [SNVR   CSVR][ 0  EV2][-SNVR  CSVR]
C
C    reconstitute the covariance matrix; round-off error shows
C    first in covariance after rotate.unrotate, due to EV1-EV2.
C
      AEE=EV1*CSVR*CSVR+EV2*SNVR*SNVR
      CNN=EV1*SNVR*SNVR+EV2*CSVR*CSVR
      BEN=(EV1-EV2)*CSVR*SNVR
C
C    DABS used in next two lines to deal with negative eigenvalues
C    that appear when BEN^2 > AEE*CNN in original covariance matrix
C
      RPAR(5)=DSQRT(DABS(AEE))
      RPAR(6)=DSQRT(DABS(CNN))
      RPAR(7)=DSIGN(1.d0,BEN)*DSQRT(DABS(BEN))
      RETURN
      END
C
      SUBROUTINE VELROT(ROTIX,RESB,GLA,GLO,VECE,VECN,VENEW,VNNEW)
C
C    routine to compute the components of a vector subject to
C    rotation described by the ROTIX matrix
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ROTIX(3,3),RESB(3),RESC(3),RHS(3)
      VAMP=DSQRT(VECN*VECN+VECE*VECE)
      IF(VAMP.EQ.0.d0)RETURN
C
C   get sines and cosines of lat and long for point B'
C
      SGLATP=RESB(3)
      CGLATP=DSQRT(RESB(1)*RESB(1)+RESB(2)*RESB(2))
      CGLONP=RESB(1)/CGLATP
      SGLONP=RESB(2)/CGLATP
C
C   get sines and cosines of lat and long for point C
C
      CLATG=DCOS(GLA)
      SLATG=DSIN(GLA)
      SLATH=CLATG*VECN/VAMP
      CLATH=DSQRT(1.d0 - SLATH*SLATH)
      CLOND=-(SLATG*SLATH)/(CLATG*CLATH)
      SLOND=VECE/(VAMP*CLATH)
      HLO=GLO+DATAN2(SLOND,CLOND)
      CLONH=DCOS(HLO)
      SLONH=DSIN(HLO)
C
C   compute location of C' by rotation of C
C
      RHS(1)=CLATH*CLONH
      RHS(2)=CLATH*SLONH
      RHS(3)=SLATH
      DO K=1,3
        RESC(K)=0.d0
        DO J=1,3
          RESC(K)=RESC(K)+RHS(J)*ROTIX(K,J)
        ENDDO
      ENDDO
C
C    get sines and cosines of lat and long for point C'
C    and then velocity components in rotated frame
C
      CHLATP=DSQRT(RESC(1)*RESC(1)+RESC(2)*RESC(2))
      CHLONP=RESC(1)/CHLATP
      SHLONP=RESC(2)/CHLATP
      VNNEW=VAMP*RESC(3)/CGLATP
      VENEW=VAMP*CHLATP*(SHLONP*CGLONP-CHLONP*SGLONP)
      RETURN
      END
C
      SUBROUTINE XROTATE(PLON,PLAT,PROT,RLON,RLAT,SLON,SLAT,
     :                  VECE,VECN,VENEW,VNNEW,IV)
C
C    routine XROTATE: computes new coordinates (SLON,SLAT) of a point 
C    (RLON,RLAT) on the surface of the sphere, given a rotation by 
C    angle PROT about pole coordinates (PLON,PLAT).  
C    Unlike routine ROTATE, which uses matrix multiply, and requires
C    that pole of rotation is located on the equator, XROTATE can be 
C    used for pole and target anywhere on the sphere.
C    Note: Input Coordinates and rotation in degrees
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C   constant factors
C
      PI=3.14159265358979323846D0
      DEGRAD=PI/180.D0
C
C     convert degrees to radians
C
      PLA=PLAT*DEGRAD
      PLO=PLON*DEGRAD
      RLA=RLAT*DEGRAD
      RLO=RLON*DEGRAD
      PRO=PROT*DEGRAD
C
C    distance between P and R (in radians)
C
      CD=DCOS(PLA)*DCOS(RLA)*DCOS(PLO-RLO) + DSIN(PLA)*DSIN(RLA)
      DELTA=DACOS(CD)
C
C    azimuth of R at P
C
      CA=1.0D0
      DENOM=DCOS(PLA)*DSIN(DELTA)
      IF(DENOM.NE.0.0)CA=(DSIN(RLA) - DSIN(PLA)*CD)/DENOM
C     WRITE(*,*)'DENOM',DENOM,'CA',CA,'PLA',PLA,'RLA',RLA
      DLON=RLO-PLO
      IF(DLON.LT.-PI)DLON=DLON+2.D0*PI
      IF(DLON.GT. PI)DLON=DLON-2.D0*PI
      IF(CA.GT.1.0D0)CA=1.0D0
      IF(CA.LT.-1.0D0)CA=-1.0D0
      IF((0.0.LE.DLON).AND.(DLON.LE.PI))THEN
        AZIM=DACOS(CA)
      ELSE
        AZIM=-DACOS(CA)
      END IF
C
C   a positive rotation (anticlockwise) is subtracted from AZIM
C
      AZIM=AZIM-PRO
C
C   compute the latitude of rotated R (now S)
C
      SINLAT=DSIN(PLA)*CD+DCOS(PLA)*DSIN(DELTA)*DCOS(AZIM)
      IF(SINLAT.GT.1.0D0)SINLAT=1.0D0
      IF(SINLAT.LT.-1.0D0)SINLAT=-1.0D0
      SLA=DASIN(SINLAT)
C
C   compute the longitude of rotated R (now S)
C
      COSLC=(CD-DSIN(PLA)*SINLAT)/(DCOS(PLA)*DCOS(SLA))
      IF(COSLC.GT.1.0D0)COSLC=1.0D0
      IF(COSLC.LT.-1.0D0)COSLC=-1.0D0
      CHLON=DACOS(COSLC)
      IF((AZIM.GT.0.D0).AND.(AZIM.LE.PI))THEN
        SLO=PLO+CHLON
      ELSE
        SLO=PLO-CHLON
      ENDIF
C
C   convert latitude and longitude of rotated point to degrees
C
      SLON=SLO/DEGRAD
      SLAT=SLA/DEGRAD
      IF(SLON.LT.-180.D0)SLON=SLON+360.D0
      IF(SLON.GT.180.D0)SLON=SLON-360.D0
C     
      RETURN
      END
      SUBROUTINE GETPOLE(ALON,ALAT,BLON,BLAT,PLON,PLAT,PROT)
C
C    routine to determine longitude of equatorial pole of rotation 
C    and angle of rotation given two points that define a great 
C    circle which is to be rotated into coincidence with the equator.  
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C   constant factors
C
      PI=3.14159265358979323846D0
      DEGRAD=PI/180.D0
C
C     convert degrees to radians and get sines and cosines
C
      ALA=ALAT*DEGRAD
      ALO=ALON*DEGRAD
      BLA=BLAT*DEGRAD
      BLO=BLON*DEGRAD
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
        PLON=ALON
        PLAT=0.D0
        PROT=90.D0
      ELSE
        PX=(-A22*SLATA+A12*SLATB)/DETM
        PY=( A21*SLATA-A11*SLATB)/DETM
C
C     compute pole coordinates and rotation from cross product pXN
C
        SINOM=DSQRT((PX*PX+PY*PY)/(1.D0+PX*PX+PY*PY))
        PROT=(DASIN(SINOM))/DEGRAD
        PLON=DATAN2(-PX,PY)/DEGRAD
        IF(PLON.LT.-180.D0)PLON=PLON+360.D0
        IF(PLON.GT.180.D0)PLON=PLON-360.D0
      ENDIF
      PLAT=0.D0
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
