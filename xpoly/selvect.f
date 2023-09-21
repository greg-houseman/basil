C
C    program selvect to read a database of GPS vectors, and select
C    and output based on criteria as follows:
C
C    arg1 is the input data file e.g., 'fname'
C    Output is written to the same fname.selv
C    if arg2 = name, a file called 'names' is expected, with a list of 
C       station names (one per line) used to choose vectors
C       if first line in 'names' is 'not ' then these names are excluded
C    if arg2 = loca, a file called 'boundary' is expected. The file should
C       contain a list of coordinate pairs interpreted as a polynomial
C       (automatically closed, last to first).  Only stations that fall
C       within the polynomial are eligible for selection.
C    arg3 and arg4 are optional parameters that may be used to over-ride
C    the default selection parameters THRESHK and ERRDIS specified below.
C    arg5 is the sequenc number of the site name on each line, e.g. 1
C    if line begins with site name, 0 if there is no site name, 8 if it 
C    follows 7 numbers on the line
C    arg6 is the number of numbers on each line, e.g. 7 if location (2) is
C    followed by velocity (2) is followed by 3 elements of covariance.
C
C    To reduce point density, all pairs of points are examined
C    to find nearest neighbours.  Points are removed in sequence.
C    At each cycle of this process the shortest connection is found
C    and the point with the largest error is then removed, the
C    other point being retained unchanged.  The operation is repeated
C    until no connecting segments < THRESHK (in km) remain.
C
C    A procedure for taking a weighted average of the pair of points
C    was also examined, but was not adopted, as individual points
C    ended up being comprised of > 10 contributing stations and the
C    weighting could not easily be applied in a rigorous manner.
C
C    The aim here is to produce a set of constraint points
C    that can be used as the basis for a basil calculation
C
C    selvect.f is derived from inplatt.f (as used for working with
C    CMT datasets) and uses the same procedure for selecting only
C    points that fall within a specified polygonal external boundary
C
C    written by: G. Houseman, July 2015; revised October 2015
C
      DOUBLE PRECISION DAT,GPSDAT,XLON,XLAT,DELTA,AZIM
      DOUBLE PRECISION XLOC,YLOC,CLON,CLAT
      PARAMETER (NBND=200,NDAT=1000)
      CHARACTER LINE*164,TLINE*128,STRINGN*128,COMMENT*128
      CHARACTER FNAME*64
      CHARACTER SITE*8,NAMES*8,STN*8,FORMT*8
      CHARACTER MODEW*4
      CHARACTER(20) THRK,ERDI,CHNSITE,CHNPAR
      DIMENSION NAMES(NBND)
      DIMENSION XLON(NBND),XLAT(NBND)
      DIMENSION GPSDAT(7,NDAT),DAT(7),STN(NDAT),NPIN(NDAT)
      DIMENSION COMMENT(NDAT),LENCOM(NDAT)
      DIMENSION DISMIN(NDAT),NEARES(NDAT),ERRMAG(NDAT)
      DIMENSION NWORD0(10),NWORD1(10)
      LLN=128
      DEG2KM=111.11111
C
C     set default thresholds.  input assumed in mm/yr, use rescale prior
C     if needed to convert input from m/yr
C     next 2 values are defaults, may be overwritten by command line args
C
      THRESHK=10.D0            ! (km) threshold for proximity
      ERRDIS=1.0               ! (mm/yr) permissible error level
      NSITE=8
      NPAR=7
C
C     arg1 is the name of the input data file, arg2 is either 'name' or 'loca'
C
C     FNAME(1:32)='                                '
C     FNAME(33:64)='                                '
      CALL GETARG(1,FNAME)
      IF(FNAME.EQ."")THEN
         WRITE(*,*)'Usage: selvect datafilename mode min-distance',
     :             ' max-err site-entry npar'
         WRITE(*,*)'where npar is number of real data on line (7); ',
     :        'site-entry is the position of the site name (8) on line'
         WRITE(*,*)"mode may take the value: 'loca' (default) or 'name'"
         WRITE(*,*)"if 'loca', file 'boundary' is required, if 'name',",
     :         " file 'names' is required"
         STOP
      ENDIF
C
C     arg3 (THRESHK) and arg4 (ERRDIS) are optional parameters that overwrite
C     above default values.  Data retained are output to the file name.selv
C
      MODE=1
      MODEW='    '
      CALL GETARG(2,MODEW)
      IF(MODEW.EQ.'name')MODE=0    ! select lines based on file name
      IF(MODEW.EQ.'loca')MODE=1    ! select lines based on location
      CALL GETARG(3,THRK)
      IF(THRK.NE."")READ(THRK,*)THRESHK
      CALL GETARG(4,ERDI)
      IF(ERDI.NE."")READ(ERDI,*)ERRDIS
      WRITE(6,10111)THRESHK,ERRDIS
10111 FORMAT('Selection parameters: THRESHK = ',F8.3,' km, ERRDIS = ',
     :       F8.3,' mm/yr')
      THRESH=THRESHK/DEG2KM ! converted to DEGREES
      CALL GETARG(5,CHNSITE)
      IF(CHNSITE.NE."")READ(CHNSITE,*)NSITE
      IF(NSITE.GT.8)NSITE=8
      CALL GETARG(6,CHNPAR)
      IF(CHNPAR.NE."")READ(CHNPAR,*)NPAR
      WRITE(6,10113)NSITE,NPAR
10113 FORMAT('Data format parameters: NSITE = ',I3,', NPAR = ',I3)
C
C    first read the external boundary file, ignoring repeat / closure
C    points
C
      IF(MODE.EQ.1)THEN
        OPEN(4,FILE='boundary',STATUS='OLD',IOSTAT=IOS)
        IF(IOS.NE.0)THEN
          WRITE(*,*)'Please provide the file "boundary"'
          STOP
        ENDIF
        NPTS=0
        DSQ=0.0
   10   CONTINUE
          READ(4,'(A128)',END=15)LINE(1:128)
          IF(LINE(1:1).NE.'#')THEN
            READ(LINE,*,IOSTAT=IOS)XLN1,XLT1
            IF(IOS.NE.0)THEN
              WRITE(*,*)'Apparent problem reading point ',NPTS,
     :                  ' of boundary file'
              STOP
            ENDIF
            IF(NPTS.GE.1)DSQ=(XLN1-XLON(NPTS))**2 +
     :                       (XLT1-XLAT(NPTS))**2
            IF((NPTS.EQ.0).OR.(DSQ.GT.1.E-6))THEN
              NPTS=NPTS+1
              IF(NPTS.GT.NBND)THEN
                WRITE(6,*)'Recompile with increased memory ',
     :                    'for boundary: NBND'
                STOP
              ENDIF                    ! within available memory
              XLON(NPTS)=XLN1
              XLAT(NPTS)=XLT1
C             WRITE(*,*)'NPTS,XLON,XLAT=',NPTS,XLON(NPTS),XLAT(NPTS),DSQ
            ENDIF                      ! distinct point
          ENDIF                        ! not a comment line
        GOTO 10
   15   CONTINUE
        CLOSE(4)
        DSQ=(XLON(NPTS)-XLON(1))**2 + (XLAT(NPTS)-XLAT(1))**2
        IF(DSQ.LE.1.E-6)NPTS=NPTS-1
        WRITE(6,10102)NPTS
10102   FORMAT('external boundary has been read: ',I5,' points')
C
C    or read the list of names
C
      ELSE IF(MODE.EQ.0)THEN
        OPEN(4,FILE='names',STATUS='OLD',IOSTAT=IOS)
        IF(IOS.NE.0)THEN
          WRITE(*,*)'Please provide the file "names"'
        ENDIF
        NMS=0
        NNEG=0
  110   CONTINUE
        NMS=NMS+1
        IF(NMS.GT.NBND)THEN
          WRITE(6,*)'Provide more memory for NAMES and recompile'
          STOP
        ENDIF
C
C   'not' at the top of the list means exclude these names
C
        NAMES(NMS)='        '
        READ(4,*,END=115)NAMES(NMS)
        IF(NAMES(NMS)(1:3).EQ.'not')THEN
          NMS=NMS-1
          NNEG=1
        ENDIF
C
C   'end' in the list of names, causes later names to be ignored
C
        IF(NAMES(NMS)(1:3).EQ.'end')GO TO 115
        GOTO 110
  115   CONTINUE
        NMS=NMS-1
        WRITE(6,10112)NMS
10112   FORMAT('list of station names has been read: ',I5,' points')
      ENDIF
C
C    initialise string and number of lines retained
C
      DO J=1,LLN
        STRINGN(J:J)=' '
      ENDDO
      DO J=1,NDAT
        COMMENT(J)=STRINGN
      ENDDO
      NACC=0
      NSET=0
      JECT=0
      WMIN=99.
      WMAX=0.0
      SCMIN=1.e37
      SCMAX=0.0
      IEXT=0
      JFND=0
        DO J=1,64
          IF(FNAME(J:J).NE.' ')JFND=J
        ENDDO
      OPEN(8,FILE=FNAME(1:JFND),STATUS='OLD')
C
C    read next line (unit 8); check it can be decoded
C
      MNQ=0
   20 CONTINUE
      LINE=STRINGN
      TLINE=STRINGN
C
C    read each line of data and check what is in it
C
      READ(8,'(A128)',END=30)LINE
      IF(LINE(1:1).EQ.'#')GO TO 20     !  skip comment lines
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
          NWORD0(NWORD)=K
        ENDIF
      ENDDO
C
C      if no content ignore line, if not parsed assume header
C
      SITE='        '
      IF(LENG.EQ.0)GO TO 20          ! blank line
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
      IF(KW1.EQ.KW2)GO TO 51
      READ(LINE(KW1:KW2),*,ERR=51)(DAT(J),J=1,NPAR)
C     WRITE(*,*)SITE(1:8),KW1,KW2,NPAR,(DAT(J),J=1,NPAR)
C
C     read whatever else is on line as text
C
      NPEND=NWORD1(NPAR)
      IF(NSITE.NE.0)NPEND=NWORD1(NPAR+1)
      LENT=LENG-NPEND
      IF(LENT.GT.0)THEN
        READ(LINE(NPEND+1:LENG),*,ERR=51)TLINE(1:LENT)
      ENDIF
      GO TO 52
C
C    line appears to not have numerical format; treat as header/comment
C
   51 CONTINUE
      WRITE(6,*)LINE(1:LENG)
      MNQ=MNQ+1
      IF(MNQ.GT.5)THEN
        WRITE(*,*)'More than 5 line errors, stopping read'
        STOP
      ENDIF
      GO TO 20
C
C    expected format: lon, lat, ue, un, se, sn, corr (mm/yr)
C    site name may be before (1), or after (NPAR+1) values
C
   52 XLOC=DAT(1)
      YLOC=DAT(2)
C
C   point selection based on list of names
C
      IF(MODE.EQ.0)THEN
        DO K=1,NMS
          IF(SITE(1:4).EQ.NAMES(K)(1:4))THEN
            IF(NNEG.EQ.1)GO TO 20
            IF(NNEG.EQ.0)GO TO 21
          ENDIF
        ENDDO
        IF(NNEG.EQ.0)GO TO 20
        IF(NNEG.EQ.1)GO TO 21
C
C    check that line from point to north pole intersects plate boundary
C
      ELSE IF(MODE.EQ.1)THEN
        INTS=0
        XA=XLON(NPTS)
        YA=XLAT(NPTS)
        DO J=1,NPTS
          XB=XLON(J)
          YB=XLAT(J)
          IF((XA.NE.XB).AND.((XLOC-XA)*(XLOC-XB).LE.0.0))THEN
C
C    meridian intersects boundary, check if north or south
C    check on FR is to catch points that fall on a segment join
C
            FR=(XLOC-XA)/(XB-XA)
            YINT=YA+FR*(YB-YA)
            IF((FR.NE.0.0).AND.(YINT.GT.YLOC))INTS=INTS+1
          END IF
          XA=XB
          YA=YB
        ENDDO
C
C    an even value of INTS means that point is outside boundary
C
        IF(MOD(INTS,2).NE.1)THEN
          WRITE(6,10104)SITE,XLOC,YLOC
10104     FORMAT(A5,2F12.5,' external, location not saved')
          IEXT=IEXT+1
          GO TO 20
        ENDIF
C
C     now check nearest approach to external boundary
C
        CALL BDIST(XLOC,YLOC,XLON,XLAT,NPTS,1,DMIN)
        IF(DMIN.LT.THRESH)THEN
          WRITE(6,10105)SITE,XLOC,YLOC
10105     FORMAT(A5,2F12.5,' too close to boundary, location not ',
     :           'saved.')
          IEXT=IEXT+1
          GO TO 20
        ENDIF
C
C     check error magnitude, reject if too large
C
        ERRG=0.5*SQRT(DAT(5)*DAT(5)+DAT(6)*DAT(6)+2.0*DAT(7)*DAT(7))
        IF(ERRG.GT.ERRDIS)THEN
          JECT=JECT+1
          WRITE(*,100)SITE,ERRG
  100     FORMAT('Rejecting record for ',A5,' with error magnitude:',
     :           F12.5,' mm/yr')
          GO TO 20
        ENDIF
      ENDIF
C
C     after excluding data values for various reasons, if we have
C     reached this point, we save the current value
C
   21 NACC=NACC+1
      IF(NACC.GT.NDAT)THEN
        WRITE(*,*)'Increase data memory (NDAT)'
        STOP
      ENDIF
      DO J=1,7
        GPSDAT(J,NACC)=DAT(J)
      ENDDO
      STN(NACC)=SITE
      COMMENT(NACC)(1:LENT)=TLINE(1:LENT)
      LENCOM(NACC)=LENT
      ERRMAG(NACC)=ERRG
      NPIN(NACC)=1
C
C    go and read the next line
C
      GO TO 20
C
C     reading of file complete; next section deals with removal of duplicates etc
C
   30 CONTINUE
      WRITE(*,101)NACC,JECT,IEXT
  101 FORMAT('Accepted data for ',I5,' points within region. ',/,I5,
     :       ' points rejected for precision, ',/,I5,' points ',
     :       'rejected as external',/,'Points rejected for location ',
     :       'might be retained by modifying boundary file and ',
     :       're-running xpoly and selvect')
      DO J=1,NACC
        NEARES(J)=0
        DISMIN(J)=181.0
      ENDDO
C
C     for each point find its nearest neighbour
C
      DO J=1,NACC
        CLON=GPSDAT(1,J)
        CLAT=GPSDAT(2,J)
        CALL CLOSER(CLON,CLAT,J,GPSDAT,NPIN,NACC,DISMIN(J),NEARES(J))
C       IF(DISMIN(J).LE.0.1)WRITE(6,103)'first   ',J,NEARES(J),DISMIN(J)
 103    FORMAT(A8,' J =',I5,'  K =',I5,'  dist =',G13.5)
      ENDDO
C
C     we will now examine points two at a time, based on the criterion of 
C     minimum separation.  At each step, we flag the point with the 
C     greater error to be henceforth ignored.  The other point is 
C     retained and for any points formerly closest to the first point, 
C     their new nearest neighbours are evaluated
C
      DO L=1,NACC
        DISM=181.0
        DO J=1,NACC
          IF((DISMIN(J).LT.DISM).AND.(NPIN(J).NE.0))THEN
            DISM=DISMIN(J)
            J1=J
            J2=NEARES(J)
          ENDIF
        ENDDO
        IF(DISM.GT.THRESH)GO TO 40      !   minimum separation retained
        E2D1=ERRMAG(J1)
        E2D2=ERRMAG(J2)
C       WGHT1=E2D2/(E2D2+E2D1)
C       WGHT2=E2D1/(E2D1+E2D2)
C
C    replace J1 with the weighted average, J2 annulled
C
C       DO K=1,7
C         GPSDAT(K,J1)=WGHT1*GPSDAT(K,J1)+WGHT2*GPSDAT(K,J2)
C       ENDDO
C       NPIN(J1)=NPIN(J1)+NPIN(J2)
C       NPIN(J2)=0
C       IF(WGHT2.GT.WGHT1)STN(J1)=STN(J2)
C
C    instead of modifying the stronger point, retain it unchanged
C
        IF(E2D2.GT.E2D1)THEN
          NPIN(J2)=0
          JCUT=J2
          JRET=J1
        ELSE
          NPIN(J1)=0
          JCUT=J1
          JRET=J2
        ENDIF
        WRITE(6,10101)STN(JCUT),STN(JRET),DISM*111.111,STN(JCUT)
10101   FORMAT('distance: ',A4,' to ',A4,' = ',F12.5,' removing ',A4)
C
C     find which is new nearest neighbour for JRET
C
        CLON=GPSDAT(1,JRET)
        CLAT=GPSDAT(2,JRET)
        CALL CLOSER(CLON,CLAT,JRET,GPSDAT,NPIN,NACC,DISMIN(JRET),
     :              NEARES(JRET))
C
C     JCUT is no more; check it isn't still paired with other points
C
        DO K=1,NACC
          IF(NEARES(K).EQ.JCUT)THEN
C
C      if it is, reassess those points
C
            CLON=GPSDAT(1,K)
            CLAT=GPSDAT(2,K)
            CALL CLOSER(CLON,CLAT,K,GPSDAT,NPIN,NACC,
     :                  DISMIN(K),NEARES(K))
C           WRITE(6,103)'second  ',K,NEARES(K),DISMIN(K)
          ENDIF
        ENDDO
C
C     repeat the point elimination process for the next nearest pair
C
      ENDDO
   40 CONTINUE
C
C     assess summary of removal operations. check all the removed 
C     points to determine min. distance to points that remain
C
      NKEPT=0
      NREM=0
C     PTDIS=180.0
      DO J=1,NACC
        IF(NPIN(J).EQ.0)THEN
          NREM=NREM+1
          CLON=GPSDAT(1,J)
          CLAT=GPSDAT(2,J)
          CALL CLOSER(CLON,CLAT,J,GPSDAT,NPIN,NACC,
     :                  DISMIN(J),NEARES(J))
C         IF(DISMIN(J).LT.PTDIS)PTDIS=DISMIN(J)
C         IF(DISMIN(J).GT.0.17)THEN
C           WRITE(*,*)'distance = ',DISMIN(J)
C           WRITE(*,*)(GPSDAT(K,J),K=1,7),STN(J),' removed'
C           WRITE(*,*)(GPSDAT(K,NEARES(J)),K=1,7),STN(NEARES(J)),
C    :                ' retained'
C         ENDIF
        ELSE
          NKEPT=NKEPT+1
        ENDIF
      ENDDO
      WRITE(6,104)NKEPT,NREM
  104 FORMAT('Points remaining: ',I5,'; Points removed:',I5,
     :       ' for proximity')
C
C     output reduced data set 
C
      OPEN(10,FILE=FNAME(1:JFND)//'.selv')
      NCORS=0
      DO J=1,NACC
        IF(NPIN(J).GT.0)THEN
C
C     finally, check magnitude of E-N corr and reset if necessary
C
          FF=GPSDAT(7,J)*GPSDAT(7,J)
          AA=GPSDAT(5,J)*GPSDAT(5,J)
          BB=GPSDAT(6,J)*GPSDAT(6,J)
          IF(FF.GT.0.5*(AA+BB))THEN
            WRITE(6,10106)STN(J),(GPSDAT(K,J),K=5,7)
            GPSDAT(7,J)=0.0
            NCORS=NCORS+1
          ENDIF
          LINE=STRINGN
          WRITE(LINE(1:66),10103)(GPSDAT(K,J),K=1,NPAR)
          IF(NSITE.NE.0)THEN
            LINE(69:76)=STN(J)
            LENTIL=78+LENCOM(J)-1
            LINE(78:LENTIL)=COMMENT(J)(1:LENCOM(J))
          ELSE
            LENTIL=69+LENCOM(J)-1
            LINE(69:LENTIL)=COMMENT(J)(1:LENCOM(J))
          ENDIF
          IF(LENTIL.GT.164)THEN
             WRITE(*,*)'Insufficient space allowed in LINE(164)'
             STOP
          ENDIF
10103     FORMAT(2F11.5,2F10.3,3F8.3)
          WRITE(FORMT(1:8),124)LENTIL
  124     FORMAT('(A',I3.3,')')
          WRITE(10,FORMT(1:8))LINE(1:LENTIL)
10106     FORMAT('Station ',A8,': E-N corr too large: ',3F13.5,
     :           ', set to zero')
        ENDIF
      ENDDO
      WRITE(6,10107)NCORS
10107 FORMAT('Warning: correlation variable has been reset on ',I4,
     :       ' sites')
C
      STOP
      END
C
C    routine to examine proximity of a point (XLOC,YLOC) to a boundary
C    (XLON(NPTS), YLAT(NPTS)), used to reject points that are too close.
C
      SUBROUTINE BDIST(XLOC,YLOC,XLON,YLAT,NPTS,ICLOSE,DMINX)
      DOUBLE PRECISION XLOC,YLOC,XLON,YLAT,AZIM,ARE2,XL,YL
      DOUBLE PRECISION DEL1,DEL2,DEL3,CS1,CS2
      DIMENSION XLON(NPTS),YLAT(NPTS)
      DMINX=180.0
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
C     perpendicular distance is 2*AREA/baseline length
C
          ELSE
            ARE2=DABS(((XLON(J)-XLOC)*(YL-YLOC)-
     :                 (XL-XLOC)*(YLAT(J)-YLOC)))
            DMIN=ARE2/DEL3
          ENDIF
C
C     minimum over all segments
C
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
      SUBROUTINE CLOSER(CLON,CLAT,J,GPSDAT,NPIN,NACC,DISMIN,NEARES)
C
C      searches a data list to find the point that is geographically
C      neares to (CLON,CLAT)
C
      DOUBLE PRECISION GPSDAT,CLON,CLAT,DLON,DLAT,DIST,AZIM
      DIMENSION GPSDAT(7,NACC),NPIN(NACC)
C
      DISMIN=181.0
      DO K=1,NACC
        IF((K.NE.J).AND.(NPIN(K).NE.0))THEN
          DLON=GPSDAT(1,K)
          IF(DABS(DLON-CLON).LT.1.D0)THEN
            DLAT=GPSDAT(2,K)
            IF(DABS(DLAT-CLAT).LT.1.D0)THEN
              CALL DISTAZ(CLON,CLAT,DLON,DLAT,DIST,AZIM)
              IF(DIST.LT.1.D0)THEN
                IF(DIST.LT.DISMIN)THEN
                  DISMIN=DIST
                  NEARES=K
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END
