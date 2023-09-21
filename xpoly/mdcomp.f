C
C    program mdcomp (model-data compare) to find nodal velocity values
C    (produced by sybil or sybilps arrow plot with verbose key on)
C    that correspond to a set of GPS measurements, to determine an
C    appropriate scaling parameter and report misfit. Raw misfit
C    (in mm/yr) and scaled misfit (relative to data covariance)
C    are both reported: for individual sites, and as rms values.
C    The weighted misfit is susceptible to failure if the e-n covariance
C    values are two large relative to the e and n variances.  I
C    recommend using the 'raw' weighting.
C
C    written by: G. Houseman, August 2015., revised Aug 2020
C
      DOUBLE PRECISION SUM2,SUMD2
      DOUBLE PRECISION QLO,QLA,XLO,XLA,DELT,AZ,EAR
C
C    double precision is used in routine DISTAZ
C
      PARAMETER (NBND=10000,NMOD=10000)
C     CHARACTER (len=180) :: STRING,STRINGN
C     CHARACTER (len=5) :: SITE
C     CHARACTER (len=80) :: FNAMEM,FNAMED,FNAMEC,FNAMET
C     CHARACTER (len=8) :: SCALEV,SCALEV0
      CHARACTER STRING*180,STRINGN*180
      CHARACTER SITE*5
      CHARACTER FNAMEM*80,FNAMED*80,FNAMEC*80,FNAMET*80
      CHARACTER SCALEV*8,SCALEV0*8
      DIMENSION STRING(NBND),SITE(NBND)
      DIMENSION XLON(NBND),XLAT(NBND),XU(NBND),XV(NBND)
      DIMENSION FIT1(NBND),FIT2(NBND)
      DIMENSION SEE(NBND),SNN(NBND),SEN(NBND)
      DIMENSION CD11(NBND),CD22(NBND),CD12(NBND)
      DIMENSION KEY(NBND)
      DIMENSION QLON(NMOD),QLAT(NMOD),QU(NMOD),QV(NMOD)
C
      DO J=1,10
        J1=(J-1)*8+1
        J2=J1+7
        FNAMET(J1:J2)='        '
      ENDDO
      SCALEV0='        '
      SCALEV=SCALEV0
      THRESH=0.0005    ! threshold for location match (degrees)
      DO K=1,180
        STRINGN(K:K)=' '
      ENDDO
C
C     first argument is the name of the GPS data file
C
      FNAMED=FNAMET
      CALL GET_COMMAND_ARGUMENT(1,FNAMED)
      IF(FNAMED.EQ.FNAMET)GO TO 10111
C
C     second argument is the name of the sybil velocity file
C
      FNAMEM=FNAMET
      CALL GET_COMMAND_ARGUMENT(2,FNAMEM)
      IF(FNAMEM.EQ.FNAMET)FNAMEM(1:8)='veln.out'
C
C     third argument is scale factor; if blank scale automatically
C   
      SCALEF=0.0
      CALL GET_COMMAND_ARGUMENT(3,SCALEV)
      IF(SCALEV.NE.SCALEV0)READ(SCALEV,*)SCALEF
C
C     new output misfit file named the same root with "_raw" or "_wgt"
C
      FNAMEC=FNAMEM
      MEND=0
      MDND=0
      MMND=0
      DO M=1,76
        IF(FNAMEC(M:M).NE.' ')MEND=M
        IF(FNAMED(M:M).NE.' ')MDND=M
        IF(FNAMEM(M:M).NE.' ')MMND=M
      ENDDO
C
C    first read the set of GPS points; expect: long, lat,
C    ve, vn, se, sn, corr, site name,  ++
C
      OPEN(8,FILE=FNAMED(1:MDND),STATUS='OLD',ERR=10111)
      NADD=0
   10 CONTINUE
      NADD=NADD+1
      IF(NADD.GT.NBND)THEN
        WRITE(6,12)NBND
   12   FORMAT('NBND =',I5,' provides insufficient memory for DATA ',
     :         'set - increase NBND and recompile mdcomp')
C       STOP
        CALL EXIT(2)
      ELSE
        STRING(NADD)=STRINGN
        READ(8,'(A180)',END=20)STRING(NADD)
C       WRITE(*,*)STRING(NADD)
      ENDIF
      READ(STRING(NADD),*)XLON(NADD),XLAT(NADD),XU(NADD),XV(NADD),
     :                    SEE(NADD),SNN(NADD),SEN(NADD),SITE(NADD)
      GO TO 10
   20 CONTINUE
      NADD=NADD-1
      CLOSE(8)
C     WRITE(6,2010)NADD
C2010 FORMAT('GPS data for ',I5,' points have been read')
C
C    now open and read model file 
C
      OPEN(4,FILE=FNAMEM(1:MMND),STATUS='OLD',ERR=10111)
      NMO=0
   30 CONTINUE
      NMO=NMO+1
      IF(NMO.GT.NMOD)THEN
        WRITE(6,32)NMOD
   32   FORMAT('NMOD =',i5,' provides insufficient memory for MODEL ',
     :         'set - increase NMOD and recompile mdcomp')
C       STOP
        CALL EXIT(2)
      ENDIF
      READ(4,*,END=40)K,QLON(NMO),QLAT(NMO),QU(NMO),QV(NMO)
      GO TO 30
   40 CONTINUE
      NMO=NMO-1
      CLOSE(4)
C     WRITE(6,2020)NMO
C2020 FORMAT(I5,' model values read')
C
C    for each data value search through the list to find 
C    a model value at (near enough) the same location
C
      MFIND=0
      DO J=1,NADD
        KEY(J)=0
      ENDDO
      DO J=1,NADD
        XLO=XLON(J)
        XLA=XLAT(J)
        DELTMIN=180.999
        DO K=1,NMO
          QLO=QLON(K)
          QLA=QLAT(K)
          CALL DISTAZ(XLO,XLA,QLO,QLA,DELT,AZ,EAR)
          IF(DELT.LT.DELTMIN)DELTMIN=DELT
          IF(DELT.LE.THRESH)THEN
            KEY(J)=K
            MFIND=MFIND+1
            GO TO 50
          ENDIF
        ENDDO
   50   CONTINUE
      ENDDO
      WRITE(6,2040)NADD,MFIND
 2040 FORMAT(' mdcomp: ',I5,' data to match; ',I5,' model points',
     ;       '  found.')
C
C    compute best fit scaling factor and rms misfit
C
      SUM1=0.0
      SUM2=0.d0
      SUMD1=0.0
      SUMD2=0.d0
      NPROB=0
      DO J=1,NADD
        K=KEY(J)
        SUM1=SUM1+XU(J)*QU(K)+XV(J)*QV(K)
        SUMD1=SUMD1+QU(K)*QU(K)+QV(K)*QV(K)
C
C    this calculation can be easily corrupted by bad values of the
C    E-N covariance.  Such values are pushed to zero here for the
C    purposes of calculating weighted scaling factor and misfit.
C
        A11=SEE(J)*SEE(J)
        A12=SIGN(1.0,SEN(J))*SEN(J)*SEN(J)
        A22=SNN(J)*SNN(J)
        IF((A12*A12).GT.(0.5*A11*A22))THEN
          SEN(J)=0.0
          A12=0.0
          NPROB=NPROB+1
        ENDIF
        DET=1.0/(A11*A22-A12*A12)
        CD11(J)= A22*DET
        CD12(J)=-A12*DET
        CD22(J)= A11*DET
        SUM2=SUM2+CD11(J)*XU(J)*QU(K)+CD22(J)*XV(J)*QV(K)+
     :            CD12(J)*(XU(J)*QV(K)+XV(J)*QU(K))
        SUMD2=SUMD2+CD11(J)*QU(K)*QU(K)+CD22(J)*QV(K)*QV(K)+
     :            CD12(J)*(QU(K)*QV(K)+QV(K)*QU(K))
      ENDDO
C
C     both scale factors are reported.  The choice of which to use
C     in basinv is set by the "misfit" variable in BASRUN
C
      IF(SCALEF.EQ.0.0)THEN
        VSCALE1=SUM1/SUMD1
        VSCALE2=SUM2/SUMD2
      ELSE
        VSCALE1=SCALEF
        VSCALE2=SCALEF
      ENDIF
      FITSM1=0.0
      FITSM2=0.0
      DRMS=0.0
      F1MIN=9.999e20
      F1MAX=0.0
      F2MIN=F1MIN
      F2MAX=0.0
      NSUM=0
      DO J=1,NADD
        K=KEY(J)
        IF(K.NE.0)THEN
          TERMX1=XU(J)-VSCALE1*QU(K)
          TERMY1=XV(J)-VSCALE1*QV(K)
          TERMX2=XU(J)-VSCALE2*QU(K)
          TERMY2=XV(J)-VSCALE2*QV(K)
          FIT1SQ=TERMX1*TERMX1+TERMY1*TERMY1
          FIT2SQ=CD11(J)*TERMX2*TERMX2 + 2.0*CD12(J)*TERMX2*TERMY2 
     :         + CD22(J)*TERMY2*TERMY2
          FIT1(J)=SQRT(FIT1SQ)
          IF(FIT2SQ.GE.0.0)THEN
            FIT2(J)=SQRT(FIT2SQ)
          ELSE
            WRITE(*,*)'problem with J =',J,' K = ',K,' FIT2SQ =',FIT2SQ
            FIT2(J)=0.0
          ENDIF
          IF(FIT1(J).LT.F1MIN)F1MIN=FIT1(J)
          IF(FIT1(J).GT.F1MAX)F1MAX=FIT1(J)
          IF(FIT2(J).LT.F2MIN)F2MIN=FIT2(J)
          IF(FIT2(J).GT.F2MAX)F2MAX=FIT2(J)
          FITSM1=FITSM1+FIT1SQ
          FITSM2=FITSM2+FIT2SQ
          DRMS=DRMS+XU(J)*XU(J)+XV(J)*XV(J)
          NSUM=NSUM+1
        ENDIF
      ENDDO
      RMSFT1=SQRT(FITSM1/FLOAT(NSUM))
      RMSFT2=SQRT(FITSM2/FLOAT(NSUM))
      DRMS=SQRT(DRMS/FLOAT(NSUM))
      WRITE(6,2050)VSCALE1,RMSFT1,DRMS,"raw"
      WRITE(6,2050)VSCALE2,RMSFT2,DRMS,"weighted"
 2050 FORMAT(' For scale factor ',G12.5,'  Misfit = ',G14.5,
     :       '  RMS data = ',G14.5,2X,A8)
      WRITE(6,2051)F1MIN,F1MAX,F2MIN,F2MAX
 2051 FORMAT(' Min and Max misfit (raw): ',2G14.5,' (weighted): '
     :       ,2G14.5)
      WRITE(*,*)'Warning: ',NPROB,' data have implausible ',
     :'covariance set to zero for weighted scale factor'
      DO J=1,NADD
        IF(KEY(J).EQ.0)THEN
        WRITE(6,2030)J,XLON(J),XLAT(J),SITE(J)(1:5)
        ENDIF
      ENDDO
 2030 FORMAT('No match for point: ',I5,2F10.3,X,A4)
C
C    write the amended data file, expected to comprise
C    lon, lat, ve (data), vn (data), ve (mod), vn (mod), misfit, site
C    model values of velocity will have had scale factor applied
C    two files are provided:  "raw" and "weighted"
C
C     WRITE(*,*)'opening file ', FNAMEC(1:MEND+4)
      OPEN(10,FILE=FNAMEC(1:MEND)//"_raw")
      OPEN(11,FILE=FNAMEC(1:MEND)//"_wgt")
      DO J=1,NADD
        K=KEY(J)
        IF(K.GT.0)THEN
          QR1=VSCALE1*QU(K)
          QR2=VSCALE1*QV(K)
          QW1=VSCALE2*QU(K)
          QW2=VSCALE2*QV(K)
        ELSE
          QR1=XU(J)
          QR2=XV(J)
          QW1=XU(J)
          QW2=XV(J)
        ENDIF
        WRITE(10,10113)XLON(J),XLAT(J),XU(J),XV(J),QR1,QR2,
     :                 FIT1(J),FIT2(J),SITE(J)
        WRITE(11,10113)XLON(J),XLAT(J),XU(J),XV(J),QW1,QW2,
     :                 FIT1(J),FIT2(J),SITE(J)
10113   FORMAT(2F11.5,4F10.4,2G14.5,X,A5)
      ENDDO
      CLOSE(10)
      CLOSE(11)
      WRITE(*,*)'Data / Model comparison is written to ',
     :FNAMEC(1:MEND)//"_raw",' and ',FNAMEC(1:MEND)//"_wgt"
      STOP
C
10111 WRITE(6,10112)
      WRITE(6,*)'looking for: ',FNAMED(1:MDND),' and ',FNAMEM(1:MMND)
10112 FORMAT('Apparent problem with command line arguments in mdcomp',
     : /,'Expected syntax: mdcomp gps-data-file sybil-vel-file')
      CALL EXIT(1)
      END
C
C    routine to make distance and azimuth calculations on a sphere
C
      SUBROUTINE DISTAZ(RLONG,RLATI,PLONG,PLATI,DELTA,AZIM,EARTH)
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
      IF(CD.GT. 1.D0)CD=1.D0
      IF(CD.LT.-1.D0)CD=-1.D0
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
