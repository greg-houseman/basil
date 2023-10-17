C*--------------------------------------------------------------------
C*    Basil / Sybil:   input.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
C       Lynn Evans   22/02/96
C
C       Routines for reading Basil input file
C
C       Routines which assign the value of the specified variables
C       (only a given set are recognised for each keyword) -
C         GETOUTPUTDATA GETREADDATA GETGEOMETRYDATA GETSERIESDATA
C         GETFORCEDATA GETDEFORMDATA GETRHEOLOGYDATA GETSOLVEDATA
C         GETFAULTDATA GETSTEPSIZEDATA GETSAVEDATA GETSTOPDATA
C         GETBCMODDATA GETLAYERDATA GETREMESHDATA GETLAGRANGEDATA
C         GETMESHDATA GETRHEOMODDATA GETTOPODATA
C       General routines - WRITEINSTR READLINE SKIP GETWORD GETVALUE
C
      SUBROUTINE INITIALREAD(NAMEBC,NMLEN,INFILE,IVARS,RVARS,
     :                       BINDIR,LNBD,OUTDIR,LNOD,NAMEDB,LNDB,
     :                       NAMEW,LNWF,NAMER,LNR,NX,NY,IMSH,IFLT,
     :                       XLEN,YLEN,AREA,QUALITY,NFP,NFPF3,
     :                       NXL1,NYL1,NPMP,NRMP,NSM,ILAG,IGRAV,
     :                       IDEN,ITEMP,ICR,IVOLD,IMREG,IVIS,NCOMP,
     :                       ITHDI,POLYFILE,IREAD,MAXNBS,JNO,KPBC,
     :                       PDIST,IERROR,IVALID,IVERB,VERSION)
C
C       Set flags for array allocation routine or read a solution
C
      INCLUDE "input.parameters"
      INCLUDE "input.data"
      INCLUDE "indices.parameters"
      INCLUDE "limits.parameters"
      COMMON/AI/LUW,LSC,LBC,LLG
      COMMON/LGDIM/NXL,NYL,NNL,NML,ILGSAVE
      COMMON/FILES/NMDUM1,LGFILE,PFILE
      COMMON/V01S/NMDUM2,NDATE,NDUM3,COMMEN
      INTEGER IVARS(64)
      REAL RVARS(64)
      CHARACTER IXY*3,IYX*1,IUT*3
      CHARACTER HNAME*12,WORD*12
      CHARACTER NMDUM1*16,NMDUM2*16,NDUM3*16
      CHARACTER BINDIR*32,OUTDIR*32,NAMEBC*32,NAMEDB*32
      CHARACTER NAMER*32,NAMEW*32
      CHARACTER NDATE*16
      CHARACTER FILER*60,FILEOUT*60
      CHARACTER POLYFILE*80,COMMEN*80
      CHARACTER LGFILE*80,PFILE*80
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER VERSION*8
      LOGICAL XIST
      DIMENSION IDAY(3)
C
      LBC=10
      LSC=6
      LUW=6
      LDIN=2
      LTMP=8
      IVALID=0
      IERROR=0
      ISET=1
      ITHDI=0
      KPBC=0
C
C    initialise strings used for filenames (and lengths)
C
      BINDIR='FD.sols                         '
      LNBD=7
      OUTDIR='FD.out                          '
      LNOD=6
      NAMEDB='debug                           '
      LNDB=5
      NAMEW='                                '
      LNWF=0
      NAMER='                                '
      LNR=0
      NDATE='                '
      DO N=1,80
        COMMEN(N:N) = ' '
      ENDDO
      IVALID=1
      OPEN(LBC,FILE=NAMEBC(1:NMLEN),STATUS='old',IOSTAT=IOS)
      IF(IOS.NE.0)THEN
        WRITE(*,100)NAMEBC
  100   FORMAT('Execution stopping due to Problem with Filename: ',A32)
        IVALID=0
      ENDIF
C
C    first check if there is an OUTPUT command, to get directory names
C
      IEND=0
      NUMLINE=0
  220 CONTINUE
        CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
        IF(IEND.EQ.1)GO TO 109
        NUMLINE=NUMLINE+1
        ITERM=0
        INDX2=1
        CALL GETWORD(INSTR,LENGTH,INDX2,WORD,LNW,MAXWORDLEN)
        IF (INDX2.NE.0) THEN
          ITERM = MATCH(SETUPTERMS,NUMSETUP,WORD)
C
C      process the output command if present
C
          IF (ITERM.EQ.INDXOUTPUT) THEN
            CALL GETOUTPUTDATA(INSTR,LENGTH,INDX2,NAMEW,LNWF,
     :                   NAMEDB,LNDB,BINDIR,LNBD,OUTDIR,LNOD,
     :                   IVERB,IERROR)
            IF (IERROR.NE.0) GO TO 600
          ENDIF
        ENDIF              ! if a non-blank line
      GO TO 220          ! go back for next line
C
C       NAMEDAT is the ascii data file for the SERIES data
C       NAMER is read data file; NAMEW is write data file; 
C       NAMEDB is checkpoint dump file
C
C      if no BIN entry provided, then output file name root
C      is set equal to the input filename
C
  109 CONTINUE
      IF(NAMEW(1:1).EQ.' ')THEN
        NAMEW=NAMEBC
        LNWF=NMLEN
      ENDIF
C
C      check that command file won't be overwritten
C
      IF((LNBD.EQ.1).AND.(BINDIR(1:1).EQ.'.'))THEN
        IF(NAMEW.EQ.NAMEBC)THEN
          WRITE(LSC,10219)BINDIR(1:LNBD),NAMEW(1:LNWF)
10219     FORMAT('Output is set to overwrite input file. ',
     :    'Suggest you change BINDIR in OUTPUT command')
          IERROR=1
          GO TO 700
        ENDIF
      ENDIF          ! if BINDIR ok
C
C       open the output file now, to receive ascii ouptut for
C       rest of this run.  Other files are managed in basil.F
C
      LUW=9
      LNWFF=LNOD+LNWF+5
      FILEOUT(1:LNWFF)=OUTDIR(1:LNOD)//'/'//NAMEW(1:LNWF)//'.out'
      INQUIRE(FILE=FILEOUT(1:LNWFF),EXIST=XIST)
      IF (XIST) THEN
        OPEN(LUW,FILE=FILEOUT(1:LNWFF),STATUS='OLD',IOSTAT=IOS)
        CLOSE(LUW,STATUS='DELETE')
      ENDIF
      OPEN(LUW,FILE=FILEOUT(1:LNWFF),STATUS='NEW',IOSTAT=IOS)
      IF (IOS.NE.0) THEN
        WRITE(LSC,*) 'Cannot open the file: ',FILEOUT(1:LNWFF)
        IERROR = 1
        GO TO 700
      ENDIF
C
C   record the host machine, date, time and input file
C
      ISTAT=HOSTNM(HNAME)
      WRITE(LSC,10095)VERSION(1:7),HNAME
      WRITE(LUW,10095)VERSION(1:7),HNAME
      CALL DATIME(IDAY,ISEC,0)
      ITHOUR=ISEC/3600
      ITMIN=ISEC/60-60*ITHOUR
      ITSEC=ISEC-3600*ITHOUR-60*ITMIN
      WRITE(NDATE(1:8),10202)IDAY(3),IDAY(2),MOD(IDAY(1),100)
      WRITE(NDATE(9:16),10203)ITHOUR,ITMIN,ITSEC
10202 FORMAT(I2.2,'/',I2.2,'/',I2.2)
10203 FORMAT(I2,':',I2.2,':',I2.2)
10095 FORMAT('basil version: ',A7,'  is running on: ',A12)
      WRITE(LUW,10101)NDATE(1:8),NDATE(9:16)
      WRITE(LSC,10101)NDATE(1:8),NDATE(9:16)
10101 FORMAT('Run commenced on ',A8,' at ',A8,/)
      WRITE(LUW,10128)JNO,NAMEBC
      WRITE(LSC,10128)JNO,NAMEBC
10128 FORMAT('Job number: ',I3,'  Reading from the input file ',A32)
      WRITE(LUW,10115)
      WRITE(LSC,10115)
10115 FORMAT('*********************************************',
     1'*********************************',/)
C
C            echo the entire input file
C
      REWIND(LBC)
      IEND=0
   50 CONTINUE
        CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
        IF (IEND.EQ.1) GO TO 51
        CALL WRITEINSTR(INSTR,LENGTH,0)
      GO TO 50
   51 WRITE(LUW,10115)
      WRITE(LSC,10115)
C
C     Now read the rest of the file from the top to get array sizes
C
      REWIND(LBC)
      IEND=0
      NUMLINE=0
 110  CONTINUE
      CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
      IF (IEND.EQ.1) GO TO 600
      NUMLINE=NUMLINE+1
      ITERM=0
      NEXTT = 1
      INDX2=1
C
C    the first word on each line is checked for command status
C
      CALL GETWORD(INSTR,LENGTH,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(SETUPTERMS,NUMSETUP,WORD)
        IF (ITERM.GT.0) THEN
C         WRITE(*,*)'Here is LINE:',INSTR(1:LENGTH)
C         WRITE(*,*)'Here is ITERM: ',ITERM,'   and WORD = ',WORD
C
C    if a MESH command present
C
          IF (ITERM.EQ.INDXMESH) THEN
            CALL GETMESHDATA(INSTR,LENGTH,INDX2,NX,NY,IMSH,IFLT,
     :                           POLYFILE,AREA,QUALITY,IERROR)
            RVARS(IELAREA)=AREA       !   2 lines added, GH, 4/5/15
            RVARS(IELQUAL)=QUALITY
C
C    if a GEOMETRY command present
C
          ELSE IF (ITERM.EQ.INDXGEOMETRY) THEN
            CALL GETGEOMETRYDATA(INSTR,LENGTH,INDX2,XLEN,
     :                        YLEN,DUM1,DUM2,DUM3,NCOMP,IGRAV,
     :                           IERROR)
C
C    if a LAGRANGE command present
C
          ELSE IF (ITERM.EQ.INDXLAGRANGE) THEN
            CALL GETLAGRANGEDATA(INSTR,LENGTH,INDX2,ILAG,NXL1,NYL1,
     :                           NPMP,NRMP,LGFILE,ILGSAVE,IERROR)
            IF (IERROR.NE.0) GO TO 600
            IF (ILAG.EQ.1.OR.ILAG.EQ.2) THEN
              REWIND(LBC)
              CALL COUNTSTRMARKERS(LBC,NSM,LSC,IERROR)
              IF (IERROR.NE.0) GO TO 700
              REWIND(LBC)
              DO I=1,NUMLINE
                CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
              ENDDO
            END IF
C
C   if a POLY command present (deprecated, use FILE in MESH command)
C
          ELSE IF (ITERM.EQ.INDXPOLY) THEN
            CALL GETPOLYDATA(INSTR,LENGTH,INDX2,POLYFILE,
     :                           IERROR)
            DO 80 K=1,80
              PFILE(K:K) = POLYFILE(K:K)
   80       CONTINUE
C
C    if a READ command present
C
          ELSE IF (ITERM.EQ.INDXREAD) THEN
            KSTART=1
            CALL GETREADDATA(INSTR,LENGTH,INDX2,NAMER,LNR,
     :                       KSTART,ISTR,IERROR)
            IF (IERROR.NE.0) GO TO 600
C
C            read the old solution from the read file
C
            IF(ISTR.GE.1)THEN
              LNRF=LNBD+LNR+1
              FILER(1:LNRF)=BINDIR(1:LNBD)//'/'//NAMER(1:LNR)
              OPEN(LDIN,FILE=FILER(1:LNRF),STATUS='old',
     :             FORM='unformatted',IOSTAT=IOS)
              IF(IOS.NE.0) THEN
                WRITE(LSC,10125)FILER(1:LNRF)
                WRITE(LUW,10125)FILER(1:LNRF)
                IERROR=1
                GO TO 700
              END IF
              IR=0
              DO 4 I = 1,KSTART
                READ(LDIN,END=7,ERR=7)
                READ(LDIN,END=7,ERR=7)(IVARS(K),K=1,64)
                READ(LDIN,END=7,ERR=7)(RVARS(K),K=1,64)
                MAXSIZE=IVARS(INROWS)
                IF (IVARS(INE)*8.GT.MAXSIZE) MAXSIZE=IVARS(INE)*8
                IF (IVARS(INELLEP).GT.MAXSIZE) MAXSIZE=IVARS(INELLEP)
                CALL SKIPSTORE(LDIN,1,IVARS,RVARS,IR,
     :                         MAXSIZE,MAXNBS,IERROR)
                WRITE(LSC,10105)IVARS(IKSTEP),RVARS(ITIMEL),ISTR
                WRITE(LUW,10105)IVARS(IKSTEP),RVARS(ITIMEL),ISTR
                IF (IERROR.NE.0) GO TO 7
                IR=IR+1
   4          CONTINUE
              NX = IVARS(INX)
              NY = IVARS(INY)
              NX1 = NX+1
              NY1 = NY+1
              IMSH = IVARS(IIMSH)
              IFLT = IVARS(IIFLT)
              NSM = IVARS(INSM)
              NPM = IVARS(INPM)
              NRM = IVARS(INRM)
              ILAG = IVARS(IILAG)
              ICR = IVARS(IICR)
              IDEN = IVARS(IIDEN)
              ITEMP = IVARS(IITEMP)
              IVOLD = IVARS(IIVOLD)
              IMREG = IVARS(IIMREG)
              NCOMP = IVARS(INCOMP)
              IGRAV = IVARS(IIGRAV)
              XLEN = RVARS(IXLEN)
              YLEN = RVARS(IYLEN)
              AREA = RVARS(IELAREA)
              CENTLNG = RVARS(IXREFM)
              YMID = RVARS(IYREFM)
              QUALITY = RVARS(IELQUAL)
              SEGMIN = RVARS(IELSEG)
              IREAD=KSTART
   7          CLOSE(LDIN)
            END IF
C
C    if a VISDENS command present
C
          ELSE IF (ITERM.EQ.INDXVISDENS) THEN
C
C           Update viscosity and density arrays if the initial values
C           of SE or RHO are changed
C
            DIFFOLD = RVARS(ITDIFF)
            ALPHAOLD = RVARS(IALPHA)
C           There does not seem to be a default TEMP so using the
C           ITEMP flag to trigger initialising TEMPT array
C           TEMPOLD = TEMP
            TEMP = 1.0
            SEOLD = SE
            RHOOLD = RHO
            GAMMAT=0.0
            BETA=0.0
            TREF=0.0
            VC=1.0
            CALL GETVISDENSDATA(INSTR,LENGTH,INDX2,VC,SE,RHO,
     :                   ITEMP,TEMP,TDIFF,ALPHA,GAMMAT,BETA,TREF,
     :                   IERROR)
            IF(GAMMAT.NE.0.0)ITHDI=1
            IF (IERROR.NE.0) GO TO 600
            REWIND(LBC)
C           Look for REG statements, set IVIS, IDEN, IMREG depending on
C            parsing of VR, DR, MR
            CALL VISVARFLAGS(LBC,SE,IVIS,IDEN,IMREG,
     :                       LUW,LSC,IDBUG,IERROR)
            IF (IERROR.NE.0) GO TO 700
            REWIND(LBC)
            DO I=1,NUMLINE
              CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
            ENDDO
            IF (ABS(SE-1.0).GT.EPS) IVIS = 1
C
C    if a LAYER command present
C
          ELSE IF (ITERM.EQ.INDXLAYER) THEN
            CALL GETLAYERDATA(INSTR,LENGTH,INDX2,ICR,
     :                      ARGAN,BRGAN,THRESH,HLENSC,
     :                      BDEPSC,RISOST,REFLEV,IERROR)
            IF (IERROR.NE.0) GO TO 600
C
C    if a FORCE command present
C
          ELSE IF (ITERM.EQ.INDXFORCE) THEN
            IDEN = 1
C
C    if a PB command present (periodic boundary)
C
          ELSEIF(ITERM.EQ.INDXPERIB)THEN
            CALL PARSEON(INSTR,IXY,PDIST,IYX,YXLIM1,YXLIM2,IBND,IUT,
     :               UTVAL1,UTVAL2,ITAP,IPR,XMIN,XMAX,YMIN,YMAX,IERR)
            IF(IXY(1:1).EQ.'X')KPBC=1
            IF(IXY(1:1).EQ.'Y')KPBC=2
            IF(ABS(PDIST).LE.2.E-5)THEN
              IERR=1
              WRITE(*,10111)PDIST
10111         FORMAT('Repeat length in PB command ',G12.5,' unaccept',
     :               'able.  Check consistent with dimension of region')
            ENDIF
          ENDIF
        ENDIF          ! processing of SETUPTERMS now complete
C
C    also examine TIMESTEPTERMS
C
        ITERM = MATCH(TIMESTEPTERMS,NUMTIME,WORD)
        IF (ITERM.GT.0) THEN
C
C          if DENSMOD or RHEOMOD statement present
C
          IF (ITERM.EQ.INDXDENSMOD) THEN
            IDEN = 1
          ELSE IF (ITERM.EQ.INDXRHEOMOD) THEN
            IVOLD = 1
          END IF
        END IF                       ! a valid TIMESTEPTERM
      END IF                         ! a non-null command found (INDX2)
      IF (IERROR.EQ.0) GO TO 110
C
C    Echo the instruction if it could not be parsed
C
 600  IF (IERROR.GT.9) CALL WRITEINSTR(INSTR,LENGTH,INDX2)
C
C    Validate the mesh data
C
      IF (XLEN.EQ.0.0.AND.YLEN.EQ.0.0.AND.IMSH.LT.3) ISET=0
      IF (NX.EQ.0.AND.NY.EQ.0) NX=10
      IF (IMSH.LT.-1.OR.IMSH.GT.3) ISET=0
      IF (IMSH.EQ.3.AND.IFLT.EQ.2) ISET=0
      IF (IFLT.LT.0.OR.IFLT.GT.4) ISET=0
      IF((IFLT.EQ.0).AND.(KPBC.GT.0))IFLT=1
      IF (ISET.NE.1) THEN
        WRITE(*,200)NAMEBC(1:NMLEN)
  200   FORMAT('MESH input not valid in file: ',A32)
        IVALID=0
      END IF
 700  CLOSE(LBC)
10105 FORMAT('Array sizes for KSTOP =',I6,'  T =',f7.3,
     1' read from unit',i3)
10125 FORMAT(' The read file could not be opened ',a60,/)
      RETURN
      END

      SUBROUTINE GETMESHDATA(INSTR,LEN,INDX2,NX,NY,IMSH,IFLT,
     :                           POLYFILE,AREA,QUALITY,IERR) 
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER TMPFILE*80
      CHARACTER POLYFILE*80
      PARAMETER(           INDXNX=1,
     :                     INDXNY=2,
     :                     INDXTYPE=3,
     :                     INDXIFLT=4,
     :                     INDXFILE=5,
     :                     INDXAREA=6,
     :                     INDXQUAL=7,
     :                     NUMPARAMETERTERMS=7)
      CHARACTER PARAMETERTERMS(1:NUMPARAMETERTERMS)*12
      DATA PARAMETERTERMS /'NX          ',
     :                     'NY          ',
     :                     'TYPE        ',
     :                     'FAULT       ',
     :                     'FILE        ',
     :                     'AREA        ',
     :                     'QUALITY     '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(PARAMETERTERMS,NUMPARAMETERTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXNX) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NX
        ELSE IF (ITERM.EQ.INDXNY) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NY
        ELSE IF (ITERM.EQ.INDXTYPE) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IMSH
        ELSE IF (ITERM.EQ.INDXIFLT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IFLT
        ELSE IF (ITERM.EQ.INDXFILE) THEN
          CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
          IF (INDX1.EQ.0) GO TO 210
          INDX2 = INDX1+1
          CALL GETWORD(INSTR,LEN,INDX2,TMPFILE,LNW,80)
          POLYFILE(1:INDX2-INDX1-1)=TMPFILE(1:INDX2-INDX1-1)
        ELSE IF (ITERM.EQ.INDXAREA) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)AREA
        ELSE IF (ITERM.EQ.INDXQUAL) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)QUALITY
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END
C
C******************************************************************************
C
      SUBROUTINE PROCESSINPUT(NAMEBC,NMLEN,BINDIR,LNBD,OUTDIR,LNOD,
     :                  NAMEDB,LNDB,NAMEW,LNWF,NAMER,LNR,IV,RV,
     :                  IFLTTIPS,IOFF,ITIMETYPE,IDENMOD,
     :                  IRHEOMOD,MFIXFLAG,NITER,IGRID,NDFORM,
     :                  RHO,TILTDEG,CENTLNG,ITIMESTEPORDER,
     :                  EX,EY,UVP,LEM,NOR,IBC,IBNGH,IBCTYP,KORNER,
     :                  QBND,POLEP,IBPOLE,TAPERF,VHB,IMAT,
     :                  VOLD,SSQ,FROT,DENS,TEMPT,IELFIX,
     :                  XYEXIT,SHAPEMIN,SEGMIN,ANGMIN,
     :                  IFBC1,IFBC2,IFEQV,ISEG,
     :                  JFBC1,JFBC2,IENDP1,IENDP2,NFPF3,
     :                  LGEM,EXLG,EYLG,STELPX,STELPY,
     :                  UXLG,UYLG,LGIBC,LGIBCF,
     :                  NE,NUP,NBP,NMP,NSEG,
     :                  GSDUDX,GSDUDY,GSDVDX,GSDVDY,INITVEL,
     :                  REMESHFLAG,NEGTIM,IERROR)
C
      INCLUDE "indices.parameters"
      INCLUDE "limits.parameters"
      INCLUDE "input.parameters"
      INCLUDE "input.data"
      CHARACTER NDATE*16,NMDUM1*16,NMDUM2*16,NDUM3*16
      CHARACTER COMMEN*80,LGFILE*80,PFILE*80
      CHARACTER BINDIR*32,OUTDIR*32,NAMEDB*32,NAMEBC*32
      CHARACTER NAMER*32,NAMEW*32
      CHARACTER FILER*60,FILEOUT*60
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER ISSEQ, ISUNF
      LOGICAL INITVEL
      LOGICAL XIST,LELEMENTFLAG,REMESHFLAG
      DOUBLE PRECISION ROTIX,PLONR,PROTR
C
      COMMON/LGDIM/NXL,NYL,NNL,NML,ILGSAVE
      COMMON/FILES/NMDUM1,LGFILE,PFILE
      COMMON/AI/LUW,LSC,LBC,LLG
      COMMON/V01S/NMDUM2,NDATE,NDUM3,COMMEN
      COMMON/BT/MSINDX,MEASUR(MAXMEAS),MSNODE(MAXMEAS),PTLOC(3,MAXMEAS)
C
      DIMENSION IV(64), RV(64)
      DIMENSION IVTMP(64), RVTMP(64)
      DIMENSION ITIMESTEPORDER(NUMTIME)
      DIMENSION NOR(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION VHB(8,NE)
      DIMENSION VOLD(8,NE)
      DIMENSION SSQ(NUP),TEMPT(NUP)
      DIMENSION DENS(7,NE)
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION EXLG(IV(INUL)),EYLG(IV(INUL))
      DIMENSION STELPX(IV(INPM)*IV(INSM)),STELPY(IV(INPM)*IV(INSM))
      DIMENSION ISEG(3,NSEG)
      DIMENSION BBOX(2,4)
      DIMENSION KORNER(10)
      DIMENSION IBC(NBP),IBCTYP(NBP*2),IBNGH(NBP*2)
      DIMENSION QBND(NBP*2)
      DIMENSION XYEXIT(4)
      DIMENSION POLEP(3,MAXPOLE)
      DIMENSION IBPOLE(NBP,2),TAPERF(NBP)
      DIMENSION ROTIX(3,3)
      DIMENSION PTALO(1),PTALA(1)
C
C     PI=3.1415926536
      RTODEG=180.0/3.1415926536
      XZERO=0.
      YZERO=0.
      LTMP=8
      LDIN=2
      LDUM=3
      LELEMENTFLAG=.FALSE.
      DO N=1,NUMTIME
        ITIMESTEPORDER(N) = 0
      ENDDO
      IEND=0
      NUMLINE=0
C
C    main loop for reading input lines begins here
C
 110  CONTINUE
      IVALID=0
      CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
      IF (IEND.EQ.1) GO TO 700
      NUMLINE=NUMLINE+1
      ITERM=0
      NEXTT = 1
      INDX2=1
      CALL GETWORD(INSTR,LENGTH,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(SETUPTERMS,NUMSETUP,WORD)
        IF (ITERM.GT.0)THEN
C         WRITE(LSC,10102)SETUPTERMS(ITERM)
10102     FORMAT(/,'Processing command: ',A12)
C
C     what follows is a series of IF blocks, one for each valid command
C     first is for the LABEL command, inserts a comment in file
C     find present end of comment, before adding new comment
C     the OUTPUT command is ignored here as it was analysed in INITIALREAD
C
        IF (ITERM.EQ.INDXLABEL) THEN
          CALL SKIP(COMMEN,80,1,'b',J)
          CALL SKIP(COMMEN,80,J,'n',J2)
          IF (J2.LT.80) CALL SKIP(COMMEN,80,J2,'b',J)
          J= J+1
          CALL SKIP(INSTR,LENGTH,INDX2,'n',INDX1)
          IF (INDX1.GT.0.AND.INDX1.LT.LINELEN) THEN
            COMMEN(J:LINELEN)=INSTR(INDX1:LINELEN-J+1)
          ELSE
            IERROR = 10
          END IF
C
C     for the READ command, reads old solution
C
        ELSE IF (ITERM.EQ.INDXREAD) THEN
          KSTART=1
          CALL GETREADDATA(INSTR,LENGTH,INDX2,NAMER,LNR,
     :                     KSTART,ISTR,IERROR)
          IF (IERROR.NE.0)THEN
            IERROR=11
            GO TO 600
          ENDIF
          IF(ISTR.GE.1)THEN
            LNRF=LNBD+LNR+1
            FILER(1:LNRF)=BINDIR(1:LNBD)//'/'//NAMER(1:LNR)
            OPEN(LDIN,FILE=FILER(1:LNRF),STATUS='old',
     :           FORM='unformatted',IOSTAT=IOS)
            IF(IOS.NE.0) THEN
              WRITE(LSC,10125)FILER(1:LNRF)
              WRITE(LUW,10125)FILER(1:LNRF)
              IERROR=12
              GO TO 700
            END IF
C
C     Read over records at the beginning of the file until
C     start of required record.
C
            IR=0
            DO 4 I = 1,KSTART-1
              READ(LDIN,END=8,ERR=8)
              READ(LDIN,END=8,ERR=8)(IVTMP(K),K=1,64)
              READ(LDIN,END=8,ERR=8)(RVTMP(K),K=1,64)
              MAXSIZE=IVTMP(INROWS)
              IF (IVTMP(INE)*8.GT.MAXSIZE) MAXSIZE=IVTMP(INE)*8
              IF (IVTMP(INELLEP).GT.MAXSIZE) MAXSIZE=IVTMP(INELLEP)
              CALL SKIPSTORE(LDIN,1,IVTMP,RVTMP,IR,
     :                       MAXSIZE,MAXNBS,IERROR)
              IF(IERROR.NE.0)THEN
                WRITE(LSC,*)'SKIPSTORE has returned IERROR =',IERROR
                WRITE(LUW,*)'SKIPSTORE has returned IERROR =',IERROR
                CLOSE(LDIN)
                IERROR=13
                GO TO 700
              ENDIF
              WRITE(LSC,10104)IVTMP(IKSTEP),RVTMP(ITIMEL),ISTR
              WRITE(LUW,10104)IVTMP(IKSTEP),RVTMP(ITIMEL),ISTR
              IR=IR+1
              GO TO 4
    8         IERROR=14
              WRITE(LSC,*)'End of file encountered on READ'
              WRITE(LUW,*)'End of file encountered on READ'
              CLOSE(LDIN)
              GO TO 700
    4       CONTINUE
C
C           Reading the input file in INITIALREAD may have changed
C           some flags e.g. VISDENS that was not in the original
C           input file.  Therefore, we should pass the values in 
C           the starting solution to READSTORE
C
            DO K=1,64
              IVTMP(K) = IV(K)
              RVTMP(K) = RV(K)
            END DO
            CALL READSTORE(LDIN,1,IVTMP,RV,NAMER,NAMEW,IR,IERROR,
     :               IVTMP(INROWS),IVTMP(INUP),IVTMP(INE),IVTMP(INBP),
     :               IVTMP(INBP2),IVTMP(INFP),IVTMP(INEL),
     :               IVTMP(INUL),IVTMP(INPM),IVTMP(INSM),IVTMP(INBL),
     :               IVTMP(INELLEP),IVTMP(INSEG),
     :               EX,EY,UVP,LEM,NOR,VHB,IMAT,VOLD,
     :               SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :               IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC1,IFBC2,
     :               IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :               UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :               IELLE,IPOLYN)
            IF(IERROR.NE.0)THEN
              WRITE(LSC,*)'READSTORE has returned IERROR =',IERROR
              WRITE(LUW,*)'READSTORE has returned IERROR =',IERROR
              CLOSE(LDIN)
              IERROR=15
              GO TO 700
            ENDIF
            IR=IR+1
            REMESHFLAG=.TRUE.
            WRITE(LSC,10104)IVTMP(IKSTEP),RVTMP(ITIMEL),ISTR
            WRITE(LUW,10104)IVTMP(IKSTEP),RVTMP(ITIMEL),ISTR
            CLOSE(LDIN)
C
C     reset time counters if ISTR=1
C
            IF(ISTR.EQ.1)THEN
              IV(IKSTEP)=0
              RV(ITIMEL)=0.0
            END IF
C
C            if read error or empty file or incompatible solution,
C            start next job
C
            IF ((IERROR.EQ.99).OR.(IERROR.EQ.98.AND.IR.LT.1)
     1                        .OR.(IERROR.EQ.97)) THEN
              GO TO 700
            ELSE
C            else using last record in file
              IERROR=0
            END IF
            IF (ISTR.LE.1)IV(IISAVE)=-10000
          END IF
C
C      The GEOMETRY command
C
        ELSE IF (ITERM.EQ.INDXGEOMETRY) THEN
C             xlen and ylen have already been set in INITIALREAD
            CALL GETGEOMETRYDATA(INSTR,LENGTH,INDX2,DUM1,DUM2,
     :                  XZERO,YZERO,CENTLNG,IV(INCOMP),IV(IIGRAV),
     :                           IERROR)
            IF (IV(IIMSH).NE.3) THEN
C
C          if not using the triangulating routine,
C          set up the finite element mesh, label and order the nodes
C
              IF(IV(IIMSH).LE.0) THEN
                CALL MESH(IV(INX),IV(INY),RV(IXLEN),RV(IYLEN),
     :                    EX,EY,LEM,IV(INUP),IV(INE),IV(INN),
     :                    IV(IIFLT),IV(IIMSH),LUW,LSC,IERROR)
               IF (IERROR.NE.0)THEN
                  IERROR=21
                  GO TO 700
                ENDIF
                CALL MPNODE(IV(INX),IV(INY),EX,EY,LEM,IV(INUP),IV(INE),
     :                      IV(INN),IV(IIFLT),IV(IIMSH),LUW,LSC)
              ELSE
                CALL MESH2(IV(INX),IV(INY),RV(IXLEN),RV(IYLEN),
     :                     EX,EY,LEM,IV(INUP),IV(INE),IV(INN),
     :                     IV(IIMSH),IV(IIFLT),LUW,LSC,IERROR)
                IF (IERROR.NE.0)THEN
                  IERROR=22
                  GO TO 700
                ENDIF
              END IF
              CALL NORDER(EX,EY,LEM,NOR,IV(INUP),IV(INE),IV(INN),
     :                    LSC,IERROR)
              DO K=1,IV(INUP)
C               NOR(K)=K
                EX(K)=EX(K)+XZERO
                EY(K)=EY(K)+YZERO
              END DO
            END IF !mesh not 3
            IF (IERROR.NE.0)THEN
              IERROR=23
              GO TO 700
            ENDIF
            IF(IV(IIFLT).EQ.2) CALL GRINIT(IV(INX),IV(INY),
     :                                     RV(IXLEN),RV(IYLEN),
     :                                     EX,EY,NOR,IENDP1,IENDP2,
     :                                     IV(INUP),IV(INN),
     :                                     LUW,LSC,IERROR)
           IF (IERROR.NE.0)THEN
              IERROR=24
              GO TO 700
            ENDIF
C
C    Establish the tables which identify boundary nodes and connections
C    VSBTAB2 should be preferable to VSBTAB but neither should be needed
C    if trimesh has been called
C
            IF (IV(IIMSH).NE.3) THEN
                    WRITE(*,*)'input.f IV(IIMSH) =',IV(IIMSH)
C           IF (IV(IIFLT).EQ.0) THEN
              CALL VSBTAB2(EX,EY,LEM,NOR,IBC,IBNGH,IBCTYP,
     :         KORNER,IV(INUP),IV(INE),IV(INBP),IV(INMP),IV(INN))
C             CALL VSBTAB(RV(IXLEN),RV(IYLEN),EX,EY,LEM,NOR,IBC,IBNGH,
C    :             XZERO,YZERO,IBCTYP,KORNER,IV(INUP),IV(INE),IV(INBP),
C    :                  IV(IIFLT),LUW,LSC,IERROR)
C             IF (IERROR.NE.0)THEN
C                IERROR=25
C               GO TO 700
C             ENDIF
C             WRITE(*,*)'IBC array after VSBTAB2'
C             CALL MITPRT(IBC,IV(INBP),1,IV(INBP),LSC)
            ENDIF
C
C          Translate mesh position, if necessary
C
C           CALL FINDBBOX(EX,EY,BBOX,NUP)
C           YMIN = BBOX(2,1)
C           YMAX = BBOX(2,4)
C           XMIN = BBOX(1,1)
C           XMAX = BBOX(1,2)
C           EPS=1e-5
C           DIFF=XZERO-XMIN
C           IF (ABS(DIFF).GT.EPS) THEN
C             DO 101 I=1,NUP
C101            EX(I) = EX(I)+DIFF
C           END IF
C           DIFF=YZERO-YMIN
C           IF (ABS(DIFF).GT.EPS) THEN
C             DO 102 I=1,NUP
C102            EY(I) = EY(I)+DIFF
C           END IF

          ELSE IF (ITERM.EQ.INDXBCOND) THEN
            CALL GETBCONDDATA(INSTR,LENGTH,INDX2,IOFF,
     :                        RV(IVELXO),RV(IVELYO),RV(IYLDSTR),
     :                        IERROR)
            IF (IERROR.NE.0)THEN
              IERROR=31
              GO TO 600
            ENDIF
C           WRITE(*,*)'Before VSBCAPPLY, IBCTYP'
C           CALL IMATPP(IBCTYP,IV(INBP)*2)
C           WRITE(*,*)'Before VSBCAPPLY, QBND'
C           CALL SMATPP(QBND,IV(INBP)*2)
C           WRITE(*,*)'Before VSBCAPPLY, UVP'
C           CALL SMATPP(UVP,IV(INUP)*2)
C
C          Read and set boundary conditions from the basil input file
C
C           IDBUG=0
            REWIND(LBC)
C           WRITE(*,*)'ISEG array properties'
C           CALL MITPRT(ISEG,IV(INSEG),1,IV(INSEG),6)
C     every 3rd entry in ISEG is the boundary label, but why is IV(INSEG) not a multiple of 3?
C
            CALL VSBCAPPLY(RV(IXLEN),RV(IYLEN),RV(IBIG),IFLTTIPS,
     :                  EX,EY,QBND,POLEP,IPOLE,RV(ISE),
     :                  RV(IARGANP),RV(IHLENSC),CENTLNG,RV(IDEFV1),
     :                  IV(INCOMP),NOR,LEM,IBC,IBNGH,IBCTYP,
     :                  IV(IIDEFTYP),ISEG,IV(INSEG),IV(INUP),
     :                  IV(INE),IV(INN),IV(INBP),IV(IIFLT),
     :                  IV(INFP),NFPF3,IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :                  LBC,LUW,LSC,IVERB,IERROR)
            IV(IIPOLE)=IPOLE
C           DO J=1,NBP
C             IF(IBPOLE(J,1).NE.0)THEN
C               IV(IIVRESET)=1
C               GO TO 31
C             ENDIF
C           ENDDO
   31       CONTINUE
            REWIND(LBC)
            IF (IERROR.NE.0)THEN
              IERROR=32
              GO TO 700
            ENDIF
            DO I=1,NUMLINE
              CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
            ENDDO
C
C      Set up fault arrays to deal with periodic boundaries if needed
C
C           WRITE(*,*)'about to call FLTBC'
C           IF((IV(IIFLT).EQ.1).OR.(IV(IIFLT).EQ.4))THEN
C             PDIST=1.0     ! these two things should come from the bc statements
C             KBXY=0
C             CALL FLTBC(EX,EY,QBND,LEM,NOR,IBC,IBNGH,IBCTYP,
C    :                 IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
C    :                 IV(IIFLT),IV(INUP),IV(INE),IV(INBP),
C    :                 IV(INFP),NFPF3,PDIST,KBXY,LUW,LSC,IERROR)
C             WRITE(*,*)'FLTBC called'
C             IF (IERROR.NE.0)THEN
C               IERROR=33
C               GO TO 700
C             ENDIF
C           ENDIF
C          in case IOFF is set, need to set TBXOFF
            CALL FINDBBOX(EX,EY,BBOX,NUP)
            YMIN = BBOX(2,1)
            YMAX = BBOX(2,4)
            XMIN = BBOX(1,1)
            XMAX = BBOX(1,2)
C           XMAX=RV(IXLEN)
C           YMAX=RV(IYLEN)
C           IF(IV(IIFLT).EQ.2) XMAX=2.0*RV(IXLEN)
            IF(IV(IIFLT).EQ.2) XMAX=2.0*XMAX
            IF(IOFF.NE.0) RV(ITBXOFF)=XMAX-XMIN
C           INITVEL=.TRUE.
C           WRITE(*,*)'QBND after VSBCAPPLY'
C           CALL SMATPP(QBND,2*IV(INBP))
C           WRITE(*,*)'IBCTYP after VSBCAPPLY'
C           CALL IMATPP(IBCTYP,2*IV(INBP))
C           WRITE(*,*)'after VSBCAPPLY, UVP'
C           CALL SMATPP(UVP,IV(INUP)*2)
C
          ELSE IF (ITERM.EQ.INDXFORCE) THEN
            CALL GETFORCEDATA(INSTR,LENGTH,INDX2,RV(IRHOG1),
     :                        IERROR)
C
C     LAYER command 
C
          ELSE IF (ITERM.EQ.INDXLAYER) THEN
            SSOLD = RV(IHLENSC)
            CALL GETLAYERDATA(INSTR,LENGTH,INDX2,IV(IICR),
     :                        RV(IARGAN),RV(IBRGAN),RV(ITHRESH),
     :                        RV(IHLENSC),RV(IBDEPSC),RV(IRISOST),
     :                        RV(IREFLEV),IERROR)
            IF (IERROR.NE.0)THEN
              IERROR=41
              GO TO 600
            ENDIF
            IF (RV(IHLENSC).NE.SSOLD) THEN
              DO I=1,IV(INUP)
                SSQ(I)=-ALOG(RV(IHLENSC))
              ENDDO
            END IF
C
C     assumes that VISDENS command already called, but may be needed for BCOND
C
            RV(IARGANP)=RV(IARGAN)*RV(IBDEPSC)*RV(IBDEPSC)*
     1                      (RV(IHLENSC)**(2.0+1.0/RV(ISE)))
C
C     TOPO command
C
          ELSE IF (ITERM.EQ.INDXTOPO) THEN
            HMIN=XZERO
            VMIN=YZERO
            HLEN=RV(IXLEN)
            VLEN=RV(IYLEN)
            CALL GETTOPODATA(INSTR,LENGTH,INDX2,FILER,ITYPE,
     :                       HMINDUM,HLENDUM,VMINDUM,VLENDUM,
     :                       IFILTER,IERROR)
            IV(ITOPO)=1
C
C          assume geometry values are longitude and
C          latitude so check range
C
            IF ((VMIN.LT.-90.0).OR.((VMIN+VLEN).GT.90.0).OR.
     :         (HMIN.LT.-180.0).OR.(HMIN.GE.RLONGMAX)) THEN
              HMAX=AMOD(HMIN+HLEN,RLONGMAX)
              WRITE(LSC,10124)HMIN,HMAX,VMIN,VMIN+VLEN
              WRITE(LUW,10124)HMIN,HMAX,VMIN,VMIN+VLEN
              IERROR=42
              GO TO 700
            END IF
C
C          Read the geographical data
C
            CALL TOPODATA(FILER,ITYPE,HMIN,VMIN,HLEN,VLEN,
     :                    EX,EY,SSQ,NOR,IV(INUP),
     :                    RV(IHLENSC),RV(IBDEPSC),
     :                    RV(IRISOST),RV(IREFLEV),IERROR)
            IF (IERROR.NE.0)THEN
              IERROR=43
              GO TO 600
            ENDIF

          ELSE IF (ITERM.EQ.INDXDEFORM) THEN
            CALL GETDEFORMDATA(INSTR,LENGTH,INDX2,IV(IIDEFTYP),
     :                         RV(IBANGL),RV(IERA),RV(IDEFV1),
     :                         TILTDEG,NDFORM,IERROR)
            IF (IERROR.NE.0)THEN
              IERROR=44
              GO TO 600
            ENDIF
            IF(IV(IIDEFTYP).NE.0.AND.IV(IIDEFTYP).LT.100)THEN
C
C            Set the initial layer thickness array SSQ
C
              IF(IV(IIDEFTYP).EQ.6) THEN
                CALL SSINIT(IV(IIDEFTYP),RV(IHLENSC),
     :                          RV(IDEFV1),RV(IDEFV2),RV(IDEFV3),
     :                          RV(IDEFV4),RV(IDEFV5),
     :                          RV(IXLEN),RV(IYLEN),
     :                          EX,EY,SSQ,NOR,IV(INUP),
     :                          LUW,LSC,IERROR)
                IF (IERROR.NE.0)THEN
                  IERROR=45
                  GO TO 600
                ENDIF
C
              ELSEIF(IV(IIDEFTYP).EQ.16) THEN
                 CALL SSINIT2(IV(IIDEFTYP),RV(IHLENSC),
     :                          RV(IDEFV1),RV(IDEFV2),RV(IDEFV3),
     :                          RV(IDEFV4),RV(IDEFV5),
     :                          RV(IXLEN),RV(IYLEN),
     :                          EX,EY,SSQ,NOR,IV(INUP),
     :                          LUW,LSC,IERROR)
               IF (IERROR.NE.0)THEN
                  IERROR=46
                  GO TO 600
                ENDIF
C
C             distort the initial mesh if necessary, DEFV contains parameters
C             defining the mesh distortion
C
              ELSE
                CALL CFORM(IV(IIDEFTYP),RV(IBANGL),
     :                   RV(IERA),RV(IXLEN),RV(IYLEN),
     :                   RV(IDEFV1),IV(INCOMP),XZERO,YZERO,
     :                   EX,EY,LEM,NOR,IV(INUP),IV(INE),IV(INN),
     :                   LUW,LSC)
              END IF !deftype 6
            END IF !deftype 0->100
            IF(IV(IIFLT).EQ.2) XMAX=2.0*RV(IXLEN)
            IF(IV(IIDEFTYP).EQ.101) CALL DFORM(NDFORM,IV(IIFCASE),
     :                                 XMIN,XMAX,YMIN,YMAX,IFLTTIPS,
     :                                 EX,EY,NOR,LEM,
     :                                 IV(INUP),IV(INE),IV(INN))
            IF(IV(IIDEFTYP).EQ.100) CALL TILT(TILTDEG,EX,EY,LEM,NOR,
     :                                   IV(INUP),IV(INE),IV(INN))
C
C      if coordinate system is rotated using the COROTATE option
C      first use two points to determine pole location and rotation matrix
C      then apply rotation to mesh and to Lagrange mesh coordinates
C      by storing two points in RV(IDEFV1 - IDEFV4), they are saved in file.
C
            IF(IV(IIDEFTYP).EQ.111)THEN
              CALL GETPOLE(RV(IDEFV1),PLONR,PROTR,ROTIX)   ! PLONR, PROTR in radians
              PTALO(1)=RV(IDEFV1)                          ! values changed by ROTATE
              PTALA(1)=RV(IDEFV2)                          ! dimension (1) to avoid compiler warning
              CALL ROTATE(ROTIX,PTALO,PTALA,1,0.0)
              IF(ABS(PTALA(1)).GT.1.e-7)THEN
                WRITE(LSC,122)PTALO,PTALA
                WRITE(LUW,122)PTALO,PTALA
  122           FORMAT('Warning: point A transformed to ',2F11.5,
     :                 '; lat should be zero !')
              ELSE
                WRITE(LSC,120)(RV(JJ),JJ=IDEFV1,IDEFV1+3)
                WRITE(LUW,120)(RV(JJ),JJ=IDEFV1,IDEFV1+3)
  120           FORMAT('New equator defined by: 2PTS',/,4F9.3,/)
              ENDIF
C
C      GETPOLE returns parameters in radians, ROTATE in degrees
C
              PLON=PLONR*RTODEG
              PROT=PROTR*RTODEG
              DPHI=PTALO(1)
              CALL ROTATE(ROTIX,EX,EY,IV(INUP),DPHI)
              IF(MOD(IV(IILAG),2).NE.0)
     :          CALL ROTATE(ROTIX,EXLG,EYLG,IV(INUL),DPHI)
              IF(IV(IILAG).GE.2)
     :          CALL ROTATE(ROTIX,STELPX,STELPY,IV(INPM),DPHI)
              WRITE(LSC,121)PLON,PLAT,-PROT,DPHI
              WRITE(LUW,121)PLON,PLAT,-PROT,DPHI
  121         FORMAT('Reverse transformation defined by: UNDO',/,
     :               4F10.4,/)
            ENDIF
C
C      spherical projection of the coordinate systems is invoked here
C      by either keyword COROTATE or SPHERICAL read by GETDEFORMDATA
C
            IF(IV(IIDEFTYP).GE.110)THEN
              PHI0=CENTLNG
              IF(IV(IIDEFTYP).EQ.111)PHI0=0.0
              CALL PROJECTXY(EX,EY,PHI0,YMID,1,IV(INUP),IERROR)
              IF (IV(IILAG).EQ.1.OR.IV(IILAG).EQ.3)
     :          CALL PROJECTXY(EXLG,EYLG,PHI0,YMID,0,IV(INUL),IERROR)
              IF (IV(IILAG).EQ.1.OR.IV(IILAG).EQ.2)
     :          CALL PROJECTXY(STELPX,STELPY,PHI0,YMID,0,
     :                         IV(INPM)*IV(INSM),IERROR)
              IF(IV(IIDEFTYP).EQ.110)THEN
                RV(IXREFM)=CENTLNG       ! save variables needed to retransform
                RV(IYREFM)=YMID
              ELSE IF(IV(IIDEFTYP).EQ.111)THEN
                RV(IXREFM)=DPHI          ! save variables needed to retransform
                RV(IYREFM)=-PROT         ! GETPOLE returns radians but
                RV(IAREFM)=PLON          ! these variables saved in degrees
              ENDIF
            END IF !deftype 110
C
C     the AREA parameter is scaled if coordinates are converted from
C     degrees to radians.  The default length scale now is radius of
C     planet, for purpose of defining HLENSC etc.
C
            IF(IV(IIDEFTYP).GE.110)THEN
              RV(IELAREA)=RV(IELAREA)*DTOR*DTOR
            END IF !deftype 110 or greater

          ELSE IF (ITERM.EQ.INDXVISDENS) THEN
C
C           Update viscosity and density arrays if the initial values
C           of SE or RHO are changed
C
            DIFFOLD = RV(ITDIFF)
            ALPHAOLD = RV(IALPHA)
C           There does not seem to be a default TEMP so using the
C           ITEMP flag to trigger initialising TEMPT array
C           TEMPOLD = TEMP
            TEMP = 1.0
            VCOLD = RV(IVC)
            SEOLD = RV(ISE)
            RHOOLD = RHO
            CALL GETVISDENSDATA(INSTR,LENGTH,INDX2,RV(IVC),
     :                   RV(ISE),RHO,
     :                   IV(IITEMP),TEMP,RV(ITDIFF),RV(IALPHA),
     :                   RV(IGAMMA),RV(IBETA),RV(ITREF),IERROR)
            IF(RV(IGAMMA).NE.0.0)THEN
              WRITE(LUW,10100)RV(IGAMMA),RV(IBETA),RV(ITREF)
              WRITE(LSC,10100)RV(IGAMMA),RV(IBETA),RV(ITREF)
10100         FORMAT('Strain weakening implemented with Gamma =',F8.3,
     :                ' Beta =',F8.3,' Tref = ',F8.3,/)
            END IF
C
C         do not allow overwriting of visdens data if read in
C         Elle poly file
C         initialise the VHB array if not already done
C         by Elle poly file
C
            IF (IV(INELLEP).GT.0) THEN
              IF (RV(IVC).NE.VCOLD) THEN
                WRITE(LUW,10106)RV(IVC),VCOLD
                WRITE(LSC,10106)RV(IVC),VCOLD
10106           FORMAT('VISCDFLT read from basil input file ',F6.2,
     :               ' does not match value from .poly file',F6.2,/)
                IERROR=51
              ELSE IF (RV(ISE).NE.SEOLD) THEN
                WRITE(LUW,10107)RV(ISE),SEOLD
                WRITE(LSC,10107)RV(ISE),SEOLD
10107           FORMAT('SE read from basil input file ',F6.2,
     :               ' does not match value from .poly file',F6.2,/)
                IERROR=52
              END IF
            ELSE IF (IV(IIVIS).NE.0) THEN
              DO I=1,IV(INE)
                DO J=1,7
                  VHB(J,I) = RV(IVC)
                ENDDO
                VHB(8,I) = RV(ISE)
              ENDDO
            END IF
            IF (IERROR.NE.0)THEN
              IERROR=53
              GO TO 700
            ENDIF
            IF (RHO.NE.RHOOLD) THEN
              DO 6 J=1,IV(INE)
                DO 17 I=1,7
                  DENS(I,J)=RHO
   17           CONTINUE
    6         CONTINUE
            END IF
            IF (IV(IITEMP).NE.0) THEN
              DO 7 J=1,IV(INUP)
                  TEMPT(J)=TEMP
    7         CONTINUE
              WRITE(*,*)'Temperature array is now initialised to', TEMP
            END IF
            IDBUG=0
            REWIND(LBC)
            CALL VISVAR(LBC,RV(ISE),RV(IHLENSC),RV(IBDEPSC),
     :                  RV(IREFLEV),RV(IRISOST),EX,EY,LEM,NOR,VHB,
     :                  SSQ,IMAT,DENS,RV(IVC),CENTLNG,RV(IDEFV1),
     :                  IV(IIDEFTYP),IV(IIVIS),IV(IIDEN),IV(IIVV),
     :                  IV(IIMREG),IV(INUP),IV(INE),LUW,LSC,IDBUG,
     :                  IERROR)
            IF (IERROR.NE.0)THEN
              IERROR=54
              GO TO 700
            ENDIF
            REWIND(LBC)
            DO I=1,NUMLINE
             CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
            ENDDO

          ELSE IF (ITERM.EQ.INDXLAGRANGE) THEN
C
C           Flags are set during INITIALREAD;           Read or
C           set up Lagrangian mesh and strain markers, if specified
C
            IF (IV(IILAG).NE.0) THEN
              XIST=.FALSE.
              IF (ILGSAVE.NE.0)
     :          INQUIRE(FILE=LGFILE,EXIST=XIST)
              IF (XIST) THEN
                IREADMESH=1
                OPEN(LLG,FILE=LGFILE,STATUS='old',IOSTAT=IOS)
                CALL READLGMESH(IV(IILAG),LGEM,EXLG,EYLG,
     :                          LGIBC,LGIBCF,
     :                          STELPX,STELPY,RV(ISTELPR),
     :                          IV(INEL),IV(INUL),IV(INBL),
     :                          IV(INSM),IV(INPM),IV(INRM),
     :                          LLG,LSC,LUW,IERROR)
                CLOSE(LLG)
              ELSE
                IREADMESH=0
                CALL LGMESH(NXL,NYL,XMIN,YMIN,XMAX,YMAX,IV(IIFCASE),
     :                    LGEM,LGIBC,LGIBCF,EXLG,EYLG,EX,EY,
     :                    IBC,IBNGH,IBCTYP,NOR,
     :                    IV(INEL),IV(INUL),NNL,NML,IV(INBL),
     :                    IV(INUP),IV(INBP),IV(IIFLT),
     :                    LUW,IERROR)
                IF (ILGSAVE.NE.0) THEN
                  OPEN(LLG,FILE=LGFILE,STATUS='new',IOSTAT=IOS)
                  CLOSE(LLG)
                ENDIF
              ENDIF ! if xist
              IF (IERROR.NE.0) GO TO 700
              IF (IV(IILAG).EQ.1.OR.IV(IILAG).EQ.2) THEN
C            markers read in READLGMESH
C            fix so it works for markers only
                IF (IREADMESH.EQ.1) THEN
C                 CALL READLGMARKERS(STELPX,STELPY,
C      :                    RV(ISTELPR),IV(INSM),IV(INPM),IV(INRM),
C      :                    LLG,LSC,LUW,IERROR)
                  CLOSE(LLG)
                ELSE
                  REWIND(LBC)
                  CALL LGMARK(LBC,STELPX,STELPY,
     :                      RV(ISTELPR),IV(INSM),IV(INPM),IV(INRM),
     :                      LUW,LSC,IERROR)
                  IF (IERROR.NE.0)THEN
                    IERROR=55
                    GO TO 700
                  ENDIF
                  REWIND(LBC)
                  DO I=1,NUMLINE
                    CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
                  ENDDO
                END IF ! if readmesh
              END IF ! if ilag is 1 or 2
            END IF ! if ilag

          ELSE IF (ITERM.EQ.INDXSOLVE) THEN
            CALL GETSOLVEDATA(INSTR,LENGTH,INDX2,RV(IAC),RV(IACNL),
     :                        RV(IWFIT),NITER,IV(IITSTOP),
     :                        IERROR)
C
C      Set up the time series data to be output
C
          ELSE IF (ITERM.EQ.INDXSERIES) THEN
            IV(IMSINDX)=0
 100        CALL GETSERIESDATA(INSTR,LENGTH,INDX2,IV(IMSINDX),
     :                         MEASUR,PTLOC,IV(INCOMP),LELEMENTFLAG,
     :                         LUW,LSC,IERROR)
            IF (IERROR.EQ.0 .AND. INSTR(LENGTH:LENGTH).EQ.'&' .AND.
     :          IV(IMSINDX).LT.MAXMEAS) THEN
              CALL READLINE(INSTR,LENGTH,LBC,LSC,IEND)
              INDX2=1
              GO TO 100
            END IF
            IF (IERROR.NE.0)THEN
              IERROR=61
              GO TO 600
            ENDIF
            IV(IIDATA)=1
C
C     set defaults for SERIES command without parameters
C
            IF(IV(IMSINDX).EQ.0)THEN
              IV(IMSINDX)=5
              MEASUR(1)=INDXTIME
              MEASUR(2)=INDXXX
                PTLOC(1,2)=0.0
                PTLOC(2,2)=0.0
              MEASUR(3)=INDXUX
                PTLOC(1,3)=0.0
                PTLOC(2,3)=0.0
              MEASUR(4)=INDXXX
                PTLOC(1,4)=RV(IXLEN)
                PTLOC(2,4)=0.0
              MEASUR(5)=INDXUX
                PTLOC(1,5)=RV(IXLEN)
                PTLOC(2,5)=0.0
            END IF
            CALL GETNODENO(MEASUR,PTLOC,MSNODE,IV(IMSINDX),
     :                     EX,EY,NOR,IV(INUP),IV(INN),
     :                     LUW,LSC,IERROR)
            IF (LELEMENTFLAG.EQV..TRUE.)
     :        CALL GETELEMENTNO(MEASUR,PTLOC,MSNODE,IV(IMSINDX),
     :                           EX,EY,NOR,IV(INUP),IV(INN),
     :                           LEM,IV(INE),LUW,LSC,IERROR)
C           IF (IERROR.NE.0) GO TO 600
          END IF ! itermparse
        ELSE
          ITERM = MATCH(TIMESTEPTERMS,NUMTIME,WORD)
          IF (ITERM.GT.0) THEN
C
C            a timestep term
C
            IVALID=1
            IF (ITERM.EQ.INDXSTOP) THEN
              CALL GETSTOPDATA(INSTR,LENGTH,INDX2,
     :                         IV(IKEXIT),RV(ITEXIT),IV(IIWRITE),
     :                         XYEXIT,IERROR)
              ITIMESTEPORDER(NEXTT) = ITERM
              NEXTT = NEXTT+1
            ELSE IF (ITERM.EQ.INDXSAVE) THEN
              CALL GETSAVEDATA(INSTR,LENGTH,INDX2,IV(IKSAVE),
     :                         RV(ITSAVE),IERROR)
              ITIMESTEPORDER(NEXTT) = ITERM
              NEXTT = NEXTT+1
            ELSE IF (ITERM.EQ.INDXSTEPSIZE) THEN
              CALL GETSTEPSIZEDATA(INSTR,LENGTH,INDX2,IV(IIDT0),
     :                        IV(IMPDEF),NEGTIM,ITIMETYPE,IERROR)
              ITIMESTEPORDER(NEXTT) = ITERM
              NEXTT = NEXTT+1
            ELSE IF (ITERM.EQ.INDXDENSMOD) THEN
              IDENMOD = 1
              ITIMESTEPORDER(NEXTT) = ITERM
              NEXTT = NEXTT+1
            ELSE IF (ITERM.EQ.INDXBCONDMOD) THEN
              CALL GETBCMODDATA(INSTR,LENGTH,INDX2,IV(IINDFIX),
     :                          RV(IOMTOT),RV(IBCV1),IV(IIBCMOD),
     :                          GSDUDX,GSDUDY,GSDVDX,GSDVDY,
     :                          IERROR)
              ITIMESTEPORDER(NEXTT) = ITERM
              NEXTT = NEXTT+1
            ELSE IF (ITERM.EQ.INDXRHEOMOD) THEN
              IRHEOMOD = 1
              ITIMESTEPORDER(NEXTT) = ITERM
              NEXTT = NEXTT+1
              CALL GETRHEOMODDATA(INSTR,LENGTH,INDX2,RV(IVISP1),
     :             IV(IRHEOTYP),IV(IIVOLD),IERROR)  
C     since RHEOMOD modifies VHB, set IV(IIVIS) to 1 to insure VHB is 
C     saved (ie.. if no modifications to VHB by REG command)
              IV(IIVIS) = 1
C     save original VHB in VOLD if needed for this RHEOTYP (MIB2002)
              IF (IV(IIVOLD).EQ.1)  THEN
                DO J=1,NE
                  DO K=1,8
                    VOLD(K,J) = VHB(K,J)
                  ENDDO
                ENDDO
              END IF
            ELSE IF (ITERM.EQ.INDXREMESH) THEN
              CALL GETREMESHDATA(INSTR,LENGTH,INDX2,MFIXFLAG,IGRID,
     :                           SHAPEMIN,RV(IELSEG),ANGMIN,IERROR)
C             RV(IELSEG)=SEGMIN
C             IF(IV(INCOMP).EQ.-1)RV(IELSEG)=SEGMIN*DTOR
C
C Hack to allow remeshing with MFIX for flt=3 (for Joerg)
C             IF (IGRID.EQ.1.AND.(IV(IIFLT).EQ.3.OR.IV(IIFLT).EQ.4))
C    :            MFIXFLAG=1
              ITIMESTEPORDER(NEXTT) = ITERM
              NEXTT = NEXTT+1
            END IF !itermparse
          END IF !iterm>0 - timestepterm
          IF (IERROR.NE.0)THEN
            IERROR=71
            GO TO 600
          ENDIF
        END IF !iterm>0 - setupterm
      END IF !indx2>0
      IF (IERROR.EQ.0) GO TO 110
C
C    Echo the instruction if it could not be parsed
C
 600  IF (IERROR.GT.9) CALL WRITEINSTR(INSTR,LENGTH,INDX2)
 700  RETURN
10104 FORMAT('Solution for KSTOP =',I6,'  T =',f7.3,
     1' read from unit',i3)
10114 FORMAT('Run commencing from the file ',A16,/)
10124 FORMAT('GEOMETRY: error in longitude, ',F6.2,' to ',F6.2,
     : ' or latitude, ',F6.2,' to ',F6.2,/)
10125 FORMAT('The read file could not be opened ',A60,/)
      END

      SUBROUTINE GETOUTPUTDATA(INSTR,LEN,INDX2,NAMEW,LNWF,
     :                   NAMEDB,LNDB,BINDIR,LNBD,OUTDIR,LNOD,
     :                   IPRINT,IERR)
      INCLUDE "limits.parameters"
      CHARACTER NAMEW*32,BINDIR*32,OUTDIR*32,NAMEDB*32
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXBIN=1,
     :                     INDXBINDIR=2,
     :                     INDXOUTDIR=3,
     :                     INDXDUMP=4,
     :                     INDXPRINT=5,
     :                     NUMOUTPUTTERMS=5)
      CHARACTER OUTPUTTERMS(1:NUMOUTPUTTERMS)*12
      DATA OUTPUTTERMS   /'BIN         ',
     :                    'BINDIR      ',
     :                    'OUTDIR      ',
     :                    'DUMP        ',
     :                    'VERBOSE     '/
C
C    read the next word on the line
C
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(OUTPUTTERMS,NUMOUTPUTTERMS,WORD)
C
C    if a valid word is found set the corresponding value
C
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXBIN) THEN
          CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
          IF (INDX1.EQ.0) GO TO 210
          INDX2 = INDX1+1
          CALL GETWORD(INSTR,LEN,INDX2,NAMEW,LNWF,32)
        ELSE IF (ITERM.EQ.INDXBINDIR) THEN
          CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
          IF (INDX1.EQ.0) GO TO 210
          INDX2 = INDX1+1
          CALL GETWORD(INSTR,LEN,INDX2,BINDIR,LNBD,32)
        ELSE IF (ITERM.EQ.INDXOUTDIR) THEN
          CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
          IF (INDX1.EQ.0) GO TO 210
          INDX2 = INDX1+1
          CALL GETWORD(INSTR,LEN,INDX2,OUTDIR,LNOD,32)
        ELSE IF (ITERM.EQ.INDXDUMP) THEN
          CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
          IF (INDX1.EQ.0) GO TO 210
          INDX2 = INDX1+1
          CALL GETWORD(INSTR,LEN,INDX2,NAMEDB,LNDB,32)
        ELSE IF (ITERM.EQ.INDXPRINT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IPRINT
        END IF
C
C     if not the end of the line, look for another word
C
        IF(INDX2.NE.0)GO TO 5
        INDX2=INDX1
      ENDIF
      RETURN
C
C      in case line can't be interpreted
C
 200  INDX2=INDX1
 210  IERR = 10
      RETURN
      END

      SUBROUTINE GETREADDATA(INSTR,LEN,INDX2,NAMER,LNR,
     :                       KSTART,ISTR,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER NAMER*32
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXBIN=1,
     :                     INDXRECORD=2,
     :                     INDXRESET=3,
     :                     INDXCONTINUE=4,
     :                     NUMREADTERMS=4)
      CHARACTER READTERMS(1:NUMREADTERMS)*12
      DATA READTERMS     /'BIN         ',
     :                    'RECORD      ',
     :                    'RESET       ',
     :                    'CONTINUE    '/
C
C   RESET is now the default setting, if a READ command is used
C    including RESET option has no specific effect.  including
C    CONTINUE option ensures that time counters not reset
C
      ISTR=1
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(READTERMS,NUMREADTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXBIN) THEN
          CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
          IF (INDX1.EQ.0) GO TO 210
          INDX2 = INDX1+1
          CALL GETWORD(INSTR,LEN,INDX2,NAMER,LNR,32)
          IF (INDX2.EQ.0) GO TO 200
C         IST=1
        ELSE IF (ITERM.EQ.INDXRECORD) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)KSTART
C       ELSE IF (ITERM.EQ.INDXRESET) THEN
C         ISTR=1
        ELSE IF (ITERM.EQ.INDXCONTINUE) THEN
          ISTR=2
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END

      SUBROUTINE GETGEOMETRYDATA(INSTR,LEN,INDX2,XLEN,YLEN,
     :                  XZERO,YZERO,CENTLNG,NCOMP,IGRAV,IERR) 
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXXLEN=1,
     :                     INDXYLEN=2,
     :                     INDXNCOMP=3,
     :                     INDXIGRAV=4,
     :                     INDXXZERO=5,
     :                     INDXYZERO=6,
     :                     INDXCENTLNG=7,
     :                     NUMGEOMETRYTERMS=7)
      CHARACTER GEOMETRYTERMS(1:NUMGEOMETRYTERMS)*12
      DATA GEOMETRYTERMS /'XLEN        ',
     :                    'YLEN        ',
     :                    'NCOMP       ',
     :                    'IGRAV       ',
     :                    'XZERO       ',
     :                    'YZERO       ',
     :                    'CENTRELONG  '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(GEOMETRYTERMS,NUMGEOMETRYTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXXLEN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)XLEN
        ELSE IF (ITERM.EQ.INDXYLEN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)YLEN
        ELSE IF (ITERM.EQ.INDXNCOMP) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NCOMP
        ELSE IF (ITERM.EQ.INDXIGRAV) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IGRAV
        ELSE IF (ITERM.EQ.INDXXZERO) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)XZERO
        ELSE IF (ITERM.EQ.INDXYZERO) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)YZERO
        ELSE IF (ITERM.EQ.INDXCENTLNG) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)CENTLNG
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END

      SUBROUTINE GETFORCEDATA(INSTR,LEN,INDX2,RHOG,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      DIMENSION RHOG(5)
      PARAMETER(           INDXRHOG=1,
     :                     NUMFORCETERMS=1)
      CHARACTER FORCETERMS(1:NUMFORCETERMS)*12
      DATA FORCETERMS    /'RHOG        '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(FORCETERMS,NUMFORCETERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXRHOG) THEN
          CALL GETARRAY(INSTR,LEN,INDX2,RHOG,5,ICNT,IERR)
          IF (IERR.NE.0) GO TO 220
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
 220  RETURN
      END

      SUBROUTINE GETDEFORMDATA(INSTR,LEN,INDX2,IDEFTYP,BANGL,ERA,
     :                   DEFV,TILTDEG,NDFORM,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      DIMENSION DEFV(5)
      PARAMETER(           INDXTYPE=1,
     :                     INDXBANGL=2,
     :                     INDXERA=3,
     :                     INDXDEFV=4,
     :                     INDXTILTDEG=5,
     :                     INDXNDFORM=6,
     :                     INDXSPHERE=7,
     :                     INDXCOROTATE=8,
     :                     NUMDEFORMTERMS=8)
C    :                     NUMDEFORMTERMS=7)
      CHARACTER DEFORMTERMS(1:NUMDEFORMTERMS)*12
      DATA DEFORMTERMS   /'TYPE        ',
     :                    'BANGL       ',
     :                    'ERA         ',
     :                    'DEFV        ',
     :                    'TILTDEG     ',
     :                    'NDFORM      ',
     :                    'SPHERICAL   ',
     :                    'COROTATE    '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(DEFORMTERMS,NUMDEFORMTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 200
        ELSE IF (ITERM.EQ.INDXTYPE) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IDEFTYP
        ELSE IF (ITERM.EQ.INDXBANGL) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)BANGL
        ELSE IF (ITERM.EQ.INDXERA) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)ERA
        ELSE IF (ITERM.EQ.INDXDEFV) THEN
          CALL GETARRAY(INSTR,LEN,INDX2,DEFV,5,ICNT,IERR)
          IF (IERR.NE.0) GO TO 220
        ELSE IF (ITERM.EQ.INDXTILTDEG) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)TILTDEG
          IDEFTYP=100
        ELSE IF (ITERM.EQ.INDXNDFORM) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NDFORM
          IDEFTYP=101
        ELSE IF (ITERM.EQ.INDXSPHERE) THEN
          IDEFTYP=110       ! change sybil if IDEFTYP changed
        ELSE IF (ITERM.EQ.INDXCOROTATE) THEN
          IDEFTYP=111
          CALL GETARRAY(INSTR,LEN,INDX2,DEFV,4,ICNT,IERR)
          IF ((IERR.NE.0).OR.(ICNT.NE.4)) GO TO 200
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
 220  RETURN
      END
C
      SUBROUTINE GETVISDENSDATA(INSTR,LEN,INDX2,VC,SE,RHO,
     :                          ITEMP,TEMP,DIFF,
     :                          ALPHA,GAMMAT,BETA,TREF,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(     INDXSE=1,INDXRHO=2,INDXTEMP=3,
     :               INDXDIFF=4,INDXALPHA=5,INDXGAM=6,
     :               INDXBET=7,INDXTREF=8,INDXVC=9,
     :               NUMVISDENSTERMS=9)
      CHARACTER VISDENSTERMS(1:NUMVISDENSTERMS)*12
      DATA VISDENSTERMS /'SE          ','RHO         ',
     :                   'TEMP        ','DIFF        ',
     :                   'ALPHA       ','GAMMA       ',
     :                   'BETA        ','TREF        ',
     :                   'VC          '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(VISDENSTERMS,NUMVISDENSTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXSE) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)SE
        ELSE IF (ITERM.EQ.INDXRHO) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)RHO
        ELSE IF (ITERM.EQ.INDXTEMP) THEN
          ITEMP=1
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)TEMP
        ELSE IF (ITERM.EQ.INDXDIFF) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)DIFF
        ELSE IF (ITERM.EQ.INDXALPHA) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)ALPHA
        ELSE IF (ITERM.EQ.INDXGAM) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)GAMMAT
        ELSE IF (ITERM.EQ.INDXBET) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)BETA
        ELSE IF (ITERM.EQ.INDXTREF) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)TREF
        ELSE IF (ITERM.EQ.INDXVC) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)VC
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
 220  RETURN
      END

      SUBROUTINE GETLAYERDATA(INSTR,LEN,INDX2,ICR,ARGAN,BRGAN,
     :                THRESH,HLENSC,BDEPSC,RISOST,REFLEV,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXTHICKNESS=1,
     :                     INDXROTATION=2,
     :                     INDXARGAN=3,
     :                     INDXBRGAN=4,
     :                     INDXTHRESH=5,
     :                     INDXHLENSC=6,
     :                     INDXBDEPSC=7,
     :                     INDXRADIUS=8,
     :                     INDXRISOST=9,
     :                     INDXREFLEV=10,
     :                     NUMLAYERTERMS=10)
      CHARACTER LAYERTERMS(1:NUMLAYERTERMS)*12
      DATA LAYERTERMS    /'THICKNESS   ',
     :                    'ROTATION    ',
     :                    'ARGAN       ',
     :                    'BRGAN       ',
     :                    'THRESH      ',
     :                    'HLENSC      ',
     :                    'BDEPSC      ',
     :                    'RADIUS      ',     ! appears to be unused
     :                    'RISOST      ',
     :                    'REFLEV      '/
      ICR = 0
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(LAYERTERMS,NUMLAYERTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXTHICKNESS) THEN
          IF (ICR.EQ.3) THEN
            ICR = 1
          ELSE
            ICR = 2
          END IF
        ELSE IF (ITERM.EQ.INDXROTATION) THEN
          IF (ICR.EQ.2) THEN
            ICR = 1
          ELSE
            ICR = 3
          END IF
        ELSE IF (ITERM.EQ.INDXARGAN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)ARGAN
        ELSE IF (ITERM.EQ.INDXBRGAN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)BRGAN
        ELSE IF (ITERM.EQ.INDXTHRESH) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)THRESH
        ELSE IF (ITERM.EQ.INDXHLENSC) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)HLENSC
        ELSE IF (ITERM.EQ.INDXBDEPSC) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)BDEPSC
        ELSE IF (ITERM.EQ.INDXRISOST) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)RISOST
        ELSE IF (ITERM.EQ.INDXREFLEV) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)REFLEV
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      IF (ICR.EQ.0) ICR = 1
 220  RETURN
      END

      SUBROUTINE GETTOPODATA(INSTR,LEN,INDX2,TOPOFILE,ITYPE,
     :                       HMIN,HLEN,VMIN,VLEN,
     :                       IFILTER,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER TOPOFILE*60

      PARAMETER(           INDXFILE=1,
     :                     INDXFTYPE=2,
     :                     INDXFILTER=3,
     :                     INDXLONGMIN=4,
     :                     INDXLONGLEN=5,
     :                     INDXLATMIN=6,
     :                     INDXLATLEN=7,
     :                     NUMTOPOTERMS=7)
      CHARACTER TOPOTERMS(1:NUMTOPOTERMS)*12
      CHARACTER TMPFILE*60
      DATA TOPOTERMS     /'FILE        ',
     :                    'TYPE        ',
     :                    'FILTER      ',
     :                    'LONGMIN     ',
     :                    'LONGLEN     ',
     :                    'LATMIN      ',
     :                    'LATLEN      '/
      DO I=1,60
        TMPFILE(I:I)=' '
        TOPOFILE(I:I)=' '
      ENDDO
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(TOPOTERMS,NUMTOPOTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXFILE) THEN
          CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
          IF (INDX1.EQ.0) GO TO 210
          INDX2 = INDX1+1
          CALL GETWORD(INSTR,LEN,INDX2,TMPFILE,LNW,60)
          TOPOFILE(1:INDX2-INDX1-1)=TMPFILE(1:INDX2-INDX1-1)
        ELSE IF (ITERM.EQ.INDXFTYPE) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)ITYPE
        ELSE IF (ITERM.EQ.INDXFILTER) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IFILTER
        ELSE IF (ITERM.EQ.INDXLONGMIN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)HMIN
        ELSE IF (ITERM.EQ.INDXLONGLEN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)HLEN
        ELSE IF (ITERM.EQ.INDXLATMIN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)VMIN
        ELSE IF (ITERM.EQ.INDXLATLEN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)VLEN
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      IF (ICR.EQ.0) ICR = 1
 220  RETURN
      END

      SUBROUTINE GETLAGRANGEDATA(INSTR,LEN,INDX2,ILAG,NXL,NYL,
     :                           NPM,NRM,LGFILE,ILGSAVE,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER LGFILE*80
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXLGMESH=1,
     :                     INDXMARKER=2,
     :                     INDXNXL=3,
     :                     INDXNYL=4,
     :                     INDXNPMP=5,
     :                     INDXNRMP=6,
     :                     INDXFILE=7,
     :                     NUMLAGRANGETERMS=7)
      CHARACTER LAGRANGETERMS(1:NUMLAGRANGETERMS)*12
      DATA LAGRANGETERMS /'LGMESH      ',
     :                    'MARKERS     ',
     :                    'NXL         ',
     :                    'NYL         ',
     :                    'NPMP        ',
     :                    'NRMP        ',
     :                    'FILE        '/
C
C    ILAG is zero if no Lagrange marks, = 1 if a Lagrange mesh added
C    = 2 if Lagrange markers used, = 3 if both are used
C
      ILAG = 0
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(LAGRANGETERMS,NUMLAGRANGETERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXLGMESH) THEN
          IF (ILAG.EQ.2) THEN
            ILAG = 1
          ELSE
            ILAG = 3
          END IF
        ELSE IF (ITERM.EQ.INDXMARKER) THEN
          IF (ILAG.EQ.3) THEN
            ILAG = 1
          ELSE
            ILAG = 2
          END IF
        ELSE IF (ITERM.EQ.INDXNXL) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NXL
        ELSE IF (ITERM.EQ.INDXNYL) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NYL
        ELSE IF (ITERM.EQ.INDXNPMP) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NPMP
        ELSE IF (ITERM.EQ.INDXNRMP) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NRMP
        ELSE IF (ITERM.EQ.INDXFILE) THEN
          CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
          IF (INDX1.EQ.0) GO TO 210
          INDX2 = INDX1+1
          CALL GETWORD(INSTR,LEN,INDX2,LGFILE,LNW,80)
          ILGSAVE=1
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      IF (ILAG.EQ.0) ILAG = 1
      RETURN
      END

      SUBROUTINE GETSOLVEDATA(INSTR,LEN,INDX2,AC,ACNL,WFIT,
     :                               NITER,ITSTOP,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXAC=1,
     :                     INDXACNL=2,
     :                     INDXWFIT=3,
     :                     INDXNITER=4,
     :                     INDXITSTOP=5,
     :                     NUMSOLVETERMS=5)
      CHARACTER SOLVETERMS(1:NUMSOLVETERMS)*12
      DATA SOLVETERMS    /'AC          ',
     :                    'ACNL        ',
     :                    'WFIT        ',
     :                    'NITER       ',
     :                    'ITSTOP      '/
   5  ISTART=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(SOLVETERMS,NUMSOLVETERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXAC) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)AC
        ELSE IF (ITERM.EQ.INDXACNL) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)ACNL
        ELSE IF (ITERM.EQ.INDXWFIT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)WFIT
        ELSE IF (ITERM.EQ.INDXNITER) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NITER
        ELSE IF (ITERM.EQ.INDXITSTOP) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)ITSTOP
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END

      SUBROUTINE GETSERIESDATA(INSTR,LEN,INDX2,MSINDX,MEASUR,
     :                            PTLOC,NCOMP,LELEFLAG,LUW,LSC,IERR)
      INCLUDE "limits.parameters"
      INCLUDE "input.parameters"
      INCLUDE "input.data"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
C
C     SERIESTERMS are in input.data.
C     NUMSERIESTERMS and indices are in input.parameters
C     refer to these two files to see how parameter names are related
C     to index numbers in the array accessed by MEASUR.  Variables
C     are of two types: those that do not require a coordinate location
C     (e.g. TIME or MIN or MAX of any field variable), and those that
C     are specific to a location. PTLOC records that (x,y) location
C
      DIMENSION PTLOC(3,MAXMEAS)
      DIMENSION MEASUR(MAXMEAS)
      LOGICAL LELEFLAG
      LELEFLAG=.FALSE.
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(SERIESTERMS,NUMSERIESTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 200
        ELSE IF ((ITERM.EQ.INDXPRES.OR.ITERM.EQ.INDXPRESMIN.OR.
     :    ITERM.EQ.INDXPRESMAX).AND.NCOMP.EQ.0) THEN
          WRITE(LSC,*)'Pressure not calculated for SERIES',
     1    ' command, NCOMP=0'
          WRITE(LUW,*)'Pressure not calculated for SERIES',
     1    ' command, NCOMP=0'
          GO TO 5
        ELSE
          MSINDX=MSINDX+1
          IF(MSINDX.GT.MAXMEAS)THEN
            WRITE(LSC,10200)MAXMEAS
            WRITE(LUW,10200)MAXMEAS
            GO TO 220
          END IF
          MEASUR(MSINDX)=ITERM
        END IF
C       IF (ITERM.EQ.INDXXX) THEN
        IF (ITERM.GE.INDXXX) THEN
          CALL GETARRAY(INSTR,LEN,INDX2,PTLOC(1,MSINDX),2,ICNT,IERR)
          IF ((IERR.NE.0).OR.(ICNT.NE.2)) GO TO 200
          IF(ITERM.GE.INDXEDXX)LELEFLAG=.TRUE.
        ENDIF
C       ELSE IF (ITERM.EQ.INDXYY) THEN
C         CALL GETARRAY(INSTR,LEN,INDX2,PTLOC(1,MSINDX),2,ICNT,IERR)
C         IF ((IERR.NE.0).OR.(ICNT.NE.2)) GO TO 200
C       ELSE IF (ITERM.EQ.INDXUX) THEN
C         CALL GETARRAY(INSTR,LEN,INDX2,PTLOC(1,MSINDX),2,ICNT,IERR)
C         IF ((IERR.NE.0).OR.(ICNT.NE.2)) GO TO 200
C       ELSE IF (ITERM.EQ.INDXUY) THEN
C         CALL GETARRAY(INSTR,LEN,INDX2,PTLOC(1,MSINDX),2,ICNT,IERR)
C         IF ((IERR.NE.0).OR.(ICNT.NE.2)) GO TO 200
C       ELSE IF (ITERM.EQ.INDXTH) THEN
C         CALL GETARRAY(INSTR,LEN,INDX2,PTLOC(1,MSINDX),2,ICNT,IERR)
C         IF ((IERR.NE.0).OR.(ICNT.NE.2)) GO TO 200
C       ELSE IF (ITERM.EQ.INDXRO) THEN
C         CALL GETARRAY(INSTR,LEN,INDX2,PTLOC(1,MSINDX),2,ICNT,IERR)
C         WRITE(*,*)'INDXRO: PTLOC =',PTLOC(1,MSINDX),PTLOC(2,MSINDX)
C         IF ((IERR.NE.0).OR.(ICNT.NE.2)) GO TO 200
C       ELSE IF (ITERM.EQ.INDXPRES) THEN
C         CALL GETARRAY(INSTR,LEN,INDX2,PTLOC(1,MSINDX),2,ICNT,IERR)
C         IF ((IERR.NE.0).OR.(ICNT.NE.2)) GO TO 200
C       ELSE IF (ITERM.GE.INDXEDXX.AND.ITERM.LE.INDXFOLT) THEN
C         CALL GETARRAY(INSTR,LEN,INDX2,PTLOC(1,MSINDX),2,ICNT,IERR)
C         IF ((IERR.NE.0).OR.(ICNT.NE.2)) GO TO 200
C         LELEFLAG=.TRUE.
C       END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
 220  CONTINUE
      RETURN
10200 FORMAT('Maximum of ',I3,' measurements allowed in SERIES command')
      END

      SUBROUTINE WRITESERIESDATA(LEM,NOR,NE,NUP,MSNODE,MEASUR,MSINDX,
     :                            PTLOC,EX,EY,TIME,LTMPDAT,IERR)
C
C     this routine reconstructs a SERIES statement for use after
C     regridding
C    
      INCLUDE "limits.parameters"
      INCLUDE "input.parameters"
      INCLUDE "input.data"
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER IWORD*12
C
C     SERIESTERMS are in input.data.
C     NUMSERIESTERMS and indices are in input.parameters
C     refer to these two files to see how parameter names are related
C     to index numbers in the array accessed by MEASUR.  Variables
C     are of two types: those that do not require a coordinate location
C     (e.g. TIME or MIN or MAX of any field variable), and those that
C     are specific to a location. PTLOC records that (x,y) location
C
      DIMENSION PTLOC(3,MAXMEAS)
      DIMENSION MEASUR(MAXMEAS)
C      local arrays
      DIMENSION TMPPTLOC(3,MAXMEAS)
C
      DO N=1,MAXMEAS
        DO I=1,3
        TMPPTLOC(I,N) = PTLOC(I,N)
        END DO
      END DO
      CALL RESETSERIES(LEM,NOR,NE,NUP,MEASUR,MSNODE,MSINDX,
     :                 TMPPTLOC,EX,EY)
      DO I=1,12
        IWORD(I:I)=' '
      END DO
      DO J=1,80
        INSTR(J:J)=' '
      END DO
      WRITE(INSTR,101) TIME,' '
      INDX1=9
      DO I=1,MSINDX
C       IF( (MEASUR(I).EQ.INDXTIME).OR.
C    :      (MEASUR(I).GE.INDXXMIN.AND.MEASUR(I).LE.INDXVMAX).OR.
C    :      (MEASUR(I).EQ.INDXTHMIN.OR.MEASUR(I).EQ.INDXTHMAX).OR.
C    :      (MEASUR(I).EQ.INDXPRESMIN.OR.MEASUR(I).EQ.INDXPRESMAX).OR.
C    :      (MEASUR(I).EQ.INDXROMIN.OR.MEASUR(I).EQ.INDXROMAX)) THEN
C
C     for parameters that don't require coordinates
C
        IF(MEASUR(I).LT.INDXXX)THEN
          IWORD(1:12)=SERIESTERMS(MEASUR(I))(1:12)
          INDX2=INDEX(IWORD,' ')
          IF (INDX2.EQ.0) INDX2=12
          INSTR(INDX1:INDX1+INDX2)=IWORD(1:INDX2)
          INDX1=INDX1+INDX2
          INSTR(INDX1:INDX1+1)=', '
          INDX1=INDX1+2
          IF (INDX1.GE.(LINELEN-14)) THEN
            INSTR(INDX1:INDX1+2)=' &'
            INDX1 = INDX1+2
            INDX1 = 1
          END IF
C
C     for parameters that require coordinates
C
        ELSE 
          IWORD(1:12)=SERIESTERMS(MEASUR(I))(1:12)
          INDX2=INDEX(IWORD,' ')
          IF (INDX2.EQ.0) INDX2=12
          INSTR(INDX1:INDX1+INDX2)=IWORD(1:INDX2)
          INDX1=INDX1+INDX2
          INSTR(INDX1:INDX1)='='
          INDX1=INDX1+1
          IF (INDX1.GE.(LINELEN-14)) THEN
            INSTR(INDX1:INDX1+2)=' &'
            INDX1 = INDX1+2
            WRITE(LTMPDAT,*)INSTR(1:INDX1)    !check file is open
            INDX1 = 1
          END IF
          WRITE(IWORD,102)TMPPTLOC(1,I)
          CALL SKIP(IWORD,12,1,'n',INDX2)
          IF (INDX2.EQ.0) INDX2=1
          CALL SKIP(IWORD,12,INDX2,'b',INDX3)
          IF (INDX3.EQ.0) INDX3=12
          LEN = INDX3-INDX2
          INSTR(INDX1:INDX1+LEN)= IWORD(INDX2:INDX3)
          INDX1=INDX1+LEN
          INSTR(INDX1:INDX1+1)=', '
          INDX1=INDX1+2
          IF (INDX1.GE.(LINELEN-14)) THEN
            INSTR(INDX1:INDX1+2)=' &'
            INDX1 = INDX1+2
            WRITE(LTMPDAT,*)INSTR(1:INDX1)
            INDX1 = 1
          END IF
          WRITE(IWORD,102)TMPPTLOC(2,I)
          CALL SKIP(IWORD,12,1,'n',INDX2)
          IF (INDX2.EQ.0) INDX2=1
          CALL SKIP(IWORD,12,INDX2,'b',INDX3)
          IF (INDX3.EQ.0) INDX3=12
          LEN = INDX3-INDX2
          INSTR(INDX1:INDX1+LEN)= IWORD(INDX2:INDX3)
          INDX1=INDX1+LEN
          INSTR(INDX1:INDX1+3)=' ; '
          INDX1=INDX1+3
          IF (INDX1.GE.(LINELEN-14)) THEN
            IF (I.NE.MSINDX) THEN
              INSTR(INDX1:INDX1+3)=' &'
              INDX1 = INDX1+3
            END IF
            WRITE(LTMPDAT,*)INSTR(1:INDX1)
            INDX1 = 1
          END IF
        END IF
      END DO
      IF (INDX1.GT.1) WRITE(LTMPDAT,*)INSTR(1:INDX1)
      RETURN
  101 FORMAT(F7.3,A1)
  102 FORMAT(G10.4)
      END

      SUBROUTINE GETBCONDDATA(INSTR,LEN,INDX2,IOFF,
     :                         VELXO,VELYO,YLDSTR,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXIOFF=1,
     :                     INDXVELXO=2,
     :                     INDXVELYO=3,
     :                     INDXYLDSTR=4,
     :                     NUMBCONDTERMS=4)
      CHARACTER BCONDTERMS(1:NUMBCONDTERMS)*12
      DATA BCONDTERMS    /'IOFF        ',
     :                    'VELXO       ',
     :                    'VELYO       ',
     :                    'YLDSTR      '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(BCONDTERMS,NUMBCONDTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXIOFF) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IOFF
        ELSE IF (ITERM.EQ.INDXVELXO) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)VELXO
        ELSE IF (ITERM.EQ.INDXVELYO) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)VELYO
        ELSE IF (ITERM.EQ.INDXYLDSTR) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)YLDSTR
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END

      SUBROUTINE GETPOLYDATA(INSTR,LEN,INDX2,POLYFILE,IERR)
      INCLUDE "limits.parameters"
      CHARACTER POLYFILE*80
      CHARACTER TMPFILE*80
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXFILE=1,
     :                     NUMPOLYTERMS=1)
      CHARACTER POLYTERMS(1:NUMPOLYTERMS)*12
      DATA POLYTERMS   /'FILE        '/

   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(POLYTERMS,NUMPOLYTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXFILE) THEN
          CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
          IF (INDX1.EQ.0) GO TO 210
          INDX2 = INDX1+1
          CALL GETWORD(INSTR,LEN,INDX2,TMPFILE,LNW,80)
          POLYFILE(1:INDX2-INDX1-1)=TMPFILE(1:INDX2-INDX1-1)
        END IF
        IF (INDX2.NE.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END

      SUBROUTINE GETSTEPSIZEDATA(INSTR,LEN,INDX2,IDT0,MPDEF,
     :                             NEGTIM,ITIMETYPE,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXTYPE=1,
     :                     INDXIDT0=2,
     :                     INDXMPDEF=3,
     :                     INDXNEGTIM=4,
     :                     NUMSTEPSIZETERMS=4)
      CHARACTER STEPSIZETERMS(1:NUMSTEPSIZETERMS)*12
      DATA STEPSIZETERMS /'TYPE        ',
     :                    'IDT0        ',
     :                    'MPDEF       ',
     :                    'NEGTIM      '/
      PARAMETER(           INDXEXPLICIT=1,
     :                     INDXRK=2,
     :                     NUMTYPETERMS=2)
      CHARACTER TYPETERMS(1:NUMTYPETERMS)*12
      DATA TYPETERMS     /'EXPLICIT    ',
     :                    'RK          '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(STEPSIZETERMS,NUMSTEPSIZETERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXTYPE) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          WORD(1:12)=INSTR(INDX1:INDX2-1)
          ITERM = MATCH(TYPETERMS,NUMTYPETERMS,WORD)
          IF (ITERM.EQ.0) GO TO 200
          ITIMETYPE=ITERM
        ELSE IF (ITERM.EQ.INDXIDT0) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IDT0
        ELSE IF (ITERM.EQ.INDXMPDEF) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)MPDEF
        ELSE IF (ITERM.EQ.INDXNEGTIM) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)NEGTIM
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END

      SUBROUTINE GETSAVEDATA(INSTR,LEN,INDX2,KSAVE,TSAVE,
     :                                         IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXKSAVE=1,
     :                     INDXTSAVE=2,
     :                     NUMSAVETERMS=2)
      CHARACTER SAVETERMS(1:NUMSAVETERMS)*12
      DATA SAVETERMS     /'KSAVE       ',
     :                    'TSAVE       '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(SAVETERMS,NUMSAVETERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXKSAVE) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)KSAVE
        ELSE IF (ITERM.EQ.INDXTSAVE) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)TSAVE
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END

      SUBROUTINE GETSTOPDATA(INSTR,LEN,INDX2,KEXIT,TEXIT,IWRITE,
     :                               XYEXIT,IERR)
      INCLUDE "limits.parameters"
      INCLUDE "indices.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXKEXIT=1,
     :                     INDXTEXIT=2,
     :                     INDXIWRITE=3,
     :                     INDXXMINEXIT=4,
     :                     INDXXMAXEXIT=5,
     :                     INDXYMINEXIT=6,
     :                     INDXYMAXEXIT=7,
     :                     INDXYEXIT=8,
     :                     NUMSTOPTERMS=8)
      CHARACTER STOPTERMS(1:NUMSTOPTERMS)*12
      DATA STOPTERMS     /'KEXIT       ',
     :                    'TEXIT       ',
     :                    'IWRITE      ',
     :                    'XMIN        ',
     :                    'XMAX        ',
     :                    'YMIN        ',
     :                    'YMAX        ',
     :                    'YEXIT       '/
      DIMENSION XYEXIT(4)
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(STOPTERMS,NUMSTOPTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXKEXIT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)KEXIT
        ELSE IF (ITERM.EQ.INDXTEXIT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)TEXIT
        ELSE IF (ITERM.EQ.INDXIWRITE) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IWRITE
        ELSE IF (ITERM.EQ.INDXXMINEXIT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)XYEXIT(IXMIN)
        ELSE IF (ITERM.EQ.INDXXMAXEXIT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)XYEXIT(IXMAX)
        ELSE IF (ITERM.EQ.INDXYMINEXIT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)XYEXIT(IYMIN)
        ELSE IF (ITERM.EQ.INDXYEXIT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)XYEXIT(IYMAX)
        ELSE IF (ITERM.EQ.INDXYEXIT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)XYEXIT(IYMAX)
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END

      SUBROUTINE GETBCMODDATA(INSTR,LEN,INDX2,INDFIX,OMTOT,
     :                        BCV,IBCMOD,DUDX,DUDY,DVDX,DVDY,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      DIMENSION BCV(5)
      PARAMETER(           INDXINDFIX=1,
     :                     INDXOMTOT=2,
     :                     INDXBCV=3,
     :                     INDXTYPE=4,
     :                     INDXDUDX=5,
     :                     INDXDVDX=6,
     :                     INDXDUDY=7,
     :                     INDXDVDY=8,
     :                     NUMBCMODTERMS=8)
      CHARACTER BCMODTERMS(1:NUMBCMODTERMS)*12
      DATA BCMODTERMS    /'INDFIX      ',
     :                    'OMTOT       ',
     :                    'BCV         ',
     :                    'TYPE        ',
     :                    'DUDX        ',
     :                    'DVDX        ',
     :                    'DUDY        ',
     :                    'DVDY        '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(BCMODTERMS,NUMBCMODTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 210
        ELSE IF (ITERM.EQ.INDXINDFIX) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)INDFIX
        ELSE IF (ITERM.EQ.INDXOMTOT) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)OMTOT
        ELSE IF (ITERM.EQ.INDXTYPE) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IBCMOD
        ELSE IF (ITERM.EQ.INDXDUDX) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)DUDX
        ELSE IF (ITERM.EQ.INDXDVDX) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)DVDX
        ELSE IF (ITERM.EQ.INDXDUDY) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)DUDY
        ELSE IF (ITERM.EQ.INDXDVDY) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)DVDY
        ELSE IF (ITERM.EQ.INDXBCV) THEN
          CALL GETARRAY(INSTR,LEN,INDX2,BCV,5,ICNT,IERR)
          IF (IERR.NE.0) GO TO 220
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
 220  RETURN
      END

      SUBROUTINE GETREMESHDATA(INSTR,LEN,INDX2,MFIXFLAG,IGRID,
     :                                  SHAPEMIN,SEGMIN,ANGMIN,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      PARAMETER(           INDXIGRID=1,
     :                     INDXSHAPEMIN=2,
     :                     INDXSEGMIN=3,
     :                     INDXANGMIN=4,
     :                     NUMREMESHTERMS=4)
      CHARACTER REMESHTERMS(1:NUMREMESHTERMS)*12
      DATA REMESHTERMS     /'IGRID       ',
     :                      'SHAPEMIN    ',
     :                      'SEGMIN      ',
     :                      'ANGMIN      '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
        ITERM = MATCH(REMESHTERMS,NUMREMESHTERMS,WORD)
        IF (ITERM.EQ.0) THEN
          GO TO 200
        ELSE IF (ITERM.EQ.INDXIGRID) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)IGRID
          IF (IGRID.EQ.3) MFIXFLAG = 1
        ELSE IF (ITERM.EQ.INDXSHAPEMIN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)SHAPEMIN
        ELSE IF (ITERM.EQ.INDXSEGMIN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)SEGMIN
        ELSE IF (ITERM.EQ.INDXANGMIN) THEN
          CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
          IF (INDX2.EQ.0) GO TO 200
          READ(INSTR(INDX1:INDX2-1),*,ERR=200)ANGMIN
        END IF
        IF (IERR.EQ.0) GO TO 5
 200    INDX2=INDX1
 210    IERR = 10
      END IF
      RETURN
      END

      SUBROUTINE GETRHEOMODDATA(INSTR,LEN,INDX2,RHEO,IRHEOTYPE,
     :IVOLD,IERR)
      INCLUDE "limits.parameters"
      CHARACTER WORD*12
      CHARACTER INSTR*(LINELEN*MAXLINES)
      DIMENSION RHEO(5)
      PARAMETER(           INDXRHEO=1,
     :                     INDXTYPE=2,
     :                     NUMRHEOTERMS=2)
      CHARACTER RHEOTERMS(1:NUMRHEOTERMS)*12
      DATA RHEOTERMS    /'RHEO        ',
     :                   'TYPE        '/
   5  INDX1=INDX2
      CALL GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAXWORDLEN)
      IF (INDX2.NE.0) THEN
         ITERM = MATCH(RHEOTERMS,NUMRHEOTERMS,WORD)
         IF (ITERM.EQ.0) THEN
            GO TO 210
         ELSE IF (ITERM.EQ.INDXTYPE) THEN
            CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
            IF (IERR.NE.0) GO TO 220
            READ(INSTR(INDX1:INDX2-1),*,ERR=200)IRHEOTYPE
       ELSE IF (ITERM.EQ.INDXRHEO) THEN
            CALL GETARRAY(INSTR,LEN,INDX2,RHEO,5,ICNT,IERR)
            IF (IERR.NE.0) GO TO 220
         END IF
         IF (IERR.EQ.0) GO TO 5
 200     INDX2=INDX1
 210     IERR = 10
      END IF
C     Depending on IRHEOTYP set IVOLD to 1 if VOLD 
C     is to be used and saved with solution.
      IF (IRHEOTYPE.EQ.1) THEN   
         IVOLD = 1
      ELSE IF (IRHEOTYPE.EQ.2) THEN   
         IVOLD = 1
      ENDIF
 220  RETURN
      END
 
      SUBROUTINE WRITEINSTR(INSTR,LENGTH,IERRINDX)
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER PNTR*80
      COMMON/AI/LUW,LSC,LBC,LLG

      IF (LENGTH.LT.1) GO TO 100
      ISTART=1
      IFINISH = 1
      IF (IERRINDX.NE.0) THEN
        DO 10 J=1,80
          PNTR(J:J)=' '
  10    CONTINUE
      END IF
      DO 20 I=LENGTH,ISTART,-1
        IF (INSTR(I:I).NE.' ') THEN
          IFINISH=I
          GO TO 80
      END IF
  20  CONTINUE

  80  ILAST=1
  70  IF (INSTR(ILAST:ILAST).EQ.'&'.OR.ILAST.EQ.IFINISH) THEN
        WRITE(LSC,*)INSTR(ISTART:ILAST)
        IF (LUW.NE.LSC) WRITE(LUW,*)INSTR(ISTART:ILAST)
        IF (IERRINDX.GT.ISTART.AND.IERRINDX.LE.ILAST) THEN
          IPTR=IERRINDX-ISTART+1
          PNTR(IPTR:IPTR)='^'
          WRITE(LSC,*)PNTR(1:IPTR)
          IF (LUW.NE.LSC) WRITE(LUW,*)PNTR(1:IPTR)
        END IF
        ISTART = ILAST+1
      END IF
      ILAST = ILAST+1
      IF (ILAST.LE.IFINISH) GO TO 70

      IF (IERRINDX.NE.0) THEN
        WRITE(LSC,*)'The last instruction could not be parsed'
        IF (LUW.NE.LSC) WRITE(LUW,*)
     :              'The last instruction could not be parsed'
      END IF
 100  RETURN
      END

      SUBROUTINE READLINE(INSTR,LENGTH,LINP,LSC,IEND)
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*MAXLINES)
C
C    read 80 character blocks while the line ends with '&'
C    ignores trailing blanks after the '&' but includes
C    the '&' (SKIP treats it as a blank)
C    up to 800 characters can be read
C
      IEND=0
      DO 10 N=1,LINELEN*MAXLINES
        INSTR(N:N) = ' '
  10  CONTINUE
      LENGTH = 0
  460 CONTINUE
      K1=1
  461 K2=K1+LINELEN-1
      IF(K2.GT.MAXLINES*LINELEN)THEN
        WRITE(LSC,*)'Available no. of continuation lines',
     :  ' (MAXLINES) exceeded'
        STOP
      END IF
      READ(LINP,'(A80)',END=700)INSTR(K1:K2)
C
C   if comment character detected, line is ignored even if 
C   continuation character occurs at end of line
C
      IF(INSTR(K1:K1).EQ.'#')RETURN
      KCONT=0
      DO 450 K = K1,K2
        IF(INSTR(K:K).EQ.'&')KCONT=K
  450 CONTINUE
C
C    this line to be continued, get next line as well
C
      IF(KCONT.NE.0)THEN
        K1=KCONT+1
        LENGTH=KCONT
        GO TO 461
      END IF
      LENGTH=LENGTH+LINELEN
      RETURN
  700 IEND = 1
      RETURN
      END

      INTEGER FUNCTION MATCH(TERMS,NUM,WORD)
      INCLUDE "limits.parameters"
      CHARACTER TERMS(1:NUM)*12
      CHARACTER WORD*12
C
C    search a set of strings
C    return the index if a match is found
C    else return 0
C

      MATCH = 0
      DO 20 N=1,NUM
        I = INDEX(TERMS(N),WORD)
        IF (I.GT.0) THEN
          MATCH = N
          GO TO 30
        END IF
  20  CONTINUE
  30  CONTINUE
      END

      SUBROUTINE SKIP(INSTR,LENGTH,K1,NNB,J1)
C
C    analyses a string to return index (J1) of next
C    blank character (NNB='b') or next non-blank
C    character (NNB='n') or equal sign (NNB='e') after K1.  
C    J1=0 is returned if suitable character is not found
C
      CHARACTER INSTR*1
      CHARACTER NNB*1
      DIMENSION INSTR(LENGTH)
C
      J1=0
      IF(K1.GT.LENGTH)RETURN
      IF(NNB.EQ.'n')THEN
        DO 10 J=K1,LENGTH
          IF((INSTR(J).NE.' ').AND.(INSTR(J).NE.',')
     :         .AND.(INSTR(J).NE.CHAR(9)).AND.(INSTR(J).NE.CHAR(13))
     :         .AND.(INSTR(J).NE.'&')) THEN
            J1=J
            GO TO 11
          END IF
   10   CONTINUE
      ELSE IF(NNB.EQ.'b')THEN
        DO 20 J=K1,LENGTH
          IF((INSTR(J).EQ.' ').OR.(INSTR(J).EQ.',').OR.
     :        (INSTR(J).EQ.CHAR(9)).OR.(INSTR(J).EQ.CHAR(13)).OR.
     :            (INSTR(J).EQ.'&').OR.(INSTR(J).EQ.'='))THEN
            J1=J
            GO TO 11
          END IF
   20   CONTINUE
      ELSE IF(NNB.EQ.'e')THEN
        DO 30 J=K1,LENGTH
          IF(INSTR(J).EQ.'=') THEN
            J1=J
            GO TO 11
          END IF
   30   CONTINUE
      END IF
   11 RETURN
      END

      SUBROUTINE GETWORD(INSTR,LEN,INDX2,WORD,LNW,MAX)
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*MAXLINES)
      CHARACTER WORD*(*)
C
C    record the index of the next non-blank character
C    find the next blank character
C    copy the non-blank characters into WORD
C
      DO 10 N=1,MAX
        WORD(N:N) = ' '
  10  CONTINUE
      LNW=0
      CALL SKIP(INSTR,LEN,INDX2,'n',INDX1)
      IF (INDX1.GT.0.AND.INDX1.LT.LEN) THEN
        CALL SKIP(INSTR,LEN,INDX1,'b',INDX2)
        IF (INDX2.NE.0.AND.INDX2-INDX1.LE.MAX) THEN
          LNW=INDX2-INDX1       ! the length of the word
          WORD(1:LNW) = INSTR(INDX1:INDX2-1)
C         WRITE(*,*)'The word is: ',WORD(1:LNW)
        END IF
      ELSE
        INDX2 = 0
      END IF
      RETURN
      END

      SUBROUTINE GETVALUE(INSTR,LEN,INDX2,INDX1)
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*MAXLINES)
C
C    read '=' and find a number or word
C    set INDX1 to the INDEX of the first non-blank after the '='
C    set INDX2 to the first blank after the value
C
      IERR = 0
      CALL SKIP(INSTR,LEN,INDX2,'e',INDX1)
      IF(INDX1.EQ.0) THEN
        INDX1 = INDX2
        INDX2 = 0
        IERR = 10
      ELSE
        INDX2=INDX1+1
        CALL SKIP(INSTR,LEN,INDX2,'n',INDX1)
        IF(INDX1.EQ.0) THEN
          INDX1 = INDX2
          INDX2 = 0
          IERR = 10
        ELSE
          CALL SKIP(INSTR,LEN,INDX1,'b',INDX2)
          IF(INDX2.EQ.0) THEN
            IERR = 10
          END IF
        END IF
      END IF
      RETURN
      END

      SUBROUTINE GETARRAY(INSTR,LEN,INDX2,RARRAY,LIMIT,ICNT,IERR)
      COMMON/AI/LUW,LSC,LBC,LLG
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*MAXLINES)
      DIMENSION RARRAY(LIMIT)

      ICNT=1
      INDX1 = INDX2
      CALL GETVALUE(INSTR,LEN,INDX2,INDX1)
      IF (INDX2.NE.0) THEN
        IF (INSTR(INDX1:INDX1).EQ.';') GO TO 700
        READ(INSTR(INDX1:INDX2-1),*,ERR=200)RARRAY(ICNT)
      END IF
  5   INDX1 = INDX2
      CALL SKIP(INSTR,LEN,INDX2,'n',INDX1)
      IF(INDX1.EQ.0) THEN
        INDX1=INDX2
        GO TO 200
      END IF
      CALL SKIP(INSTR,LEN,INDX1,'b',INDX2)
      IF(INDX2.EQ.0) GO TO 200
      IF (INSTR(INDX1:INDX1).EQ.';') GO TO 700
      ICNT=ICNT+1
      IF (ICNT.GT.LIMIT) THEN
        IERR=12
        GO TO 300
      ELSE
        READ(INSTR(INDX1:INDX2-1),*,ERR=200)RARRAY(ICNT)
      END IF
      IF (IERR.EQ.0) GO TO 5
      GO TO 700
 200  IERR=10
      INDX2 = INDX1
 300  WRITE(LUW,*)'An error occurred while reading an array -'
      WRITE(LSC,*)'An error occurred while reading an array -'
      IF (IERR.EQ.10) THEN
        WRITE(LUW,*) 'check punctuation'
        WRITE(LSC,*) 'check punctuation'
      ELSE IF (IERR.EQ.12) THEN
        WRITE(LUW,10200) LIMIT
        WRITE(LSC,10200) LIMIT
10200 FORMAT('number of entries exceeds array size (',I2,')')
      END IF
 700  RETURN
      END
      SUBROUTINE FINDMAXNBS(LEM,NE,NN,MAXNBS,IERR)
      DIMENSION LEM(6,NE)
      INTEGER NUMNBS(NN)

      DO I=1,NN
        NUMNBS(I) = 0
      END DO
      DO I=1,NE
        DO J=1,3
          N=LEM(J,I)
          NUMNBS(N)=NUMNBS(N)+1
        END DO
      END DO
      MAXNBS=0
      DO I=1,NN
        IF (NUMNBS(I).GT.MAXNBS) MAXNBS=NUMNBS(I)
      END DO
      RETURN
      END
