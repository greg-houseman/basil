C*--------------------------------------------------------------------
C*    Basil / Sybil:   basil.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
      SUBROUTINE NEXTFILENAME(NAMEBC,NMLEN,INFILE,IEND,IVALID)
C
C       Get the next input filename, from line number INFILE.
C       Set the IEND flag if the end of basil.in is reached.
C       Set the IVALID flag if the listed input file exists
C
      CHARACTER NAMEBC*32
      LOGICAL XIST
      LINP=5
      NMLEN=0
      IVALID=1
      OPEN(LINP,FILE='basil.in',STATUS='old',IOSTAT=IOS)
      DO I=1,INFILE
        READ(LINP,10098,END=3)NAMEBC
      ENDDO
      CLOSE(LINP)
C
C    check validity of name
C
      CALL SKIP(NAMEBC,32,1,'n',I)
      IF (I.EQ.0) THEN
        IVALID=0
        RETURN
      ELSE
        CALL SKIP(NAMEBC,32,I,'b',J)
        NMLEN=J-I
      END IF
      INQUIRE(FILE=NAMEBC(1:NMLEN), EXIST=XIST)
      IF((NMLEN.LE.0).OR.(.NOT.XIST))THEN
        IVALID=0
        WRITE(*,100)NAMEBC
  100   FORMAT('Filename not found: ',A32)
      ENDIF
      RETURN
C
C     if end of file encountered
C
   3  IEND=1
      IVALID=0
      CLOSE(LINP)
      RETURN
      WRITE(*,50003)
50003 FORMAT('Processing completed for all input files')
10098 FORMAT(A32)
      END

      SUBROUTINE INPUTEXISTS(FNAME,NMLEN,IVALID)
C
C      Sets the IVALID flag to 0 if the file is not found
C
      CHARACTER FNAME*32
      LOGICAL XIST
      IVALID=1
      INQUIRE(FILE=FNAME(1:NMLEN), EXIST=XIST)
      IF (.NOT.XIST) THEN
C       WRITE(*,100)FNAME
  100   FORMAT('File not found: ',A32)
        IVALID=0
      END IF
      RETURN
      END

      SUBROUTINE ARRAYDIM(IV,NUPVAL,NBPVAL,NEVAL,NNVAL,NROWSVAL,
     :                    NSEGVAL,NFPVAL,NFPF3VAL,NYVAL,NULVAL,
     :                    NBLVAL,NELVAL,NSMVAL,NPMVAL)
C
C      sets the passed in scalars to the values of the array of
C      integers.
C      Allows the index values to be listed in the one file -
C      indices.parameters
C
      INCLUDE "indices.parameters"
      DIMENSION IV(64)

      NUPVAL = IV(INUP)
      NBPVAL = IV(INBP)
      NEVAL = IV(INE)
      NNVAL = IV(INN)
      NMPVAL = IV(INMP)
      NROWSVAL = IV(INROWS)
      NSEGVAL = IV(INSEG)
      NFPVAL = IV(INFP)
      NFPF3VAL=IV(INFPF3)
      NYVAL = IV(INY)
      NULVAL = IV(INUL)
      NBLVAL = IV(INBL)
      NELVAL = IV(INEL)
      NSMVAL = IV(INSM)
      NPMVAL = IV(INPM)
C     WRITE(6,*)'ARRAYDIM: NUPVAL,NBPVAL,NEVAL,NNVAL,NROWSVAL',
C    :           NUPVAL,NBPVAL,NEVAL,NNVAL,NROWSVAL
      RETURN
      END

      SUBROUTINE FLAGS(IV,ILAG,IFLT,IMSH,IGRAV,ICR,ITEMP,IDEN,
     :                 IVIS,IVOLD,IMREG,ITHDI,NCOMP,NFPF3,KPBC)
      INCLUDE "indices.parameters"
      DIMENSION IV(64)

      ILAG = IV(IILAG)
      IFLT = IV(IIFLT)
      IMSH = IV(IIMSH)
      IGRAV = IV(IIGRAV)
      ICR = IV(IICR)
      ITEMP = IV(IITEMP)
      IDEN = IV(IIDEN)
      IVIS = IV(IIVIS)
      IVOLD = IV(IIVOLD)
      IMREG = IV(IIMREG)
      ITHDI = IV(IITHDI)
      NCOMP = IV(INCOMP)
      NFPF3 = IV(INFPF3)
      KPBC = IV(IKPBC)
      RETURN
      END

      SUBROUTINE STKARRAYDIM(NKVAL,NK2VAL,NK4VAL)
      COMMON/STKDIM/KBW,K2BW,NK,NK2,K4BW,NK4
      SAVE /STKDIM/
      NKVAL = NK
      NK2VAL = NK2
      NK4VAL = NK4
      RETURN
      END

      SUBROUTINE RUNBASIL(IERROR,INTV,RLV,NAMEBC,NMLEN,NAMEW,LNWF,
     :               NAMER,LNR,BINDIR,LNBD,OUTDIR,LNOD,NAMEDB,LNDB,
     :               IELFIX,IENDP1,IENDP2,LGEM,LGIBC,LGIBCF,
     :               IELLE,IPOLYN,TEMPT,NUP,IREAD,
     :               EXO,EYO,EXLG,EYLG,UXLG,UYLG,NUL,
     :               STELPX,STELPY,ISTELP,
     :               UVELPX,UVELPY,NLAG,IVERB)
C
C    the following parameters, previously passed here, now enter 
C    via the C_DATA_TO_F routine as used when remeshing is invoked
C
C               SSQ,FROT,DENS,VHB,THDISS,THDINT,IMAT,VOLD,
C    :          EXTS,EYTS,SSQTS,FROTTS,UVPTS,SSQTS,FROTTS )
C
C  *****  BASIL  *****
C
C  This program sets up the finite element problem for deformation
C  about a fault embedded in a non-linear viscous medium.
C
C  developed by Terence Barr and Greg Houseman
C  based on finite element code originally developed by Greg Houseman
C    for deformation of a thin viscous sheet
C
C  Main program sets up common blocks and matrix sizes and
C  initialises matrices to zero.  The physical constants
C  are input and the main steps in the solution are given.
C
C  Basic flow of the program:
C
C  Initial conditions:
C    1) set up geometry and boundary conditions (on fault and external)
C  For each time step (2-5):
C    2) set up matrices
C    3) solve the matrices (using conjugate gradient method)
C    4) adjust viscosity if rheology is non-linear
C         and go to step 2, iterating until convergence
C    5) calculate new position for mesh and lagrangian markers
C         and crustal thickness if necessary
C
C  The convention used to set up the various geometries and
C  boundary conditions for the various problems is:
C
C    1) create a mesh with a fault always on a constant X-coordinate
C    2) set up the arrays that define the fault
C    3) if desired, rotate the mesh 90 degrees and deform the mesh
C         about the fault tip (must remain a rectangle)
C    4) apply external boundary conditions
C    5) do any extreme deformations to the mesh - such as tilt it for the
C         tiltblock case, etc - i.e. any deformations which cause the mesh
C         to no longer be a rectangle oriented in the x and y directions.
C
C    The array sizes are set by the parameters in the following file
C
      INTERFACE
        SUBROUTINE REMESH(EX,EY,IBCTYP,IBC,IBNGH,SEGLIM,ANGLIM,
     :                  PTLOC,ISEG,LEM,NOR,INTV,IVNEW,IFLT,NFP,
     :                  NFPF3,IFBC1,IFEQV,PDIST,MSNODE,MEASUR,
     :                  LUW,LSC,IERR)
        INCLUDE "indices.parameters"
        INCLUDE "limits.parameters"
        INTEGER INTV(64),IVNEW(64)
        DIMENSION MEASUR(MAXMEAS), MSNODE(MAXMEAS), PTLOC(3,MAXMEAS)
        REAL, DIMENSION(:), POINTER :: EX, EY
        INTEGER, DIMENSION(:,:), POINTER :: ISEG, LEM
        INTEGER, DIMENSION(:), POINTER :: NOR,IBC,IFBC1,IFEQV
        INTEGER, DIMENSION(:), POINTER :: IBCTYP,IBNGH
C       INTEGER, ALLOCATABLE :: IPROTECTED(:)
        LOGICAL, ALLOCATABLE :: FOUND(:)
        END SUBROUTINE REMESH
        SUBROUTINE C_DATA_TO_F(EX,EY,EXTS,EYTS,SSQ,FROT,SSQTS,FROTTS,
     :                   QBND,QLOAD,UVP,UVPTS,DENS,VHB,VOLD,
     :                   THDINT,THDISS,TAPERF,STK,STK2,STK3,STK4,
     :                   ISEG,LEM,IMAT,IBC,
     :                   IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :                   IBPOLE,IBCTYP,IBNGH,NOR,IK,IK2,IK4,INTV)
        REAL, DIMENSION(:), POINTER :: EX, EY, EXTS, EYTS
        REAL, DIMENSION(:), POINTER :: SSQ, FROT, SSQTS, FROTTS
        REAL, DIMENSION(:), POINTER :: QBND, QLOAD, UVP, UVPTS
        REAL, DIMENSION(:,:), POINTER :: DENS, VHB, VOLD, THDINT, THDISS
        REAL, DIMENSION(:),POINTER :: TAPERF
        DOUBLE PRECISION, DIMENSION(:), POINTER :: STK, STK2, STK3, STK4
        INTEGER, DIMENSION(:,:), POINTER :: ISEG, LEM
        INTEGER, DIMENSION(:), POINTER :: IBC,NOR,IMAT,IBCTYP,IBNGH
        INTEGER, DIMENSION(:), POINTER :: IFBC1, IFBC2, IFEQV
        INTEGER, DIMENSION(:), POINTER :: JFBC1, JFBC2
        INTEGER, DIMENSION(:), POINTER :: IK, IK2, IK4
        INTEGER, DIMENSION(:,:),POINTER :: IBPOLE
        INTEGER INTV(64)
        END SUBROUTINE C_DATA_TO_F
        SUBROUTINE GETREALPTR(CXREF,CXARRAY,NSIZE)
        INTEGER NSIZE
        REAL, DIMENSION(NSIZE), TARGET :: CXARRAY
        REAL, DIMENSION(:), POINTER :: CXREF
        END SUBROUTINE

        SUBROUTINE GETDBLPTR(CXREF,CXARRAY,NSIZE)
        INTEGER NSIZE
        DOUBLE PRECISION, DIMENSION(NSIZE), TARGET :: CXARRAY
        DOUBLE PRECISION, DIMENSION(:), POINTER :: CXREF
        END SUBROUTINE

        SUBROUTINE GETINTPTR(CIREF,CIARRAY,NSIZE)
        INTEGER NSIZE
        INTEGER, DIMENSION(NSIZE), TARGET :: CIARRAY
        INTEGER, DIMENSION(:), POINTER :: CIREF
        END SUBROUTINE
      END INTERFACE
C
      INCLUDE 'input.parameters'
      INCLUDE 'indices.parameters'
      INCLUDE 'limits.parameters'
C     DOUBLE PRECISION THDISS,THDINT,DTDISS
      DOUBLE PRECISION PNI,PLI
C
C   all common blocks used in program are given here
C
      COMMON/V01S/NDUM2,NDATE,NDUM3,COMMEN
      COMMON/FILES/NDUM1,LGFILE,PFILE
      COMMON/AI/LUW,LSC,LBC,LLG
      COMMON/BT/MSINDX,MEASUR(MAXMEAS),MSNODE(MAXMEAS),PTLOC(3,MAXMEAS)
      COMMON/STKDIM/KBW,K2BW,NK,NK2,K4BW,NK4
      SAVE /STKDIM/
      COMMON/LGDIM/NXL,NYL,NNL,NML,ILGSAVE
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      DIMENSION INTV(64)
      DIMENSION RLV(64)
      DIMENSION ISEGC(3,INTV(INSEG))
      REAL, DIMENSION(:),POINTER :: EX, EXTS, EY, EYTS
      REAL, DIMENSION(:),POINTER :: SSQ, FROT, SSQTS, FROTTS
      REAL, DIMENSION(:),POINTER :: QBND, QLOAD, UVP, UVPTS
      REAL, DIMENSION(:,:),POINTER :: DENS, VHB, VOLD, THDINT, THDISS
      REAL, DIMENSION(:),POINTER :: TAPERF
      DOUBLE PRECISION, DIMENSION(:),POINTER :: STK, STK2, STK3, STK4
C     INTEGER, DIMENSION(:,:),POINTER :: ISEG, LEM, IBNGH
      INTEGER, DIMENSION(:,:),POINTER :: ISEG, LEM
      INTEGER, DIMENSION(:),POINTER :: IBC, NOR, IMAT, IBCTYP, IBNGH
      INTEGER, DIMENSION(:),POINTER :: IFBC1, IFBC2, IFEQV, JFBC1, JFBC2
      INTEGER, DIMENSION(:),POINTER :: IK, IK2, IK4
      INTEGER, DIMENSION(:,:),POINTER :: IBPOLE
C     INTEGER, ALLOCATABLE :: IPROTECTED(:)
      DIMENSION IVNEW(64)
      INTEGER*4 IDAY(3)
      CHARACTER NDATE*16,NDUM1*16,NDUM2*16,NDUM3*16
      CHARACTER NAMEDAT*20
C     CHARACTER SERDAT*20
      CHARACTER COMMEN*80
      CHARACTER BINDIR*32,OUTDIR*32,NAMEDB*32,NAMEBC*32
      CHARACTER NAMER*32,NAMEW*32
      CHARACTER FILED*60,FILEW*60,FILEDAT*60
C     CHARACTER FILESER*60
      CHARACTER WFORM*10
      CHARACTER LGFILE*80,PFILE*80
      LOGICAL INITVEL
      LOGICAL XIST
      LOGICAL REMESHFLAG
      DIMENSION EXO(NUP),EYO(NUP),TEMPT(NUP)
      DIMENSION IELFIX(NUP)
      DIMENSION EXLG(NUL),EYLG(NUL),UXLG(NUL),UYLG(NUL)
      DIMENSION UVELPX(NLAG),UVELPY(NLAG)
      DIMENSION PNI(42),PLI(21)
      DIMENSION IFLTTIPS(MAXTIPS*2)
      DIMENSION ITIMESTEPORDER(NUMTIME)
      DIMENSION PVALUE(MAXMEAS),SMIN(3),SMAX(3)
      DIMENSION XMIN(2),XMAX(2),YMIN(2),YMAX(2)
      DIMENSION XYEXIT(4)
      DIMENSION KORNER(10)
      DIMENSION POLEP(3,MAXPOLE)
C
      LSC=6
      LOUT=4
      LDUM=3
C     LSERDAT=8
      LUW=9
      LBC=10
      LDAT=11
      LLG=12
C
C     XYEXIT is an array of 4 variables which may be used to allow 
C     the calculation to stop at a defined length condition instead
C     of a time condition
C
      XYEXIT(IXMIN) = -99999.9
      XYEXIT(IXMAX) = 99999.9
      XYEXIT(IYMIN) = -99999.9
      XYEXIT(IYMAX) = 99999.9
      SHAPELIM = 1.01
      SEGLIM = 1.0E-3
      ANGLIM=5.0
      REMESHFLAG=.FALSE.
C
C    The following will contain filenames depending on input that
C    is read by the routine INITIALREAD -> GETOUTPUTDATA
C
      DO J=1,6
        JB=(J-1)*10+1
        FILED(JB:JB+9)='          '
        FILEW(JB:JB+9)='          '
        FILEDAT(JB:JB+9)='          '
C       FILESER(JB:JB+9)='          '
      ENDDO
C
      CALL C_DATA_TO_F(EX,EY,EXTS,EYTS,SSQ,FROT,SSQTS,FROTTS,
     :                 QBND,QLOAD,UVP,UVPTS,DENS,VHB,VOLD,
     :                 THDINT,THDISS,TAPERF,STK,STK2,STK3,STK4,
     :                 ISEG,LEM,IMAT,IBC,IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :                 IBPOLE,IBCTYP,IBNGH,NOR,IK,IK2,IK4,INTV)
C
C    initialise variables
C
      CALL INITVARS(INTV,RLV,IOFF,ITIMETYPE,IDENMOD,IRHEOMOD,
     :              MFIXFLAG,NITER,IGRID,NDFORM,RHO,TILTDEG,TTSV,
     :              CENTLNG,INITVEL,NEGTIM)
      CALL INITARRAYS(LEM,NOR,DENS,TEMPT,VHB,THDISS,THDINT,IMAT,VOLD,
     :                QBND,POLEP,IBPOLE,TAPERF,UVP,QLOAD,INTV,
     :                IBC,IBNGH,IBCTYP,SSQ,FROT,IELFIX,IFLTTIPS,
     :                INTV(INUP),INTV(INE),
     :                INTV(INROWS),INTV(INBP),
     :                RHO,RLV(ISE),RLV(IHLENSC))
      IF (INTV(IILAG).EQ.1.OR.INTV(IILAG).EQ.3)
     :  CALL  INITLGARRAYS(EXLG,EYLG,UXLG,UYLG,INTV(INUL))
      IF (INTV(IILAG).EQ.1.OR.INTV(IILAG).EQ.2)
     :  CALL  INITMARKERARRAYS(UVELPX,UVELPY,INTV(INSM),INTV(INPM))
      IERROR=0
C
C    establish Fault arrays, depending on action of trimesh
C
      NFPF3=INTV(INFPF3)
      NFP=INTV(INFP)
      IFLT=INTV(IIFLT)
      IF((IFLT.GT.0).AND.(INTV(IIMSH).EQ.3).AND.(IREAD.EQ.0))THEN
        CALL FLTINIT(EX,EY,LEM,NOR,IBC,IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :               INTV(INUP),INTV(INE),NFP,NFPF3,
     :               INTV(INBP),INTV(INSEG),LUW,LSC,IERROR)
        IF (IERROR.NE.0)THEN
          WRITE(LSC,*)'FLTINIT has returned IERROR =',IERROR
          WRITE(LUW,*)'FLTINIT has returned IERROR =',IERROR
          GO TO 555
        ENDIF
      ENDIF
C
C    time and date now recorded when *.out file is set in INITIALREAD
C    but can use this and next call to DATIME to get elapsed time
C
      CALL DATIME(IDAY,ISEC,0)
C
C    Open and read the input file for this job
C
      OPEN(LBC,FILE=NAMEBC(1:NMLEN),STATUS='old',IOSTAT=IOS)
      CALL PROCESSINPUT(NAMEBC,NMLEN,BINDIR,LNBD,OUTDIR,LNOD,
     :                  NAMEDB,LNDB,NAMEW,LNWF,NAMER,LNR,INTV,RLV,
     :                  IFLTTIPS,IOFF,ITIMETYPE,IDENMOD,
     :                  IRHEOMOD,MFIXFLAG,NITER,IGRID,NDFORM,
     :                  RHO,TILTDEG,CENTLNG,ITIMESTEPORDER,
     :                  EX,EY,UVP,LEM,NOR,IBC,IBNGH,IBCTYP,KORNER,
     :                  QBND,POLEP,IBPOLE,TAPERF,VHB,IMAT,VOLD,
     :                  SSQ,FROT,DENS,TEMPT,IELFIX,
     :                  XYEXIT,SHAPELIM,SEGLIM,ANGLIM,
     :                  IFBC1,IFBC2,IFEQV,ISEG,
     :                  JFBC1,JFBC2,IENDP1,IENDP2,NFPF3,
     :                  LGEM,EXLG,EYLG,STELPX,STELPY,
     :                  UXLG,UYLG,LGIBC,LGIBCF,
     :                  INTV(INE),INTV(INUP),INTV(INBP),
     :                  INTV(INMP),INTV(INSEG),
     :                  GSDUDX,GSDUDY,GSDVDX,GSDVDY,
     :                  INITVEL,REMESHFLAG,NEGTIM,IERROR)
      CLOSE(LBC)
      IF (IERROR.NE.0)THEN
        WRITE(LSC,*)'PROCESSINPUT has returned IERROR =',IERROR
        WRITE(LUW,*)'PROCESSINPUT has returned IERROR =',IERROR
        GO TO 555
      ENDIF
C     WRITE(*,*)'IBNGH in two parts, NBP =',INTV(INBP)
C     CALL IMATPP(IBNGH,INTV(INBP))
C     CALL IMATPP(IBNGH(INTV(INBP)+1),INTV(INBP))
      WRITE(LSC,10218)NAMEBC
      WRITE(LUW,10218)NAMEBC
10218 FORMAT('Processing completed successfully for input file: ',A32)
C      WRITE(*,*)'IBCTYP after PROCESSINPUT'
C      CALL IMATPP(IBCTYP,INTV(INBP)*2)
C      WRITE(*,*)'QBND after PROCESSINPUT'
C      CALL SMATPP(QBND,INTV(INBP)*2)
C      WRITE(*,*)'UVP after PROCESSINPUT'
C      CALL SMATPP(UVP,INTV(INUP)*2)
C     CALL MITPRT(IBCTYP,INTV(INBP),2,2*INTV(INBP),6)
C
C    Check any inconsistencies in input values
C
      IF (IFLT.EQ.0.AND.(IGRID.EQ.1.OR.IGRID.EQ.2)) THEN
        IGRID=0
        WRITE(LUW,10113)
        WRITE(LSC,10113)
      ENDIF
      IF(NEGTIM*NEGTIM.NE.1)NEGTIM=1
C
C    set up values which depend on input data
C
      DFLTSSQ=-ALOG(RLV(IHLENSC))
      SEGLIM=RLV(IELSEG)
C
C    if SE set from poly input, check RV(ISE)
C
C     WRITE(*,*)'VHB array after initialisation or reading'
C     CALL SMATPP(VHB,8*INTV(INE))
      CALL SETIVV(INTV(IIVV),INTV(IIVIS),RLV(ISE),VHB,INTV(INE))
C
C    print the input values
C
      WRITE(LUW,10115)
      WRITE(LSC,10115)
      CALL PRINTINP(INTV,RLV,IOFF,ITIMETYPE,IDENMOD,IRHEOMOD,
     :                  MFIXFLAG,NITER,IGRID,NDFORM,NEGTIM,NFPF3,
     :                  RHO,TILTDEG,TTSV,SHAPELIM,SEGLIM,ANGLIM)
C
C   ANGLIM assumes radians, input in degrees.  SEGLIM input in
C   non-d length units, square is used, in radians if NCOMP=-1
C
      ANGLIM=ANGLIM*DTOR
      IF(INTV(INCOMP).EQ.-1)SEGLIM=SEGLIM*DTOR
      SEGLIM=SEGLIM*SEGLIM
C
C   The setting of IVV has been moved out of the ist.eq.0 loop
C   as IVV is not read from a previous solution
C
      IF(INTV(IIVV).EQ.0)WRITE(LUW,10114)
      IF(INTV(IIVV).EQ.1)WRITE(LUW,10116)
      IF(INTV(IIVV).EQ.2)WRITE(LUW,10117)RLV(ISE)
      IF(INTV(IIVV).EQ.3)WRITE(LUW,10118)RLV(ISE)
      IF(INTV(IIVV).EQ.4)WRITE(LUW,10119)
c
c    OMDOT is the rotation rate of the indenter.  enter OMTOT
c    the total rotation (in degrees) in 0.5 dtu.
c
      OMDOT=RLV(IOMTOT)*3.14159265/(180.0*0.5)
      IF(OMDOT.NE.0.0)WRITE(LUW,10105)OMDOT
C
C    THRESH is used to determine sudden change in Argand number
C     given in terms of factor*original crustal thickness
C   [next line already done in input.f - some tidying needed ?]
      RLV(IARGANP)=RLV(IARGAN)*RLV(IBDEPSC)*RLV(IBDEPSC)*
     1                      (RLV(IHLENSC)**(2.0+1.0/RLV(ISE)))
      RLV(IBRGANP)=RLV(IBRGAN)*RLV(IBDEPSC)*RLV(IBDEPSC)*
     1                      (RLV(IHLENSC)**(2.0+1.0/RLV(ISE)))
      IF(RLV(ITHRESH).GT.0.0)THRESHLOG=ALOG(REAL(RLV(ITHRESH)))
     1                           -ALOG(REAL(RLV(IHLENSC)))

      CALL LNCORD(PNI,PLI)
C
C     create the file used for debug; delete old version
C
      LNDBF=LNDB+LNBD+1
      FILED(1:LNDBF)=BINDIR(1:LNBD)//'/'//NAMEDB(1:LNDB)
      INQUIRE(FILE=FILED(1:LNDBF),EXIST=XIST)
      IF (XIST) THEN
        OPEN(LDUM,FILE=FILED(1:LNDBF),STATUS='OLD',FORM='UNFORMATTED',
     :       IOSTAT=IOS)
        CLOSE(LDUM,STATUS='DELETE')
      ENDIF
      OPEN(LDUM,FILE=FILED(1:LNDBF),STATUS='NEW',FORM='UNFORMATTED',
     :     IOSTAT=IOS)
      CLOSE(LDUM)
C
C     create the file used for binary output; delete old version
C
      LNWFF=LNBD+LNWF+1
      FILEW(1:LNWFF)=BINDIR(1:LNBD)//'/'//NAMEW(1:LNWF)
      INQUIRE(FILE=FILEW(1:LNWFF),EXIST=XIST)
      IF (XIST) THEN
        OPEN(LOUT,FILE=FILEW(1:LNWFF),STATUS='OLD',FORM='UNFORMATTED',
     :       IOSTAT=IOS)
        IF(IOS.EQ.0)CLOSE(LOUT,STATUS='DELETE')
C       CALL SYSTEM('ls FD.sols/MODEL')
        WRITE(LSC,10135)FILEW(1:LNWFF)
        WRITE(LUW,10135)FILEW(1:LNWFF)
10135   FORMAT('File deleted: ',A30)
      ENDIF
      OPEN(LOUT,FILE=FILEW(1:LNWFF),STATUS='NEW',FORM='UNFORMATTED',
     :     IOSTAT=IOS)
      CLOSE(LOUT)
C
C      these files are used for series data
C
      IF (INTV(IIDATA).NE.0) THEN
        LNDT=LNOD+LNWF+5
        FILEDAT=OUTDIR(1:LNOD)//'/'//NAMEW(1:LNWF)//'.dat'
        OPEN(LDAT,file=FILEDAT,STATUS='old',IOSTAT=IOS)
        IF(IOS.EQ.0)CLOSE(LDAT,STATUS='delete',IOSTAT=IOS)
        OPEN(LDAT,file=FILEDAT,status='new')
        CLOSE(LDAT)
C       IF (IGRID.EQ.4) THEN
C         LNDT=LNOD+LNWF+5
C         FILESER(1:LNDT)=OUTDIR(1:LNOD)//'/'//NAMEW(1:LNWF)//'.ser'
C         OPEN(LSERDAT,file=FILESER,STATUS='old',IOSTAT=IOS)
C         IF(IOS.EQ.0)CLOSE(LSERDAT,STATUS='delete',IOSTAT=IOS)
C         OPEN(LSERDAT,file=FILESER,status='new')
C         CLOSE(LSERDAT)
C       END IF
      END IF
C
C   Save the original grid to enable regridding
C
      IF (IFLT.NE.0.AND.(IGRID.EQ.1.OR.IGRID.EQ.2)) THEN
        CALL MTCOPY(EX,EXO,INTV(INUP),0.0,1.0)
        CALL MTCOPY(EY,EYO,INTV(INUP),0.0,1.0)
      END IF
C
C       save the starting condition
C
      WRITE(LUW,10115)
      WRITE(LSC,10115)
      OPEN(LDUM,FILE=FILED,STATUS='OLD',FORM='UNFORMATTED')
      CALL WRITESTORE(LDUM,0,INTV,RLV,NAMER,NAMEW,0,IERROR,
     :          INTV(INROWS),INTV(INUP),INTV(INE),INTV(INBP),
     :          INTV(INBP2),INTV(INFP),INTV(INEL),
     :          INTV(INUL),INTV(INPM),INTV(INSM),INTV(INBL),
     :          INTV(INELLEP),INTV(INSEG),
     :          EX,EY,UVP,LEM,NOR,VHB,THDINT,
     :          IMAT,VOLD,SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :          IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC1,IFBC2,
     :          IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :          UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :          IELLE,IPOLYN)
      CLOSE(LDUM)
      WRITE(LUW,10125)FILED(1:60)
      WRITE(LSC,10125)FILED(1:60)
C
C    if INDFIX = 1 or 3 the velocity on all fixed velocity segments
C    of the boundary is zeroed.  If INDFIX = 2 or 3 the stress on all
C    fixed stress segments of the boundary is zeroed.
C
C     IF(INTV(IIBCMOD).EQ.6) THEN
      IF((INTV(IINDFIX).EQ.1).OR.(INTV(IINDFIX).EQ.3))THEN
        CALL INDEND(QBND,IBCTYP,INTV(INBP2),0)
        WRITE(LUW,10122)
        WRITE(LSC,10122)
      END IF
C     ELSE IF(INTV(IIBCMOD).EQ.7)THEN
      IF((INTV(IINDFIX).EQ.2).OR.(INTV(IINDFIX).EQ.3))THEN
        CALL INDEND(QBND,IBCTYP,INTV(INBP2),1)
        WRITE(LUW,10123)
        WRITE(LSC,10123)
      END IF
C
C    check that IDT0 is >= 1/TSAVE
      IF (RLV(ITSAVE).GT.0.AND.
     :             INTV(IIDT0).LT.NINT(1.0/RLV(ITSAVE))) THEN
        RLV(ITSAVE) = 1.0/INTV(IIDT0)
        WRITE(LUW,10128)RLV(ITSAVE),INTV(IIDT0)
        WRITE(LSC,10128)RLV(ITSAVE),INTV(IIDT0)
      ENDIF

      IF (INTV(IKEXIT).LT.0) GO TO 555

      IW=0
      IR=0
      ISET=0
      SIGMA=0
      KRMINT=10
      KREMESH=0
      IREM=0    ! changed to 1 on first remesh
C
C     -------------------------------------------------start timestep
C     commence the main loop
C
   99 CONTINUE
      WRITE(LUW,10115)
      WRITE(LSC,10115)
C     WRITE(*,*)'Commencing main solver loop at statement 99'
C
C    Use a 2 step Runge-Kutta predictor-corrector algorithm for timestepping
C    First store the current mesh and finite deformation arrays so that they
C    are preserved until the 2nd step of the R-K iteration.  First step
C    computes velocity estimate at time level n, 2nd step computes velocity
C    estimate from time level n+1, then takes average of two estimates in 
C    order to do time integration of mesh and crustal thickness. Duplicate
C    arrays are also used for remeshing, if required, so copies made here
C    prior to possible remeshing
C
        CALL MTCOPY(EX,EXTS,INTV(INUP),0.0,1.0)
        CALL MTCOPY(EY,EYTS,INTV(INUP),0.0,1.0)
        IF (ISSQACTIVE.EQ.1) CALL MTCOPY(SSQ,SSQTS,INTV(INUP),0.0,1.0)
        IF (IROTACTIVE.EQ.1) CALL MTCOPY(FROT,FROTTS,INTV(INUP),0.0,1.0)
C
C    recondition the finite element mesh if required
C
      IF (IGRID.EQ.3) THEN
        CALL MFIX(JNONUM,EX,EY,RLV(IVC),VHB,LEM,NOR,SSQ,FROT,
C    :            EX1,EY1,NOMOV,SSQ1,FROT1,IELF1,
     :            IELFIX,IBC,INTV(IIVIS),
     :            INTV(IITHDI),THDINT,THDISS,
     :            INTV(INUP),INTV(INE),INTV(INBP),PNI,PLI,LSC,IERROR)
        IF (IERROR.NE.0) GO TO 555
      ELSE IF((IGRID.EQ.4).AND.(REMESHFLAG.EQV..TRUE.).AND.
     :    (INTV(INSEG).GT.0)) THEN
        KREMESH=INTV(IKSTEP)
        WRITE(LUW,10129)RLV(ITIMEL)
        WRITE(LSC,10129)RLV(ITIMEL)
C
C       save the solution before the remeshing
C       this could probably be dispensed with, unless debugging needed
C
        OPEN(LDUM,FILE=FILED,STATUS='OLD',FORM='UNFORMATTED')
        CALL WRITESTORE(LDUM,0,INTV,RLV,NAMER,NAMEW,1,IERROR,
     :          INTV(INROWS),INTV(INUP),INTV(INE),INTV(INBP),
     :          INTV(INBP2),INTV(INFP),INTV(INEL),
     :          INTV(INUL),INTV(INPM),INTV(INSM),INTV(INBL),
     :          INTV(INELLEP),INTV(INSEG),
     :          EX,EY,UVP,LEM,NOR,VHB,THDINT,
     :          IMAT,VOLD,SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :          IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC1,IFBC2,
     :          IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :          UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :          IELLE,IPOLYN)
        CLOSE(LDUM)
C
C     remeshing grid
C
        CALL REMESH(EX,EY,IBCTYP,IBC,IBNGH,SEGLIM,ANGLIM,
     :              PTLOC,ISEG,LEM,NOR,INTV,IVNEW,IFLT,NFP,
     :              NFPF3,IFBC1,IFEQV,PDIST,MSNODE,MEASUR,
     :              LUW,LSC,IERROR)
        IF (IERROR.NE.0) THEN
          WRITE(*,*)'Remesh error: terminating current job'
          GO TO 555
        END IF
        IREM=1
C
C     this routine redirects the pointers to the new arrays
C     produced in the remeshing operation
C
        CALL C_DATA_TO_F(EX,EY,EXTS,EYTS,SSQ,FROT,SSQTS,FROTTS,
     :                   QBND,QLOAD,UVP,UVPTS,DENS,VHB,VOLD,
     :                   THDINT,THDISS,TAPERF,STK,STK2,STK3,STK4,
     :                   ISEG,LEM,IMAT,IBC,
     :                   IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :                   IBPOLE,IBCTYP,IBNGH,
     :                   NOR,IK,IK2,IK4,INTV)
C
C    re-establish Fault arrays after remeshing - if fault present
C
      NFPF3=INTV(INFPF3)
      NFP=INTV(INFP)
      IF((IFLT.GT.0).AND.(INTV(IIMSH).EQ.3).AND.(IREAD.EQ.0))THEN
        CALL FLTINIT(EX,EY,LEM,NOR,IBC,IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :               INTV(INUP),INTV(INE),NFP,NFPF3,
     :               INTV(INBP),INTV(INSEG),LUW,LSC,IERROR)
        IF (IERROR.NE.0)THEN
          WRITE(LSC,*)'FLTINIT has returned IERROR =',IERROR
          WRITE(LUW,*)'FLTINIT has returned IERROR =',IERROR
          GO TO 555
        ENDIF
C
C    reset fault arrays that act on periodic boundaries
C
        IF(NFP.GT.NFPF3)THEN
          CALL PERIBC(EX,EY,QBND,LEM,NOR,IBC,IBNGH,IBCTYP,
     :                IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :                IFLT,INTV(INUP),INTV(INE),INTV(INBP),
     :                NFP,NFPF3,INTV(IKPBC),PDIST,IREM,LUW,LSC,IERROR)
C         WRITE(*,*)'PERIBC called PDIST,NFP,NFPF3=',PDIST,NFP,NFPF3
          IF (IERROR.NE.0)THEN
            WRITE(LSC,*)'PERIBC has returned IERROR =',IERROR
            WRITE(LUW,*)'PERIBC has returned IERROR =',IERROR
            GO TO 555
          ENDIF
        ENDIF
C
C       WRITE(*,*)'IFEQV array follows'
C       CALL IMATPP(IFEQV,NFP)
C       CALL MITPRT(IFEQV,2,NFP/2,NFP,6)
      ENDIF
C
C   reset series arrays
C
        IF (INTV(IMSINDX).GT.0) THEN
          CALL GETNODENO(MEASUR,PTLOC,MSNODE,INTV(IMSINDX),
     :                     EX,EY,NOR,INTV(INUP),INTV(INN),
     :                     LUW,LSC,IERROR)
          CALL GETELEMENTNO(MEASUR,PTLOC,MSNODE,INTV(IMSINDX),
     :                       EX,EY,NOR,INTV(INUP),INTV(INN),
     :                       LEM,INTV(INE),LUW,LSC,IERROR)

#ifndef GFORTRAN
C==      the call below for SunOS Fortran and g77
          OPEN(LDAT,FILE=FILEDAT,STATUS='old',
     1                ACCESS='append',IOSTAT=IOS)
#else
C==      the call below for gfortran. g77 does not recognize POSITION
          OPEN(LDAT,FILE=FILEDAT,STATUS='old',
     1                POSITION='append',IOSTAT=IOS)
#endif
          CALL WRITESERIESDATA(LEM,NOR,INTV(INE),INTV(INUP),
     :                MSNODE,MEASUR,INTV(IMSINDX),
     :                PTLOC,EX,EY,RLV(ITIMEL),LDAT,IERROR)
          CLOSE(LDAT)
        END IF
C
C       save the solution from the remeshing
C
        OPEN(LDUM,FILE=FILED,STATUS='OLD',FORM='UNFORMATTED')
        CALL WRITESTORE(LDUM,0,INTV,RLV,NAMER,NAMEW,1,IERROR,
     :          INTV(INROWS),INTV(INUP),INTV(INE),INTV(INBP),
     :          INTV(INBP2),INTV(INFP),INTV(INEL),
     :          INTV(INUL),INTV(INPM),INTV(INSM),INTV(INBL),
     :          INTV(INELLEP),INTV(INSEG),
     :          EX,EY,UVP,LEM,NOR,VHB,THDINT,
     :          IMAT,VOLD,SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :          IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC1,IFBC2,
     :          IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :          UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :          IELLE,IPOLYN)
        CLOSE(LDUM)
        WRITE(LUW,10132)FILED
        WRITE(LSC,10132)FILED
        REMESHFLAG = .FALSE.
      END IF
C
C    2 step Runge-Kutta predictor-corrector algorithm used in timestep
C    loop. Can use single step explicit method by choosing TYPE=EXPLICIT
C    in the input file. The above arrays will still be used, but there 
C    will only be one pass per timestep.
C
      DO 199 IRK=1,ITIMETYPE
C
C    Calculate the load vector
C         IGRAV = +-1 for X direction
C         IGRAV = +-2 for Y direction
C         IGRAV = 3 for arbitrary direction
C         IGRAV = 4 for z direction
C     A.  using gradients of crustal thickness if IGRAV=4
C
        IF(INTV(IIGRAV).EQ.4.AND.ISSQACTIVE.EQ.1)THEN
          IF((RLV(IARGAN).NE.0.0).OR.(RLV(IBRGAN).NE.0.0))
     :         CALL PLOAD(ISET,EX,EY,LEM,NOR,IELFIX,PNI,PLI,
     :                 SSQ,QLOAD,RLV(IARGANP),RLV(IBRGANP),
     :                 THRESHLOG,INTV(INROWS),INTV(INUP),
     :                 INTV(INE),INTV(INCOMP),LSC)
C     B. using body force distribution if IGRAV!=4
C     RHOG(1) is the buoyancy coeffiecient F, RHOG(2) is the background
C     relative density in the solution region.  RHOG(3) is the relative
C     density of a special region, and RHOG(4) is a parameter used to
C     define the special region in subroutine DENSIT
C
        ELSE IF(INTV(IIGRAV).NE.4.AND.INTV(IIGRAV).NE.0)THEN
          IF(IDENMOD.GT.0)
     :      CALL DENSIT(DENS,EX,EY,NOR,LEM,INTV(INE),INTV(INUP),
     :                  RLV(IRHOG2),RLV(IRHOG3),RLV(IRHOG4),
     :                  LUW,LSC,IDBUG)
          IF(INTV(IIDEN).GT.0)
     :      CALL BFORCE(INTV(IIGRAV),RLV(IRHOG1),INTV(INCOMP),EX,EY,
     :                  LEM,NOR,PNI,PLI,DENS,QLOAD,
     :                  INTV(INUP),INTV(INE),INTV(INROWS))
C       WRITE(*,*)'QLOAD after BFORCE NROWS...',INTV(INROWS)
C       CALL SMATPP(QLOAD,INTV(INROWS))
        END IF
C
C     The following group of subroutines modify the boundary condition
C     arrays for specific purposes:
C     IF (MB)  set by BCONDMOD
C
C     MB1: adjust the applied normal stress boundary conditions where
C     a fixed external pressure is applied
C
C    n.b. clash of purpose between use of IBCMOD=14 by LITHBC and MENDBC
C
        IF((INTV(IIBCMOD).EQ.1).OR.(INTV(IIBCMOD).EQ.3).OR.
     :     ((INTV(IIBCMOD).GE.11).AND.(INTV(IIBCMOD).LE.14)))
     :      CALL LITHBC(RLV(IHLENSC),EX,EY,NOR,IBC,IBNGH,IBCTYP,
     :                  QBND,RLV(IARGANP),INTV(INUP),
     :                  INTV(INBP),INTV(INN),INTV(IIBCMOD),
     :                  INTV(INCOMP),RLV(IXREFM),RLV(IBCV1),LUW)
C
C     apply plate rotations to set QBND on boundary segments as required
C     this function now completed by VSBSET called from LOADBC in cginit.f
C
C       IF(INTV(IIVRESET).NE.0)
C    :    CALL VSBSET(QBND,EX,EY,NOR,IBC,IBCTYP,POLEP,
C    :                IBPOLE,TAPERF,CENTLNG,INTV(INCOMP),
C    :                INTV(INBP),INTV(INUP),LUW,LSC)
C      
C     MB3: adjust the boundary conditions as the slab goes past some x-coordinate
C     BCV(1) is 0.0 to annul this routine, +ve to change boundary to applied
C      traction, -ve to change boundary to frictional condition
C     BCV(2) is the X-coordinate at which the fixed velocity condition is
C     changed. BCV(3) is the new value of traction, or friction to be applied
C
C    n.b. clash of purpose between use of IBCMOD=14 by LITHBC and MENDBC
C
        IF((INTV(IIBCMOD).EQ.4).OR.(INTV(IIBCMOD).EQ.5).OR.
     :     (INTV(IIBCMOD).EQ.8).OR.(INTV(IIBCMOD).EQ.14).OR.
     :     (INTV(IIBCMOD).EQ.15).OR.(INTV(IIBCMOD).EQ.21))
     :      CALL MENDBC(INTV(IIBCMOD),RLV(IBCV1),RLV(ITIMEL),
     :                  RLV(IXLEN),RLV(IYLEN),RLV(IBIG),EX,EY,
     :                  NOR,IBC,IBCTYP,QBND,INTV(INUP),INTV(INBP),
     :                  LSC,LUW)

C
C      Adjust boundary velocities on xmin and xmax boundaries
C
        IF(INTV(IIBCMOD).EQ.9)
     :      CALL UPDBNDVEL(QBND,EX,EY,IBC,NOR,
     :                   INTV(INBP),INTV(INUP),RLV(IBIG))
C
C      Adjust boundary velocities based on general shear strain rate
C
        IF(INTV(IIBCMOD).EQ.10)
     :      CALL GENSHBC(INTV(INUP),INTV(INBP),RLV(IBIG),
     :                   QBND,EX,EY,IBC,IBCTYP,NOR,
     :                   RLV(IVELXO), RLV(IVELYO),RLV(ITBXOFF),
     :                   GSDUDX,GSDUDY,GSDVDX,GSDVDY)
C       IF (INITVEL.EQV..TRUE.) THEN
C           WRITE(*,*)'INITVELOCITY is called'
C           CALL INITVELOCITY(QBND,UVP,EX,EY,IBC,NOR,
C    :                        INTV(INBP),INTV(INUP),
C    :                        RLV(IBIG),INTV(IIDEFTYP))
C           INITVEL = .FALSE.
C       ENDIF
C
C     Adjust Viscosity between timesteps if required
C      type 1 or 2 for Eulerian dependence on position
C
        IF((IRHEOMOD.GT.0).AND.(IRHEOTYP.LE.2)) THEN
          CALL RHEOL(EX,EY,NOR,LEM,VHB,VOLD,IMAT,RLV(IVISP1),
     :         INTV,IRHEOTYP,INTV(INE),INTV(INUP),INTV(INN),IRK,IRHEOXC)
        END IF
C
C      solve for the velocity and pressure field
C
C      WRITE(*,*)'IBCTYP prior CGRUN'
C      CALL IMATPP(IBCTYP,INTV(INBP)*2)
C      WRITE(*,*)'QBND prior CGRUN'
C      CALL SMATPP(QBND,INTV(INBP)*2)
C      WRITE(*,*)'UVP prior CGRUN'
C      CALL SMATPP(UVP,INTV(INUP)*2)
        CALL CGRUN(INTV,RLV,UVP,EX,EY,LEM,NOR,VHB,THDISS,THDINT,
     :             SSQ,QLOAD,QBND,IBC,IBCTYP,IBNGH,IFBC1,IFBC2,IFEQV,
     :             JFBC1,JFBC2,PNI,PLI,STK,STK2,STK3,STK4,
     :             IK,IK2,IK4,VISCMIN,VISCMAX,ED2IMIN,ED2IMAX,
     :             THDIMIN,THDIMAX,POLEP,INTV(INE),INTV(INBP),INTV(INN),
     :             INTV(INUP),INTV(INFP),INTV(INROWS),IREAD,
C    :             KBW,K2BW,K4BW,NKP,NK2P,NK4P,
     :             IREM,IT,RMSU,RMSV,NITER,IOFF,NFPF3,IVERB,IERROR)
C
        IF (IERROR.NE.0) THEN
          IF ((IT.GE.INTV(IITSTOP)).OR.(IERROR.EQ.2)) THEN
C
C       save the solution from the last iteration
C
            OPEN(LDUM,FILE=FILED,STATUS='OLD',FORM='UNFORMATTED')
            CALL WRITESTORE(LDUM,0,INTV,RLV,NAMER,NAMEW,1,IERROR,
     :          INTV(INROWS),INTV(INUP),INTV(INE),INTV(INBP),
     :          INTV(INBP2),INTV(INFP),INTV(INEL),
     :          INTV(INUL),INTV(INPM),INTV(INSM),INTV(INBL),
     :          INTV(INELLEP),INTV(INSEG),
     :          EX,EY,UVP,LEM,NOR,VHB,THDINT,
     :          IMAT,VOLD,SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :          IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC1,IFBC2,
     :          IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :          UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :          IELLE,IPOLYN)
            CLOSE(LDUM)
            WRITE(LUW,10124)NAMEDB
          END IF
          IF(IERROR.EQ.2)THEN
            WRITE(LSC,10126)
            WRITE(LUW,10126)
            GO TO 100
          END IF
          GO TO 100
        END IF
C
C    Following section only for first step of R-K cycle
C     output solution parameters, check exit condition,
C     update time, report solution parameters to .out file
C
        IF(IRK.EQ.1)THEN
          WRITE(LSC,10103)
          WRITE(LUW,10103)
          WRITE(LSC,10111)INTV(IKSTEP),RLV(ITIMEL),IT,RLV(IAREA),
     :                    RLV(IVOLUM),RLV(IVISCD),RLV(IDVMX),
     :                    RMSU,RMSV,ISET
          WRITE(LUW,10111)INTV(IKSTEP),RLV(ITIMEL),IT,RLV(IAREA),
     :                    RLV(IVOLUM),RLV(IVISCD),RLV(IDVMX),
     :                    RMSU,RMSV,ISET
          CALL UVPLIM(UVVAL,UVP,EX,EY,NOR,
     :                  UMAX,VMAX,UMIN,VMIN,INTV(INROWS),INTV(INUP),
     :                  LUW,LSC,1)
          CALL XYUVP(PVALUE,EX,EY,UVP,NOR,THDINT,
     :            INTV(INUP),INTV(INN),INTV(INCOMP),INTV(INROWS),
     :            FROT,SSQ,SMIN,SMAX,UMAX,VMAX,UMIN,VMIN,
     :            XMIN,XMAX,YMIN,YMAX,
     :            ED2IMIN,ED2IMAX,VISCMIN,VISCMAX,
     :            THDIMIN,THDIMAX,WKINTMIN,WKINTMAX,
     :            RLV(ITIMEL),PTLOC,VHB,RLV,LEM,IELFIX,
     :            INTV(INE),INTV(INFP),INTV(IIVV),INTV(IICR),
     :            INTV(IMSINDX),MEASUR,MSNODE,LUW,LSC,1)
          IF (INTV(IIDATA).NE.0) THEN
            IF (INTV(IKSTEP).EQ.0) THEN
              OPEN(LDAT,FILE=FILEDAT,STATUS='old',IOSTAT=IOS)
              CALL WRITESERIESDATA(LEM,NOR,INTV(INE),INTV(INUP),
     :                MSNODE,MEASUR,INTV(IMSINDX),
     :                PTLOC,EX,EY,RLV(ITIMEL),LDAT,IERROR)
            ELSE
#ifndef GFORTRAN
C==      the call below for SunOS Fortran and g77
              OPEN(LDAT,FILE=FILEDAT,STATUS='old',
     1                ACCESS='append',IOSTAT=IOS)
#else
C==      the call below for gfortran. g77 does not recognize POSITION
              OPEN(LDAT,FILE=FILEDAT,STATUS='old',
     1                POSITION='append',IOSTAT=IOS)
#endif
            ENDIF
            WRITE(WFORM,10002)INTV(IMSINDX)
10002       FORMAT('(',I3,'G15.7)')
            WRITE(LDAT,WFORM)(PVALUE(J),J=1,INTV(IMSINDX))
            CLOSE(LDAT)
          END IF
C
C      save the solution and mesh in the dumpfile every step
C
          OPEN(LDUM,FILE=FILED,STATUS='old',FORM='unformatted')
          CALL WRITESTORE(LDUM,0,INTV,RLV,NAMER,NAMEW,1,IERROR,
     :          INTV(INROWS),INTV(INUP),INTV(INE),INTV(INBP),
     :          INTV(INBP2),INTV(INFP),INTV(INEL),
     :          INTV(INUL),INTV(INPM),INTV(INSM),INTV(INBL),
     :          INTV(INELLEP),INTV(INSEG),
     :          EX,EY,UVP,LEM,NOR,VHB,THDINT,
     :          IMAT,VOLD,SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :          IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC1,IFBC2,
     :          IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :          UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :          IELLE,IPOLYN)
          CLOSE(LDUM)
C
C      save the solution in output file if  increment in KSTEP is > ksave,
C        or if the increment in TIMEL is > tsave.
C
          IF(((INTV(IKSTEP)-INTV(IISAVE)).GE.INTV(IKSAVE)).OR.
     1         ((RLV(ITIMEL)-TTSV).GE.(RLV(ITSAVE)*0.999)))THEN
            TTSV=RLV(ITIMEL)
            INTV(IISAVE)=INTV(IKSTEP)
            IF (INTV(IISAVE).EQ.0) THEN
              OPEN(LOUT,FILE=FILEW(1:LNWFF),STATUS='old',
     :             FORM='unformatted',IOSTAT=IOS)
            ELSE
#ifndef GFORTRAN
C==      the call below for SunOS Fortran and g77
              OPEN(LOUT,FILE=FILEW(1:LNWFF),STATUS='old',
     :             FORM='unformatted',ACCESS='append',IOSTAT=IOS)
#else
C==      the call below for gfortran. g77 does not recognize POSITION
              OPEN(LOUT,FILE=FILEW(1:LNWFF),STATUS='old',
     :             FORM='unformatted',POSITION='append',IOSTAT=IOS)
#endif
            ENDIF
            CALL WRITESTORE(LOUT,1,INTV,RLV,NAMER,NAMEW,1,IERROR,
     :          INTV(INROWS),INTV(INUP),INTV(INE),INTV(INBP),
     :          INTV(INBP2),INTV(INFP),INTV(INEL),
     :          INTV(INUL),INTV(INPM),INTV(INSM),INTV(INBL),
     :          INTV(INELLEP),INTV(INSEG),
     :          EX,EY,UVP,LEM,NOR,VHB,THDINT,
     :          IMAT,VOLD,SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :          IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC1,IFBC2,
     :          IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :          UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :          IELLE,IPOLYN)
            CLOSE(LOUT)
            WRITE(LSC,10121)FILEW(1:LNWFF)
            WRITE(LUW,10121)FILEW(1:LNWFF)
            IF (ILGSAVE.NE.0) THEN
              OPEN(LLG,FILE=LGFILE,STATUS='old',IOSTAT=IOS)
              CALL WRITELGMESH( INTV(IILAG),LGEM,EXLG,EYLG,
     :                          LGIBC,LGIBCF,
     :                          STELPX,STELPY,RLV(ISTELPR),
     :                          INTV(INEL),INTV(INUL),INTV(INBL),
     :                          INTV(INSM),INTV(INPM),INTV(INRM),
     :                          LLG,LSC,LUW,IERROR)
              CLOSE(LLG)
            ENDIF
C            IF (INTV(IMSINDX).GT.0.) THEN
CC    save SERIES data mesh coordinates
C              IF (INTV(IISAVE).EQ.0) THEN
C                OPEN(LSERDAT,FILE=FILESER,STATUS='old',IOSTAT=IOS)
C              ELSE
C#ifndef GFORTRAN
CC==      the call below for SunOS Fortran and g77
C                OPEN(LSERDAT,FILE=FILESER,STATUS='old',
C     1                ACCESS='append',IOSTAT=IOS)
C#else
CC==      the call below for gfortran. g77 does not recognize POSITION
C                OPEN(LSERDAT,FILE=FILESER,STATUS='old',
C     1                POSITION='append',IOSTAT=IOS)
C#endif
C              ENDIF
C              CALL WRITESERIESDATA(LEM,NOR,INTV(INE),INTV(INUP),
C     :                MSNODE,MEASUR,INTV(IMSINDX),
C     :                PTLOC,EX,EY,RLV(ITIMEL),LSERDAT,IERROR)
C              CLOSE(LSERDAT)
C            ENDIF
            IW=IW+1
            IF(IW.GE.INTV(IIWRITE))THEN
              WRITE(LUW,13001) INTV(IIWRITE)
              WRITE(LSC,13001) INTV(IIWRITE)
              GO TO 100
            END IF
          END IF
C
C      exit if necessary
C
          IF((INTV(IKSTEP).GE.INTV(IKEXIT)).OR.
     :                 (RLV(ITIMEL).GE.0.9999*RLV(ITEXIT)).OR.
     :                   (XMIN(1).LE.XYEXIT(IXMIN)).OR.
     :                   (XMAX(1).GE.XYEXIT(IXMAX)).OR.
     :                   (YMIN(2).LE.XYEXIT(IYMIN)).OR.
     :                   (YMAX(2).GE.XYEXIT(IYMAX)))THEN
            WRITE(LSC,13002) INTV(IKSTEP),RLV(ITIMEL),
     :                       INTV(IKEXIT),RLV(ITEXIT)
            WRITE(LUW,13002) INTV(IKSTEP),RLV(ITIMEL),
     :                       INTV(IKEXIT),RLV(ITEXIT)
            GO TO 100
          END IF
          CALL DATIME(IDAY,ISEC,1)
C
C      Calculate timestep using velocity estimate at time level n
C
          CALL TSTEPS(RLV(ITIMEL),INTV(IIDT0),RLV(ITSAVE),TTSV,
     :                INTV(IMPDEF),RLV(IDVMX),SIGMA,DT)
          IF(DT.LE.1.E-4)THEN
            WRITE(LSC,10127)
            WRITE(LUW,10127)
            GO TO 100
          END IF
C
C       Integrate the thermal dissipation through the first half timestep
C
          DTDISS=DT
          IF(ITIMETYPE.EQ.2)DTDISS=0.5*DT
          IF(INTV(IITHDI).GT.0)
     :    CALL WORKED(THDISS,THDINT,INTV(INE)*7,DTDISS,
     :         RLV(IBETA),WKINTMIN,WKINTMAX,LUW,LSC)
C
C       1st step: store the current solution estimate for time level n
C
          CALL MTCOPY(UVP,UVPTS,INTV(INROWS),0.0,1.0)
C
C       2nd step: estimate the current growth rate of the velocity field for
C         use in setting the timestep (next time round)
C
        ELSE IF(IRK.EQ.2)THEN
          CALL UVPLIM(UVVALP,UVP,EX,EY,NOR,
     :                  UMAX,VMAX,UMIN,VMIN,INTV(INROWS),INTV(INUP),
     :                  LUW,LSC,0)
          CALL XYUVP(PVALUE,EX,EY,UVP,NOR,THDINT,
     :            INTV(INUP),INTV(INN),INTV(INCOMP),INTV(INROWS),
     :            FROT,SSQ,SMIN,SMAX,UMAX,VMAX,UMIN,VMIN,
     :            XMIN,XMAX,YMIN,YMAX,
     :            ED2IMIN,ED2IMAX,VISCMIN,VISCMAX,
     :            THDIMIN,THDIMAX,WKINTMIN,WKINTMAX,
     :            RLV(ITIMEL),PTLOC,VHB,RLV,LEM,IELFIX,
     :            INTV(INE),INTV(INFP),INTV(IIVV),INTV(IICR),
     :            INTV(IMSINDX),MEASUR,MSNODE,LUW,LSC,1)
          SDT=ALOG(UVVALP/UVVAL)
          SIGMA=ABS(SDT)/DT
          WRITE(LSC,*)'Current estimate of velocity growth rate =',SIGMA
C
C        calculate average of solution estimate from
C          time levels n and n+1 for 2nd order Runge-Kutta integration
C
          BRK=(1.0+SDT/3.0)/(2.0+SDT)
          ARK=1.0-BRK
          CALL MTCOPY(UVP,UVPTS,INTV(INROWS),ARK,BRK)
C
C       Integrate the thermal dissipation through the second half timestep
C        (THDISS will have been updated after second pass through CGRUN)
C
          IF(INTV(IITHDI).GT.0)
     :    CALL WORKED(THDISS,THDINT,INTV(INE)*7,DTDISS,
     :         RLV(IBETA),WKINTMIN,WKINTMAX,LUW,LSC)
C
        END IF
C
C      Calculate new finite deformation parameters related to div
C       and curl of current velocity field.
C
        DTP=DT*NEGTIM
        IF(INTV(IICR).NE.0)
     :      CALL CRUST(DTP,SSQTS,FROTTS,SSQ,FROT,UVPTS,EX,EY,
     :                 PNI,PLI,LEM,NOR,
     :                 INTV(INUP),INTV(INE),INTV(INROWS),INTV(INCOMP))
C 
C      deform the lagrangian mesh - doesn't affect any calculations 
C      this mesh is just for display purposes 
C
        IF (INTV(IILAG).NE.0.AND.IRK.EQ.ITIMETYPE)
     :    CALL LGDEF(DTP,UVPTS,UXLG,UYLG,UVELPX,UVELPY,
     :                 EXTS,EYTS,EXLG,EYLG,STELPX,STELPY,ISTELP,
     :                 LEM,NOR,
     :                 INTV(INUP),INTV(INE),INTV(INROWS),INTV(INUL),
     :                 NNL,NML,INTV(INSM),INTV(INPM),
     :                 INTV(IILAG),IFLT)
C
C      Deform the finite element mesh using the current velocity
C
        CALL DFSTEP(DTP,EXTS,EYTS,EX,EY,UVPTS,NOR,LEM,
     :              INTV(INUP),INTV(INE),INTV(INN),INTV(INCOMP),
     :              IOFF,RLV(ITBXOFF),RLV(ITBYOFF),
     :              RLV(IVELXO),RLV(IVELYO),
     :              SEGMIN,LSEGMIN,SHAPEMIN,LSHAPEMIN,
     :              ANGMIN,LANGMIN)
C
C      report the mesh shape parameters and decide if REMESHING is needed
C
        IF(IRK.EQ.1)THEN
          WRITE(LUW,10133)SEGMIN,SHAPEMIN,ANGMIN*180.0/PI
          WRITE(LSC,10133)SEGMIN,SHAPEMIN,ANGMIN*180.0/PI
C
C      if one of the mesh limits has been exceeded after the mesh interval
C      or two of the mesh limits have been exceeded then remesh
C
          IF((SEGMIN.LT.SEGLIM).OR.(SHAPEMIN.LT.SHAPELIM).OR.
     :              (ANGMIN.LT.ANGLIM)) THEN
            IF( ((SEGMIN.LT.SEGLIM).AND.(ANGMIN.LT.ANGLIM)).OR.
     :          (INTV(IKSTEP).GE.KREMESH+KRMINT)) THEN
              REMESHFLAG=.TRUE.
            END IF
          END IF
        END IF
C
C      Exit this loop after the 2nd Runge-Kutta step
C
  199 CONTINUE
C
C    Regridding even if fault is locked because of discontinuous
C    boundary conditions
C    Should the fault regridding be moved to the start of the loop?
      IF(IGRID.EQ.1) CALL FLTGRD(RLV(IHLENSC),RLV(ITBXOFF),RLV(ITBYOFF),
     :                  EX,EY,EXO,EYO,SSQ,FROT,
     :                  NOR,LEM,IFBC1,IFEQV,
C    :                  EX1,EY1,SSQ1,FROT1,NOMOV,JNOAL2,
     :                  INTV(INUP),INTV(INE),INTV(INN),INTV(INFP),
     :                  LUW,LSC)
      IF(IGRID.EQ.2) CALL REGRID(RLV(IHLENSC),EX,EY,EXO,EYO,SSQ,FROT,
     :                  NOR,LEM,IFBC1,IFEQV,IBCTYP,IENDP1,IENDP2,
C    :                  EX1,EY1,SSQ1,FROT1,NOMOV,JNOAL2,
     :                  INTV(INUP),INTV(INE),INTV(INN),INTV(INBP),
     :                  INTV(INFP))
C
C    Update the time counters
C
      INTV(IKSTEP)=INTV(IKSTEP)+1
      RLV(ITIMEL)=RLV(ITIMEL) + DT
C
C    ---------------------------------------------------end timestep
C     Return for next timestep
      GO TO 99
  100 CONTINUE
C
C     Save the final solution if it wasn't already saved above
C
      IF(INTV(IKSTEP).GT.INTV(IISAVE))THEN
#ifndef GFORTRAN
C==      the call below for SunOS Fortran and g77
        OPEN(LOUT,FILE=FILEW,STATUS='old',FORM='unformatted',
     1                ACCESS='append')
#else
C==      the call below for gfortran. g77 does not recognize POSITION
        OPEN(LOUT,FILE=FILEW,STATUS='old',FORM='unformatted',
     1                POSITION='append')
#endif
        CALL WRITESTORE(LOUT,1,INTV,RLV,NAMER,NAMEW,1,IERROR,
     :          INTV(INROWS),INTV(INUP),INTV(INE),INTV(INBP),
     :          INTV(INBP2),INTV(INFP),INTV(INEL),
     :          INTV(INUL),INTV(INPM),INTV(INSM),INTV(INBL),
     :          INTV(INELLEP),INTV(INSEG),
     :          EX,EY,UVP,LEM,NOR,VHB,THDINT,
     :          IMAT,VOLD,SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :          IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC1,IFBC2,
     :          IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :          UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :          IELLE,IPOLYN)
        CLOSE(LOUT)
        WRITE(LSC,10121)FILEW(1:LNWFF)
        WRITE(LUW,10121)FILEW(1:LNWFF)
        IF (ILGSAVE.NE.0) THEN
          OPEN(LLG,FILE=LGFILE,STATUS='old',IOSTAT=IOS)
          CALL WRITELGMESH( INTV(IILAG),LGEM,EXLG,EYLG,
     :                          LGIBC,LGIBCF,
     :                          STELPX,STELPY,RLV(ISTELPR),
     :                          INTV(INEL),INTV(INUL),INTV(INBL),
     :                          INTV(INSM),INTV(INPM),INTV(INRM),
     :                          LLG,LSC,LUW,IERROR)
          CLOSE(LLG)
        ENDIF
C
C    save SERIES data if required; different OPEN statements
C    for SunOS Fortran and g77
C
C        IF (INTV(IMSINDX).GT.0.) THEN
C#ifndef GFORTRAN
C          OPEN(LSERDAT,FILE=FILESER,STATUS='old',
C     1                ACCESS='append',IOSTAT=IOS)
C#else
C          OPEN(LSERDAT,FILE=FILESER,STATUS='old',
C     1                POSITION='append',IOSTAT=IOS)
C#endif
C          CALL WRITESERIESDATA(LEM,NOR,INTV(INE),INTV(INUP),
C     :          MSNODE,MEASUR,INTV(IMSINDX),
C     :          PTLOC,EX,EY,RLV(ITIMEL),LSERDAT,IERROR)
C          CLOSE(LSERDAT)
C        ENDIF
       END IF
C
C     print the solution vector
C
      IF(INTV(IIPRINT).NE.0)THEN
        WRITE(LUW,10106)RLV(ITIMEL),INTV(IKSTEP)
        WRITE(LUW,10107)
         CALL MATPRT(UVP(1),INTV(INUP),2,INTV(INUP2),LUW)
        IF(INTV(INCOMP).NE.0)THEN
          WRITE(LUW,10120)
          CALL MATPRT(UVP(INTV(INUP2)+INTV(INFP)+1),INTV(INN),
     :                1,INTV(INN),LUW)
        END IF
        IF(INTV(IKSTEP).EQ.0)GO TO 54
        WRITE(LUW,10109)
        CALL MATPRT(SSQ(1),INTV(INUP),1,INTV(INUP),LUW)
   54   IF(INTV(IIVV).EQ.0)GO TO 55
      END IF
C
C     get elapsed time
C
   55 CALL DATIME(IDAY,ISEC,0)
      ITHOUR=ISEC/3600
      ITMIN=ISEC/60-60*ITHOUR
      ITSEC=ISEC-3600*ITHOUR-60*ITMIN
      WRITE(LUW,12101)IDAY(3),IDAY(2),MOD(IDAY(1),100),
     :                ITHOUR,ITMIN,ITSEC
      WRITE(LSC,12101)IDAY(3),IDAY(2),MOD(IDAY(1),100),
     :                ITHOUR,ITMIN,ITSEC
      WRITE(LUW,10115)
      WRITE(LSC,10115)
C      do not overwrite an error setting
      IF (IERROR.EQ.0) THEN
C      do not report a small timestep error if input file forces dt=0
        IF(DT.LE.1.E-4.AND.
     :       (INTV(IKSTEP).NE.0.AND.RLV(ITEXIT).NE.0.0))
     :   IERROR=1
      END IF
  555 RETURN
C
10098 FORMAT(A16)
12101 FORMAT('Job finished on ',I2.2,'/',I2.2,'/',I2.2,' at ',
     : I2,':',I2.2,':',I2.2,/)
10100 FORMAT(I2.2,'/',I2.2,'/',I2.2)
10101 FORMAT(I2,':',I2.2,':',I2.2)
C10102 format(' execution required ',I6,' seconds')
10103 FORMAT('KSTP    TIME ITER   AREA       VOLUME    ',
     1'VISDIS     DVMX   URMS    VRMS    ISET')
10105 format(' rotation rate of indenter is',g12.5,/)
10106 format(' execution stopped at TIMEL =',g14.5,
     1'    after',i6,'  timesteps',//)
10107 format(' the velocity vector follows ....',//)
10108 format(' the lithosphere has been set to a uniform ',
     1'non-dimensional thickness of 1.0/',f4.1,/)
10109 format(' the crustal thickness matrix follows ..',//)
10111 FORMAT(I4,F8.4,I5,G12.5,2G10.3,F8.3,2F8.5,I5)
10113 format(' RESETTING IGRID to 0 (IGRID=1,2 only valid with a fault',
     1//)
10115 FORMAT('*********************************************',
     1'*********************************',/)
10114 format(' viscosity is assumed constant',/)
10116 format(' viscosity is assumed newtonian but spatially varying',/)
10117 format(' viscosity is a function of strain rate',
     1' i.e., fluid is non-newtonian, with stress exponent',
     2f6.2,/)
10118 format(' viscosity is spatially varying (set from input file)',
     1' and is non-newtonian, with stress exponent',f6.2,/)
10119 format(' viscosity is spatially varying and is non-newtonian, ',
     1' with stress exponent spatially varying (set from input file)',/)
10120 FORMAT(' the pressure vector follows ....',//)
10121 format('The solution for this step has been written to: ',A32)
10122 format('Velocity has been set to zero on all fixed velocity',
     1       ' external boundaries.',/)
10123 format('Traction has been set to zero on all fixed traction',
     1       ' external boundaries.',/)
10124 format(' Last calculation is saved in the file ',a16,/)
10125 FORMAT('      to file ',A40,/)
10126 FORMAT('Calculation terminating due to mesh distortion')
10127 FORMAT('Calculation terminating as timestep is now < 0.0002')
10128 FORMAT(' WARNING: TSAVE set to',F7.4,' as IDT0 is',I5)
10129 FORMAT('REMESHING: Time is ',F9.4)
C10130 FORMAT('Mesh: SEGMIN = ',G12.5,' at N = ',I5,'   SEGLIM = ',G12.5,
C     :  /,'    SHAPEMIN = ',G12.5,' at N = ',I5,' SHAPELIM = ',G12.5,
C     :  /,'      ANGMIN = ',G12.5,' at N = ',I5,'   ANGLIM = ',G12.5)
C10131 FORMAT('       IGRID = ',I4,'  REMESHFLAG = ',L3,'  NSEG = ',I7)
10132 FORMAT('Remeshed solution now saved to file: ',A32,/)
10133 FORMAT('Mesh condition: SEGMIN = ',G12.5,
     :       ' SHAPEMIN =',G12.5,' ANGMIN =',G12.5)
10202 FORMAT("Problem deallocating memory for protected nodes.")
13001 FORMAT("Record limit exceeded, job stopping. IWRITE >",I4,/)
13002 FORMAT("KSTEP =",I6,", TIMEL =",F9.4,' KEXIT =',I6,
     1 ' TEXIT =',F9.4,/)
C
      END

      SUBROUTINE CALCMESHXY(NX,NY,IMSH,XLEN,YLEN)
C
C    NY and NX determine the number of elements spanning the mesh
C    in Y and X direction respectively in a rectangular area of dimension
C    YLEN*XLEN.  An approximately hexagonal (isotropic) mesh results
C        if NX = (sqr3)*NY*(XLEN/YLEN)/2 approximately.
C           (NX and NY should be even numbers >=4 if IMSH=0)
C
C    IMSH=0  or -1 - triangular mesh with hexagonal symmetry (-1 has
C      original configuration, 0 has row order switched)
C    IMSH=1 or 2 - triangular mesh obtained by subdividing rectangles
C      (orientation of subdividing lines is opposite from 1 to 2)
C
      IF (IMSH.LE.0) THEN
        IF (NX.EQ.0) THEN
          IF (MOD(NY,2).NE.0) NY=NY+1
          IF (NY.LT.4) NY=4
          NX = SQRT(3.0)*NY*(XLEN/YLEN)/2
          IF (MOD(NX,2).NE.0) NX=NX+1
          IF (NX.LT.4) NX=4
        ELSE
          IF (MOD(NX,2).NE.0) NX=NX+1
          IF (NX.LT.4) NX=4
          NY = NX*2/SQRT(3.0)*YLEN/XLEN
          IF (MOD(NY,2).NE.0) NY=NY+1
          IF (NY.LT.4) NY=4
        END IF
      ELSE IF (IMSH.EQ.1.OR.IMSH.EQ.2) THEN
        IF (NX.EQ.0) NX = XLEN*NY/YLEN
        IF (NY.EQ.0) NY = YLEN*NX/XLEN
      END IF
      RETURN
      END

      SUBROUTINE INITDIMENSIONS(INTV,RLV,
     :                          NUP,NN,NMP,NBP,NE,NFP,NFPF3,NSEG,
     :                          NX,NY,IMSH,IFLT,MAXNBS,NELLEP,
     :                          XLEN,YLEN,EX,EY,VC,SE,
     :                          NXLVAL,NYLVAL,NPMP,NRMP,NSM,
     :                          ITHDI,ILAG,NUL,IVIS,IPOLY,IDEN,ITEMP,
     :                          ICR,IVOLD,IMREG,NCOMP,KPBC,PDIST)
C    :                          KBW,K2BW,K4BW,NK,NK2,NK4)
      INCLUDE "indices.parameters"
      COMMON/LGDIM/NXL,NYL,NNL,NML,ILGSAVE
      DIMENSION INTV(64),RLV(64)
      DIMENSION EX(NUP),EY(NUP)
      EPS=1E-4
C
C      check whether nx or ny needs to be calculated
C
      IF (NX.EQ.0.OR.NY.EQ.0)
     :  CALL CALCMESHXY(NX,NY,IMSH,XLEN,YLEN)
      RLV(IXLEN) = XLEN
      RLV(IYLEN) = YLEN
      IFLTTYP=IFLT
      IF (IFLT.EQ.4) IFLTTYP=3
      IFAC=MAX(1,IFLTTYP)
      IADD=MAX(0,1-IFLTTYP)
      IMUL=MIN(1,IFLTTYP)
      IMSHQ = 1-MIN(1,IMSH)
      INTV(INX)=NX
      INTV(INY)=NY
      INTV(IIMSH)=IMSH
      INTV(IIFLT)=IFLT
      INTV(IIVIS)=IVIS
      INTV(IIVV)=0
      INTV(IIPOLY)=IPOLY
      INTV(INCOMP)=NCOMP
      INTV(IKPBC)=KPBC
C
      IF (IMSH.EQ.3) THEN
        IF (IFLT.EQ.3.OR.IFLT.EQ.4) IFAC=1
        INTV(INE)=NE
        INTV(INN)=NN
        INTV(INMP)=NMP
        INTV(INBP)=NBP
        INTV(INBP2)=2*INTV(INBP)
        INTV(INSEG)=NSEG
C
C    if periodic boundary conditions invoked, find how many points
C    are used. Note that this method relies on domain being rectangular
C
        IF(KPBC.NE.0)THEN
          YMIN = EY(1)
          YMAX = EY(1)
          XMIN = EX(1)
          XMAX = EX(1)
          DO 9 I=2,NUP
            IF (EY(I).GT.YMAX) YMAX = EY(I)
            IF (EY(I).LT.YMIN) YMIN = EY(I)
            IF (EX(I).GT.XMAX) XMAX = EX(I)
            IF (EX(I).LT.XMIN) XMIN = EX(I)
   9      CONTINUE
          NXL=0
          NYL=0
          NXLMIN=0
          NXLMAX=0
          DO 10 N = 1,NUP
            IF (EY(N).EQ.YMIN) NXL=NXL+1
            IF (EX(N).EQ.XMIN) NYL=NYL+1
  10      CONTINUE
C         IF(KPBC.EQ.1)THEN
C           INTV(INFP)=NYL*2
C         ELSE IF(KPNC.EQ.2)THEN
C           INTV(INFP)=NXL*2
C         ENDIF
        ENDIF
C       IF(IFLT.GE.3)THEN
C         INTV(INFP)= NFP + INTV(INFP)
C       ENDIF
C  if 'PB' in input, periodic bnd pts have been added to NFP in trimesh 
        INTV(INFP)=NFP
        INTV(INFPF3)=NFPF3
        INTV(INUP)=NUP
C
C    For the meshes in which triangles are formed by subdivision
C      of rectangles
C
      ELSE
        IF(IMSH.GE.1)THEN
          INTV(INE)=2*NX*NY*IFAC
          INTV(INN)=(NY+1)*(NX+1)*IFAC
          INTV(INMP)=(3*NY*NX+NX+NY)*IFAC
          INTV(INBP)=4*(NX+NY)*IFAC
C
C    For the mesh with hexagonal symmetry (original configuration)
C
        ELSE IF (IMSH.EQ.-1)THEN
          INTV(INE)=(2*NX*NY+NY)*IFAC
          INTV(INN)=((NY+1)*(NX+1)+NY/2+1)*IFAC
          INTV(INMP)=(3*NY*NX+NX+5*NY/2+1)*IFAC
          INTV(INBP)=(4*(NX+NY)+4)*IFAC
C
C    For the mesh with hexagonal symmetry (modified configuration)
C
        ELSE IF(IMSH.EQ.0) THEN
          INTV(INE)=(2*NX*NY+NY)*IFAC
          INTV(INN)=((NY+1)*(NX+1)+NY/2)*IFAC
          INTV(INMP)=(3*NY*NX+NX+5*NY/2)*IFAC
          INTV(INBP)=4*(NX+NY)*IFAC
        END IF
C
        INTV(INBP2)=2*INTV(INBP)
        IF(KPBC.EQ.1)INTV(INFP)=2*(NY*2+1)*IMUL
        IF(KPBC.EQ.2)INTV(INFP)=2*(NX*2+1)*IMUL
        INTV(INUP)=INTV(INN)+INTV(INMP)
      END IF
C
C   space is allocated for pressure on nodes, even when NCOMP<1
C   but matrix solver does not use this part of vector
C
      INTV(INUP2)=2*INTV(INUP)
      INTV(INROWS)=INTV(INUP2)+INTV(INFP)+INTV(INN)
C     IF(INTV(INCOMP).GE.1)INTV(INROWS)=INTV(INROWS)+INTV(INN)
C
C     initialise the lagrangian common block and array dimensions
C
      NXL = NXLVAL
      NYL = NYLVAL
      INTV(IILAG)=ILAG
      INTV(INEL)=2*(NXL*NYL)*IFAC
      NML=(3*NYL*NXL+NXL+NYL)*IFAC
      NNL=((NYL+1)*(NXL+1))*IFAC
      NUL=NNL+NML
      INTV(INUL)=NUL
      INTV(INBL)=4*(NYL+NXL)*IFAC
      INTV(INSM)=NSM
      INTV(INPM)=NPMP
      INTV(INRM)=NRMP
C
C     for stiffness matrix array dimensions STK, STK2, STK4
C
      CALL INITSTKDIM(INTV,IMSH,MAXNBS)
C
      INTV(INELLEP)=NELLEP
      INTV(IIDEN)=IDEN
      INTV(IITEMP)=ITEMP
      INTV(IICR)=ICR
      INTV(IIVOLD)=IVOLD
      INTV(IIMREG)=IMREG

      INTV(IIBCMOD)=0
      INTV(IIDT0)=25
      INTV(IKEXIT)=100
      INTV(IIPRINT)=0
      INTV(IITSTOP)=100
      INTV(IKSAVE)=4
      INTV(IISAVE)=-10000
      INTV(IIWRITE)=10
      INTV(IINDFIX)=0
      INTV(IJLV)=0
      INTV(ILFLAT)=0
      INTV(IIFCASE)=3
      INTV(IMPDEF)=5
      INTV(IIDEFTYP) = 0
C
C    Initialise flags specifying which data will be saved
C
      INTV(IIDATA)=0
      INTV(IMSINDX)=0
      INTV(IRHEOTYP)=0
      INTV(IITHDI)=ITHDI
C
C    initialise time counters
C
      INTV(IKSTEP)=0

      RLV(ISTELPR)=0.05
      RLV(IBIG)=1.e+20
      RLV(IVELXO)=0.0
      RLV(IVELYO)=0.0
      RLV(ITBXOFF)=0.0
      RLV(ITBYOFF)=0.0
      RLV(IARGAN)=0.0
      RLV(IHLENSC)=50.0
      RLV(IBDEPSC)=0.35
      RLV(IRISOST)=7.75/100.0
      RLV(IREFLEV)=100.0
      RLV(ITDIFF)=1.0
      RLV(IALPHA)=3.0E-5
      RLV(IGAMMA)=0.0
      RLV(ITEXIT)=0.0
      RLV(IWFIT)=1.0
      RLV(IAC)=5.0e-6
      RLV(IACNL)=RLV(IAC)*10.0
      RLV(ITSAVE)=0.08
      RLV(IBANGL)=0.0
      RLV(IERA)=0.0
      RLV(IBRGAN)=0.0
      RLV(ITHRESH)=10.0
      RLV(IAREA)=0.0
      RLV(IVOLUM)=0.0
      RLV(IVISCD)=0.0
      RLV(IDVMX)=0.0
      RLV(IAREFM)=0.0
      RLV(IXREFM)=0.0
      RLV(IYREFM)=0.0
      RLV(IYLDSTR)=0.0
      RLV(IOMTOT)=0.0
C
C   DEFV is used to define the mesh distortion
C
      RLV(IDEFV1)=0.0
      RLV(IDEFV2)=0.0
      RLV(IDEFV3)=0.0
      RLV(IDEFV4)=0.0
      RLV(IDEFV5)=0.0
      RLV(IBCV1)=0.0
      RLV(IBCV2)=0.0
      RLV(IBCV3)=0.0
      RLV(IBCV4)=0.0
      RLV(IBCV5)=0.0
      RLV(IRHOG1)=0.0
      RLV(IRHOG2)=0.0
      RLV(IRHOG3)=0.0
      RLV(IRHOG4)=0.0
      RLV(IRHOG5)=0.0
C
C   IVISP is now used for RHEOMOD parameters MIB 0202
C
      RLV(IVISP1)=0.0
      RLV(IVISP2)=0.0
      RLV(IVISP3)=0.0
      RLV(IVISP4)=0.0
      RLV(IVISP5)=0.0
C
C    initialise time counters
C
      RLV(ITIMEL)=0.0
      IF((IMSH.EQ.3).AND.(NELLEP.GT.0))THEN
        RLV(IVC) = VC
        RLV(ISE) = SE
      END IF
C
C    check VC and SE are sensible
C
      IF(ABS(RLV(ISE)).LT.EPS) RLV(ISE)=1.0
      IF(ABS(RLV(IVC)).LT.EPS) RLV(IVC)=1.0
C
C    the next 3 parameters were only relevant to the maspar version
C
      INTV(IMROWS)=4
      INTV(INXMAX)=((NX+IMSHQ)*2+1)*IFAC
      INTV(INYMAX)=NY*2+1
      RETURN
      END
C
      SUBROUTINE INITSTKDIM(INTV,IMSH,MAXNBS)
      INCLUDE "indices.parameters"
      PARAMETER (KBWP=32, K4BWP=38)
      COMMON/STKDIM/KBW,K2BW,NK,NK2,K4BW,NK4
      SAVE /STKDIM/
      DIMENSION INTV(64)
C
C   if mesh generated by triangle
C
      IF (IMSH.EQ.3.AND.MAXNBS.NE.0) THEN
        KBW = 2*(MAXNBS*3+1)
        K4BW = KBW+MAXNBS+1
C
C   regular (original) mesh
C
      ELSE
        KBW = KBWP
        K4BW = K4BWP
      ENDIF
      K2BW = 2
      NK = INTV(INUP2)*KBW
      NK2 = INTV(INUP2)*K2BW
      NK4 = INTV(INN)*K4BW
C     WRITE(*,*)'INITSTKDIM: NK =',NK,' NK2 =',NK2,' NK4 =',NK4
      RETURN
      END

      SUBROUTINE INITREMESHDIM(INTV,RLV,
     :                          NUP,NN,NMP,NBP,NE,NFP,NFPF3,NSEG,
     :                          NX,NY,IMSH,IFLT,MAXNBS,NELLEP,
     :                          XLEN,YLEN,EY,
     :                          NXLVAL,NYLVAL,NPMP,NRMP,NSM)
      INCLUDE "indices.parameters"
      PARAMETER (KBWP=32, K2BWP=2, K4BWP=38)
C       Meshing from poly file - allow for 12 neighbouring elements
      PARAMETER (KBWPOLYP=62, K2BWPOLYP=2, K4BWPOLYP=74)
      COMMON/STKDIM/KBW,K2BW,NK,NK2,K4BW,NK4
      SAVE /STKDIM/
      COMMON/LGDIM/NXL,NYL,NNL,NML,ILGSAVE
      DIMENSION INTV(64),RLV(64)
      DIMENSION EY(NUP)
      IF (IFLT.EQ.4) IFLTTYP=3
      ELSE IFLTTYP=IFLT
      IFAC=MAX(1,IFLTTYP)
      IADD=MAX(0,1-IFLTTYP)
      IMUL=MIN(1,IFLTTYP)
      IMSHQ = 1-MIN(1,IMSH)
      INTV(INX)=NX
      INTV(INY)=NY
      INTV(IIMSH)=IMSH
      INTV(IIFLT)=IFLT
      INTV(INFP)=NFP
      INTV(INFPF3)=NFPF3
C
        IF (IFLT.EQ.3.OR.IFLT.EQ.4) IFAC=1
        INTV(INE)=NE
        INTV(INN)=NN
        INTV(INMP)=NMP
        INTV(INBP)=NBP
        INTV(INBP2)=2*INTV(INBP)
        INTV(INSEG)=NSEG
C
        IF (IFLT.EQ.1) THEN
          YMIN = EY(1)
          YMAX = EY(1)
          DO 9 I=2,NUP
            IF (EY(I).GT.YMAX) YMAX = EY(I)
            IF (EY(I).LT.YMIN) YMIN = EY(I)
   9      CONTINUE
          NXL=0
          DO 10 N = 1,NUP
            IF (EY(N).EQ.YMAX) NXL=NXL+1
            IF (EY(N).EQ.YMIN) NXL=NXL+1
  10      CONTINUE
          NYL = NBP - NXL + 4
C
          INTV(INFP)=NYL*IMUL
        END IF
C
      INTV(INUP)=NUP
      INTV(INUP2)=2*INTV(INUP)
      INTV(INROWS)=INTV(INUP2)+INTV(INFP)+INTV(INN)
C
C     initialise the lagrangian common block and array dimensions
C
C
C     initialise the common block for STK array dimensions
C
      IF (IMSH.EQ.3.AND.MAXNBS.NE.0) THEN
        KBW = 2*(MAXNBS*3+1)
        K2BW = K2BWPOLYP
        K4BW = KBW+MAXNBS+1
      ELSE
        KBW = KBWP
        K2BW = K2BWP
        K4BW = K4BWP
      ENDIF
      NK = INTV(INUP2)*KBW
      NK2 = INTV(INUP2)*K2BW
      NK4 = INTV(INN)*K4BW

C
C   May need to adjust RLV values for BCONDMOD etc
C
      RETURN
      END

      SUBROUTINE INITARRAYS(LEM,NOR,DENS,TEMPT,VHB,THDISS,THDINT,
     :                      IMAT,VOLD,QBND,POLEP,IBPOLE,TAPERF,UVP,
     :                      QLOAD,IV,IBC,IBNGH,IBCTYP,SSQ,FROT,
     :                      IELFIX,IFLTTIPS,NUP,NE,NROWS,NBP,RHO,
     :                      SE,HLENSC)
      INCLUDE 'indices.parameters'
      INCLUDE 'limits.parameters'
      COMMON/BT/MSINDX,MEASUR(MAXMEAS),MSNODE(MAXMEAS),PTLOC(3,MAXMEAS)
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      DIMENSION LEM(6,NE),NOR(NUP)
      DIMENSION UVP(NROWS)
      DIMENSION QLOAD(NROWS)
      DIMENSION VHB(8,NE)
      DIMENSION THDISS(7,NE),THDINT(7,NE)
      DIMENSION VOLD(8,NE)
      DIMENSION IMAT(NE)
      DIMENSION SSQ(NUP)
      DIMENSION FROT(NUP)
      DIMENSION DENS(7,NE)
      DIMENSION TEMPT(NUP)
      DIMENSION IELFIX(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION IBC(NBP),IBNGH(NBP*2),IBCTYP(NBP*2)
      DIMENSION POLEP(3,MAXPOLE)
      DIMENSION IBPOLE(NBP*2),TAPERF(NBP)
      DIMENSION IV(64)
      DIMENSION IFLTTIPS(MAXTIPS*2)
C
C     initialise the matrices
      IF (IV(IIMSH).NE.3) THEN
        DO 2 J=1,NE
          DO 12 I=1,6
            LEM(I,J)=0
   12     CONTINUE
    2   CONTINUE
        DO 10 I=1,NBP
          IBC(I)=0
          IBNGH(I)=0
          IBNGH(I+NBP)=0
   10   CONTINUE
      END IF
C       polep array created in basil.F (runbasil)
      DO K=1,MAXPOLE
        DO J=1,3
          POLEP(J,K)=0.0
        ENDDO
      ENDDO
C
C    Initialise viscosity/exponent array to background constant
C    and material property array to default (lithosphere).
C    VHB array only used if IVIS>0
C    VHB array initialised when the array is allocated
C    Viscosity values may be read from a poly file so
C    VHB array is already initialised if IMSH=3
C
C     IF (IV(IIVIS).GT.0) THEN
C     DO 3 J=1,NE
C       DO 13 I=1,7
C         IF (IV(IIMSH).NE.3) VHB(I,J)=1.0
C  13   CONTINUE
C   3 CONTINUE
C     END IF

      IF(IV(IITHDI).GT.0)THEN
      DO 203 J=1,NE
        DO 213 I=1,7
          THDISS(I,J)=0.0
          THDINT(I,J)=0.0
  213   CONTINUE
  203 CONTINUE
      END IF

      IF (IV(IIMREG).NE.0) THEN
        DO 4 J=1,NE
          IMAT(J) = 1
    4   CONTINUE
      ENDIF

      IF (IV(IIDEN).NE.0) THEN
        DO 5 J=1,NE
          DO 15 I=1,7
            DENS(I,J)=RHO
   15     CONTINUE
    5   CONTINUE
      ENDIF

      IF (IV(IIVOLD).NE.0) THEN
        DO 6 J=1,NE
          DO 16 I=1,7
            VOLD(I,J)=SE
   16     CONTINUE
    6   CONTINUE
      ENDIF

C  This should be done if not LE mesh 3  ???
      IF (IV(IIMSH).NE.3.AND.IV(IIVIS).GT.0) THEN
        DO 7 J=1,NE
          VHB(8,J)=SE
    7   CONTINUE
      END IF

Cle      DO 17 J=1,NBP
Cle        QBND(J)=0.0
Cle        QBND(J+NBP)=0.0
Cle        IBPOLE(J)=0
Cle        IBPOLE(J+NBP)=0
Cle        TAPERF(J)=0.0
Cle   17 CONTINUE

Cle      DO 18 I=1,NROWS
Cle        UVP(I)=0.0
Cle        QLOAD(I)=0.0
Cle   18 CONTINUE

Cle      IF (IV(IIMSH).NE.3) THEN
Cle        DO 19 J=1,NBP*2
Cle          IBCTYP(J)=-1
Cle   19   CONTINUE
Cle      END IF
 
      DO 9 I=1,NUP
        NOR(I)=I
C   arbitrary velocity field to start (to help VISK)
C       UVP(I)=1.0/float(I)
    9 CONTINUE
C
C    initial temperature
C
Cle      IF (IV(IITEMP).NE.0) THEN
Cle        DO 20 I=1,NUP
Cle          TEMPT(I)=1.0
Cle   20   CONTINUE
Cle      END IF
C
C    initial crustal thickness and rotation
C
      IF (ISSQACTIVE.EQ.1) THEN
        DO 21 I=1,NUP
Cle          IELFIX(I)=0
          SSQ(I)=DFLTSSQ
   21   CONTINUE
      END IF
C
C    initial temperature, crustal thickness and rotation
C
Cle      IF (IROTACTIVE.EQ.1) THEN
Cle        DO I=1,NUP
Cle          FROT(I)=0.0
Cle        ENDDO
Cle      END IF

      DO 117 I=1,MAXTIPS*2
        IFLTTIPS(I) = 0
 117  CONTINUE
C
C    initialise arrays used if saving to a series data file
C
      DO 118 I=1,MAXMEAS
        MEASUR(I) = 0
        MSNODE(I) = 0
 118  CONTINUE
      RETURN
      END

      SUBROUTINE INITREMESHARRAYS(IV,NOR,DENS,VHB,
     :                            RHO,SE,VHB_VAL)
      INCLUDE 'indices.parameters'
      INCLUDE 'limits.parameters'
      DIMENSION IV(64)
      DIMENSION NOR(IV(INUP))
      DIMENSION DENS(7,IV(INE))
      DIMENSION VHB(8,IV(INE))
C
C     initialise the matrices
      NE = IV(INE)
      NROWS = IV(INROWS)
      NUP = IV(INUP)
C      IF (IV(IIDEN).NE.0) THEN
C        DO 5 J=1,NE
C          DO 15 I=1,7
C            DENS(I,J)=RHO
C   15     CONTINUE
C    5   CONTINUE
C      ENDIF
      IF (IV(IIVV).GE.1) THEN
        DO J=1,NE
C         DO  I=1,7
C           VHB(I,J)=VHB_VAL
C         END DO
          VHB(8,J)=SE
        END DO
      ENDIF

      RETURN
      END

      SUBROUTINE INITCOMMENTS(TXT,LEN)
C     COMMON/V01S/NDUM2,NDATE,NDUM3,COMMEN
C     CHARACTER NDUM2*16,NDUM3*16,NDATE*16
      CHARACTER COMMEN*80
      CHARACTER TXT*80

      DO 10 I=1,80
        COMMEN(I:I)=' '
  10  CONTINUE
      IF (LEN.GT.0) COMMEN(1:LEN) = TXT(1:LEN)
      RETURN
      END

      SUBROUTINE INITLGARRAYS(EXLG,EYLG,UXLG,UYLG,NUL)
      DIMENSION EXLG(NUL),EYLG(NUL),UXLG(NUL),UYLG(NUL)
C
C    lagrangian arrays
C
      DO 50 I=1,NUL
        EXLG(I)=0.0
        EYLG(I)=0.0
        UXLG(I)=0.0
        UYLG(I)=0.0
   50 CONTINUE
      RETURN
      END

      SUBROUTINE INITMARKERARRAYS(UVELPX,UVELPY,NSM,NPM)
      DIMENSION UVELPX(NPM,NSM),UVELPY(NPM,NSM)
C
C    strain marker arrays
C
      DO  55 JL=1,NSM
      DO  60 IL=1,NPM
        UVELPX(IL,JL)=0.0
        UVELPY(IL,JL)=0.0
   60 CONTINUE
   55 CONTINUE
      RETURN
      END

      SUBROUTINE INITVARS(INTV,RLV,IOFF,ITIMETYPE,IDENMOD,
     :                IRHEOMOD,MFIXFLAG,NITER,IGRID,NDFORM,RHO,TILTDEG,
     :                TTSV,CENTLNG,INITVEL,NEGTIM)
      INCLUDE 'indices.parameters'
      INCLUDE 'limits.parameters'
      COMMON/BT/MSINDX,MEASUR(MAXMEAS),MSNODE(MAXMEAS),PTLOC(3,MAXMEAS)
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      DIMENSION INTV(64)
      DIMENSION RLV(64)
      LOGICAL INITVEL
C
C    default values for physical parameters
C
      IOFF=0
      NDFORM=0
      TILTDEG=0.0
      CENTLNG=0.0
      IF (INTV(IIFLT).EQ.4) IGRID=2
      ELSE IGRID=(INTV(IIFLT)+1)/2
      IDENMOD=0
      IRHEOMOD=0
      TTSV=0.0
      NEGTIM=1
      RHO=1.0
      MFIXFLAG = 0
      NITER = INTV(INROWS)*3
C   initialise to explicit timestepping - 1=explicit, 2=RK
      ITIMETYPE = 1
C
C    initialise flag controlling initialisation of velocity field
C
      INITVEL=.FALSE.
C
C    Set flags controlling use of SSQ, IELFIX, FROT
C
      DFLTSSQ=-ALOG(RLV(IHLENSC))
      ISSQACTIVE = 0
      IROTACTIVE = 0
      IF (INTV(IIGRAV).EQ.4 .OR. INTV(IICR).EQ.1 .OR. INTV(IICR).EQ.2
     : .OR. INTV(IJLV).NE.0 )  ISSQACTIVE = 1
      IF (INTV(IICR).EQ.1 .OR. INTV(IICR).EQ.3)
     :  IROTACTIVE = 1
C
      RETURN
      END

      SUBROUTINE PRINTINP(IV,RV,IOFF,ITIMETYPE,IDENMOD,IRHEOMOD,
     :                  MFIXFLAG,NITER,IGRID,NDFORM,NEGTIM,NFPF3,
     :                  RHO,TILTDEG,TTSV,SHAPELIM,SEGLIM,ANGLIM)
      INCLUDE 'indices.parameters'
      COMMON/AI/LUW,LSC,LBC,LLG
      COMMON/LGDIM/NXL,NYL,NNL,NML,ILGSAVE
      DIMENSION IV(64)
      DIMENSION RV(64)
C     DIMENSION DEFV(5), BCV(5), RHOG(5), VISP(5)
C
C     record the runtime values of variables which may have been
C     altered by the input file
C
      WRITE(LUW,10110)
      WRITE(LSC,10110)
10110 FORMAT('Input Variables:',/)
      IF (IV(IKPBC).NE.0) THEN
        WRITE(LSC,10134)IV(INFP)
        WRITE(LUW,10134)IV(INFP)
10134   FORMAT('Periodic boundary is set, KPBC = ',I4)
      ENDIF
      IF ((IV(IIFLT).NE.0).AND.(NFPF3.GT.0)) THEN
        WRITE(LSC,10133)NFPF3
        WRITE(LUW,10133)NFPF3
10133   FORMAT('One or more internal faults is set, NFPF3 = ',I4)
      END IF
      WRITE(LUW,10112)IV(INX),IV(INY),IV(IIMSH),IV(IIFLT),
     :                IV(INE),IV(INN)
      WRITE(LSC,10112)IV(INX),IV(INY),IV(IIMSH),IV(IIFLT),
     :                IV(INE),IV(INN)
10112 FORMAT('NX        = ',I7,'  NY        = ',I7,
     :      '  MESH      = ',I7,/,
     :      'FAULT     = ',I7,'  NE        = ',I7,
     :      '  NN        = ',I7)
      WRITE(LUW,10114)IV(INMP),IV(INUP),IV(INBP)
      WRITE(LSC,10114)IV(INMP),IV(INUP),IV(INBP)
10114 FORMAT('NMP       = ',I7,'  NUP       = ',I7,
     :      '  NBP       = ',I7)
      WRITE(LUW,10115)IV(INFP),IV(INFPF3),IV(IKPBC)
      WRITE(LSC,10115)IV(INFP),IV(INFPF3),IV(IKPBC)
10115 FORMAT('NFP       = ',I7,'  NFPF3     = ',I7,
     :      '  KPBC      = ',I7)
      IF (IV(IILAG).EQ.1.OR.IV(IILAG).EQ.3) THEN
        WRITE(LUW,10116)NXL,NYL,IV(INEL),IV(INUL),IV(INBL)
        WRITE(LSC,10116)NXL,NYL,IV(INEL),IV(INUL),IV(INBL)
10116   FORMAT('NXL       = ',I7,'  NYL       = ',I7,
     :      '  NEL       = ',I7,/,
     :      'NUL       = ',I7,'  NBL       = ',I7)
      END IF
      IF (IV(IILAG).EQ.1.OR.IV(IILAG).EQ.2) THEN
        WRITE(LUW,10118)IV(INSM),IV(INPM),IV(INRM)
        WRITE(LSC,10118)IV(INSM),IV(INPM),IV(INRM)
10118   FORMAT('NSM       = ',I7,'  NPM       = ',I7,
     :      '  NRM       = ',I7)
      END IF
      WRITE(LUW,10200)RV(IXLEN),RV(IYLEN),IV(INCOMP),IV(IIGRAV),
     :          RV(IHLENSC),RV(IBDEPSC),RV(IRISOST),RV(IREFLEV)
      WRITE(LSC,10200)RV(IXLEN),RV(IYLEN),IV(INCOMP),IV(IIGRAV),
     :          RV(IHLENSC),RV(IBDEPSC),RV(IRISOST),RV(IREFLEV)
10200 FORMAT('XLEN      = ',F7.3,'  YLEN      = ',F7.3,
     :      '  NCOMP     = ',I7,/,
     :      T22,'IGRAV     = ',I7,'  HLENSC    = ',F7.2,/,
     :      T22,'BDEPSC    = ',F7.2,'  RISOST    = ',F7.2,/,
     :      T22,'REFLEV    = ',F7.2)
      WRITE(LUW,10201)RV(IAC),RV(IACNL),RV(IWFIT),NITER,IV(IITSTOP)
      WRITE(LSC,10201)RV(IAC),RV(IACNL),RV(IWFIT),NITER,IV(IITSTOP)
10201 FORMAT('AC        = ',E7.1,'  ACNL      = ',E7.1,
     :       '  WFIT      = ',F7.4,/,
     1       T22,'NITER     = ',I7,'  ITSTOP    = ',I7)
      WRITE(LUW,10202)RV(IARGAN),RV(IBRGAN),RV(ITHRESH),RHO,
     :      RV(IRHOG1),RV(IRHOG2),RV(IRHOG3),RV(IRHOG4),RV(IRHOG5)
      WRITE(LSC,10202)RV(IARGAN),RV(IBRGAN),RV(ITHRESH),RHO,
     :      RV(IRHOG1),RV(IRHOG2),RV(IRHOG3),RV(IRHOG4),RV(IRHOG5)
10202 FORMAT('ARGAN     = ',F7.2,'  BRGAN     = ',F7.2,
     :       '  THRESH    = ',G15.8,/,
     :       T22,'RHO       = ',F7.4,/T22,'RHOG = ',5F7.3)
      write(LUW,10203)IOFF,RV(IVELXO),RV(IVELYO),RV(IYLDSTR)
      write(LSC,10203)IOFF,RV(IVELXO),RV(IVELYO),RV(IYLDSTR)
10203 FORMAT('IOFF      = ',I7,'  VELXO     = ',F7.2,
     :       '  VELYO     = ',F7.2,/,
     :       T22,'YLDSTR    = ',F7.2)
      WRITE(LUW,10204)IV(IIDEFTYP),RV(IBANGL),RV(IERA),TILTDEG,
     :                NDFORM,RV(IDEFV1),RV(IDEFV2),RV(IDEFV3),
     :                RV(IDEFV4),RV(IDEFV5)
      WRITE(LSC,10204)IV(IIDEFTYP),RV(IBANGL),RV(IERA),TILTDEG,
     :                NDFORM,RV(IDEFV1),RV(IDEFV2),RV(IDEFV3),
     :                RV(IDEFV4),RV(IDEFV5)
10204 FORMAT('IDEFTYP   = ',I7,'  BANGL     = ',F7.3,
     :       '  ERA       = ',F7.3,/,
     :       T22,'TILTDEG   = ',F7.3,'  NDFORM    = ',I7/,
     :       T22,'DEFV = ',5F7.2)
      WRITE(LUW,10205)RV(ISE),RV(IVC),IV(IIVV),
     :                RV(IGAMMA),RV(IBETA)
      WRITE(LSC,10205)RV(ISE),RV(IVC),IV(IIVV),
     :                RV(IGAMMA),RV(IBETA)
10205 FORMAT('SE        = ',F7.2,'  VC        = ',F7.2,
     :       '  IVV       = ',I7/,
     :       'GAMMA     = ',F7.2,'  BETA      = ',F7.2)
      WRITE(LUW,10217)SHAPELIM,SEGLIM,ANGLIM
      WRITE(LSC,10217)SHAPELIM,SEGLIM,ANGLIM
10217 FORMAT('SHAPELIM  = ',F7.3,'  SEGLIM    = ',F7.3,
     :       '  ANGLIM    = ',F7.3)
      WRITE(LUW,10206)IV(IILAG)
      WRITE(LSC,10206)IV(IILAG)
10206 FORMAT('ILAG      = ',I7)
      WRITE(LUW,10207)IV(IKEXIT),RV(ITEXIT),IV(IIWRITE)
      WRITE(LSC,10207)IV(IKEXIT),RV(ITEXIT),IV(IIWRITE)
10207 FORMAT('KEXIT     = ',I7,'  TEXIT     = ',F7.3,
     :       '  IWRITE    = ',I7)
      WRITE(LUW,10208)RV(ITSAVE),IV(IKSAVE),NEGTIM
      WRITE(LSC,10208)RV(ITSAVE),IV(IKSAVE),NEGTIM
10208 FORMAT('TSAVE     = ',F7.4,'  KSAVE     = ',I7,
     :       '  NEGTIM    = ',I7)
      WRITE(LUW,10209)IV(IIDATA)
      WRITE(LSC,10209)IV(IIDATA)
10209 FORMAT('IDATA     = ',I7)
      WRITE(LUW,10210)ITIMETYPE,IV(IIDT0),IV(IMPDEF)
      WRITE(LSC,10210)ITIMETYPE,IV(IIDT0),IV(IMPDEF)
10210 FORMAT('ITIMETYPE = ',I7,'  IDT0      = ',I7,
     :       '  MPDEF     = ',I7)
      WRITE(LUW,10211)IV(IICR),IV(IIDEN),IV(IIMREG)
      WRITE(LSC,10211)IV(IICR),IV(IIDEN),IV(IIMREG)
10211 FORMAT('ICR       = ',I7,'  IDEN      = ',I7,
     :       '  IMREG     = ',I7)
      WRITE(LUW,10212)MFIXFLAG,IGRID
      WRITE(LSC,10212)MFIXFLAG,IGRID
10212 FORMAT('MFIXFLAG  = ',I7,'  IGRID     = ',I7)
      WRITE(LUW,10213)IDENMOD
      WRITE(LSC,10213)IDENMOD
10213 FORMAT('IDENMOD   = ',I7)
      WRITE(LUW,10214)IV(IIBCMOD),IV(IINDFIX),RV(IOMTOT),
     :             RV(IBCV1),RV(IBCV2),RV(IBCV3),RV(IBCV4),RV(IBCV5)
      WRITE(LSC,10214)IV(IIBCMOD),IV(IINDFIX),RV(IOMTOT),
     :             RV(IBCV1),RV(IBCV2),RV(IBCV3),RV(IBCV4),RV(IBCV5)
      WRITE(LUW,10215)IRHEOMOD,RV(IVISP1),RV(IVISP2)
      WRITE(LUW,10216)RV(IVISP3), RV(IVISP4),RV(IVISP5)
      WRITE(LSC,10215)IRHEOMOD,RV(IVISP1),RV(IVISP2)
      WRITE(LSC,10216)RV(IVISP3), RV(IVISP4),RV(IVISP5)
10215 FORMAT('IRHEOMOD  = ',I7,' RHEO(1) = ',F7.2,' RHEO(2) = ',F7.2)
10216 FORMAT(' RHEO(3) = ',F7.2,' RHEO(4) = ',F7.2,' RHEO(5) = ',F7.2)
10214 FORMAT('IBCMOD    = ',I7,'  INDFIX    = ',I7,
     :       '  OMTOT     = ',F7.2/,
     :       T22,'BCV = ',5F7.3)
      WRITE(LUW,'(/)')
      WRITE(LSC,'(/)')

      RETURN
      END

      SUBROUTINE MTCOPY(ARFROM,ARTO,IDIM,ARK,BRK)
C
C     routine copies ARFROM to ARTO, so that values are preserved
C     for future use.  Also permits a linear combination of the
C     2 arrays to be used
C
      DIMENSION ARFROM(IDIM), ARTO(IDIM)
      IF(ARK.EQ.0.0)THEN
        DO J=1,IDIM
          ARTO(J)=ARFROM(J)
        ENDDO
      ELSE
        DO J=1,IDIM
          ARTO(J)=ARK*ARTO(J)+BRK*ARFROM(J)
        ENDDO
      END IF
      RETURN
      END

      SUBROUTINE MTICOPY(IARFROM,IARTO,IDIM,ARK,BRK)
C
C     routine copies IARFROM to IARTO, so that values are preserved
C     for future use.  Also permits a linear combination of the
C     2 arrays to be used
C
      DIMENSION IARFROM(IDIM), IARTO(IDIM)
      IF(ARK.EQ.0.0)THEN
        DO J=1,IDIM
          IARTO(J)=IARFROM(J)
        ENDDO
      ELSE
        DO J=1,IDIM
          IARTO(J)=ARK*IARTO(J)+BRK*IARFROM(J)
        ENDDO
      END IF
      RETURN
      END
C
      SUBROUTINE FILADD (K,NEXT)
C===Position the record pointer of IUNIT K at EOF and return the
C   next-record pointer, NEXT.
      CHARACTER BUF*128
      INTEGER K,NEXT
      INTEGER ITMP(64)
      NEXT=1
      READ(K,END=4) BUF
      READ(K,END=4) 
      DO 1 NEXT=2,999999999
        READ(K,END=3)(ITMP(I),I=1,64)
        GOTO 1
C---To be conformable to the Fortran 77 Standard, ANSI X3.9-1978, and
C   proposed X3J3/S8.115 (Fortran 90), Subroutine FILEADD should not be
C   called when the file associated with Unit number K does not exist.
C   Reason: e.g. there the implementation of a null file (file
C   containing no records) is file-system (operating system) dependent.
3       IF (NEXT.GT.1)BACKSPACE(K)
        RETURN
1     CONTINUE
4     RETURN
      END

      SUBROUTINE GETVISCDFLT(VHB_VAL,VHB,NE)
      DIMENSION VHB(8,NE)
      VHB_VAL=VHB(1,1)
      END

      SUBROUTINE GETREALPTR(CXREF,CXARRAY,NSIZE)
      INTEGER NSIZE
      REAL, DIMENSION(NSIZE), TARGET :: CXARRAY
      REAL, DIMENSION(:), POINTER :: CXREF
      CXREF => CXARRAY
      END SUBROUTINE

      SUBROUTINE GETINTPTR(CIREF,CIARRAY,NSIZE)
      INTEGER NSIZE
      INTEGER, DIMENSION(NSIZE), TARGET :: CIARRAY
      INTEGER, DIMENSION(:), POINTER :: CIREF
      CIREF => CIARRAY
      END SUBROUTINE

      SUBROUTINE GETDBLPTR(CXREF,CXARRAY,NSIZE)
      INTEGER NSIZE
      DOUBLE PRECISION, DIMENSION(NSIZE), TARGET :: CXARRAY
      DOUBLE PRECISION, DIMENSION(:), POINTER :: CXREF
      CXREF => CXARRAY
      END SUBROUTINE
