C*--------------------------------------------------------------------
C*    Basil / Sybil:   r8read.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      PROGRAM r8tor4
C   this Program is used to convert a REAL*8 solution from 
C   basil to a REAL*4 solution for use by sybil 
C   This file must be compiled with the -r8 option and the
C   file containing the r4write procedure must be compiled
C   without -r8. (Records written by unformatted WRITE statements
C   begin and end with an integer specifying the record length.
C   If compiled with -r8, READ and WRITE assume an 8-byte integer.
C   If compiled without -r8, READ and WRITE assume an 4-byte integer.)
C 
C   The cg.parameters files should match that of the solution to be
C   converted
C
      INCLUDE 'cg.parameters' 
C
C    the parameters below match the appropriate index into the
C    array of 64 integers written by WRITESTORE (kdvaux.f)
C
      PARAMETER (INDXNX=1,INDXNY=2,INDXMSH=3)
      PARAMETER (INDXCR=5,INDXVIS=6,INDXDEN=7,INDXLAG=8,INDXFLT=9)
      PARAMETER (INDXNUP=14,INDXNFP=37)

      COMMON/STR/ NAMER,NAMEW,NDATE,COMMEN

      COMMON/REALS/ RVARS(64)

      COMMON/RARRAY/ EX(NUPP),EY(NUPP),UVP(NROWSP),
     4 QBND(NBP2P),SSQ(NUPP),
     5 FROT(NUPP),VHB(8,NEP),DENS(7,NEP),EXLG(NULP),EYLG(NULP),
     6 STELPX(NPMP,NSMP),STELPY(NPMP,NSMP),UXLG(NULP),UYLG(NULP)

      COMMON/INTS/ IVARS(64)

      COMMON/IARRAY/ NOR(NUPP),IBC(NBPP),IBNGH(NBP2P),
     5 IBCTYP(NBP2P),LEM(6,NEP),
     6 IELFIX(NUPP),LGEM(6,NELP),LGIBC(NBLP),LGIBCF(NBLP),
     7 IFBC(NFPP),IFBC2(NFPP),IFEQV(NFPP),JFBC1(NFPP),JFBC2(NFPP)

      CHARACTER*12 FNAME,GNAME
      CHARACTER*16 NAMER,NAMEW,NDATE
      CHARACTER*80 COMMEN
      REAL*4 RVARS,
     4 EX(NUPP),EY(NUPP),UVP(NROWSP),QBND(NBP2P),SSQ(NUPP),
     5 FROT(NUPP),VHB(8,NEP),DENS(7,NEP),EXLG(NULP),EYLG(NULP),
     6 STELPX(NPMP,NSMP),STELPY(NPMP,NSMP),UXLG(NULP),UYLG(NULP)
      INTEGER*4 IVARS,
     5 NOR(NUPP),IBC(NBPP),IBNGH(NBP2P),IBCTYP(NBP2P),LEM(6,NEP),
     6 IELFIX(NUPP),LGEM(6,NELP),LGIBC(NBLP),LGIBCF(NBLP),
     7 IFBC(NFPP),IFBC2(NFPP),IFEQV(NFPP),JFBC1(NFPP),JFBC2(NFPP)
C
C   4 byte variables above, 8 byte variables below
C
      REAL RVARS8(64),
     4 EX8(NUPP),EY8(NUPP),UVP8(NROWSP),QBND8(NBP2P),SSQ8(NUPP),
     5 FROT8(NUPP),VHB8(8,NEP),DENS8(7,NEP),EXLG8(NULP),EYLG8(NULP),
     6 STELPX8(NPMP,NSMP),STELPY8(NPMP,NSMP),UXLG8(NULP),UYLG8(NULP)
      INTEGER IVARS8(64),
     5 NOR8(NUPP),IBC8(NBPP),IBNGH8(NBP2P),IBCTYP8(NBP2P),LEM8(6,NEP),
     6 IELFIX8(NUPP),LGEM8(6,NELP),LGIBC8(NBLP),LGIBCF8(NBLP),
     7 IFBC8(NFPP),IFBC28(NFPP),IFEQV8(NFPP),JFBC18(NFPP),JFBC28(NFPP)
C
C    Input file names
C
      WRITE(6,*)'Enter name of file to be reformatted from r8 to r4'
      READ(5,*)FNAME
      JUNIT=3
      OPEN(JUNIT,FILE='FD.sols/'//FNAME,STATUS='old',
     1 FORM='unformatted')
      NCOUNT=0
C
C      read 8-byte variables from unit J
C
 1000 NCOUNT=NCOUNT+1
      READ(JUNIT,END=6,ERR=5)NAMEW,NDATE,NAMER,COMMEN
C    KEXIT and BIG are written but not read
      READ(JUNIT,END=5,ERR=5) IVARS8
      READ(JUNIT,END=5,ERR=5) RVARS8
      
      DO 1 I=1,64
      IVARS(I)=IVARS8(I)
  1   RVARS(I)=RVARS8(I)
      
C    check array sizes compatible with current program before reading
C    arrays
C
      IF((NXP.NE.IVARS(INDXNX)).OR.(NYP.NE.IVARS(INDXNY)).OR.
     1  (IMSHP.NE.IVARS(INDXMSH)).OR.(INFLT.NE.IVARS(INDXFLT)))THEN
        WRITE(LUW,*)
     1    'Solution and program array sizes are not compatible'
        WRITE(LUW,*)'NXP=',NXP,' NYP=',NYP,' IMSHP=',IMSHP,
     1                                           ' INFLT=',INFLT
        WRITE(LUW,*)'NX=',IVARS(INDXNX),' NY=',IVARS(INDXNY),
     1            ' IMSH=',IVARS(INDXMSH),' IFLT=',IVARS(INDXFLT)
      STOP
      END IF

C
C     (b) mandatory mesh / solution / boundary condition blocks
C
      READ(JUNIT,END=5,ERR=5)((LEM8(K,N),K=1,6),N=1,NEP),
     1 (NOR8(K),K=1,NUPP),(IBC8(K),K=1,NBPP),(IBNGH8(K),K=1,NBP2P),
     2 (IBCTYP8(K),K=1,NBP2P)
      READ(JUNIT,END=5,ERR=5)(EX8(K),K=1,NUPP),(EY8(K),K=1,NUPP),
     1 (UVP8(K),K=1,NROWSP),(QBND8(K),K=1,NBP2P)
      DO 20 I=1,NUPP
      NOR(I)=NOR8(I)
      EX(I)=EX8(I)
   20 EY(I)=EY8(I)

      DO 30 J=1,NEP
      DO 30 I=1,6
   30 LEM(I,J)=LEM8(I,J)

      DO 40 I=1,NROWSP
   40 UVP(I)=UVP8(I)

      DO 50 I=1,NBP2P
      QBND(I)=QBND8(I)
      IBCTYP(I)=IBCTYP8(I)
   50 IBNGH(I)=IBNGH8(I)

      DO 60 I=1,NBP
   60 IBC(I)=IBC8(I)
C
C     (c) optional crustal thickness / rotation block
C
      IF(IVARS(INDXCR).NE.0)THEN
      READ(JUNIT,END=5,ERR=5)(IELFIX8(N),N=1,NUPP)
      READ(JUNIT,END=5,ERR=5)(SSQ8(K),K=1,NUPP),(FROT8(K),K=1,NUPP)
      DO 70 I=1,NUPP
      SSQ(I)=SSQ8(I)
      FROT(I)=FROT8(I)
   70 IELFIX(I)=IELFIX8(I)
      END IF
C
C     (d) optional viscosity arrays (if IVIS .ne. 0)
C
      IF(IVARS(INDXVIS).NE.0)THEN
      READ(JUNIT,END=5,ERR=5)((VHB8(K,N),K=1,8),N=1,NEP)
      DO 80 J=1,NEP
      DO 80 I=1,8
   80 VHB(I,J)=VHB8(I,J)
      END IF
C
C     (e) optional density distribution arrays (if IDEN .ne. 0)
C
      IF(IVARS(INDXDEN).NE.0)THEN
      READ(JUNIT,END=5,ERR=5)((DENS8(K,N),K=1,7),N=1,NEP)
      DO 90 J=1,NEP
      DO 90 I=1,7
   90 DENS(I,J)=DENS8(I,J)
      END IF
C
C     (f) optional Lagrangian arrays (if ILAG .ne. 0)
C
      IF(IVARS(INDXLAG).NE.0)THEN
      READ(JUNIT,END=5,ERR=5)((LGEM8(K,N),K=1,6),N=1,NELP),
     1 (LGIBC8(K),K=1,NBLP),(LGIBCF8(K),K=1,NBLP)
      READ(JUNIT,END=5,ERR=5)(EXLG8(K),K=1,NULP),(EYLG8(K),K=1,NULP),
     1 (UXLG8(K),K=1,NULP),(UYLG8(K),K=1,NULP),
     2((STELPX8(K,N),K=1,NPMP),N=1,NSMP),
     3((STELPY8(K,N),K=1,NPMP),N=1,NSMP)
      DO 100 J=1,NELP
      DO 100 I=1,6
  100 LGEM(I,J)=LGEM8(I,J)
      DO 110 I=1,NBLP
      LGIBC(I)=LGIBC8(I)
  110 LGIBCF(I)=LGIBCF8(I)
      DO 120 I=1,NULP
      EXLG(I)=EXLG8(I)
      EYLG(I)=EYLG8(I)
      UXLG(I)=EXLG8(I)
  120 UYLG(I)=EYLG8(I)
      DO 130 J=1,NSMP
      DO 130 I=1,NPMP
      STELPX(I,J)=STELPX8(I,J)
  130 STELPY(I,J)=STELPY8(I,J)
      END IF
C
C     (g) optional fault arrays block (if IFLT. ne. 0)
C
      IF(IVARS(INDXFLT).NE.0)THEN
      READ(JUNIT,END=5,ERR=5)(IFBC8(N),N=1,NFPP),(IFBC28(N),N=1,NFPP),
     1 (IFEQV8(K),K=1,NFPP),(JFBC18(K),K=1,NFPP),
     2 (JFBC28(K),K=1,NFPP)
      DO 140 I=1,NFPP
      IFBC(I)=IFBC8(I)
      IFBC2(I)=IFBC28(I)
      IFEQV(I)=IFEQV8(I)
      JFBC1(I)=JFBC18(I)
  140 JFBC2(I)=JFBC28(I)
      END IF
      CLOSE(JUNIT)
C
C    Open up output file
C
      IF(NCOUNT.EQ.1)THEN
      WRITE(6,*)'Enter new file name (default is to overwrite old)'
      READ(5,'(A12)')GNAME
      IF(GNAME(1:4).EQ.'    ')GNAME=FNAME
      OPEN(JUNIT,FILE='FD.sols/'//GNAME,FORM='UNFORMATTED')
      END IF

C   write the 4-byte solution
      CALL WRITERECORD(JUNIT)

      WRITE(6,201)NCOUNT,IVARS(INDXNUP),IVARS(INDXNFP)
      IF(NCOUNT.GT.50)STOP
      GO TO 1000
c
    5 WRITE(6,*)'Error encountered on Read/Write operation'
      STOP
    6 WRITE(6,*)'End of file encountered on Read operation'
      STOP
 200  FORMAT('Solution ',I4,' with NX1 =',I4,'  NY1 =', I4,
     1'  IFLT =',I4,'  has been read')
 201  FORMAT('Solution ',I4,' with NUP =',I6,'  NFP =', I4,
     1'  has been written')
      END

