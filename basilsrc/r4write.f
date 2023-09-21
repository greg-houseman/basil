C*--------------------------------------------------------------------
C*    Basil / Sybil:   r4write.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
      SUBROUTINE WRITERECORD(JUNIT)

      INCLUDE 'cg.parameters'
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

      COMMON/IARRAY/ NOR(NUPP),IBC(NBPP),IBNGH(NBP2P),IBCTYP(NBP2P),
     6 LEM(6,NEP),IELFIX(NUPP),LGEM(6,NELP),LGIBC(NBLP),LGIBCF(NBLP),
     7 IFBC(NFPP),IFBC2(NFPP),IFEQV(NFPP),JFBC1(NFPP),JFBC2(NFPP)
      CHARACTER*16 NAMER,NAMEW,NDATE
      CHARACTER*80 COMMEN
C
C     write 4-byte variables to unit I
C
      WRITE(JUNIT,ERR=5)NAMEW,NDATE,NAMER,COMMEN
C    KEXIT and BIG are written
      WRITE(JUNIT,ERR=5)IVARS
      WRITE(JUNIT,ERR=5)RVARS
C
C     (b) mandatory mesh / solution / boundary condition blocks
C
      WRITE(JUNIT,ERR=5)((LEM(K,N),K=1,6),N=1,NEP),
     1 (NOR(K),K=1,NUPP),(IBC(K),K=1,NBPP),(IBNGH(K),K=1,NBP2P),
     2 (IBCTYP(K),K=1,NBP2P)
      WRITE(JUNIT,ERR=5)(EX(K),K=1,NUPP),(EY(K),K=1,NUPP),
     1 (UVP(K),K=1,NROWSP),(QBND(K),K=1,NBP2P)
C
C     (c) optional crustal thickness / rotation block
C
      IF(IVARS(INDXCR).NE.0)THEN
      WRITE(JUNIT,ERR=5)(IELFIX(N),N=1,NUPP)
      WRITE(JUNIT,ERR=5)(SSQ(K),K=1,NUPP),(FROT(K),K=1,NUPP)
      END IF

C     (d) optional viscosity arrays (if IVIS .ne. 0)
C
      IF(IVARS(INDXVIS).NE.0)THEN
      WRITE(JUNIT,ERR=5)((VHB(K,N),K=1,8),N=1,NEP)
      END IF
C
C     (e) optional density distribution arrays (if IDEN .ne. 0)
C
      IF(IVARS(INDXDEN).NE.0)THEN
      WRITE(JUNIT,ERR=5)((DENS(K,N),K=1,7),N=1,NEP)
      END IF
C
C     (f) optional Lagrangian arrays (if ILAG .ne. 0)
C
      IF(IVARS(INDXLAG).NE.0)THEN
      WRITE(JUNIT,ERR=5)((LGEM(K,N),K=1,6),N=1,NELP),
     1 (LGIBC(K),K=1,NBLP),(LGIBCF(K),K=1,NBLP)
      WRITE(JUNIT,ERR=5)(EXLG(K),K=1,NULP),(EYLG(K),K=1,NULP),
     1 (UXLG(K),K=1,NULP),(UYLG(K),K=1,NULP),
     2((STELPX(K,N),K=1,NPMP),N=1,NSMP),
     3((STELPY(K,N),K=1,NPMP),N=1,NSMP)
      END IF
C
C     (g) optional fault arrays block (if IFLT. ne. 0)
C
      IF(IVARS(INDXFLT).NE.0)THEN
      WRITE(JUNIT,ERR=5)(IFBC(N),N=1,NFPP),(IFBC2(N),N=1,NFPP),
     1 (IFEQV(K),K=1,NFPP),(JFBC1(K),K=1,NFPP),
     2 (JFBC2(K),K=1,NFPP)
      END IF
      RETURN
c
    5 WRITE(6,*)'Error encountered on Write operation'
      RETURN
      END
