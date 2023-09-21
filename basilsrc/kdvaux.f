C*--------------------------------------------------------------------
C*    Basil / Sybil:   kdvaux.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE CATSTR(WORD1,WORD2,WORD3)
C
C    This subroutine removes the blank spaces
C    before concatanating two strings
C
      CHARACTER WORD1*60,WORD2*60,WORD3*60
      LW1=LEN(WORD1)
      LW2=LEN(WORD2)
      J1=0
      J2=0
      DO 20010 I=1,LW1
       IF(WORD1(I:I).EQ.' ') GO TO 20020
       J1=J1+1
20010 CONTINUE
20020 DO 20030 I=1,LW2
       IF(WORD2(I:I).EQ.' ') GO TO 20040
       J2=J2+1
20030 CONTINUE
20040 WORD3=WORD1(1:J1)//WORD2(1:J2)
      RETURN
      END
C
      SUBROUTINE INDEND(QBND,IBCTYP,NBP2P,ITYPE)
C
C     Routine to zero the boundary conditions:
C       zero velocity if ITYPE=0
C       zero traction if ITYPE=1
C
      DIMENSION QBND(NBP2P),IBCTYP(NBP2P)
C
C   The following loop sets all boundary velocity components to zero.
C
      IF(ITYPE.EQ.0)THEN
        DO 10 I=1,NBP2P
          IF(IBCTYP(I).EQ.0)QBND(I)=0.0
C         IBCTYP(I)=0
   10   CONTINUE
C
C   The following loop sets all boundary traction components to zero.
C
      ELSE IF(ITYPE.EQ.1)THEN
        DO 20 I=1,NBP2P
          IF(IBCTYP(I).EQ.1)QBND(I)=0.0
C         IBCTYP(I)=1
   20   CONTINUE
      END IF
C
C     Output the boundary condition matrices (for debugging)
C
C     WRITE(LUW,*)' Boundary condition array IBCTYP(NBP,2)'
C     CALL MITPRT(IBCTYP,NBP2P/2,2,NBP2P,LUW)
C     WRITE(LUW,*)' Boundary condition array QBND(NBP,2)'
C     CALL MATPRT(QBND,NBP2P/2,2,NBP2P,LUW)
C
      RETURN
      END

      LOGICAL FUNCTION INARRAY(IVEC,ISIZE,N)
      DIMENSION IVEC(ISIZE)
      INARRAY = .FALSE.
      DO 10 I=1,ISIZE
        IF (IVEC(I).EQ.N) THEN
          INARRAY = .TRUE.
          GO TO 12
        END IF
 10   CONTINUE
 12   END
C
      SUBROUTINE UPDBNDVEL(QBND,EX,EY,IBC,NOR,NBP,NUP,BIG)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION NOR(NUP)
      DIMENSION IBC(NBP)
      DIMENSION BBOX(2,4)
      INTEGER PNTONSEGMENT

      EPS = 1E-5
C       find the box dimensions
      CALL FINDBBOX(EX,EY,BBOX,NUP)
C     XMIN = EX(1)
C     XMAX = EX(1)
C     YMIN = EY(1)
C     YMAX = EY(1)
C     DO 10 I=2,NUP
C       IF (EX(I).GT.XMAX) XMAX = EX(I)
C       IF (EX(I).LT.XMIN) XMIN = EX(I)
C       IF (EY(I).GT.YMAX) YMAX = EY(I)
C       IF (EY(I).LT.YMIN) YMIN = EY(I)
C 10  CONTINUE
C       find the velocity in the y-direction (set in input)
      YMIN = BBOX(2,1)
      YMAX = BBOX(2,3)
      XMIN = BBOX(1,1)
      XMAX = BBOX(1,2)
      DO 15 I=1,NBP
        NODE=IBC(I)
        NNODE=NOR(NODE)
        XN=EX(NNODE)
        YN=EY(NNODE)
        IF(ABS(YN-YMIN).LT.EPS) THEN
          VVAL = QBND(I+NBP)
          GOTO 100
        ENDIF
  15  CONTINUE
C       calculate the velocity in the x-direction
 100  UVAL = VVAL * (XMAX-XMIN)/(YMAX-YMIN)
C     print *,uval,XMAX,XMIN,YMAX,YMIN
C       set the boundary velocities for side boundaries
      DO 20 I=1,NBP
        NODE=IBC(I)
        NNODE=NOR(NODE)
        XN=EX(NNODE)
        YN=EY(NNODE)
C       IF(ABS(XN-XMIN).LT.EPS) THEN
        IF (PNTONSEGMENT(XN,YN,BBOX(1,1),BBOX(2,1),
     :                          BBOX(1,4),BBOX(2,4)).NE.0) THEN
          QBND(I)=-UVAL
C       ELSE IF(ABS(XN-XMAX).LT.EPS) THEN
        ELSE IF (PNTONSEGMENT(XN,YN,BBOX(1,2),BBOX(2,2),
     :                          BBOX(1,3),BBOX(2,3)).NE.0) THEN
          QBND(I)=UVAL
C     ELSE IF(ABS(YN-YMAX).GT.EPS.AND.ABS(YN-YMIN).GT.EPS) THEN
C     print *,NNODE
        ENDIF
  20  CONTINUE
      RETURN
      END

      SUBROUTINE INITVELOCITY(QBND,UVP,EX,EY,IBC,NOR,NBP,NUP,
     :                        BIG,IDEFTYP)
C
C    SUBROUTINE TO BE DELETED IF NO LONGER NEEDED ??
C    this routine initialises velocity array, purpose and method
C    need explanation here, and possibly redesign.
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION UVP(NUP*2)
      DIMENSION NOR(NUP)
      DIMENSION IBC(NBP)
      DIMENSION BBOX(2,4)
C
      EPS = 1E-6
      IRECT = 1
      U1 = 0
      U2 = 0
      U3 = 0
      U4 = 0
      V1 = 0
      V2 = 0
      V3 = 0
      V4 = 0
C       find the box dimensions
      CALL FINDBBOX(EX,EY,BBOX,NUP)
      IF (IDEFTYP.NE.0) THEN ! basil1.2.7
        IRECT = 0
C       U1 is velocity at point (XMIN,YMIN) then anticlockwise for 2,3,4
        DO 10 I=1,NBP
          NODE=IBC(I)
          NNODE=NOR(NODE)
          XN=EX(NNODE)
          YN=EY(NNODE)
          IF(ABS(YN-BBOX(2,1)).LT.EPS.AND.ABS(XN-BBOX(1,1)).LT.EPS) THEN
            U1 = QBND(I)
            V1 = QBND(I+NBP)
            IRECT = 1
          ENDIF
          IF(ABS(YN-BBOX(2,2)).LT.EPS.AND.ABS(BBOX(1,2)-XN).LT.EPS) THEN
            U2 = QBND(I)
            V2 = QBND(I+NBP)
            IRECT = 1
          ENDIF
          IF(ABS(BBOX(2,3)-YN).LT.EPS.AND.ABS(BBOX(1,3)-XN).LT.EPS) THEN
            U3 = QBND(I)
            V3 = QBND(I+NBP)
            IRECT = 1
          ENDIF
          IF(ABS(BBOX(2,4)-YN).LT.EPS.AND.ABS(XN-BBOX(1,4)).LT.EPS) THEN
            U4 = QBND(I)
            V4 = QBND(I+NBP)
            IRECT = 1
          ENDIF
  10    CONTINUE
      ENDIF
  15  IF ((IDEFTYP.NE.0).OR.(IRECT.EQ.0).OR.
     :      (U1.EQ.0.AND.U2.EQ.0.AND.U3.EQ.0.AND.U4.EQ.0)) THEN
C         arbitrary velocity field to start (to help VISK)
C       DO 20 I=1,NUP
C         UVP(I)=1.0/float(I)
C  20   CONTINUE
      ELSE
C         calculate the node velocities
        XLEN = BBOX(1,2)-BBOX(1,1)
        YLEN = BBOX(2,4)-BBOX(2,1)
        XLENxYLEN = XLEN * YLEN
        DO 25 I=1,NUP
          NNODE=NOR(I)
          XN=EX(NNODE)-BBOX(1,1)
          YN=EY(NNODE)-BBOX(2,1)
          UVP(I) = (XN*YN*U3 + XN*(YLEN-YN)*U2 +
     :           YN*(XLEN-XN)*U4 +(XLEN-XN)*(YLEN-YN)*U1)/XLENxYLEN
          UVP(I+NUP) = (XN*YN*V3 + XN*(YLEN-YN)*V2 +
     :           YN*(XLEN-XN)*V4 +(XLEN-XN)*(YLEN-YN)*V1)/XLENxYLEN
  25    CONTINUE
      ENDIF
      RETURN
      END
C
      SUBROUTINE MATPRT(Z,NX,NY,NK,LUW)
C
C   this routine, and its double precision version (following) is used
C   for debugging to write array contents when needed
C
      DIMENSION Z(NK)
      DIMENSION L(10)
C
C            FIND NUMBER OF BLOCKS
C
      JB=10
      NB=NX/JB+1
      MM=MOD(NX,JB)
      IF (MM.EQ.0) NB=NB-1
C
C            CYCLE OVER BLOCKS
C
      DO I=1,NB
        JA=(I-1)*JB + 1
        JC=I*JB
        IF (I.EQ.NB) JC=NX
        JJ=JC-JA+1
C
        DO J=1,JJ
          L(J)=JA+J-1
        ENDDO
        WRITE (LUW,101) (L(J),J=1,JJ)
  101   FORMAT (7X,10(' J=',I5,3X))
C
C            CYCLE OVER ROWS
C
        DO K1=1,NY
          WRITE (LUW,102) K1, (Z((K1-1)*NX+J),J=JA,JC)
        ENDDO
  102   FORMAT (1X,'K=',I4,1X,10G14.5)
        WRITE(LUW,103)
      ENDDO       ! on I
  103 FORMAT(/)
C
      RETURN
      END
C
      SUBROUTINE DMATPRT(Z,NX,NY,NK,LUW)
C
C   this routine, and its single precision version (preceding) is used
C   for debugging to write array contents when needed
C
      DOUBLE PRECISION Z
      DIMENSION Z(NK)
      DIMENSION L(10)
C
C            FIND NUMBER OF BLOCKS
C
      JB=10
      NB=NX/JB+1
      MM=MOD(NX,JB)
      IF (MM.EQ.0) NB=NB-1
C
C            CYCLE OVER BLOCKS
C
      DO I=1,NB
        JA=(I-1)*JB + 1
        JC=I*JB
        IF (I.EQ.NB) JC=NX
        JJ=JC-JA+1
C
        DO J=1,JJ
          L(J)=JA+J-1
        ENDDO
        WRITE (LUW,101) (L(J),J=1,JJ)
  101   FORMAT (7X,10(' J=',I5,3X))
C
C            CYCLE OVER ROWS
C
        DO K1=1,NY
          WRITE (LUW,102) K1, (Z((K1-1)*NX+J),J=JA,JC)
        ENDDO
  102   FORMAT (1X,'K=',I4,1X,10G14.5)
        WRITE(LUW,103)
      ENDDO       ! on I
  103 FORMAT(/)
C
      RETURN
      END

      SUBROUTINE SMATPP(ARR,NK)
C
C    subroutine to provide diagnostic data for a single precision array 
C               for use in debugging
C
      DIMENSION ARR(NK)
C
      AMIN=ARR(1)
      JMIN=1
      AMAX=ARR(1)
      JMAX=1
      NZERO=0
      NSUM=0
      ASUM=0.0
      VAR=0.0
      AVE=0.0
C
C    find max, min, average of ARR
C
      DO J=1,NK
        IF(ISNAN(ARR(J)))THEN
          WRITE(*,*)J,'th element of array is a NaN'
          STOP
        ENDIF
        IF(ARR(J).EQ.0)THEN
          NZERO=NZERO+1
        ELSE
          ASUM=ASUM+ARR(J)
          NSUM=NSUM+1
        ENDIF
        IF(ARR(J).GT.AMAX)THEN
          AMAX=ARR(J)
          JMAX=J
        ELSEIF(ARR(J).LT.AMIN)THEN
          AMIN=ARR(J)
          JMIN=J
        ENDIF
      ENDDO
      IF(NSUM.GT.0)AVE=ASUM/FLOAT(NSUM)
C
C    get rms variation
C
      DO J=1,NK
        DEV=ARR(J)-AVE
        DEV=DEV*DEV
        VAR=VAR+DEV
      ENDDO
      ARMS=SQRT(VAR/FLOAT(NK))
C
C    report array facts
C
      WRITE(*,10001)NK,NSUM,AMIN,AMAX,ARMS,AVE
10001 FORMAT('no. of entries = ',I6,'  no. non-zero = ',I6,
     :       '  min = ',G12.5,'  max = ',G12.5,/,
     :       '  rms = ',G12.5,'  average = ',G12.5)
      RETURN
      END

      SUBROUTINE DMATPP(ARR,NK,LUW)
C
C    subroutine to provide diagnostic data for a double precision
C            array for use in debugging
C
      DOUBLE PRECISION ARR,VAR,DEV
      DIMENSION ARR(NK)
C
      AMIN=ARR(1)
      JMIN=1
      AMAX=ARR(1)
      JMAX=1
      NZERO=0
      NSUM=0
      ASUM=0.0
      VAR=0.d0
      AVE=0.d0
C
C    find max, min, average of ARR
C
      DO J=1,NK
        IF(ISNAN(ARR(J)))THEN
          WRITE(*,*)J,'th element of array is a NaN'
          STOP
        ENDIF
        IF(ARR(J).EQ.0)THEN
          NZERO=NZERO+1
        ELSE
          ASUM=ASUM+ARR(J)
          NSUM=NSUM+1
        ENDIF
        IF(ARR(J).GT.AMAX)THEN
          AMAX=ARR(J)
          JMAX=J
        ELSEIF(ARR(J).LT.AMIN)THEN
          AMIN=ARR(J)
          JMIN=J
        ENDIF
      ENDDO
      IF(NSUM.GT.0)AVE=ASUM/FLOAT(NSUM)
C
C    get rms variation
C
      DO J=1,NK
        DEV=ARR(J)-AVE
        DEV=DEV*DEV
        VAR=VAR+DEV
      ENDDO
      ARMS=DSQRT(VAR/DFLOAT(NK))
C
C    report array facts
C
      WRITE(*,10001)NK,NSUM,AMIN,AMAX,ARMS,AVE
10001 FORMAT('no. of entries = ',I6,'  no. non-zero = ',I6,
     :       '  min = ',G12.5,'  max = ',G12.5,/,
     :       '  rms = ',G12.5,'  average = ',G12.5)
      RETURN
      END

      SUBROUTINE MITPP(IARR,NK)
C
C    subroutine to provide diagnostic data for an integer array for use in debugging
C
      DIMENSION IARR(NK)
C
      AMIN=FLOAT(IARR(1))
      JMIN=1
      AMAX=FLOAT(IARR(1))
      JMAX=1
      NZERO=0
      NSUM=0
      ASUM=0.0
      VAR=0.0
      AVE=0.0
C
C    find max, mix, average of ARR
C
      DO J=1,NK
C       IF(ISNAN(ARR(J)))THEN
C         WRITE(*,*)J,'th element of array is a NaN'
C         STOP
C       ENDIF
        IF(IARR(J).EQ.0)THEN
          NZERO=NZERO+1
        ELSE
          ASUM=ASUM+FLOAT(IARR(J))
          NSUM=NSUM+1
        ENDIF
        IF(FLOAT(IARR(J)).GT.AMAX)THEN
          AMAX=FLOAT(IARR(J))
          JMAX=J
        ELSEIF(FLOAT(IARR(J)).LT.AMIN)THEN
          AMIN=FLOAT(IARR(J))
          JMIN=J
        ENDIF
      ENDDO
      IF(NSUM.GT.0)AVE=ASUM/FLOAT(NSUM)
C
C    get rms variation
C
      DO J=1,NK
        DEV=FLOAT(IARR(J))-AVE
        DEV=DEV*DEV
        VAR=VAR+DEV
      ENDDO
      ARMS=SQRT(VAR/FLOAT(NK))
C
C    report array facts
C
      WRITE(*,10001)NK,NSUM,AVE,AMIN,AMAX,ARMS
10001 FORMAT('no. of entries = ',I6,'  no. non-zero = ',I6,
     :       '  average = ',G12.5,'  min = ',G12.5,'  max = ',
     :       G12.5,'  rms = ',G12.5)
      RETURN
      END

      SUBROUTINE  MATREAD(Z,NX,NY,NK,LUW)
C
      DIMENSION Z(NK)
      DIMENSION L(10)
C
C            FIND NUMBER OF BLOCKS
C
      JB=10
      NB=NX/JB+1
      MM=MOD(NX,JB)
      IF (MM.EQ.0) NB=NB-1
C
C            CYCLE OVER BLOCKS
C
      DO I=1,NB
        JA=(I-1)*JB + 1
        JC=I*JB
        IF (I.EQ.NB) JC=NX
        JJ=JC-JA+1
C
        DO J=1,JJ
          L(J)=JA+J-1
        ENDDO
        READ (LUW,101) (L(J),J=1,JJ)
  101   FORMAT (7X,10(3X,I5,3X))
C
C            CYCLE OVER ROWS
C
        DO K1=1,NY
          READ (LUW,102,ERR=104) K2,(Z((K1-1)*NX+J),J=JA,JC)
        ENDDO
  102   FORMAT (3X,I4,1X,10G14.5)
        READ(LUW,103)
      ENDDO         !  on I
  103 FORMAT(/)
      RETURN
  104 WRITE(LSC,*)K2
      STOP
C
      RETURN
      END
      SUBROUTINE  MITPRT(JZ,NX,NY,NK,LUW)
C
      DIMENSION JZ(NK)
      DIMENSION L(10)
C
C            FIND NUMBER OF BLOCKS
C
      JB=10
      NB=NX/JB+1
      MM=MOD(NX,JB)
      IF (MM.EQ.0) NB=NB-1
C
C            CYCLE OVER BLOCKS
C
      DO I=1,NB
        JA=(I-1)*JB + 1
        JC=I*JB
        IF (I.EQ.NB) JC=NX
        JJ=JC-JA+1
C
        DO J=1,JJ
          L(J)=JA+J-1
        ENDDO
        WRITE (LUW,101) (L(J),J=1,JJ)
  101   FORMAT (10X,10(' J=',I5,3X))
C
C            CYCLE OVER ROWS
C
        DO K1=1,NY
          WRITE (LUW,102) K1, (JZ((K1-1)*NX+J),J=JA,JC)
        ENDDO
  102   FORMAT (1X,'K=',I4,1X,10I11)
      WRITE(LUW,103)
      ENDDO          ! on I
  103 FORMAT(/)
C
      RETURN
      END
      SUBROUTINE  MITREAD(JZ,NX,NY,NK,LUW)
C
      DIMENSION JZ(NK)
      DIMENSION L(10)
C
C            FIND NUMBER OF BLOCKS
C
      JB=10
      NB=NX/JB+1
      MM=MOD(NX,JB)
      IF (MM.EQ.0) NB=NB-1
C
C            CYCLE OVER BLOCKS
C
      DO I=1,NB
        JA=(I-1)*JB + 1
        JC=I*JB
        IF (I.EQ.NB) JC=NX
        JJ=JC-JA+1
C
        DO J=1,JJ
          L(J)=JA+J-1
        ENDDO
        READ (LUW,101) (L(J),J=1,JJ)
  101   FORMAT (10X,10(3X,I5,3X))
C
C            CYCLE OVER ROWS
C
        DO K1=1,NY
          READ (LUW,102,ERR=104) K2,(JZ((K1-1)*NX+J),J=JA,JC)
        ENDDO
  102   FORMAT (3X,I4,1X,10I11)
        READ(LUW,103)
      ENDDO
  103 FORMAT(/)
      RETURN
  104 WRITE(LSC,*)K2
      STOP
C
      RETURN
      END

      SUBROUTINE READSTORE(JUNIT,JPOS,INTV,RLV,NAMER,NAMEW,
     :               IRPT,JPEND,NROWS,NUP,NE,NBP,NBP2,NFP,NEL,
     :               NUL,NPM,NSM,NBL,NELLEP,NSEG,
     :               EX,EY,UVP,LEM,NOR,VHB,IMAT,VOLD,
     :               SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :               IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC,IFBC2,
     :               IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :               UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :               IELLE,IPOLYN)
      INCLUDE 'indices.parameters'
      INCLUDE 'limits.parameters'
C     DOUBLE PRECISION UVP
      COMMON/V01S/NMDUM2,NDATE,NMDUM3,COMMEN
      COMMON/AI/LUW,LSC,LBC,LLG
      COMMON/BT/MSINDX,MEASUR(MAXMEAS),MSNODE(MAXMEAS),PTLOC(3,MAXMEAS)
      CHARACTER NMDUM2*16,NMDUM3*16,NDATE*16
      CHARACTER NAMER*32,NAMEW*32
      CHARACTER COMMEN*80,BUFFER*80
      DIMENSION INTV(64), RLV(64)
      DIMENSION UVP(NROWS)
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION LEM(6,NE),NOR(NUP)
      DIMENSION VHB(8,NE)
      DIMENSION VOLD(8,NE)
      DIMENSION IMAT(NE)
      DIMENSION SSQ(NUP),FROT(NUP),TEMPT(NUP)
      DIMENSION DENS(7,NE),IELFIX(NUP)
      DIMENSION QBND(NBP2)
      DIMENSION POLEP(3,MAXPOLE),IBPOLE(NBP2),TAPERF(NBP)
      DIMENSION IBC(NBP),IBNGH(NBP2),IBCTYP(NBP2)
      DIMENSION IFBC(NFP),IFBC2(NFP),IFEQV(NFP)
      DIMENSION JFBC1(NFP),JFBC2(NFP)
      DIMENSION LGEM(6,NEL),EXLG(NUL),EYLG(NUL)
      DIMENSION UXLG(NUL),UYLG(NUL)
      DIMENSION STELPX(NPM,NSM)
      DIMENSION STELPY(NPM,NSM)
      DIMENSION LGIBC(NBL),LGIBCF(NBL)
      DIMENSION IELLE(NELLEP),IPOLYN(NE)
      DIMENSION ISEG(3,NSEG)
C     DIMENSION DEFV(5), BCV(5), RHOG(5), VISP(5), DUMMY(14)
C
      IF (IRPT.EQ.0) THEN
        DO I=1,80
          BUFFER(I:I) = ' '
        ENDDO
        BUFFER(1:8) = 'Reading '
        ISTART = 9
        LEN = ISTART
      END IF
C
c      read common from logical unit j
c
      IF(JPOS.EQ.0) REWIND JUNIT
C
C    memory has been already been allocated for the mesh 
C    the following values are set from current mesh specifcations for
C    comparison with mesh from the input file
C
      NXP = INTV(INX)
      NYP = INTV(INY)
      IMSHP = INTV(IIMSH)
      INFLT = INTV(IIFLT)
      NSMP = INTV(INSM)
      NPMP = INTV(INPM)
      NRMP = INTV(INRM)
      NULP = INTV(INUL)
C
C    If eof then file is empty or last record has been read
C    Let the calling routine determine whether it is an error
C
C    The variables in the first record of the file are ignored. They will be reset
C
      READ(JUNIT,END=6,ERR=5)
      READ(JUNIT,END=5,ERR=5)(INTV(K),K=1,64)
C
C    check array sizes compatible with current program before reading
C    arrays
C
      IF((NXP.NE.INTV(INX)).OR.(NYP.NE.INTV(INY)).OR.
     :          (IMSHP.NE.INTV(IIMSH)).OR.(INFLT.NE.INTV(IIFLT))) THEN
      WRITE(LUW,*)'Solution and program array sizes are not compatible'
      WRITE(LUW,*)'NXP=',NXP,' NYP=',NYP,' IMSHP=',IMSHP,' INFLT=',INFLT
      WRITE(LUW,*)'NX=',INTV(INX),' NY=',INTV(INY),' IMSH=',INTV(IIMSH),
     :            ' IFLT=',INTV(IIFLT)
      WRITE(LSC,*)'Solution and program array sizes are not compatible'
      WRITE(LSC,*)'NXP=',NXP,' NYP=',NYP,' IMSHP=',IMSHP,' INFLT=',INFLT
      WRITE(LSC,*)'NX=',INTV(INX),' NY=',INTV(INY),' IMSH=',INTV(IIMSH),
     :            ' IFLT=',INTV(IIFLT)
      GO TO 7
      END IF
C
C    basil does not overwrite the value of BIG
C    it is read and used by sybil
C
      TMPBIG = RLV(IBIG)
      READ(JUNIT,END=5,ERR=5)(RLV(K),K=1,64)
C
C    check value of BIG compatible with current program before reading
C    arrays
C
      IF (TMPBIG.NE.RLV(IBIG)) THEN
      WRITE(LUW,*)'Solution and program values for BIG differ'
      WRITE(LUW,*)'Solution BIG=',TMPBIG,' BIG=',RLV(IBIG)
      GO TO 7
      END IF
      NX1=INTV(INX)+1
      NY1=INTV(INY)+1
C
C     (b) mandatory mesh / solution / boundary condition blocks
C
      READ(JUNIT,END=5,ERR=5)((LEM(K,N),K=1,6),N=1,INTV(INE)),
     :                       (NOR(K),K=1,INTV(INUP)),
     :                       (IBC(K),K=1,INTV(INBP)),
     :                       (IBNGH(K),K=1,INTV(INBP2)),
     :                       (IBCTYP(K),K=1,INTV(INBP2))
      READ(JUNIT,END=5,ERR=5)(EX(K),K=1,INTV(INUP)),
     :                       (EY(K),K=1,INTV(INUP)),
     :                       (UVP(K),K=1,INTV(INROWS)),
     :                       (QBND(K),K=1,INTV(INBP2))
C
C     (c) optional crustal thickness / rotation block
C
      IF(INTV(IICR).NE.0)THEN
C
C    three separate cases for consistency with old solutions
C
        IF (INTV(IICR).EQ.1) THEN
          READ(JUNIT,END=5,ERR=5)(IELFIX(N),N=1,INTV(INUP))
          READ(JUNIT,END=5,ERR=5)(SSQ(K),K=1,INTV(INUP)),
     :                           (FROT(K),K=1,INTV(INUP))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+19) = 'thickness rotation '
            LEN = LEN+19
          END IF
        ELSE IF (INTV(IICR).EQ.2) THEN
          READ(JUNIT,END=5,ERR=5)(IELFIX(N),N=1,INTV(INUP))
          READ(JUNIT,END=5,ERR=5)(SSQ(K),K=1,INTV(INUP))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+10) = 'thickness '
            LEN = LEN+10
          END IF
        ELSE IF (INTV(IICR).EQ.3) THEN
          READ(JUNIT,END=5,ERR=5)(FROT(K),K=1,INTV(INUP))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+9) = 'rotation '
            LEN = LEN+9
          END IF
        ENDIF
      END IF

C     (d) optional viscosity arrays (if IVIS .ne. 0)
C
      IF(INTV(IIVIS).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)((VHB(K,N),K=1,8),N=1,INTV(INE))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+10) = 'viscosity '
          LEN = LEN+10
        END IF
      END IF
C
C     (e) optional density distribution arrays (if IDEN,ITEMP .ne. 0)
C
      IF(INTV(IIDEN).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)((DENS(K,N),K=1,7),N=1,INTV(INE))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+8) = 'density '
          LEN = LEN+8
        END IF
      END IF
      IF(INTV(IITEMP).NE.0)THEN
         READ(JUNIT,END=5,ERR=5)(TEMPT(N),N=1,INTV(INUP))
         IF (IRPT.EQ.0) THEN
           BUFFER(LEN:LEN+12) = 'temperature '
           LEN = LEN+12
         END IF
       END IF
C
C     (f) optional Lagrangian arrays (if ILAG .ne. 0)
C
      IF(INTV(IILAG).NE.0)THEN
C
C    check array sizes compatible with current program before reading
C    arrays
C
        IF((NSMP.NE.INTV(INSM)).OR.(NPMP.NE.INTV(INPM)).OR.
     :         (NULP.NE.INTV(INUL)).OR.(NRMP.NE.INTV(INRM))) THEN
          WRITE(LUW,*)'Solution and program lgmesh sizes',
     :                '  are not compatible'
          WRITE(LUW,*)'NULP=',NULP,' NSMP=',NSMP,' NPMP=',NPMP,
     :                ' NRMP=',NRMP
          WRITE(LUW,*)'NUL=',INTV(INUL),' NSM=',INTV(INSM),
     :                ' NPM=',INTV(INPM),' NRM=',INTV(INRM)
          GO TO 7
        END IF
C
C    three separate cases for consistency with old solutions
C
        IF (INTV(IILAG).EQ.1) THEN
          READ(JUNIT,END=5,ERR=5)((LGEM(K,N),K=1,6),N=1,INTV(INEL)),
     :                           (LGIBC(K),K=1,INTV(INBL)),
     :                           (LGIBCF(K),K=1,INTV(INBL))
          READ(JUNIT,END=5,ERR=5)(EXLG(K),K=1,INTV(INUL)),
     :                           (EYLG(K),K=1,INTV(INUL)),
     :                           (UXLG(K),K=1,INTV(INUL)),
     :                           (UYLG(K),K=1,INTV(INUL)),
     :                   ((STELPX(K,N),K=1,INTV(INPM)),N=1,INTV(INSM)),
     :                   ((STELPY(K,N),K=1,INTV(INPM)),N=1,INTV(INSM))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+14) = 'Lgmesh marker '
            LEN = LEN+14
          END IF
        ELSE IF (INTV(IILAG).EQ.2) THEN
          READ(JUNIT,END=5,ERR=5)
     :               ((STELPX(K,N),K=1,INTV(INPM)),N=1,INTV(INSM)),
     :               ((STELPY(K,N),K=1,INTV(INPM)),N=1,INTV(INSM))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+7) = 'marker '
            LEN = LEN+7
          END IF
        ELSE IF (INTV(IILAG).EQ.3) THEN
          READ(JUNIT,END=5,ERR=5)((LGEM(K,N),K=1,6),N=1,INTV(INEL)),
     :                           (LGIBC(K),K=1,INTV(INBL)),
     :                           (LGIBCF(K),K=1,INTV(INBL))
          READ(JUNIT,END=5,ERR=5)(EXLG(K),K=1,INTV(INUL)),
     :                           (EYLG(K),K=1,INTV(INUL)),
     :                           (UXLG(K),K=1,INTV(INUL)),
     :                           (UYLG(K),K=1,INTV(INUL))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+7) = 'Lgmesh '
            LEN = LEN+7
          END IF
        END IF
      END IF
C
C     (g) optional fault arrays block (if IFLT. ne. 0)
C
      IF(INTV(IIFLT).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)(IFBC(N),N=1,INTV(INFP)),
     :                         (IFBC2(N),N=1,INTV(INFP)),
     1                         (IFEQV(K),K=1,INTV(INFP)),
     :                         (JFBC1(K),K=1,INTV(INFP)),
     2                         (JFBC2(K),K=1,INTV(INFP))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+6) = 'fault '
          LEN = LEN+6
        END IF
      END IF
C
C     (h) optional series data arrays block (if MSINDX. ne. 0)
C
      IF(INTV(IMSINDX).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)(MEASUR(N),N=1,INTV(IMSINDX)),
     1                         (MSNODE(N),N=1,INTV(IMSINDX))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+6) = 'series '
          LEN = LEN+7
        END IF
      END IF
C
C     (i) optional elle node number array
C
      IF(INTV(INELLEP).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)(IELLE(N),N=1,INTV(INELLEP))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+4) = 'elle '
          LEN = LEN+5
        END IF
      END IF
C
C     (j) optional polygon numbers for triangles
C
      IF(INTV(IIPOLY).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)(IPOLYN(N),N=1,INTV(INE))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+7) = 'polygon '
          LEN = LEN+8
        END IF
      END IF
C
C     (k) optional starting viscosity array (if IIVOLD .ne. 0)
C
      IF(INTV(IIVOLD).NE.0)THEN
         READ(JUNIT,END=5,ERR=5)((VOLD(K,N),K=1,8),N=1,INTV(INE))
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+8) = 'viscold '
            LEN = LEN+8
         END IF
      END IF
C
C     (l) optional material property array (if IIMREG .ne. 0)
C
      IF(INTV(IIMREG).NE.0)THEN
         READ(JUNIT,END=5,ERR=5)(IMAT(N),N=1,INTV(INE))
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+9) = 'material '
            LEN = LEN+9
         END IF
      END IF
C
C     (m) optional segment array block (if NSEG . gt. 0)
C
      IF(INTV(INSEG).GT.0)THEN
         READ(JUNIT,END=5,ERR=5)((ISEG(K,N),K=1,3),N=1,INTV(INSEG))
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+9) = 'segment '
            LEN = LEN+9
         END IF
      END IF
C
C     (n) optional IBPOLE/TAPERF array block (if IVRESET .gt. 0)
C
      IF(INTV(IIVRESET).GT.0)THEN
         READ(JUNIT,END=5,ERR=5)(IBPOLE(K),K=1,2*NBP)
         READ(JUNIT,END=5,ERR=5)(TAPERF(K),K=1,NBP)
         IF(INTV(IIPOLE).GT.0)READ(JUNIT,END=5,ERR=5)
     :     ((POLEP(J,K),J=1,3),K=1,INTV(IIPOLE))
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+5) = 'pole '
            LEN = LEN+5
         END IF
      END IF

      IF (IRPT.EQ.0.AND.LEN.GT.ISTART) THEN
        BUFFER(LEN:LEN+6) = 'arrays'
        LEN = LEN+6
        WRITE(LSC,'(A80/)')BUFFER
        WRITE(LUW,'(A80/)')BUFFER
      END IF
C     WRITE(*,*)'QBND array after reading'
C     CALL SMATPP(QBND,NBP*2)
C     WRITE(*,*)'IBCTYP after reading'
C     CALL MITPP(IBCTYP,NBP*2)
      RETURN
  5   WRITE(LUW,*)'Problem encountered with READ operation in STORE'
      WRITE(LSC,*)'Problem encountered with READ operation in STORE'
      JPEND=99
      RETURN
  6   JPEND=98
      RETURN
  7   JPEND=97
      RETURN
      END

      SUBROUTINE WRITESTORE(JUNIT,JPOS,INTV,RLV,NAMER,NAMEW,
     :               IRPT,JPEND,NROWS,NUP,NE,NBP,NBP2,NFP,NEL,
     :               NUL,NPM,NSM,NBL,NELLEP,NSEG,
     :               EX,EY,UVP,LEM,NOR,VHB,THDINT,IMAT,
     :               VOLD,SSQ,FROT,DENS,TEMPT,POLEP,IBPOLE,TAPERF,
     :               IELFIX,QBND,IBC,IBNGH,IBCTYP,IFBC,IFBC2,
     :               IFEQV,JFBC1,JFBC2,ISEG,LGEM,EXLG,EYLG,
     :               UXLG,UYLG,STELPX,STELPY,LGIBC,LGIBCF,
     :               IELLE,IPOLYN)

C
C   integer variables passed in as some compilers won't allow
C   an array index in a dimension statement
C
      INCLUDE 'indices.parameters'
      INCLUDE 'limits.parameters'
C     DOUBLE PRECISION UVP
      COMMON/V01S/NMDUM2,NDATE,NMDUM3,COMMEN
      COMMON/AI/LUW,LSC,LBC,LLG
      COMMON/BT/MSINDX,MEASUR(MAXMEAS),MSNODE(MAXMEAS),PTLOC(3,MAXMEAS)
      CHARACTER NMDUM2*16,NMDUM3*16,NDATE*16
      CHARACTER NAMER*32,NAMEW*32
      CHARACTER COMMEN*80,BUFFER*80
      DIMENSION INTV(64), RLV(64)
      DIMENSION UVP(NROWS)
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION LEM(6,NE),NOR(NUP)
      DIMENSION VHB(8,NE)
      DIMENSION THDINT(7,NE)
      DIMENSION VOLD(8,NE)
      DIMENSION IMAT(NE)
      DIMENSION SSQ(NUP),FROT(NUP)
      DIMENSION DENS(7,NE),IELFIX(NUP)
      DIMENSION QBND(NBP2)
      DIMENSION POLEP(3,MAXPOLE),IBPOLE(NBP2),TAPERF(NBP)
      DIMENSION IBC(NBP),IBNGH(NBP2),IBCTYP(NBP2)
      DIMENSION IFBC(NFP),IFBC2(NFP),IFEQV(NFP)
      DIMENSION JFBC1(NFP),JFBC2(NFP)
      DIMENSION LGEM(6,NEL),EXLG(NUL),EYLG(NUL)
      DIMENSION UXLG(NUL),UYLG(NUL)
      DIMENSION STELPX(NPM,NSM)
      DIMENSION STELPY(NPM,NSM)
      DIMENSION LGIBC(NBL),LGIBCF(NBL)
      DIMENSION IELLE(NELLEP),IPOLYN(NE)
      DIMENSION ISEG(3,NSEG)
C
      IF (IRPT.EQ.0) THEN
        DO I=1,80
          BUFFER(I:I) = ' '
        ENDDO
        BUFFER(1:7) = 'Saving '
        ISTART = 8
        LEN = ISTART
      END IF
C
c      write common from logical unit j
c
      NSPACE = 0
      IF(JPOS.EQ.0) REWIND JUNIT
C
C     only 16 bytes of names saved in order not to break backward compatibility
C
      WRITE(JUNIT,ERR=5)NAMEW(1:16),NDATE,NAMER(1:16),COMMEN
C
C     every record begins with the header variables
C
      WRITE(JUNIT,ERR=5)(INTV(K),K=1,64)
      WRITE(JUNIT,ERR=5)(RLV(K),K=1,64)
C
C     (b) mandatory mesh / solution / boundary condition blocks
C
      WRITE(JUNIT,ERR=5)((LEM(K,N),K=1,6),N=1,INTV(INE)),
     :                  (NOR(K),K=1,INTV(INUP)),
     :                  (IBC(K),K=1,INTV(INBP)),
     :                  (IBNGH(K),K=1,INTV(INBP2)),
     :                  (IBCTYP(K),K=1,INTV(INBP2))
      WRITE(JUNIT,ERR=5)(EX(K),K=1,INTV(INUP)),
     :                  (EY(K),K=1,INTV(INUP)),
     :                  (UVP(K),K=1,INTV(INROWS)),
     :                  (QBND(K),K=1,INTV(INBP2))
C
C     (c) optional crustal thickness / rotation block
C
      IF(INTV(IICR).NE.0)THEN
C
C    three separate cases for consistency with old solutions
C
        IF (INTV(IICR).EQ.1) THEN
          WRITE(JUNIT,ERR=5)(IELFIX(N),N=1,INTV(INUP))
          WRITE(JUNIT,ERR=5)(SSQ(K),K=1,INTV(INUP)),
     :                      (FROT(K),K=1,INTV(INUP))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+19) = 'thickness rotation '
            LEN = LEN+19
          END IF
        ELSE IF (INTV(IICR).EQ.2) THEN
          WRITE(JUNIT,ERR=5)(IELFIX(N),N=1,INTV(INUP))
          WRITE(JUNIT,ERR=5)(SSQ(K),K=1,INTV(INUP))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+10) = 'thickness '
            LEN = LEN+10
          END IF
        ELSE IF (INTV(IICR).EQ.3) THEN
          WRITE(JUNIT,ERR=5)(FROT(K),K=1,INTV(INUP))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+9) = 'rotation '
            LEN = LEN+9
          END IF
        END IF
      END IF

C     (d) optional viscosity arrays (if IVIS .ne. 0)
C
      IF(INTV(IIVIS).NE.0)THEN
        IF(INTV(IITHDI).NE.0)THEN
          GAMMAT=RLV(IGAMMA)
          TREF=RLV(ITREF)
C
C    array manipulation here to output the current effective strength coeficients
C    while preserving original VHB array and avoiding declaring yet another array
C
          IF(TREF.EQ.0.0)THEN
            DO N=1,INTV(INE)
              DO K=1,7
                VHB(K,N)=VHB(K,N)*EXP(-GAMMAT*THDINT(K,N))
              ENDDO
            ENDDO
          ELSE      !  if TREF.ne.0
            DO N=1,INTV(INE)
              DO K=1,7
                WRKTRM=GAMMAT*THDINT(K,N)
                BLOG=LOG(VHB(K,N))
                DNUMER=BLOG+(TREF*BLOG-1.0)*WRKTRM
                DENOMR=1.0+TREF*WRKTRM
                VHB(K,N)=EXP(DNUMER/DENOMR)
              ENDDO
            ENDDO
          ENDIF
          WRITE(JUNIT,ERR=5)((VHB(K,N),K=1,8),N=1,INTV(INE))
          IF(TREF.EQ.0.0)THEN
            DO N=1,INTV(INE)
              DO K=1,7
                VHB(K,N)=VHB(K,N)*EXP(GAMMAT*THDINT(K,N))
              ENDDO
            ENDDO
          ELSE      !  if TREF.ne.0
            DO N=1,INTV(INE)
              DO K=1,7
                WRKTRM=GAMMAT*THDINT(K,N)
                BLOG=LOG(VHB(K,N))
                DNUMER=WRKTRM+(TREF*WRKTRM+1.0)*BLOG
                DENOMR=1.0+TREF*WRKTRM
                VHB(K,N)=EXP(DNUMER/DENOMR)
              ENDDO
            ENDDO
          ENDIF
C
        ELSE        !  if IITHDI==0
          WRITE(JUNIT,ERR=5)((VHB(K,N),K=1,8),N=1,INTV(INE))
        END IF      !  IITHDI
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+10) = 'viscosity '
          LEN = LEN+10
        END IF
      END IF        !  IIVIS
C
C     (e) optional density distribution arrays (if IDEN,ITEMP .ne. 0)
C
      IF(INTV(IIDEN).NE.0)THEN
        WRITE(JUNIT,ERR=5)((DENS(K,N),K=1,7),N=1,INTV(INE))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+8) = 'density '
          LEN = LEN+8
         END IF
       END IF
       IF(INTV(IITEMP).NE.0)THEN
         WRITE(JUNIT,ERR=5)(TEMPT(N),N=1,INTV(INUP))
         IF (IRPT.EQ.0) THEN
           BUFFER(LEN:LEN+12) = 'temperature '
           LEN = LEN+12
         END IF
       END IF
C
C     (f) optional Lagrangian arrays (if ILAG .ne. 0)
C
      IF(INTV(IILAG).NE.0)THEN
C
C    three separate cases for consistency with old solutions
C
        IF (INTV(IILAG).EQ.1) THEN
          WRITE(JUNIT,ERR=5)((LGEM(K,N),K=1,6),N=1,INTV(INEL)),
     :                      (LGIBC(K),K=1,INTV(INBL)),
     :                      (LGIBCF(K),K=1,INTV(INBL))
          WRITE(JUNIT,ERR=5)(EXLG(K),K=1,INTV(INUL)),
     :                      (EYLG(K),K=1,INTV(INUL)),
     :                      (UXLG(K),K=1,INTV(INUL)),
     :                      (UYLG(K),K=1,INTV(INUL)),
     :                   ((STELPX(K,N),K=1,INTV(INPM)),N=1,INTV(INSM)),
     :                   ((STELPY(K,N),K=1,INTV(INPM)),N=1,INTV(INSM))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+14) = 'Lgmesh marker '
            LEN = LEN+14
          END IF
        ELSE IF (INTV(IILAG).EQ.2) THEN
          WRITE(JUNIT,ERR=5)
     :         ((STELPX(K,N),K=1,INTV(INPM)),N=1,INTV(INSM)),
     :         ((STELPY(K,N),K=1,INTV(INPM)),N=1,INTV(INSM))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+7) = 'marker '
            LEN = LEN+7
          END IF
        ELSE IF (INTV(IILAG).EQ.3) THEN
          WRITE(JUNIT,ERR=5)((LGEM(K,N),K=1,6),N=1,INTV(INEL)),
     :                      (LGIBC(K),K=1,INTV(INBL)),
     :                      (LGIBCF(K),K=1,INTV(INBL))
          WRITE(JUNIT,ERR=5)(EXLG(K),K=1,INTV(INUL)),
     :                      (EYLG(K),K=1,INTV(INUL)),
     :                      (UXLG(K),K=1,INTV(INUL)),
     :                      (UYLG(K),K=1,INTV(INUL))
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+7) = 'Lgmesh '
            LEN = LEN+7
          END IF
        END IF
      END IF
C
C     (g) optional fault arrays block (if IFLT. ne. 0)
C
      IF(INTV(IIFLT).NE.0)THEN
        WRITE(JUNIT,ERR=5)(IFBC(N),N=1,INTV(INFP)),
     :                    (IFBC2(N),N=1,INTV(INFP)),
     :                    (IFEQV(K),K=1,INTV(INFP)),
     :                    (JFBC1(K),K=1,INTV(INFP)),
     :                    (JFBC2(K),K=1,INTV(INFP))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+6) = 'fault '
          LEN = LEN+6
        END IF
      END IF
C
C     (h) optional series data arrays block (if MSINDX. ne. 0)
C
      IF(INTV(IMSINDX).NE.0)THEN
        WRITE(JUNIT,ERR=5)(MEASUR(N),N=1,INTV(IMSINDX)),
     :                    (MSNODE(N),N=1,INTV(IMSINDX))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+6) = 'series '
          LEN = LEN+7
        END IF
      END IF
C
C     (i) optional node numbers for Elle files
C
      IF(INTV(INELLEP).NE.0)THEN
        WRITE(JUNIT,ERR=5)(IELLE(N),N=1,INTV(INELLEP))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+4) = 'elle '
          LEN = LEN+5
        END IF
      END IF
C
C     (j) optional polygon numbers for triangles
C
      IF(INTV(IIPOLY).NE.0)THEN
        WRITE(JUNIT,ERR=5)(IPOLYN(N),N=1,INTV(INE))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+7) = 'polygon '
          LEN = LEN+8
        END IF
      END IF
C
C     (k) optional starting viscosity array (if IIVOLD .ne. 0)
C
      IF(INTV(IIVOLD).NE.0)THEN
         WRITE(JUNIT,ERR=5)((VOLD(K,N),K=1,8),N=1,INTV(INE))
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+8) = 'viscold '
            LEN = LEN+8
         END IF
      END IF
C
C     (l) optional material property array (if IIMREG .ne. 0)
C
      IF(INTV(IIMREG).NE.0)THEN
         WRITE(JUNIT,ERR=5)(IMAT(N),N=1,INTV(INE))
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+9) = 'material '
            LEN = LEN+9
         END IF
      END IF
C
C     (m) optional segment array block (if NSEG . gt. 0)
C
      IF(INTV(INSEG).GT.0)THEN
         WRITE(JUNIT,ERR=5)((ISEG(K,N),K=1,3),N=1,INTV(INSEG))
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+8) = 'segment '
            LEN = LEN+8
         END IF
      END IF
C
C     (n) optional IBPOLE/TAPERF array block (if IVRESET .gt. 0)
C
      IF(INTV(IIVRESET).GT.0)THEN
         WRITE(JUNIT,ERR=5)(IBPOLE(K),K=1,2*INTV(INBP))
         WRITE(JUNIT,ERR=5)(TAPERF(K),K=1,INTV(INBP))
         IF(INTV(IIPOLE).GT.0)WRITE(JUNIT,ERR=5)
     :     ((POLEP(J,K),J=1,3),K=1,INTV(IIPOLE))
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+5) = 'pole '
            LEN = LEN+5
         END IF
      END IF
C
      IF(IRPT.EQ.0)THEN
        BUFFER(LEN:LEN+6) = 'arrays'
        LEN = LEN+6
        WRITE(LSC,'(A80)')BUFFER(1:80)
        WRITE(LUW,'(A80)')BUFFER(1:80)
      END IF
      RETURN
c
  5   WRITE(LUW,*)'Problem encountered with WRITE operation in STORE'
      WRITE(LSC,*)'Problem encountered with WRITE operation in STORE'
      JPEND=99
      RETURN
      END
C
C******************************************************************************
C
      SUBROUTINE WRITELGMESH(ILAG,LGEM,EXLG,EYLG,
     :                          LGIBC,LGIBCF,
     :                          STELPX,STELPY,RAD,
     :                          NEL,NUL,NBL,NSM,NPM,NRM,
     :                          LLG,LSC,LUW,IERROR)
      DIMENSION STELPX(NPM*NSM),STELPY(NPM*NSM)
      DIMENSION EXLG(NUL),EYLG(NUL)
      DIMENSION LGIBC(NBL),LGIBCF(NBL)
      DIMENSION LGEM(NEL)
C
      WRITE(LLG,*)NUL,NEL,NBL,NSM,NPM,NRM,RAD
      CALL MITPRT(LGEM,NEL,6,6*NEL,LLG)
      CALL MITPRT(LGIBC,NBL,1,NBL,LLG)
      CALL MITPRT(LGIBCF,NBL,1,NBL,LLG)
      CALL MATPRT(EXLG,NUL,1,NUL,LLG)
      CALL MATPRT(EYLG,NUL,1,NUL,LLG)
      IF (ILAG.EQ.1.OR.ILAG.EQ.2) THEN
        CALL MATPRT(STELPX,NPM*NSM,1,NPM*NSM,LLG)
        CALL MATPRT(STELPY,NPM*NSM,1,NPM*NSM,LLG)
      ENDIF
C
      RETURN
      END
C
C******************************************************************************
C
      SUBROUTINE READLGMESH(ILAG,LGEM,EXLG,EYLG,
     :                          LGIBC,LGIBCF,
     :                          STELPX,STELPY,RAD,
     :                          NEL,NUL,NBL,NSM,NPM,NRM,
     :                          LLG,LSC,LUW,IERROR)
      DIMENSION STELPX(NPM*NSM),STELPY(NPM*NSM)
      DIMENSION EXLG(NUL),EYLG(NUL)
      DIMENSION LGIBC(NBL),LGIBCF(NBL)
      DIMENSION LGEM(NEL)
C
      READ(LLG,*)NUL,NEL,NBL,NSM,NPM,NRM,RAD
      CALL MITREAD(LGEM,NEL,6,6*NEL,LLG)
      CALL MITREAD(LGIBC,NBL,1,NBL,LLG)
      CALL MITREAD(LGIBCF,NBL,1,NBL,LLG)
      CALL MATREAD(EXLG,NUL,1,NUL,LLG)
      CALL MATREAD(EYLG,NUL,1,NUL,LLG)
      IF (ILAG.EQ.1.OR.ILAG.EQ.2) THEN
        CALL MATREAD(STELPX,NPM*NSM,1,NPM*NSM,LLG)
        CALL MATREAD(STELPY,NPM*NSM,1,NPM*NSM,LLG)
      ENDIF
C
      RETURN
      END

      SUBROUTINE SKIPSTORE(JUNIT,JPOS,INTV,RLV,IRPT,
     :                     MAXSIZE,MAXNBS,JPEND)
      INCLUDE 'indices.parameters'
      INCLUDE 'limits.parameters'
      COMMON/AI/LUW,LSC,LBC,LLG
      CHARACTER BUFFER*80
      DIMENSION INTV(64), RLV(64)
      DIMENSION ITMP(MAXSIZE)
      DIMENSION LEMTMP(MAXSIZE)
C
      IF (IRPT.EQ.0) THEN
        DO I=1,80
          BUFFER(I:I) = ' '
        ENDDO
        BUFFER(1:8) = 'Reading '
        ISTART = 9
        LEN = ISTART
      END IF
      NXP = INTV(INX)
      NYP = INTV(INY)
      IMSHP = INTV(IIMSH)
      INFLT = INTV(IIFLT)
      NSMP = INTV(INSM)
      NPMP = INTV(INPM)
      NRMP = INTV(INRM)
      NULP = INTV(INUL)
C
C    If eof then file is empty or last record has been read
C    Let the calling routine determine whether it is an error
C
C    The first comments, IntegerVars and RealVars have been read
C
C     (b) mandatory mesh / solution / boundary condition blocks
C
      READ(JUNIT,END=5,ERR=5)(LEMTMP(K),K=1,INTV(INE)*6),
     :                       (ITMP(K),K=1,INTV(INUP)),
     :                       (ITMP(K),K=1,INTV(INBP)),
     :                       (ITMP(K),K=1,INTV(INBP2)),
     :                       (ITMP(K),K=1,INTV(INBP2))
      CALL FINDMAXNBS(LEMTMP,INTV(INE),INTV(INN),MAXNBS,IERROR)
      READ(JUNIT,END=5,ERR=5)
C
C     (c) optional crustal thickness / rotation block
C
      IF(INTV(IICR).NE.0)THEN
C
C    three separate cases for consistency with old solutions
C
        IF (INTV(IICR).EQ.1) THEN
          READ(JUNIT,END=5,ERR=5)
          READ(JUNIT,END=5,ERR=5)
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+19) = 'thickness rotation '
            LEN = LEN+19
          END IF
        ELSE IF (INTV(IICR).EQ.2) THEN
          READ(JUNIT,END=5,ERR=5)
          READ(JUNIT,END=5,ERR=5)
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+10) = 'thickness '
            LEN = LEN+10
          END IF
        ELSE IF (INTV(IICR).EQ.3) THEN
          READ(JUNIT,END=5,ERR=5)
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+9) = 'rotation '
            LEN = LEN+9
          END IF
        ENDIF
      END IF

C     (d) optional viscosity arrays (if IVIS .ne. 0)
C
      IF(INTV(IIVIS).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+10) = 'viscosity '
          LEN = LEN+10
        END IF
      END IF
C
C     (e) optional density distribution arrays (if IDEN,ITEMP .ne. 0)
C
      IF(INTV(IIDEN).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+8) = 'density '
          LEN = LEN+8
        END IF
      END IF
      IF(INTV(IITEMP).NE.0)THEN
         READ(JUNIT,END=5,ERR=5)
         IF (IRPT.EQ.0) THEN
           BUFFER(LEN:LEN+12) = 'temperature '
           LEN = LEN+12
         END IF
       END IF
C
C     (f) optional Lagrangian arrays (if ILAG .ne. 0)
C
      IF(INTV(IILAG).NE.0)THEN
C
C    check array sizes compatible with current program before reading
C    arrays
C
        IF((NSMP.NE.INTV(INSM)).OR.(NPMP.NE.INTV(INPM)).OR.
     :         (NULP.NE.INTV(INUL)).OR.(NRMP.NE.INTV(INRM))) THEN
          WRITE(LUW,*)'Solution and program lgmesh sizes',
     :                '  are not compatible'
          WRITE(LUW,*)'NULP=',NULP,' NSMP=',NSMP,' NPMP=',NPMP,
     :                ' NRMP=',NRMP
          WRITE(LUW,*)'NUL=',INTV(INUL),' NSM=',INTV(INSM),
     :                ' NPM=',INTV(INPM),' NRM=',INTV(INRM)
          GO TO 7
        END IF
C
C    three separate cases for consistency with old solutions
C
        IF (INTV(IILAG).EQ.1) THEN
          READ(JUNIT,END=5,ERR=5)
          READ(JUNIT,END=5,ERR=5)
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+14) = 'Lgmesh marker '
            LEN = LEN+14
          END IF
        ELSE IF (INTV(IILAG).EQ.2) THEN
          READ(JUNIT,END=5,ERR=5)
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+7) = 'marker '
            LEN = LEN+7
          END IF
        ELSE IF (INTV(IILAG).EQ.3) THEN
          READ(JUNIT,END=5,ERR=5)
          READ(JUNIT,END=5,ERR=5)
          IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+7) = 'Lgmesh '
            LEN = LEN+7
          END IF
        END IF
      END IF
C
C     (g) optional fault arrays block (if IFLT. ne. 0)
C
      IF(INTV(IIFLT).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+6) = 'fault '
          LEN = LEN+6
        END IF
      END IF
C
C     (h) optional series data arrays block (if MSINDX. ne. 0)
C
      IF(INTV(IMSINDX).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)
     1                         (ITMP(N),N=1,INTV(IMSINDX))
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+6) = 'series '
          LEN = LEN+7
        END IF
      END IF
C
C     (i) optional elle node number array
C
      IF(INTV(INELLEP).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+4) = 'elle '
          LEN = LEN+5
        END IF
      END IF
C
C     (j) optional polygon numbers for triangles
C
      IF(INTV(IIPOLY).NE.0)THEN
        READ(JUNIT,END=5,ERR=5)
        IF (IRPT.EQ.0) THEN
          BUFFER(LEN:LEN+7) = 'polygon '
          LEN = LEN+8
        END IF
      END IF
C
C     (k) optional starting viscosity array (if IIVOLD .ne. 0)
C
      IF(INTV(IIVOLD).NE.0)THEN
         READ(JUNIT,END=5,ERR=5)
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+8) = 'viscold '
            LEN = LEN+8
         END IF
      END IF
C
C     (l) optional material property array (if IIMREG .ne. 0)
C
      IF(INTV(IIMREG).NE.0)THEN
         READ(JUNIT,END=5,ERR=5)
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+9) = 'material '
            LEN = LEN+9
         END IF
      END IF
C
C     (m) optional segment array block (if NSEG . gt. 0)
C
      IF(INTV(INSEG).GT.0)THEN
         READ(JUNIT,END=5,ERR=5)
         IF (IRPT.EQ.0) THEN
            BUFFER(LEN:LEN+9) = 'segment '
            LEN = LEN+9
         END IF
      END IF

      IF (IRPT.EQ.0.AND.LEN.GT.ISTART) THEN
        BUFFER(LEN:LEN+6) = 'arrays'
        LEN = LEN+6
        WRITE(LSC,'(A80/)')BUFFER
        WRITE(LUW,'(A80/)')BUFFER
      END IF
      RETURN
  5   WRITE(LUW,*)'Problem encountered with SKIP operation in STORE'
      WRITE(LSC,*)'Problem encountered with SKIP operation in STORE'
      JPEND=99
      RETURN
  6   JPEND=98
      RETURN
  7   JPEND=97
      RETURN
      END

