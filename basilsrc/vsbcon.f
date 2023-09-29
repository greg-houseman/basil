C*--------------------------------------------------------------------
C*    Basil / Sybil:   vsbcon.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE VSBTAB(XLEN,YLEN,EX,EY,LEM,NOR,IBC,IBNGH,
     :      XZERO,YZERO,IBCTYP,KORNER,NUP,NE,NBP,IFLT,LUW,LSC,IERR)
C
C    This routine establishes two tables: IBC(NBPP), which is
C    a list of node numbers on the external boundary, and IBNGH(NBPP,2)
C    which is a set of pointers to node numbers to the adjacent
C    vertex nodes which are also on the boundary. IBNGH(J) is the vertex
C    neighbour of node J, in a clockwise direction, while IBNGH(J+NBP) is
C    the neighbouring vertex in an anticlockwise direction.
C    These tables are primarily used by VSBCON to set boundary condition
C    arrays.
C    The criterion used here to identify a boundary node is that it
C    lies on X=0 or X=XLEN, or on Y=0 or Y=YLEN.
C
C    The array IBCTYP recognises 3 basic types of boundary condition
C    IBCTYP=0 is prescribed velocity, IBCTYP=1 is prescribed traction
C    IBCTYP=2 is prescribed friction in a resistive type boundary cond.
C    IBCTYP of 11 is the default value for a locked fault
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION IBC(NBP),IBNGH(NBP*2),IBCTYP(NBP*2)
      DIMENSION KORNER(10),ICL(3)
C
C     All elements are searched to locate a boundary segment
C     on an external boundary.  
C
      EPS = 1.0E-4
      NBPREG=0
      ICORN=0
      XB2=XZERO+XLEN
      YB2=YZERO+YLEN
C
      DO K=1,3
        ICL(K)=0
      ENDDO
      DO K=1,10
        KORNER(K)=0
      ENDDO
C
C    loop on elements
C
      DO 400 N=1,NE
        ICN=0
C
C    loop on midpoint nodes
C
        DO J1=1,3
          J2=MOD(J1,3)+1
          J4=J2+3
          K1=LEM(J1,N)
          K2=LEM(J2,N)
          K4=LEM(J4,N)
          LK4=NOR(K4)
          XM=EX(LK4)
          YM=EY(LK4)
          XT1=XM-XZERO
          XT2=XM-XB2
          YT1=YM-YZERO
          YT2=YM-YB2
          IF((ABS(XT1).LT.EPS).OR.(ABS(XT2).LT.EPS).OR.
     1       (ABS(YT1).LT.EPS).OR.(ABS(YT2).LT.EPS).OR.
     2       (IFLT.NE.0.AND.ABS(XT2-XLEN).LT.EPS))THEN
C
C     add midpoint node to table
C
            ICN=ICN+1
            ICL(ICN)=J1
            NBPREG=NBPREG+1
            IBC(NBPREG)=K4
            IBNGH(NBPREG)=K1
            IBNGH(NBPREG+NBP)=K2
C
C     check whether vertex nodes already in table
C
            ICK1=0
            ICK2=0
            DO I=1,NBPREG
              IF(IBC(I).EQ.K1)ICK1=I
              IF(IBC(I).EQ.K2)ICK2=I
            ENDDO
C
C     add vertex nodes to table if not already present
C     identify neighbouring vertex nodes on boundary
C
            IF(ICK1.EQ.0)THEN
              NBPREG=NBPREG+1
              IBC(NBPREG)=K1
              IBNGH(NBPREG+NBP)=K2
            ELSE
              IBNGH(ICK1+NBP)=K2
            END IF               ! on ICK1
            IF(ICK2.EQ.0)THEN
              NBPREG=NBPREG+1
              IBC(NBPREG)=K2
              IBNGH(NBPREG)=K1
            ELSE
              IBNGH(ICK2)=K1
            END IF               ! on ICK2
C
C     store node numbers of external corners (a corner is a node shared
C     by two boundary segments)
C
            IF(ICN.EQ.2)THEN
              ICORN=ICORN+1
              IF(ICORN.GT.10)THEN
                WRITE(LUW,*)'Insufficient storage for KORNER in VSBTAB'
                STOP
              END IF
              INDEX=ICL(1)*ICL(2)
              IF(INDEX.EQ.2)JCN=2
              IF(INDEX.EQ.6)JCN=3
              IF(INDEX.EQ.3)JCN=1
              KORNER(ICORN)=LEM(JCN,N)
            END IF      ! on value of ICN
          END IF        ! on values of XT1,XT2,YT1,YT2
        ENDDO           ! on J1
  400 CONTINUE          ! on N
C     WRITE(LUW,*)'Identified corner nodes are: '
C     WRITE(LUW,1010)(KORNER(K),K=1,10)
C1010 FORMAT(10I6)
      IF(NBPREG.NE.NBP)THEN
        WRITE(LUW,10100)NBPREG,NBP
        WRITE(LSC,10100)NBPREG,NBP
10100   FORMAT('VSBTAB: Incorrect no. of entries in IBC:  ',
     :         'NBPREG = ',I8,' NBP = ',I8)
        IERR = 1
      END IF
 
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE VSBTAB2(EX,EY,LEM,NOR,IBC,IBNGH,
     :                   IBCTYP,KORNER,NUP,NE,NBP,NMP,NN)
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION IBC(NBP),IBNGH(2*NBP),IBCTYP(2*NBP)
      DIMENSION KORNER(10)
C    mpe and mpn are local arrays
      INTEGER MPE(NMP,2),MPN(NMP)

      EPS=1.E-6
C
C***********************************************************
C  LOCATE NODES ON BOUNDARY AND STORE NODE NUMBERS IN IBC
C***********************************************************
C
C  Step 1: Record numbers of elements adjacent to each midpoint
C          node
C 
C  Search element mid-points
C
      DO N=1,NMP
        MPN(N)=0
        MPE(N,1)=0
        MPE(N,2)=0
      ENDDO
C
      DO 300 I=1,NE
        DO K=4,6
          ITEMP=NOR(LEM(K,I))-NN
          IF((ITEMP.LE.0).OR.(ITEMP.GT.NMP) ) THEN
            WRITE(6,10500)K,I,ITEMP
10500       FORMAT('VSBTAB2 problem, stopping; K, I, ITEMP =',3I7)
            STOP
          END IF
          IF(MPE(ITEMP,1).EQ.0) THEN
            MPE(ITEMP,1)=I
            MPN(ITEMP)=K
          ELSE
            MPE(ITEMP,2)=I
          END IF        !  on value of MPE
        ENDDO           !  on K
  300 CONTINUE          !  on I
C
C  Midpoints with only one adjacent element are on boundary
C  and element number in first entry is boundary adjacent
C
C  Step 2: Store node numbers in IBC and adjacent vertex
C          nodes in IBNGH
      J=1
      DO 400 I=1,NMP
C
C  Store midpoint node in IBC
C
        IF(MPE(I,2).EQ.0) THEN
      IBC(J)=LEM(MPN(I),MPE(I,1))
      IF(MPN(I).EQ.4) THEN
        ICVN=3
      ELSE
        ICVN=MPN(I)-4
      ENDIF
          IAVN=MPN(I)-3
C
C  Store adjacent vertex nodes in IBNGH
C
      IBNGH(J)=LEM(ICVN,MPE(I,1))
      IBNGH(J+NBP)=LEM(IAVN,MPE(I,1))
          J=J+1
C
C  Store clockwise vertex node in IBC if not already stored
C  AND anti-clockwise adjacent vertex
C
          NIA=0
          DO 401 K=1,J-1
            IF (LEM(ICVN,MPE(I,1)).EQ.IBC(K)) THEN
C           Node is already in array
            IBNGH(K+NBP)=LEM(IAVN,MPE(I,1))
            NIA=1
            ENDIF
  401     CONTINUE

          IF(NIA.EQ.0) THEN
            IBC(J)=LEM(ICVN,MPE(I,1))
            IBNGH(J+NBP)=LEM(IAVN,MPE(I,1))
            J=J+1
          ENDIF
C
C  Store anti-clockwise vertex node in IBC if not already stored
C  AND clockwise adjacent vertex
C
          NIA=0
          DO 402 K=1,J-1
            IF (LEM(IAVN,MPE(I,1)).EQ.IBC(K)) THEN
C           Node already in array
            IBNGH(K)=LEM(ICVN,MPE(I,1))
            NIA=1
            ENDIF
  402     CONTINUE

          IF(NIA.EQ.0) THEN
            IBC(J)=LEM(IAVN,MPE(I,1))
            IBNGH(J)=LEM(ICVN,MPE(I,1))
            J=J+1
      END IF

        ENDIF
  400 CONTINUE
      IF(NBP.NE.J-1) THEN
      WRITE(6,*) 'ERROR: NBP =',NBP,' .NE. boundary nodes found =',J-1
C     STOP
      ENDIF
C
C***************************************************
C  FINDING CORNER NODES
C***************************************************
C
C  Step 1: Find max and min x values
C
      XMAX=0
      XMIN=1
      DO 600 I=1,NBP
        IF(EX(NOR(IBC(I))).GT.XMAX) THEN
          XMAX=EX(NOR(IBC(I)))
        END IF
        IF(EX(NOR(IBC(I))).LT.XMIN) THEN
          XMIN=EX(NOR(IBC(I)))
        END IF
  600 CONTINUE
C
C  Find min and max y for both min and max x
C
      YMINM=1
      YMAXM=0
      YMINX=1
      YMAXX=0
      DO 601 I=1,NBP
        IF((EX(NOR(IBC(I)))-XMIN.LT.EPS).AND.
     :    (EY(NOR(IBC(I)))).LT.YMINM) THEN
          YMINM=EY(NOR(IBC(I)))
          KORNER(1)=IBC(I)
        END IF
        IF((EX(NOR(IBC(I)))-XMIN.LT.EPS).AND.
     :    (EY(NOR(IBC(I)))).GT.YMAXM) THEN
          YMAXM=EY(NOR(IBC(I)))
          KORNER(4)=IBC(I)
        END IF
        IF((XMAX-EX(NOR(IBC(I))).LT.EPS).AND.
     :    (EY(NOR(IBC(I)))).GT.YMAXX) THEN
          YMAXX=EY(NOR(IBC(I)))
          KORNER(3)=IBC(I)
        END IF
        IF((XMAX-EX(NOR(IBC(I))).LT.EPS).AND.
     :    (EY(NOR(IBC(I)))).LT.YMINX) THEN
          YMINX=EY(NOR(IBC(I)))
          KORNER(2)=IBC(I)
        END IF
  601 CONTINUE
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE VSBCAPPLY(XLEN,YLEN,BIG,IFLTTIPS,EX,EY,QBND,
     :                  POLEP,IPOLE,VSE,ARGAN,HLENSC,
     :                  CENTLNG,RPT2,NCOMP,NOR,LEM,IBC,IBNGH,IBCTYP,
     :                  IDEFTYP,ISEG,NSEG,NUP,NE,NN,NBP,IFLT,NFP,
     :                  NFPF3,IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :                  LBC,LUW,LSC,IDBUG,IERR)
C
C    this routine is called to apply the boundary conditions using
C    keywords ON, BETWEEN, POLE, CREF, PB
C
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*2)
      CHARACTER FRM160*6
      CHARACTER IXY*3,IYX*1,IUT*3
      DOUBLE PRECISION ROTIX(3,3)
      DOUBLE PRECISION PLON8,PROT8
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION IBC(NBP),IBNGH(NBP*2),IBCTYP(NBP*2)
      DIMENSION IFBC1(NFP),IFBC2(NFP),IFEQV(NFP),JFBC1(NFP),JFBC2(NFP)
      DIMENSION IFLTTIPS(MAXTIPS*2)
      DIMENSION ISEG(3,NSEG),POLEP(3,MAXPOLE)
      DIMENSION RPT2(4)
C
      XMID=CENTLNG
      LENGTH=LINELEN*2
      WRITE(FRM160,10160)LENGTH
10160 FORMAT('(A',I3.3,')')
C
C   The following line adjusting the boundary tractions relative to Ar 
C   was used in an earlier series of numerical experiments.
C   number.  This may be an unhelpful complication.  If the velocities
C   are too large or too small, adjusting the background viscosity
C   coefficient may be the best way to bring the solution back into
C   values that are of order 1, without entailing unnecessary scaling
C   complexity.
C
C     CT0=ARGAN/(HLENSC**(2))
C     WRITE(LSC,10033)CT0
C     WRITE(LUW,10033)CT0
C0033 FORMAT('Please note that boundary tractions are scaled by ',
C    :       "Ar'.(L/D)^2 = ",G12.5)
      CT0=1.0
C
C     Zero the QBND array prior to adding in boundary condition constraints
C
      NBP2=2*NBP
      DO J=1,NBP2
        QBND(J)=0.0
      ENDDO
C
      IPOLE=0
      IF (NSEG.GT.0)
     :  CALL STORETIPS(ISEG,NSEG,IFLTTIPS,IBCTYP,IBC,NBP,NOR,NUP)
C
      CALL FINDMLIM(EX,EY,NUP,XMIN,XMAX,YMIN,YMAX)
      XPT0=0.5*(XMIN+XMAX)
      YPT0=0.5*(YMIN+YMAX)
C
C    if spherical shell coordinate rotation is invoked pass
C    rotation operators to PARSEBE via VSBCBE for poles and cref
C
      IF(IDEFTYP.EQ.111)THEN
        CALL GETPOLE(RPT2,PLON8,PROT8,ROTIX) ! find equatorial pole for rotation
        PTALO=RPT2(1)
        PTALA=RPT2(2)
        CALL ROTATE(ROTIX,PTALO,PTALA,1,0.0) ! forward rotation to get DPHI
        DPHI=PTALO
        XMID=DPHI
      ENDIF
C
C    Read boundary instructions, one at a time
C
      WRITE(LUW,1010)
      WRITE(LSC,1010)
 1010 FORMAT('Boundary Conditions:')
  460 CONTINUE
      READ(LBC,FRM160,END=700)INSTR
      CALL SKIP(INSTR,LENGTH,1,'n',J1)
      IF((J1.EQ.0).OR.(J1.GE.(LINELEN*2)))GO TO 460
C     WRITE(*,*)'VSBCAPPLY: INSTR(1:60)=',INSTR(1:60)
C
C   'ON' statements are processed here
C
      IF(INSTR(J1:J1+1).EQ.'ON')THEN
        CALL VSBCON(INSTR,IPR,XLEN,YLEN,BIG,IFLTTIPS,EX,EY,QBND,
     :              ARGAN,HLENSC,CENTLNG,CT0,IDEFTYP,NOR,LEM,IBC,
     :              IBNGH,IBCTYP,ISEG,NSEG,NUP,NE,NN,NBP,LBC,
     :              LUW,LSC,IDBUG,IERR)
        IF (IERR.NE.0)GO TO 50
C
C   PB periodic boundary KBXY (direction) and PDIST (wavelength)
C
      ELSEIF(INSTR(J1:J1+1).EQ.'PB')THEN
        CALL PARSEON(INSTR,IXY,PDIST,IYX,YXLIM1,YXLIM2,IBND,IUT,
     :               UTVAL1,UTVAL2,ITAP,IPR,XMIN,XMAX,YMIN,YMAX,IERR)
        KPBC=1
        IF(IXY(1:1).EQ.'Y')KPBC=2
        CALL PERIBC(EX,EY,QBND,LEM,NOR,IBC,IBNGH,IBCTYP,
     :              IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :              IFLT,NUP,NE,NBP,NFP,NFPF3,KPBC,
     :              PDIST,0,LUW,LSC,IERR)
        IF (IERR.NE.0)GO TO 50
C
C   'BETWEEN, POLE, CREF' statements are processed by this routine
C
      ELSE
        CALL VSBCBE(INSTR,IPR,BIG,EX,EY,QBND,ARGAN,HLENSC,
     :              CENTLNG,CT0,ROTIX,DPHI,NOR,LEM,IBC,IBNGH,IBCTYP,
     :              NCOMP,NUP,NE,NN,NBP,XPT0,YPT0,POLEP,
     :              IPOLE,IDEFTYP,LBC,LUW,LSC,IDBUG,IERR)
        IF (IERR.NE.0) GO TO 50
      END IF
      GO TO 460
C
C    Finished reading boundary condition instructions
C
  700 CONTINUE
C
C     Check all conditions set; only one warning message output
C     even if multiple nodes unset
C
      CALL RESETTIPS(IFLTTIPS,IBCTYP,IBC,NBP,NOR,NUP)  ! purpose ?
      NUMCNR=0
      IWARNX=0
      IWARNY=0
      DO J1=1,NBP
        J2=J1+NBP
C       IF (IBCTYP(J1).EQ.10) NUMCNR=NUMCNR+1
        IF(IBCTYP(J1).LT.0)THEN
          IWARNX=IWARNX+1
          IF(IWARNX.EQ.1)NJNX=IBC(J1)
          IBCTYP(J1)=1
          QBND(J1)=0.0
C         WRITE(*,*)'IWARNX,J1,IBC,EX,EY =',
C    :      IWARNX,J1,IBC(J1),EX(IBC(J1)),EY(IBC(J1))
        ENDIF
        IF(IBCTYP(J2).LT.0)THEN
          IWARNY=IWARNY+1
          IF(IWARNY.EQ.1)NJNY=IBC(J1)
          IBCTYP(J2)=1
          QBND(J2)=0.0
        ENDIF
      ENDDO
C
C    if spherical shell coordinate rotation is invoked
C    undo the coordinate transformation for reporting
C    ROTIX is modified for inverse rotation here
C
      IF(IDEFTYP.EQ.111)THEN
        CALL MATSET(PLON8,-PROT8,ROTIX)
      ENDIF
      IERR = 0
      IF(IWARNX.NE.0)THEN
        XPTU=EX(NJNX)
        YPTU=EY(NJNX)
        IF(IDEFTYP.GE.110)THEN
          CALL PROJECTDEG(XPTU,YPTU,XMID,YMID,1,NCOMP,IERR)
          IF(IDEFTYP.EQ.111)CALL ROTATE(ROTIX,XPTU,YPTU,1,0.0)
        ENDIF
        WRITE(LSC,10100)'X',IWARNX,NJNX,EX(NJNX),EY(NJNX),XPTU,YPTU
      ENDIF
      IF(IWARNY.NE.0)THEN
        XPTU=EX(NJNY)
        YPTU=EY(NJNY)
        IF(IDEFTYP.GE.110)THEN
          CALL PROJECTDEG(XPTU,YPTU,XMID,YMID,1,NCOMP,IERR)
          IF(IDEFTYP.EQ.111)CALL ROTATE(ROTIX,XPTU,YPTU,1,0.0)
        ENDIF
        WRITE(LSC,10100)'Y',IWARNY,NJNY,EX(NJNY),EY(NJNY),XPTU,YPTU
      ENDIF
10100 FORMAT('Nodes unset, ',A1,': ',I5,' set to zero stress, including'
     : ,' node ',I5,/,' at x = ',F8.3,' y = ',F8.3,' long = ',F8.3,
     :  ' lat = ',F8.3)
C
      IF (NUMCNR.GT.4) 
     :  WRITE(LSC,*)'Found extra corners ',NUMCNR
      WRITE(LUW,10103)
      WRITE(LSC,10103)
10103 FORMAT(/)
C
C     Output the boundary condition matrices (if debugging only)
C
   50 CONTINUE
      IF(IERR.NE.0)RETURN
      IDOUT=11
      IF(IDOUT.EQ.0)RETURN
        OPEN(IDOUT,FILE='BCS')
        WRITE(IDOUT,10101)
10101   FORMAT(' ',/,'Array data from VSBCON:',/,
     :  '   J      X        Y      IBC  NGH1  NGH2 TYPX   QBND_X',
     :  '    TYPY    QBND_Y')
        DO 600 J=1,NBP
          QEX=EX(NOR(IBC(J)))
          QEY=EY(NOR(IBC(J)))
          IF(IDEFTYP.GE.110)THEN
            CALL PROJECTDEG(QEX,QEY,XMID,YMID,1,NCOMP,IERR)
          ENDIF
          JP=J+NBP
          UU=QBND(J)
          IF(IBCTYP(J).EQ.0.0)UU=UU
          VV=QBND(JP)
          IF(IBCTYP(JP).EQ.0.0)VV=VV
          WRITE(IDOUT,10102)J,QEX,QEY,IBC(J),
     1               IBNGH(J),IBNGH(JP),IBCTYP(J),UU,IBCTYP(JP),VV
10102     FORMAT(I5,2F9.3,3I6,I4,G14.5,I4,G14.5)
  600   CONTINUE
        CLOSE(IDOUT)
C
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE VELROT(ROTIX,RLON,RLAT,VECE,VECN)
C
C    routine VELROT: reprojects velocity components (VECE,VECN)
C    on rotation of coordinate system define by matrix ROTIX
C    given a rotation by angle PROT about a pole at (0, PLON)
C    this routine follows the algorithm of COROTATE.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*4 RLON,RLAT,VECE,VECN
      DIMENSION ROTIX(3,3),RHS(3),RES(3)
C
C   constant factors
C
      PI=3.14159265358979323846D0
      DEGRAD=PI/180.D0
C
C   evaluate the right hand side using the long/lat of the point
C   and effect the rotation
C
      GLA=RLAT*DEGRAD
      GLO=RLON*DEGRAD
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
C   rotate the vector whose east and north components are (VECE,VECN)
C
      VAMP=SQRT(VECN*VECN+VECE*VECE)
      IF(VAMP.EQ.0.d0)RETURN
C
C   get sines and cosines of lat and long for point G'
C
      SGLATP=RES(3)
      CGLATP=DSQRT(RES(1)*RES(1)+RES(2)*RES(2))
      CGLONP=RES(1)/CGLATP
      SGLONP=RES(2)/CGLATP
C
C   get sines and cosines of lat and long for point H
C
      SLATH=CLATG*VECN/VAMP
      CLATH=DSQRT(1.d0 - SLATH*SLATH)
      CLOND=-(SLATG*SLATH)/(CLATG*CLATH)
      SLOND=VECE/(VAMP*CLATH)
      HLO=GLO+DATAN2(SLOND,CLOND)
      CLONH=DCOS(HLO)
      SLONH=DSIN(HLO)
C
C   compute location of H' by rotation of H
C
      RHS(1)=CLATH*CLONH
      RHS(2)=CLATH*SLONH
      RHS(3)=SLATH
      DO K=1,3
        RES(K)=0.d0
        DO J=1,3
          RES(K)=RES(K)+RHS(J)*ROTIX(K,J)
        ENDDO
      ENDDO
C
C    get sines and cosines of lat and long for point H'
C    and then velocity components in rotated frame
C
      CHLATP=DSQRT(RES(1)*RES(1)+RES(2)*RES(2))
      CHLONP=RES(1)/CHLATP
      SHLONP=RES(2)/CHLATP
      VECN=VAMP*RES(3)/CGLATP
      VECE=VAMP*CHLATP*(SHLONP*CGLONP-CHLONP*SGLONP)
C     WRITE(*,100)'V: ',VAMP,VECE,VECN
C 100 FORMAT(A4,5F12.4)
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE VSBCON(INSTR,IPR,XLEN,YLEN,BIG,IFLTTIPS,EX,EY,QBND,
     :                  ARGAN,HLENSC,CENTLNG,CT0,IDEFTYP,NOR,LEM,IBC,
     :                  IBNGH,IBCTYP,ISEG,NSEG,NUP,NE,NN,NBP,LBC,
     :                  LUW,LSC,IDBUG,IERR)
C
C    This routine tests each point on the IBC table using the boundary
C    condition statements found in the file on unit LBC
C
C    Types of boundary condition:
C        ITYP=0 if prescribed velocity condition is used
C        ITYP=1 if prescribed traction condition is used
C        ITYP=2 if prescribed friction coefficient is used
C        ITYP=3 at point of discontinuity in applied traction
C        ITYP=4 at point of discontinuity in applied friction
C        and for the fault:
C
C         IBCTYP     FAULT  CORNER   EXT BC
C         ------  --------  ------   ------
C	       10       either   yes      vel
C	       11       locked    no       -
C	       12       locked   yes      stress
C	       13     unlocked    no       -
C	       14     unlocked   yes      stress
C
C        Ensure that velocity takes precedence over traction condition
C         at a discontinuity by placing the velocity condition after
C         the traction condition in the bc file.
C        (A discontinuity between applied friction and traction is not
C         rigorously dealt with in present version of program)
C       
C        For the case with a fault, the external boundary must come
C        first in the bc file and then the fault. At the present the 
C        code only considers the cases where the corner nodes (i.e. the
C        nodes shared by the fault and the external boundary) only have
C        an applied velocity or stress on the external boundary (equiv.
C        to ITYPE 0 or 1 only).
C
C     For applied traction boundary condition, QBND holds p.Dl.T
C     where Dl is length of boundary segment, T is traction component
C     (TX or TY); p = 1/6 for vertex nodes, p = 2/3 for midpoint nodes.
C     This is an exact representation of the traction boundary
C     integral if traction varies linearly along element boundary
C     For applied velocity condition, QBND holds velocity component
C
      INCLUDE "limits.parameters"
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION IBC(NBP),IBNGH(NBP*2),IBCTYP(NBP*2)
      CHARACTER IYX*1
      CHARACTER IXY*3
      CHARACTER IUT*3
      CHARACTER INSTR*(LINELEN*2)
      DIMENSION IFLTTIPS(MAXTIPS*2)
      DIMENSION BBOX(2,4)
      DIMENSION ISEG(3,NSEG)
      LOGICAL VALID, VALIDA, VALIDB
      INTEGER PNTONSEGMENT

      EPSR=5.E-6
      EPSR=1.E-5  !basil1.2.7
      PMID=2.0/3.0
      PVERT=1.0/6.0
      IYX=' '
      IXY='   '
      IUT=IXY
      CALL FINDBBOX(EX,EY,BBOX,NUP)
      YMIN = BBOX(2,1)
      YMAX = BBOX(2,4)
      XMIN = BBOX(1,1)
      XMAX = BBOX(1,2)
      IRECT = 0
      IF (BBOX(1,1).GT.BBOX(1,4)-EPSR.AND.BBOX(1,1).LT.BBOX(1,4)+EPSR)
     :  IRECT = 1
      YMID=(YMAX-YMIN)*0.5
C
C   only 'ON' statements are recognised by this routine
C
      CALL PARSEON(INSTR,IXY,XYVAL,IYX,YXLIM1,YXLIM2,IBND,
     1 IUT,UTVAL1,UTVAL2,ITAP,IPR,XMIN,XMAX,YMIN,YMAX,IERR)
      IF (IERR.NE.0) THEN
        WRITE(*,*)'Error detected in parsing ON statement'
        WRITE(*,*)INSTR
        RETURN
      END IF
      IF(IPR.EQ.0)RETURN     ! Not a valid ON statement
C
C     IUT is the parameter to be set, i.e., UX, TX, FY
C     IXY is the key for where it is to be set, e.g. X, Y, BND
C     XYVAL is the value of IXY to activate setting.
C     IYX is the secondary coordinate used to limit application
C     YXLIM1,YXLIM2 are the limits on the secondary coordinate
C     UTVAL1, UTVAL2 are the values to be set
C     ITAP determines the application of a taper
C
C     WRITE(*,*)'VSBCON: ON statement read, UTVAL,IUT,NBP,IXY = ',
C    :          UTVAL,IUT,NBP,IXY
C
      JSET=0
C
C     For each node on the boundary
C
      DO 500 JX=1,NBP
        JY=JX+NBP
        NODE=IBC(JX)
        NODEA=IBNGH(JX)
        NODEB=IBNGH(JY)
        NNODE=NOR(NODE)
        NNODEA=NOR(NODEA)
        NNODEB=NOR(NODEB)
        XN=EX(NNODE)
        YN=EY(NNODE)
        XNA=EX(NNODEA)
        YNA=EY(NNODEA)
        XNB=EX(NNODEB)
        YNB=EY(NNODEB)
        NODSET=0         ! to be changed if BC to be set
C
C     If setting on constant X boundary
C
        IF(IXY(1:1).EQ.'X')THEN
C
C VALIDA and VALIDB used in setting BND below, but why needed ?
C 
          VALID = .FALSE.
          VALIDA = .FALSE.
          VALIDB = .FALSE.
          IF (IRECT.EQ.1) THEN
          VALID= (ABS(XN-XYVAL).LT.EPSR).AND.
     :            (YN.GT.YXLIM1-EPSR).AND.(YN.LT.YXLIM2+EPSR)
          VALIDA= (ABS(XNA-XYVAL).LT.EPSR).AND.
     :            (YNA.GT.YXLIM1-EPSR).AND.(YNA.LT.YXLIM2+EPSR)
          VALIDB= (ABS(XNB-XYVAL).LT.EPSR).AND.
     :            (YNB.GT.YXLIM1-EPSR).AND.(YNB.LT.YXLIM2+EPSR)
C         ELSE   !  when exactly would this ELSE block be used ?
C           IF (XYVAL.EQ.XMIN) THEN
C           VALID = (PNTONSEGMENT(XN,YN,BBOX(1,1),BBOX(2,1),
C    :                           BBOX(1,4),BBOX(2,4)).EQ.1)
C           VALIDA = (PNTONSEGMENT(XNA,YNA,BBOX(1,1),BBOX(2,1),
C    :                           BBOX(1,4),BBOX(2,4)).EQ.1)
C           VALIDB = (PNTONSEGMENT(XNB,YNB,BBOX(1,1),BBOX(2,1),
C    :                           BBOX(1,4),BBOX(2,4)).EQ.1)
C           ELSE
C           VALID = (PNTONSEGMENT(XN,YN,BBOX(1,2),BBOX(2,2),
C    :                           BBOX(1,3),BBOX(2,3)).EQ.1)
C           VALIDA = (PNTONSEGMENT(XNA,YNA,BBOX(1,2),BBOX(2,2),
C    :                           BBOX(1,3),BBOX(2,3)).EQ.1)
C           VALIDB = (PNTONSEGMENT(XNB,YNB,BBOX(1,2),BBOX(2,2),
C    :                           BBOX(1,3),BBOX(2,3)).EQ.1)
C           ENDIF
          ENDIF
          IF(ABS(XN-XYVAL).LT.EPSR)NODSET=1
C
C     If setting on constant Y boundary
C
        ELSE IF(IXY(1:1).EQ.'Y')THEN
           IF(ABS(YN-XYVAL).LT.EPSR)NODSET=1
C
C     If setting on a specified boundary segment number
C
        ELSE IF(IXY(1:1).EQ.'B')THEN
          DO I=1,NSEG
            IF (ISEG(3,I).EQ.IBND)THEN     ! careful don't double-dip on stress
              IF((ISEG(1,I).EQ.NNODE.OR.ISEG(2,I).EQ.NNODE).OR.    !vertex
     :           (ISEG(1,I).EQ.NNODEA.AND.ISEG(2,I).EQ.NNODEB).OR. !midpoint
     :           (ISEG(1,I).EQ.NNODEB.AND.ISEG(2,I).EQ.NNODEA))
     :        NODSET=1
            ENDIF
          ENDDO
        ENDIF
C
C     use secondary coordinate to restrict application
C
        IF(IYX.NE.' ')THEN
          XPTU=XN
          YPTU=YN
          IF(IDEFTYP.GE.110)
     :      CALL PROJECTDEG(XPTU,YPTU,CENTLNG,YMID,1,NCOMP,IERR)
          YXL1=YXLIM1-EPSR   ! to ensure end points captured
          YXL2=YXLIM2+EPSR
          IF(IYX.EQ.'X')THEN
            IF((XPTU.LT.YXL1).OR.(XPTU.GT.YXL2))NODSET=0
            IF((XNA.LT.YXL1).OR.(XNA.GT.YXL2))XNA=XN
            IF((XNB.LT.YXL1).OR.(XNB.GT.YXL2))XNB=XN
          ELSEIF(IYX.EQ.'Y')THEN
            IF((YPTU.LT.YXL1).OR.(YPTU.GT.YXL2))NODSET=0
            IF((YNA.LT.YXL1).OR.(YNA.GT.YXL2))YNA=XN
            IF((YNB.LT.YXL1).OR.(YNB.GT.YXL2))YNB=XN
          ENDIF
        ENDIF
C
C    If NODSET is set, we now set a boundary condition on this node
C    first set conditions that don't depend on segment lengths
C
        IF(NODSET.EQ.1)THEN
          JSET=JSET+1
C
C    Apply a taper if the FOR block is defined and ITAP != 0
C
          UTVAL=UTVAL1
          IF(ITAP.NE.0)THEN
            IF(IYX.EQ.'X')XYN=XN
            IF(IYX.EQ.'Y')XYN=YN
            CALL TAPER(ITAP,YXLIM1,YXLIM2,UTVAL1,
     :                              UTVAL2,XYN,UTVAL)
          ENDIF
C
C    this one doesn't set a boundary condition, but labels the boundary
C    VALIDA and VALIDB defined above - but to what purpose ?
C
          IF(IUT.EQ.'BND')THEN
            DO I=1,NSEG
              IF ((ISEG(1,I).EQ.NNODE).AND.
     :            ((ISEG(2,I).EQ.NNODEA.AND.VALIDA).OR.
     :              (ISEG(2,I).EQ.NNODEB.AND.VALIDB))) THEN
                ISEG(3,I)=INT(UTVAL+EPSR)
              ELSE IF ((ISEG(2,I).EQ.NNODE).AND.
     :            ((ISEG(1,I).EQ.NNODEA.AND.VALIDA).OR.
     :              (ISEG(1,I).EQ.NNODEB.AND.VALIDB))) THEN
                ISEG(3,I)=INT(UTVAL+EPSR)
              ENDIF
            ENDDO
C
C   set fixed displacement rate conditions (UX, UY or UB (both))
C
          ELSE IF(IUT.EQ.'UX')THEN
C           IF (IXY(2:2).EQ.' ') THEN   ! test doesn't work with ON BND ....
              IBCTYP(JX)=0
              QBND(JX)=UTVAL
C           ELSE            !   for a fault surface (+ or -)
C     n.b. if these blocks activated move setting of XNA,XNB below
C             XDIFF=(XNA+XNB)/2.0-XN
C             IF ((IXY(2:2).EQ.'+'.AND.XDIFF.GT.0).OR.
C    1                    (IXY(2:2).EQ.'-'.AND.XDIFF.LT.0)) THEN
C                   IF(IBCTYP(JX).LT.0) IBCTYP(JX)=0
C                   QBND(JX)=UTVAL
C             END IF
C           END IF
          ELSE IF(IUT.EQ.'UY')THEN     ! set Y component of velocity
            IBCTYP(JY)=0
            QBND(JY)=UTVAL
          ELSE IF(IUT.EQ.'UB')THEN     ! set both velocity components to same value
            IBCTYP(JX)=0
            QBND(JX)=UTVAL
            IBCTYP(JY)=0
            QBND(JY)=UTVAL
C
C    in relation to nodes on faults: in general pairs of nodes are listed
C    in sequence in IBC, we set the conditions on both nodes in each pair
C    the same even though we only need one pair.  In any case they should 
C    be the same for continuity of traction to hold.
C    locks a fault (default anyway, but may be useful in some cases)
C
          ELSE IF(IUT.EQ.'FL')THEN
            IBCTYP(JX)=11
            IBCTYP(JY)=11              
            QBND(JX)=0.0
            QBND(JY)=0.0
C
C    unlocks a fault for zero shear stress, continuous nornmal stress
C
          ELSE IF(IUT.EQ.'FU')THEN
            IBCTYP(JX)=11
            IBCTYP(JY)=21              
            QBND(JX)=0.0
            QBND(JY)=0.0
C
C    set friction coefficient on a fault (or unlock it if zero)
C    setting of conditions including segment length dependence 
C    is managed in SEMFLT for FT and FN
C
          ELSE IF(IUT.EQ.'FN')THEN             ! normal
            IBCTYP(JX)=21
            QBND(JX)=UTVAL
          ELSE IF(IUT.EQ.'FT')THEN             ! tangential
            IBCTYP(JY)=21
            QBND(JY)=UTVAL
C
C    all other conditions require segment lengths so these are first determined
C    Setting segment lengths here is not robust to geometry change if time-stepping
C    Solution would be to move actual setting to QLOADBC called by CGRUN
C
          ELSE
            XDIFA=XN-XNA
            XDIFB=XNB-XN
            YDIFA=YN-YNA
            YDIFB=YNB-YN
            DLENA=SQRT(XDIFA*XDIFA+YDIFA*YDIFA)
            DLENB=SQRT(XDIFB*XDIFB+YDIFB*YDIFB)
C
C   the following conditionals seem to deal with corners; can simplify ??
C
            IF(NNODE.LE.NN)THEN                   ! if a vertex node
              PC=PVERT
              XBIF=0.0
              YBIF=0.0
              DLEN=0.0
              IF(IXY(1:1).EQ.'X')THEN
                IF (IRECT.EQ.1) THEN               ! if region is rectangular
                  IF(ABS(XNA-XYVAL).LT.EPSR)THEN
                    DLEN=DLEN+DLENA
                    XBIF=XBIF+XDIFA
                    YBIF=YBIF+YDIFA
                  ENDIF
                  IF(ABS(XNB-XYVAL).LT.EPSR)THEN
                    DLEN=DLEN+DLENB
                    XBIF=XBIF+XDIFB
                    YBIF=YBIF+YDIFB
                  ENDIF
                ELSE IF (XYVAL.EQ.XMIN) THEN        ! for irregular shape
                  IF (PNTONSEGMENT(XNA,YNA,BBOX(1,1),BBOX(2,1),
     :                           BBOX(1,4),BBOX(2,4)).EQ.1)THEN
                    DLEN=DLEN+DLENA
                    XBIF=XBIF+XDIFA
                    YBIF=YBIF+YDIFA
                  ENDIF
                  IF (PNTONSEGMENT(XNB,YNB,BBOX(1,1),BBOX(2,1),
     :                           BBOX(1,4),BBOX(2,4)).EQ.1)THEN
                    DLEN=DLEN+DLENB
                    XBIF=XBIF+XDIFB
                    YBIF=YBIF+YDIFB
                  ENDIF
                ELSE IF (XYVAL.EQ.XMAX) THEN
                  IF (PNTONSEGMENT(XNA,YNA,BBOX(1,2),BBOX(2,2),
     :                           BBOX(1,3),BBOX(2,3)).EQ.1)THEN
                    DLEN=DLEN+DLENA
                    XBIF=XBIF+XDIFA
                    YBIF=YBIF+YDIFA
                  ENDIF
                  IF (PNTONSEGMENT(XNB,YNB,BBOX(1,1),BBOX(2,1),
     :                           BBOX(1,3),BBOX(2,3)).EQ.1)THEN
                    DLEN=DLEN+DLENB
                    XBIF=XBIF+XDIFB
                    YBIF=YBIF+YDIFB
                  ENDIF
                END IF
              ELSE IF(IXY(1:1).EQ.'Y')THEN
                IF(ABS(YNA-XYVAL).LT.EPSR)THEN
                  DLEN=DLEN+DLENA
                  XBIF=XBIF+XDIFA
                  YBIF=YBIF+YDIFA
                ENDIF
                IF(ABS(YNB-XYVAL).LT.EPSR)THEN
                  DLEN=DLEN+DLENB
                  XBIF=XBIF+XDIFB
                  YBIF=YBIF+YDIFB
                ENDIF
              END IF
            ELSE                      ! if a midpoint node (no corners)
              PC=PMID
              XBIF=XNB-XNA
              YBIF=YNB-YNA
              DLEN=SQRT(XBIF*XBIF+YBIF*YBIF)
            ENDIF                     ! setting of segment lengths
C
C    set fixed traction conditions (adding to QBND, so that corners
C    may get contributions from both boundaries)
C
            IF(IUT.EQ.'TX')THEN
              IF(IBCTYP(JX).NE.1)QBND(JX)=0.0
              IBCTYP(JX)=1
              QBND(JX)=QBND(JX)+UTVAL*PC*DLEN
            ELSE IF(IUT.EQ.'TY')THEN
              IF(IBCTYP(JY).NE.1)QBND(JY)=0.0
              IBCTYP(JY)=1
              QBND(JY)=QBND(JY)+UTVAL*PC*DLEN
C
C    for TN: normal traction, magnitude is relative to the body
C    forces as scaled by Argand number (as used in MENDBC)
C
            ELSE IF(IUT.EQ.'TN')THEN
              RMAG=PC*CT0
              IF(IBCTYP(JX).NE.1)QBND(JX)=0.0
              IBCTYP(JX)=1
              QBND(JX)=QBND(JX)+UTVAL*YBIF*RMAG
              IF(IBCTYP(JY).NE.1)QBND(JY)=0.0
              IBCTYP(JY)=1
              QBND(JY)=QBND(JY)-UTVAL*XBIF*RMAG
C
C    for TT: tangential traction, magnitude is relative to the body
C    forces as scaled by Argand number (+ve is anti-clockwise)
C
            ELSE IF(IUT.EQ.'TT')THEN
              RMAG=PC*CT0
              IF(IBCTYP(JX).NE.1)QBND(JX)=0.0
              IBCTYP(JX)=1
              QBND(JX)=QBND(JX)+UTVAL*XBIF*RMAG
              IF(IBCTYP(JY).NE.1)QBND(JY)=0.0
              IBCTYP(JY)=1
              QBND(JY)=QBND(JY)+UTVAL*YBIF*RMAG
C
C    set friction coefficient on external boundary (QBND incremented, as above)
C
            ELSE IF(IUT.EQ.'FX')THEN
              IBCTYP(JX)=2
              QBND(JX)=QBND(JX)+UTVAL*PC*DLEN
            ELSE IF(IUT.EQ.'FY')THEN
              IBCTYP(JY)=2
              QBND(JY)=QBND(JY)+UTVAL*PC*DLEN
            ENDIF            ! end of traction-related BC options
          ENDIF              ! conditions that require segment lengths
        ENDIF                ! BC to be set; NODSET == 1
  500 CONTINUE               ! end of main loop on boundary nodes
C
C    print message to confirm what is set
C
      IF(IYX.EQ.' ')THEN
        IF (IXY.EQ.'BND') THEN
          WRITE(LUW,10006)JSET,IXY,IBND,IUT,UTVAL1
          WRITE(LSC,10006)JSET,IXY,IBND,IUT,UTVAL1
        ELSE
          WRITE(LUW,10001)JSET,IXY,XYVAL,IUT,UTVAL1
          WRITE(LSC,10001)JSET,IXY,XYVAL,IUT,UTVAL1 
        END IF
10001   FORMAT(I4,' nodes set ON ',A3,' = ',F8.3,' : ',
     1   A3,' = ',F8.3)
10006   FORMAT(I4,' nodes set ON ',A3,' = ',I4,' : ',
     1   A3,' = ',F8.3)
      ELSE
        IF (IXY.EQ.'BND') THEN
          WRITE(LUW,10005)JSET,IXY,IBND,YXLIM1,IYX,YXLIM2,IUT,UTVAL1
          WRITE(LSC,10005)JSET,IXY,IBND,YXLIM1,IYX,YXLIM2,IUT,UTVAL1
        ELSE
          WRITE(LUW,10002)JSET,IXY,XYVAL,YXLIM1,IYX,YXLIM2,IUT,UTVAL1
          WRITE(LSC,10002)JSET,IXY,XYVAL,YXLIM1,IYX,YXLIM2,IUT,UTVAL1
        END IF
10002   FORMAT(I4,' nodes set ON ',A3,' = ',F8.3,' For ',
     1   F8.3,' <= ',A1,' <= ',F8.3,'  :  ',A2,' = ',F8.3)
10005   FORMAT(I4,' nodes set ON ',A3,' = ',I4,' For ',
     1   F8.3,' <= ',A1,' <= ',F8.3,'  :  ',A2,' = ',F8.3)
      END IF
      IF(ITAP.NE.0)THEN
        WRITE(LUW,10003)UTVAL2,ITAP
        WRITE(LSC,10003)UTVAL2,ITAP
10003 FORMAT('   tapered to ',F8.3,' using taper function ',I1)
      END IF
C
      RETURN
      END
C
C     store flt tip node number  ! not sure what to do with this block
C      should it be in VSBCON ?  Depends how IFLTTIPS is used
C
C              DO K=1,MAXTIPS
C                IF (IFLTTIPS(K).EQ.0)THEN
C                  IFLTTIPS(K) = NNODE 
C                  GO TO 1010
C                ENDIF
C                    only save one node if internal fault
C                IF (EX(IFLTTIPS(K)).EQ.XN.AND.EY(IFLTTIPS(K)).EQ.YN)
C    :             GO TO 1010
C              ENDDO
C              IERR = 1
C              WRITE(LUW,*)'Too many fault tips - max=',MAXTIPS
C              WRITE(LSC,*)'Too many fault tips - max=',MAXTIPS
C              RETURN
C1010          CONTINUE
C
C****************************************************************************************
C
      SUBROUTINE VSBCBE(INSTR,IPR,BIG,EX,EY,QBND,ARGAN,HLENSC,
     :                  CENTLNG,CT0,ROTIX,DPHI,NOR,LEM,IBC,IBNGH,
     :                  IBCTYP,NCOMP,NUP,NE,NN,NBP,XPT0,YPT0,
     :                  POLEP,IPOLE,IDEFTYP,
     :                  LBC,LUW,LSC,IDBUG,IERR)
C
C    This routine tests each point on the IBC table using the boundary
C    condition BETWEEN statements found in the file on unit LBC
C
C    Types of boundary condition:
C        ITYP=0 if prescribed velocity condition is used
C        ITYP=1 if prescribed traction condition is used
C        ITYP=2 if prescribed friction coefficient is used
C        ITYP=3 at point of discontinuity in applied traction
C        ITYP=4 at point of discontinuity in applied friction
C        and for the fault:
C
C        Ensure that velocity takes precedence over traction condition
C         at a discontinuity by placing the velocity condition after
C         the traction condition in the bc file.
C        (A discontinuity between applied friction and traction is not
C         rigorously dealt with in present version of program)
C       
C        For the case with a fault, the external boundary must come
C        first in the bc file and then the fault. At the present the 
C        code only cosiders the cases where the corner nodes (i.e. the
C        nodes shared by the fault and the external boundary) only have
C        an aplied velocity or stress on the external boundary (equiv.
C        to ITYPE 0 or 1 only).
C
C     For applied traction boundary condition, QBND holds p.Dl.T
C     where Dl is length of boundary segment, T is traction component
C     (TX or TY); p = 1/6 for vertex nodes, p = 2/3 for midpoint nodes.
C     This is an exact representation of the traction boundary
C     integral if traction varies linearly along element boundary
C     For applied velocity condition, QBND holds velocity component
C
      INCLUDE "limits.parameters"
      CHARACTER IYX*1
      CHARACTER IXY*3,IUT*3,COMM*3
      CHARACTER INSTR*(LINELEN*2)
C
      DOUBLE PRECISION ROTIX(3,3)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION IBC(NBP),IBNGH(NBP*2),IBCTYP(NBP*2)
      DIMENSION POLEP(3,MAXPOLE),POLET(3),RPT2(4)
C
      TOOPI=2.0*PI      ! PI and DTOR are in limits.parameters
      EPSR=1.E-6
      PMID=2.0/3.0
      PVERT=1.0/6.0
C
C    Read boundary instructions, one at a time
C
      IERR=0
      CALL PARSEBE(INSTR,COMM,IXY,XYVAL,IYX,YXLIM1,YXLIM2,
     : IUT,UTVAL1,UTVAL2,ITAP,IPR,XMIN,XMAX,YMIN,YMAX,
     : XPT1,YPT1,XPT2,YPT2,XPT0,YPT0,POLEP,ROTIX,DPHI,
     : CENTLNG,IPOLE,IP1,IP2,IDEFTYP,IERR)
      IF (IERR.NE.0) GO TO 50
C
C   only 'BETWEEN', 'POLE', 'CREF' statements are recognised
C   by this routine
C
      IF(IPR.LT.4)RETURN  ! if POLE or CREF statement
C
C    set up the limits on azimuth relative to XPT0, YPT0 (CREF)
C
      XDIF=XPT1-XPT0
      YDIF=YPT1-YPT0
      AZ1=ATAN2(YDIF,XDIF)
      XDIF=XPT2-XPT0
      YDIF=YPT2-YPT0
      AZ2=ATAN2(YDIF,XDIF)
      IF(IPR.EQ.4)UTVAL=UTVAL1
C
C    if using a rotation pole, set pole parameters for fixed velocity
C    will be overwritten inside loop if ITAP != 0
C
      IF(IPR.EQ.5)THEN
        DO K=1,3
          POLET(K)=POLEP(K,IP1)
        ENDDO
      END IF
C
C     For each node on the boundary
C
C     write(*,*)'VSBCBE, NN NBP',NN,NBP
C     WRITE(*,*)'VSBCBE IBNGH in two parts, NBP =',NBP
C     CALL IMATPP(IBNGH,NBP)
C     CALL IMATPP(IBNGH(NBP+1),NBP)
C     CALL MITPRT(IBNGH,NBP,2,2*NBP,6)
      DO 500 JX=1,NBP
        IF((IBCTYP(JX).LT.10).OR.(IBCTYP(JX).GT.500))THEN     !don't hit fault nodes
          JY=JX+NBP
          NODE=IBC(JX)
          NNODE=NOR(NODE)
          NNODE=IBC(JX)
          XN=EX(NNODE)
          YN=EY(NNODE)
          AZ=ATAN2((YN-YPT0),(XN-XPT0))   ! relative to CREF
C
C     Check if boundary node lies within azimuth defined by
C      (XPT1,YPT1) and (XPT2,YPT2) relative to (XPT0,YPT0)
C     This comparison in projected coordinates if applicable
C
          IF(((AZ2.GE.AZ1).AND.(AZ.GE.AZ1).AND.(AZ.LE.AZ2)).OR.
     :       ((AZ2.LT.AZ1).AND.((AZ.GE.AZ1).OR.(AZ.LE.AZ2))))THEN
C
C     For traction conditions get boundary segment length
C     NODEA is next node clockwise, NODEB is next node anticlockwise
C     For problems on the sphere, distances are already projected
C
            IF(IUT(1:1).EQ.'T'.OR.IUT(1:1).EQ.'F')THEN
              NODEA=IBNGH(JX)
              NODEB=IBNGH(JY)
              NNODEA=NOR(NODEA)
              NNODEB=NOR(NODEB)
              XNA=EX(NNODEA)
              YNA=EY(NNODEA)
              XNB=EX(NNODEB)
              YNB=EY(NNODEB)
C
C     doesn't take account of curvature on sphere, but as a means 
C     of segregating groups of points on a boundary is adequate
C
              AZA=ATAN2((YNA-YPT0),(XNA-XPT0))
              AZB=ATAN2((YNB-YPT0),(XNB-XPT0))
C
C     get segment lengths for traction: For midpoint nodes
C
              IF(NNODE.GT.NN)THEN     ! if a vertex node
                PC=PMID
                XBIF=XNB-XNA
                YBIF=YNB-YNA
                DLEN=SQRT(XBIF*XBIF+YBIF*YBIF)
C
C     For vertex nodes; segment lengths in dimensionless distance units
C
              ELSE                    ! if a midpoint node
                PC=PVERT
                XDIFA=XN-XNA
                XDIFB=XNB-XN
                YDIFA=YN-YNA
                YDIFB=YNB-YN
                DLENA=SQRT(XDIFA*XDIFA+YDIFA*YDIFA)
                DLENB=SQRT(XDIFB*XDIFB+YDIFB*YDIFB)
                DLEN=0.0
                XBIF=0.0
                YBIF=0.0
                IF(((AZ2.GE.AZ1).AND.(AZA.GE.AZ1).AND.(AZA.LE.AZ2)).OR.
     :             ((AZ2.LT.AZ1).AND.((AZA.GE.AZ1).OR.(AZA.LE.AZ2))))
     :            THEN
                  DLEN=DLEN+DLENA
                  XBIF=XBIF+XDIFA
                  YBIF=YBIF+YDIFA
                ENDIF
                IF(((AZ2.GE.AZ1).AND.(AZB.GE.AZ1).AND.(AZB.LE.AZ2)).OR.
     :             ((AZ2.LT.AZ1).AND.((AZB.GE.AZ1).OR.(AZB.LE.AZ2))))
     :            THEN
                  DLEN=DLEN+DLENB
                  XBIF=XBIF+XDIFB
                  YBIF=YBIF+YDIFB
                ENDIF                  ! if in azimuth range
              END IF                   ! if a vertex node / midpoint node
            END IF ! if traction or friction condition
C
C     apply the taper based on azimuth of point relative to CREF
C
            CSRAT=1.0
            IF(ITAP.NE.0)THEN
              S1=AZ1
              S2=AZ2
              SZ=AZ
C             IF((AZ1.GE.AZ2).AND.(AZ.LE.AZ2))THEN
C               S1=AZ1-TOOPI
C             ELSE IF((AZ1.GE.AZ2).AND.(AZ.GE.AZ1))THEN
C               S2=AZ2+TOOPI
              IF(AZ1.GT.AZ2)THEN
                S2=AZ2+TOOPI
                IF(AZ.LT.AZ2)SZ=SZ+TOOPI
              ENDIF                     ! crossing azimuth discontinuity
              IF(IPR.EQ.4)THEN          ! Cartesian taper
                CALL TAPER(ITAP,S1,S2,UTVAL1,UTVAL2,SZ,UTVAL)
              ELSE IF(IPR.EQ.5)THEN     ! spherical taper
                IF(IP1.NE.IP2)THEN
                  CALL PTAPER(ITAP,S1,S2,SZ,CSRAT,
     :                  POLEP(1,IP1),POLEP(1,IP2),POLET)
                ENDIF                    ! endpoints differ
C               WRITE(*,*)'vsbcon: S1,S2,SZ,JX,CSRAT=',S1,S2,SZ,JX,CSRAT
              END IF                     ! spherical/cartesian taper
            END IF                       ! if taper applied
C
C   displacement rates derived from a pole of rotation evaluated by
C   VSBSET using Pole numbers and possible taper values saved here
C   POLET(1,2,3) are long, lat and rate for pole
C
            IF(IPR.EQ.5)THEN             ! rotation on sphere
              IBCTYP(JX)=500+IP1
              IBCTYP(JY)=500+IP2
              QBND(JX)=CSRAT
              QBND(JY)=0.0
C
C   set fixed displacement rate conditions
C
            ELSE                         ! cartesian conditions
              IF(IUT.EQ.'UX')THEN
                IBCTYP(JX)=0
                QBND(JX)=UTVAL
              ELSE IF(IUT.EQ.'UY')THEN
                IBCTYP(JY)=0
                QBND(JY)=UTVAL
              ELSE IF(IUT.EQ.'UB')THEN
                IBCTYP(JX)=0
                QBND(JX)=UTVAL
                IBCTYP(JY)=0
                QBND(JY)=UTVAL
C
C    traction boundary conditions are added, assuming QBND was zeroed
C    just prior to reading boundary condition statements (input.f)
C    [but QBND values are re-set now if overwriting a non-traction]
C
C          integral{Pnx.ds} = integral{P(-dy)}
C          integral{Pny.ds} = integral{P dx}
C   for vertices, int{P dx}+ int{P dx'} = int{P (dx+dx')}
C
              ELSE IF(IUT.EQ.'TX')THEN
                RMAG=PC*CT0
                IF(IBCTYP(JX).NE.1)QBND(JX)=0.0
                IBCTYP(JX)=1
                QBND(JX)=QBND(JX)+UTVAL*DLEN*RMAG
              ELSE IF(IUT.EQ.'TY')THEN
                RMAG=PC*CT0
                IF(IBCTYP(JY).NE.1)QBND(JY)=0.0
                IBCTYP(JY)=1
                QBND(JY)=QBND(JY)+UTVAL*DLEN*RMAG
C
C    for TN: normal traction, magnitude is relative to the body
C    forces as scaled by Argand number (as used in MENDBC)
C
              ELSE IF(IUT.EQ.'TN')THEN
                RMAG=PC*CT0
                IF(IBCTYP(JX).NE.1)QBND(JX)=0.0
                IBCTYP(JX)=1
                QBND(JX)=QBND(JX)+UTVAL*YBIF*RMAG
                IF(IBCTYP(JY).NE.1)QBND(JY)=0.0
                IBCTYP(JY)=1
                QBND(JY)=QBND(JY)-UTVAL*XBIF*RMAG
C
C    for TT: tangential traction, magnitude is relative to the body
C    forces as scaled by Argand number (+ve is anti-clockwise)
C
              ELSE IF(IUT.EQ.'TT')THEN
                RMAG=PC*CT0
                IF(IBCTYP(JX).NE.1)QBND(JX)=0.0
                IBCTYP(JX)=1
                QBND(JX)=QBND(JX)+UTVAL*XBIF*RMAG
                IF(IBCTYP(JY).NE.1)QBND(JY)=0.0
                IBCTYP(JY)=1
                QBND(JY)=QBND(JY)+UTVAL*YBIF*RMAG
C
C    set friction coefficient on boundary (QBND incremented, as above)
C
              ELSE IF(IUT.EQ.'FX')THEN
                IBCTYP(JX)=2
                QBND(JX)=QBND(JX)+UTVAL*PC*DLEN
              ELSE IF(IUT.EQ.'FY')THEN
                IBCTYP(JY)=2
                QBND(JY)=QBND(JY)+UTVAL*PC*DLEN
              END IF                 !(if value of IUT: UX, UY, etc)
            END IF                   ! set poles or UX, TX etc
          END IF                     !(if point in right azimuth range)
        ENDIF                        ! if not a fault node
  500 CONTINUE                       ! loop on JX
      RETURN
C
  50  WRITE(*,*)'Error detected in parsing BETWEEN statement'
      WRITE(*,*)INSTR
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE VSBSET(XLOAD,QBND,EX,EY,IBC,IBCTYP,
     :                  POLEP,BIG,CENTLNG,NBP,NROWS,NUP)
C
C    This is applied to reset velocity boundary conditions that are based on
C    a plate tectonic rotation pole.  It uses arrays that are set in VSBCBE
C    when the initial boundary conditions are read.  This routine should reset
C    the QBND values only for boundary nodes on which a pole is applied
C    the pole number is held in IBCTYP in case a tapered pole is required
C    on a boundary segment two pole numbers are held in IBCTYP(JX), IBCTYP(JY),
C    the proportion of which is defined in QBND(JX).  Actual boundary
C    velocity values are set in QLOAD when VSBSET is called from LOADBC.
C
      INCLUDE "limits.parameters"
      DOUBLE PRECISION XLOAD,BIG
      DIMENSION XLOAD(NROWS)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION IBC(NBP),IBCTYP(NBP*2)
      DIMENSION POLEP(3,MAXPOLE),POLET(3)
      DATA IFRST/0/
      SAVE IFRST
C
      RTOD=1.0/DTOR
      S1=0.0
      S2=1.0
      SS=0.0
      NC=0
C
C     For each node on the boundary
C
      DO 500 JX=1,NBP
        JY=JX+NBP
C
C     check if pole rotation applied to this point
C
        IP1=IBCTYP(JX)-500
        IP2=IBCTYP(JY)-500
        IF(IP1.GT.0)THEN
          IF(IP2.EQ.IP1)THEN         ! use only pole1
            DO K=1,3
              POLET(K)=POLEP(K,IP1)
            ENDDO
          ELSE                     ! use pole tapered between pole1 and pole2
            TAPER=QBND(JX)
            CALL PTAPER(-1,S1,S2,SS,TAPER,
     :                  POLEP(1,IP1),POLEP(1,IP2),POLET)
          ENDIF
          NODE=IBC(JX)
          XN=EX(NODE)
          YN=EY(NODE)
C
C   compute and set displacement rates derived from the pole of rotation
C
          IF(POLET(3).NE.0.0)THEN
C
C     get (lat,long) from projected coordinates (radians -> degrees)
C     if in rotated coordinate system CENTLNG should be zero.
C
            BPLON=CENTLNG+RTOD*(XN/SIN(YN))
            BPLAT=RTOD*YN - 90.0
            CALL VELNE(POLET(2),POLET(1),BPLAT,BPLON,UYV,UXV)  ! (lat, long)
            XLOAD(NODE)=BIG*UXV*POLET(3)
            XLOAD(NODE+NUP)=BIG*UYV*POLET(3)
            NC=NC+1
          ENDIF    ! non-zero rotation rate
        ENDIF      ! IP1 != 0, pole applied
  500 CONTINUE
      IF((NC.GT.0).AND.(IFRST.EQ.0))THEN
        WRITE(*,10001)NC
        IFRST=1
      ENDIF
10001 FORMAT('Velocity boundary conditions have been reset on',I6,
     :       ' nodes')
C
      RETURN
      END
C
C****************************************************************************
C
      SUBROUTINE PARSEBE(INSTR,COMM,IXY,XYVAL,IYX,YXLIM1,YXLIM2,
     : IUT,UTVAL1,UTVAL2,ITAP,IPR,XMIN,XMAX,YMIN,YMAX,XPT1,YPT1,
     : XPT2,YPT2,XPT0,YPT0,POLEP,ROTIX,DPHI,CENTLNG,IPOLE,IP1,IP2,
     : IDEFTYP,IERR)
C
C     This subroutine parses the boundary condition statements
C     Used by VSBCBE.  IPR is set to 4 (specified traction or velocity)
C     or 5 (velocity defined by pole) for 'BETWEEN', 3 for 'POLE' command,
C     2 for 'CREF' if parsing is succesful, else IPR = 0 is returned
C     on return from POLE command POLEP = (long, lat, rate)
C
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*2),PNTR*(LINELEN*2)
      CHARACTER IYX*1
      CHARACTER IXY*3,IUT*3,COMM*3
      CHARACTER FRM160*6
      DOUBLE PRECISION ROTIX(3,3)
      DIMENSION POLEP(3,MAXPOLE)
      COMMON/AI/LUW,LSC,LBC,LLG
      IPR=0
      IYX=' '
      ITAP=0
      LENGTH=LINELEN*2
      WRITE(FRM160,10160)LENGTH
10160 FORMAT('(A',I3.3,')')
      CENTL=CENTLNG
      IF(IDEFTYP.EQ.111)CENTL=0.0
C
C     check for valid keywords
C
      CALL SKIP(INSTR,LENGTH,1,'n',J1)
      IF(J1.EQ.0)RETURN
      IF(J1.GE.(LINELEN*2))RETURN
      IF((INSTR(J1:J1+6).NE.'BETWEEN').AND.
     :  (INSTR(J1:J1+3).NE.'CREF').AND.(INSTR(J1:J1+3).NE.'POLE'))RETURN
C
C     A boundary condition command is now recognised
C     Problems in parsing the string go to 200 below
C
      IF(INSTR(J1:J1+3).EQ.'CREF')THEN
        IPR=2
        K1=J1+4
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)XPT0D
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)YPT0D
C
C     If coordinate rotation / spherical projection, CREF in radians
C
        XPT0=XPT0D
        YPT0=YPT0D
        IF(IDEFTYP.EQ.111)THEN
          CALL ROTATE(ROTIX,XPT0,YPT0,1,DPHI)
C         XPT0D=XPT0
C         YPT0D=YPT0
        ENDIF
        IF(IDEFTYP.GE.110)
     :    CALL PROJECTXY(XPT0,YPT0,CENTL,YMID,0,1,IERR)
        WRITE(LSC,10134)XPT0,YPT0,XPT0D,YPT0D
        WRITE(LUW,10134)XPT0,YPT0,XPT0D,YPT0D
10134   FORMAT('Reference point CREF: ',
     :         '(',F8.3,',',F8.3,'), in degrees: (',
     :         F8.3,',',F8.3,')')
        RETURN
      ELSE IF(INSTR(J1:J1+3).EQ.'POLE')THEN
        IPR=3
        IPOLE=IPOLE+1
        IF(IPOLE.GT.MAXPOLE) THEN
          WRITE(LSC,10136)MAXPOLE
          WRITE(LUW,10136)MAXPOLE
10136   FORMAT('Attempting to input more poles than space ',
     :    'permits: ',I3,' in VSBCBE')
        GO TO 200
        END IF
        K1=J1+4
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)PLONG
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)PLAT
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.':')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)RRATE
C
C     If coordinate rotation in effect POLE coords in degrees
C
        IF(IDEFTYP.EQ.111)THEN
          CALL ROTATE(ROTIX,PLONG,PLAT,1,DPHI)
        ENDIF
        WRITE(LSC,10135)IPOLE,PLONG,PLAT,RRATE
        WRITE(LUW,10135)IPOLE,PLONG,PLAT,RRATE
        POLEP(1,IPOLE)=PLONG
        POLEP(2,IPOLE)=PLAT
        POLEP(3,IPOLE)=RRATE
10135   FORMAT('Pole ',I3,' parameters: Long =',F8.3,'  Lat =',F8.3,
     :  '  Rate =',F8.3)
        RETURN
C
C    for the BETWEEN command
C
      ELSE IF(INSTR(J1:J1+6).EQ.'BETWEEN')THEN
        IPR=4
C
C    read the coordinates of the first point
C
        K1=J1+7
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)XPT1D
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)YPT1D
C
C    read 'AND'
C
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        IF(INSTR(J1:J2-1).NE.'AND')GO TO 200
C
C    read the coordinates of the second point
C
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)XPT2D
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)YPT2D
C
C     read ':' and then IUT
C
      K1=J2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.':')GO TO 200
  100 K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      CALL SKIP(INSTR,LENGTH,J1,'b',J2)
      IF(J2.EQ.0)GO TO 200
      IUT=INSTR(J1:J1+1)
C
C    identify pole parameters for a particular segment
C
      IF(IUT.EQ.'PO')THEN
        IPR=5
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)IP1
C
C     look for 'TO' clause if present
C
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)THEN
          IP2=IP1
          GO TO 150
        END IF
        IF(INSTR(J1:J1+1).NE.'TO')GO TO 200
        K1=J1+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        CALL SKIP(INSTR,LENGTH,J1,'b',J2)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J2-1),*)IP2
C
C     read ':' and then 'TP' clause if present
C
        K1=J2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 150
        IF(INSTR(J1:J1).NE.':')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 150
        IF(INSTR(J1:J1+1).NE.'TP')GO TO 200
C
C     read '=' and then ITAP
C
        K1=J1+2
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J1.EQ.0)GO TO 200
        IF(INSTR(J1:J1).NE.'=')GO TO 200
        K1=J1+1
        CALL SKIP(INSTR,LENGTH,K1,'n',J1)
        IF(J2.EQ.0)GO TO 200
        READ(INSTR(J1:J1),*)ITAP
        GO TO 150
      END IF
C
C    identify specific velocity or traction parameters for a segment
C
      IF((IUT.NE.'UX').AND.(IUT.NE.'UY').AND.
     1(IUT.NE.'TX').AND.(IUT.NE.'TY').AND.
     2(IUT.NE.'FL').AND.(IUT.NE.'FU').AND.
     3(IUT.NE.'FX').AND.(IUT.NE.'FY').AND.
     4(IUT.NE.'FN').AND.(IUT.NE.'UB').AND.
     5(IUT.NE.'TN').AND.(IUT.NE.'TT'))GO TO 200
C
C     read '=' and then UTVAL1
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.'=')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      CALL SKIP(INSTR,LENGTH,J1,'b',J2)
      IF(J2.EQ.0)GO TO 200
      READ(INSTR(J1:J2-1),*)UTVAL1
C
C     look for 'TO' clause if present
C
      K1=J2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)THEN
        UTVAL2=UTVAL1
        GO TO 150
      END IF
      IF(INSTR(J1:J1+1).NE.'TO')GO TO 200
C
C     read UTVAL2
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      CALL SKIP(INSTR,LENGTH,J1,'b',J2)
      IF(J2.EQ.0)GO TO 200
      READ(INSTR(J1:J2-1),*)UTVAL2
C
C     read ':' and then 'TP' clause if present
C
      K1=J2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 150
      IF(INSTR(J1:J1).NE.':')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 150
      IF(INSTR(J1:J1+1).NE.'TP')GO TO 200
C
C     read '=' and then ITAP
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.'=')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J2.EQ.0)GO TO 200
      READ(INSTR(J1:J1),*)ITAP
C
C     reading successfully completed
C
  150 CONTINUE
C
C     If coordinate rotation in effect locations are in radians
C      while rotation poles are defined in (rotated) degrees
C
      XPT1=XPT1D
      YPT1=YPT1D
      XPT2=XPT2D
      YPT2=YPT2D
      IF(IDEFTYP.EQ.111)THEN
        CALL ROTATE(ROTIX,XPT1,YPT1,1,DPHI)
C       XPT1D=XPT1
C       YPT1D=YPT1
        CALL ROTATE(ROTIX,XPT2,YPT2,1,DPHI)
C       XPT2D=XPT2
C       YPT2D=YPT2
      ENDIF
      IF(IDEFTYP.GE.110)THEN
        CALL PROJECTXY(XPT1,YPT1,CENTL,YMID,0,1,IERR)
        CALL PROJECTXY(XPT2,YPT2,CENTL,YMID,0,1,IERR)
      ENDIF
      IF(IPR.EQ.4)THEN
        IF(IDEFTYP.GE.110)THEN
          WRITE(LUW,10005)XPT1,YPT1,XPT2,YPT2,IUT,UTVAL1
          WRITE(LSC,10005)XPT1,YPT1,XPT2,YPT2,IUT,UTVAL1
10005     FORMAT('Between(x,y) (',F8.3,',',F8.3,') and (',F8.3,',',
     :            F8.3,') ',A2,' = ',F8.3)
        ENDIF
        WRITE(LUW,10015)XPT1D,YPT1D,XPT2D,YPT2D,IUT,UTVAL1
        WRITE(LSC,10015)XPT1D,YPT1D,XPT2D,YPT2D,IUT,UTVAL1
10015 FORMAT('Between(deg) (',F8.3,',',F8.3,') and (',F8.3,',',
     :       F8.3,') ',A2,' = ',F8.3)
        IF(ITAP.NE.0)THEN
          WRITE(LUW,10003)UTVAL2,ITAP
          WRITE(LSC,10003)UTVAL2,ITAP
10003     FORMAT('   tapered to ',F8.3,' using taper function ',I1)
        END IF
      ELSE IF(IPR.EQ.5)THEN
        IF(DEFTYP.GE.110)THEN
          WRITE(LUW,10006)XPT1,YPT1,XPT2,YPT2,IP1
          WRITE(LSC,10006)XPT1,YPT1,XPT2,YPT2,IP1
10006     FORMAT('Between(x,y) (',F8.3,',',F8.3,') and (',F8.3,',',
     :            F8.3,') use POLE ',I3)
        ENDIF
        WRITE(LUW,10016)XPT1D,YPT1D,XPT2D,YPT2D,IP1
        WRITE(LSC,10016)XPT1D,YPT1D,XPT2D,YPT2D,IP1
10016 FORMAT('Between(deg) (',F8.3,',',F8.3,') and (',F8.3,',',
     :       F8.3,') use POLE ',I3)
        IF(ITAP.NE.0)THEN
          WRITE(LUW,10007)IP2,ITAP
          WRITE(LSC,10007)IP2,ITAP
10007     FORMAT('   tapered to pole',I3,' using taper function ',I1)
        END IF
      END IF
      END IF
      RETURN
C
C     problem in parsing the instruction
C
  200 CONTINUE
      DO J=1,LINELEN*2
        PNTR(J:J)=' '
      ENDDO
      PNTR(K1:K1)='^'
      WRITE(LUW,10004)
      WRITE(LSC,10004)
10004 FORMAT('The following instruction could not be ',
     1'successfully parsed:')
      WRITE(LUW,FRM160)INSTR
      WRITE(LUW,FRM160)PNTR
      WRITE(LSC,FRM160)INSTR
      WRITE(LSC,FRM160)PNTR
      IERR = 1
      RETURN
      END
C
C***************************************************************************************
C
      SUBROUTINE PARSEVALUE(VAL,XMIN,XMAX,YMIN,YMAX,INSTR,J1,J2,ISET)
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*2)

      ISET = 1
C     check for keyword rather than number
      IF (INSTR(J1:J1).EQ.'X') THEN
        IF (INSTR(J1:J2-1).EQ.'XMIN') THEN
          VAL=XMIN
        ELSE IF (INSTR(J1:J2-1).EQ.'XMAX') THEN
          VAL=XMAX
        ELSE
          ISET = 0
        END IF
      ELSE IF (INSTR(J1:J1).EQ.'Y') THEN
        IF (INSTR(J1:J2-1).EQ.'YMIN') THEN
          VAL=YMIN
        ELSE IF (INSTR(J1:J2-1).EQ.'YMAX') THEN
          VAL=YMAX
        ELSE
          ISET = 0
        END IF
      ELSE 
        READ(INSTR(J1:J2-1),*)VAL
      END IF
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE PARSEON(INSTR,IXY,XYVAL,IYX,YXLIM1,YXLIM2,IBND,
     1 IUT,UTVAL1,UTVAL2,ITAP,IPR,XMIN,XMAX,YMIN,YMAX,IERR)
C
C     This subroutine parses the boundary condition statements
C     Used by VSBCON.  IPR = 1 is returned if parsing of 'ON'
C     command is successful, else IPR = 0 is returned
C
      INCLUDE "limits.parameters"
      CHARACTER INSTR*(LINELEN*2),PNTR*(LINELEN*2)
      CHARACTER IYX*1,IXY*3,IUT*3,CMD*2
      CHARACTER FRM160*6
      COMMON/AI/LUW,LSC,LBC,LLG
      IPR=0
      IYX=' '
      ITAP=0
      IBND=0
      YXLIM1=-99999.9
      YXLIM2= 99999.9
      LENGTH=(LINELEN*2)
      WRITE(FRM160,10160)LENGTH
10160 FORMAT('(A',I3.3,')')
C
C     check for keyword 'ON'
C
      CALL SKIP(INSTR,LENGTH,1,'n',J1)
      IF(J1.EQ.0)RETURN
      IF(J1.GE.(LINELEN*2))RETURN
      CMD=INSTR(J1:J1+1)
      IF((CMD.NE.'ON').AND.(CMD.NE.'PB'))RETURN
C
C     A boundary condition command is now recognised
C     Problems in parsing the string go to 200 below
C
C     find the next non-zero character and set the IXY value
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IXY=INSTR(J1:J1+1)
      IF((IXY(1:1).NE.'X').AND.(IXY(1:1).NE.'Y').AND.(IXY(1:1).NE.'B'))
     1    GO TO 200
      IF((IXY(2:2).NE.'N').AND.(IXY(2:2).NE.'+')
     1    .AND.(IXY(2:2).NE.'-').AND.(IXY(2:2).NE.' '))GO TO 200
      IF(IXY(2:2).EQ.'N')THEN
        J1=J1+1
        IF(INSTR(J1+1:J1+1).NE.'D')GO TO 200
        IXY(3:3)=INSTR(J1+1:J1+1)
      END IF
C
C     read '=' and then XYVAL
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.'=')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      CALL SKIP(INSTR,LENGTH,J1,'b',J2)
      IF(J2.EQ.0)GO TO 200
      CALL PARSEVALUE(XYVAL,XMIN,XMAX,YMIN,YMAX,INSTR,J1,J2,ISET)
      IF (IXY(1:2).EQ.'BN') IBND=INT(XYVAL)
      IF(CMD.EQ.'PB')RETURN
      IF(ISET.EQ.0)GO TO 200
C
C     look for 'FOR' clause if present
C
      K1=J2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).EQ.':')GO TO 100
      IF(INSTR(J1:J1+2).NE.'FOR')GO TO 200
C
C     find the next non-zero character and set the IYX value
C
      K1=J1+3
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IYX=INSTR(J1:J1)
      IF((IXY(1:1).EQ.'X').AND.(IYX.NE.'Y'))GO TO 200
      IF((IXY(1:1).EQ.'Y').AND.(IYX.NE.'X'))GO TO 200
C
C     read '=' and then YXLIM1
C
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.'=')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      CALL SKIP(INSTR,LENGTH,J1,'b',J2)
      IF(J2.EQ.0)GO TO 200
      CALL PARSEVALUE(YXLIM1,XMIN,XMAX,YMIN,YMAX,INSTR,J1,J2,ISET)
      IF (ISET.EQ.0) GO TO 200
C     READ(INSTR(J1:J2-1),*)YXLIM1
C
C     look for 'TO' clause if present
C
      K1=J2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).EQ.':')THEN
        YXLIM2=YXLIM1
        GO TO 100
      END IF
      IF(INSTR(J1:J1+1).NE.'TO')GO TO 200
C
C     read YXLIM2
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      CALL SKIP(INSTR,LENGTH,J1,'b',J2)
      IF(J2.EQ.0)GO TO 200
      CALL PARSEVALUE(YXLIM2,XMIN,XMAX,YMIN,YMAX,INSTR,J1,J2,ISET)
      IF (ISET.EQ.0) GO TO 200
C     READ(INSTR(J1:J2-1),*)YXLIM2
C
C     read ':' and then IUT
C
      K1=J2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.':')GO TO 200
  100 K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IUT=INSTR(J1:J1+1)
      IF((IUT.NE.'UX').AND.(IUT.NE.'UY').AND.
     :(IUT.NE.'TX').AND.(IUT.NE.'TY').AND.(IUT.NE.'TN').AND.
     :(IUT.NE.'TT').AND.(IUT.NE.'FL').AND.(IUT.NE.'FU').AND.
     :(IUT.NE.'FT').AND.(IUT.NE.'FN').AND.(IUT.NE.'FX').AND.
     :(IUT.NE.'FY').AND.(IUT.NE.'UB').AND.(IUT.NE.'BN'))GO TO 200
      IF(IUT.EQ.'BN')THEN
        J1=J1+1
        IF(INSTR(J1+1:J1+1).NE.'D')GO TO 200
        IUT(3:3)=INSTR(J1+1:J1+1)
      ELSE IF(IUT.EQ.'FL')THEN      ! locked fault does not require parameters
        UTVAL1=0.0
        GO TO 150
      END IF
C
C     read '=' and then UTVAL1
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.'=')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      CALL SKIP(INSTR,LENGTH,J1,'b',J2)
      IF(J2.EQ.0)GO TO 200
      READ(INSTR(J1:J2-1),*,ERR=200)UTVAL1
C
C     look for 'TO' clause if present
C
      K1=J2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)THEN
        UTVAL2=UTVAL1
        GO TO 150
      END IF
      IF(INSTR(J1:J1+1).NE.'TO')GO TO 200
C
C     read UTVAL2
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      CALL SKIP(INSTR,LENGTH,J1,'b',J2)
      IF(J2.EQ.0)GO TO 200
      READ(INSTR(J1:J2-1),*,ERR=200)UTVAL2
C
C     read ':' and then 'TP' clause if present
C
      K1=J2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 150
      IF(INSTR(J1:J1).NE.':')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 150
      IF(INSTR(J1:J1+1).NE.'TP')GO TO 200
C
C     read '=' and then ITAP
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.'=')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      CALL SKIP(INSTR,LENGTH,J1,'b',J2)
      IF(J2.EQ.0)GO TO 200
      READ(INSTR(J1:J2-1),*)ITAP
C
C     read ':' and then 'BND' clause if present
C
      K1=J2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 150
      IF(INSTR(J1:J1).NE.':')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 150
      IF(INSTR(J1:J1+1).NE.'TP')GO TO 200
C
C     read '=' and then IBND
C
      K1=J1+2
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      IF(INSTR(J1:J1).NE.'=')GO TO 200
      K1=J1+1
      CALL SKIP(INSTR,LENGTH,K1,'n',J1)
      IF(J1.EQ.0)GO TO 200
      READ(INSTR(J1:J1),*)IBND
C
C     reading successfully completed
C
  150 CONTINUE
      IPR=1
      RETURN
C
C     problem in parsing the instruction
C
  200 CONTINUE
      DO J=1,LINELEN*2
        PNTR(J:J)=' '
      ENDDO
      PNTR(K1:K1)='^'
      WRITE(LUW,10004)
      WRITE(LSC,10004)
10004 FORMAT('The following instruction could not be ',
     1'successfully parsed:')
      WRITE(LUW,FRM160)INSTR
      WRITE(LUW,FRM160)PNTR
      WRITE(LSC,FRM160)INSTR
      WRITE(LSC,FRM160)PNTR
      IERR = 1
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE LOOK(INSTR,LENG,SUBJ,JL)
C
C   Locates the substring SUBJ in the string INSTR
C
      CHARACTER INSTR*1
      CHARACTER SUBJ*2
      DIMENSION INSTR(LENG)
      JL=0
      IF(LENG.LE.1)RETURN
      DO 10 J=1,LENG-1
        IF(INSTR(J).EQ.SUBJ(1:1))THEN
          IF(INSTR(J+1).EQ. SUBJ(2:2))THEN
            JL=J
            GO TO 20
          END IF
        END IF
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE TAPER(ITAP,YXLIM1,YXLIM2,UTVAL1,UTVAL2,
     1 YXN,UTVAL)
C
C     applies a taper between limits YXLIM1, YXLIM2, to 
C     return a value UTVAL between UTVAL1 and UTVAL2
C
C      ITAP = 0 is constant (no taper)
C      ITAP = 1 is linear ramp
C      ITAP = 2 is sin**2 taper (zero gradient at both ends)
C      ITAP = 3 is cos taper (zero gradient at YXLIM1)
C      ITAP = 4 is sin taper (zero gradient at YXLIM2)
C      ITAP = 5 is parabolic taper, zero gradient at XYLIM2)
C       other values of ITAP may be included later
C
      SAVE PION2
      DATA PION2/1.570796327/
C
      UTVAL=UTVAL1
      IF(YXLIM2.EQ.YXLIM1)RETURN
      RAT=((YXN-YXLIM1)/(YXLIM2-YXLIM1))
      CSRAT=0.0
      IF(ITAP.EQ.1)THEN
        CSRAT=RAT
      ELSE IF(ITAP.EQ.2)THEN
        CSRAT=SIN(RAT*PION2)
        CSRAT=CSRAT*CSRAT
      ELSE IF(ITAP.EQ.3)THEN
        CSRAT=1.0-COS(RAT*PION2)
      ELSE IF(ITAP.EQ.4)THEN
        CSRAT=SIN(RAT*PION2)
      ELSE IF(ITAP.EQ.5)THEN
        CSRAT=(2.0-RAT)*RAT
      END IF
C
C    set the value for return
C
      UTVAL=UTVAL1+(UTVAL2-UTVAL1)*CSRAT
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE PTAPER(ITAP,S1,S2,SS,CSRAT,POLE1,POLE2,POLE)
C
C     applies a taper to rotation pole vectors to return
C     pole parameters based on a  Cartesian interpolation
C     between pole vectors POLE1 and POLE2, each vector
C     a triple (long, lat, rotation rate)
C     Note: returns POLE1 and CSRAT=1 when SS=S1
C     Note: returns POLE2 and CSRAT=0 when SS=S2
C
C      ITAP < 0 use CSRAT from input line
C      ITAP = 0 is constant (no taper)
C      ITAP = 1 is linear ramp
C      ITAP = 2 is sin**2 taper (zero gradient at both ends)
C      ITAP = 3 is cos taper (zero gradient at S1)
C      ITAP = 4 is sin taper (zero gradient at S2)
C      ITAP = 5 is parabolic taper, zero gradient at S2)
C       other values of ITAP may be included later
C
      DIMENSION POLE1(3),POLE2(3),POLE(3)
      SAVE PION2,DTOR
      DATA PION2,DTOR/1.570796327,0.0174532925/
C
C    defaults allow for zero pole etc
C
      DO K=1,3
        POLE(K)=POLE1(K)
      ENDDO
C
C    calculate the taper factor and set CSRAT
C
      IF(ITAP.GE.0)THEN
        CSRAT=0.0
        IF(S2.EQ.S1)GOTO 10
        RAT=(SS-S1)/(S2-S1)
        IF(ITAP.EQ.1)THEN
          CSRAT=RAT
        ELSE IF(ITAP.EQ.2)THEN
          CSRAT=SIN(RAT*PION2)
          CSRAT=CSRAT*CSRAT
        ELSE IF(ITAP.EQ.3)THEN
          CSRAT=1.0-COS(RAT*PION2)
        ELSE IF(ITAP.EQ.4)THEN
          CSRAT=SIN(RAT*PION2)
        ELSE IF(ITAP.EQ.5)THEN
          CSRAT=(2.0-RAT)*RAT
        END IF
      ENDIF
C
C    convert POLE1 and POLE2 to Cartesian vectors
C
      XLON1=POLE1(1)*DTOR
      XLAT1=POLE1(2)*DTOR
      XLON2=POLE2(1)*DTOR
      XLAT2=POLE2(2)*DTOR
      CALL SPHCAR(XLAT1,XLON1,X1,Y1,Z1)
      CALL SPHCAR(XLAT2,XLON2,X2,Y2,Z2)
      X1=X1*POLE1(3)
      Y1=Y1*POLE1(3)
      Z1=Z1*POLE1(3)
      X2=X2*POLE2(3)
      Y2=Y2*POLE2(3)
      Z2=Z2*POLE2(3)
C
C    calculate the interpolated pole vector
C
      XX=X1+(X2-X1)*CSRAT
      YY=Y1+(Y2-Y1)*CSRAT
      ZZ=Z1+(Z2-Z1)*CSRAT
   10 CSRAT=1.0-CSRAT
C
C   get magnitude and pole coordinates of interpolated pole
C
      POLE(3)=SQRT(XX*XX+YY*YY+ZZ*ZZ)
      IF(POLE(3).EQ.0.0)THEN
        POLE(1)=0.0
        POLE(2)=0.0
      ELSE
        ZZ=ZZ/POLE(3)
        RTOD=1.0/DTOR
        POLE(2)=RTOD*ASIN(ZZ)
        POLE(1)=RTOD*ATAN2(YY,XX)
      ENDIF
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE FINDMLIM(EX,EY,NUP,XMIN,XMAX,YMIN,YMAX)
C
C    This routine is used to find X and Y limits of finite element mesh
C
      DIMENSION EX(NUP),EY(NUP)
      XMIN=EX(1)
      YMIN=EY(1)
      XMAX=EX(1)
      YMAX=EY(1)
      DO JE=2,NUP
        IF(EX(JE).LT.XMIN)XMIN=EX(JE)
        IF(EY(JE).LT.YMIN)YMIN=EY(JE)
        IF(EX(JE).GT.XMAX)XMAX=EX(JE)
        IF(EY(JE).GT.YMAX)YMAX=EY(JE)
      ENDDO
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE FINDBBOX(EX,EY,BBOX,NUP)
      DIMENSION EX(NUP), EY(NUP), BBOX(2,4)

C    finds the corners of the bounding box (assuming y values
C    for top and bottom boundaries are constant). Stores the
C    coordinates in BBOX starting at the bottom left and in
C    anticlockwise order

      EPS = 1.0E-6
      EPS = 1.0E-5  !basil1.2.7
      BBOX(2,1) = EY(1)
      BBOX(2,2) = EY(1)
      BBOX(2,3) = EY(1)
      BBOX(2,4) = EY(1)
      DO 10 I=2,NUP
        IF (EY(I).GT.BBOX(2,3)) BBOX(2,3) = EY(I)
        IF (EY(I).LT.BBOX(2,2)) BBOX(2,2) = EY(I)
  10  CONTINUE
      BBOX(2,1) = BBOX(2,2)
      BBOX(2,4) = BBOX(2,3)
      YMIN = BBOX(2,2)+EPS
      YMAX = BBOX(2,3)-EPS
      ISETMIN1 = 0
      ISETMIN2 = 0
      ISETMAX1 = 0
      ISETMAX2 = 0
      DO 20 I=1,NUP
        IF (EY(I).LT.YMIN) THEN
          IF (ISETMIN1.EQ.0) THEN
            BBOX(1,1) = EX(I)
            ISETMIN1 = 1
          ELSE IF (EX(I).LT.BBOX(1,1)) THEN
            BBOX(1,1) = EX(I)
          ENDIF
          IF (ISETMAX1.EQ.0) THEN
            BBOX(1,2) = EX(I)
            ISETMAX1 = 1
          ELSE IF (EX(I).GT.BBOX(1,2)) THEN
            BBOX(1,2) = EX(I)
          ENDIF
        ELSE IF (EY(I).GT.YMAX) THEN
          IF (ISETMIN2.EQ.0) THEN
            BBOX(1,4) = EX(I)
            ISETMIN2 = 1
          ELSE IF (EX(I).LT.BBOX(1,4)) THEN
            BBOX(1,4) = EX(I)
          ENDIF
          IF (ISETMAX2.EQ.0) THEN
            BBOX(1,3) = EX(I)
            ISETMAX2 = 1
          ELSE IF (EX(I).GT.BBOX(1,3)) THEN
            BBOX(1,3) = EX(I)
          ENDIF
        END IF
  20  CONTINUE
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE STORETIPS(ISEG,NSEG,IFLTTIPS,IBCTYP,IBC,NBP,NOR,NUP)
      INCLUDE "limits.parameters"
      DIMENSION ISEG(3,NSEG)
      DIMENSION IFLTTIPS(MAXTIPS*2)
      DIMENSION IBC(NBP)
      DIMENSION IBCTYP(NBP*2)
      DIMENSION NOR(NUP)

      J=1
      DO 25 I=1,NBP
        IF (IBCTYP(I).EQ.15) THEN
          N=NOR(IBC(I))
          IFLTTIPS(J)=N
          DO 20 K=1,NSEG
            IF ((ISEG(1,K).EQ.N.OR.ISEG(2,K).EQ.N).AND.
     :        (ISEG(3,K).GE.FLTMIN.AND.ISEG(3,K).LE.FLTMAX)) THEN
              IFLTTIPS(J+MAXTIPS)= ISEG(3,K)
              J=J+1
C             IBCTYP(I)=11
C             IBCTYP(I+NBP)=11
              GO TO 25
            END IF
  20      CONTINUE
        END IF
  25  CONTINUE
      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE RESETTIPS(IFLTTIPS,IBCTYP,IBC,NBP,NOR,NUP)
      INCLUDE "limits.parameters"
      DIMENSION IFLTTIPS(MAXTIPS*2)
      DIMENSION IBC(NBP)
      DIMENSION IBCTYP(NBP*2)
      DIMENSION NOR(NUP)

      MAX=MAXTIPS
      DO 25 I=1,NBP
          N=NOR(IBC(I))
          DO 20 K=1,MAX
            IF (IFLTTIPS(K).EQ.N) THEN
              IBCTYP(I)=11
              IBCTYP(I+NBP)=11
              GO TO 25
            END IF
  20      CONTINUE
  25  CONTINUE
      RETURN
      END
C
C   next 3 subroutines provided by Richard Gordon
C   July 2005 (VELNE, SPHCAR, GETNED)
C
C****************************************************************************************
C
C   Output from Public domain Ratfor, version 1.0
      SUBROUTINE VELNE(POLAT,POLON,SITLAT,SITLON,VNORTH,VEAST)
      CON=3.1415926/180.
      omega=1.
      CALL SPHCAR(con*polat,con*polon,p1,p2,p3)
      CALL SPHCAR(con*sitlat,con*sitlon,a1,a2,a3)
      r = 1
      vx = omega*r*(p2*a3-a2*p3)
      vy = omega*r*(p3*a1-a3*p1)
      vz = omega*r*(p1*a2-a1*p2)
      CALL GETNED(con*sitlat,con*sitlon,vx,vy,vz,vnorth,veast,vdown)
      RETURN
      END
C
C****************************************************************************************
C
C Output from Public domain Ratfor, version 1.0
      SUBROUTINE SPHCAR(LAT,LON,X1,X2,X3)
      real lat, lon, x1, x2, x3
      x1 = cos(lat) * cos(lon)
      x2 = cos(lat) * sin(lon)
      x3 = sin(lat)
      RETURN
      END
C
C****************************************************************************************
C
C Output from Public domain Ratfor, version 1.0
      SUBROUTINE GETNED(SITLAT,SITLON,X,Y,Z,NORTH,EAST,DOWN)
      real north
      north = -sin(sitlat)*cos(sitlon)*x-sin(sitlat)*sin(sitlon)*y+cos(s
     *itlat)*z
      east = -sin(sitlon)*x+cos(sitlon)*y
      down = -cos(sitlat)*cos(sitlon)*x-cos(sitlat)*sin(sitlon)*y-sin(si
     *tlat)*z
      RETURN
      END
