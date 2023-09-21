C*--------------------------------------------------------------------
C*    Basil / Sybil:   fltbcon.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE PERIBC(EX,EY,QBND,LEM,NOR,IBC,IBNGH,IBCTYP,
     :                 IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :                 IFLT,NUP,NE,NBP,NFP,NFPF3,KPBC,
     :                 PDIST,IREM,LUW,LSC,IERR)
C
C   this routine replaces FLTBC, and adds entries to the Fault arrays 
C   for periodic external boundaries.  The setting of the same arrays
C   for internal fault boundaries is now done partly in createafault
C   (trimesh.c) and partly in FLTINIT called from RUNBASIL
C
      INCLUDE "limits.parameters"
      CHARACTER IXY*3
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION QBND(NBP*2)
      DIMENSION IBC(NBP),IBNGH(NBP*2),IBCTYP(NBP*2)  ! IBNGH not set here; needed ?
      DIMENSION IFBC1(NFP),IFBC2(NFP),IFEQV(NFP)
      DIMENSION JFBC1(NFP),JFBC2(NFP)
C
C  NFPF3 = the number of nodes on internal faults
C  NFP   = number of fault nodes + periodic boundary nodes
C  IFBC1 = the node numbers that are to be treated as faults
C          one side of a fault is assigned negative values
C  IFBC2 = the position in the array IBC of the node given in IFBC1
C  IFEQV = gives the position in the array IFBC1 of the equivalent node 
C          on the other side of the fault (e.g. if J=IFEQV(I),
C          then IFBC1(I) and IFBC1(J) are across the fault from each other)
C  JFBC1 = the element number that contains the node in IFBC1
C          if the node is a midpoint node, 0 if it is a vertex node.
C  JFBC2 = the position (i.e. 4-6) in the element JFBC1 of the
C          node in IFBC1. If JFBC1=0 then JFBC2=0
C
      EPS=1.E-5
      IF(NFP.LE.NFPF3)THEN
        WRITE(LUW,10099)
        WRITE(LSC,10099)
        RETURN
      ENDIF
C
C    find matching points on opposite boundaries by looking at all pairs
C    of points in IBC; the criterion is that the first coordinate differs
C    by PDIST and the second coordinate must be the same (within EPS)
C    N.B. complication with paired fault nodes needs further thought
C    based on assignment of IPRES below
C
      INDXF=NFPF3
      IF(KPBC.EQ.1)THEN              ! periodicity in the X-direction
        DO I=1,NBP
          II=IBC(I)
          IPRES=0
          IF(NFPF3.GT.0)THEN
            DO K=1,NFPF3
              IF((IFBC1(K).EQ.II).OR.(IFBC1(K).EQ.-II))THEN
                IPRES=IFBC1(K)
C               WRITE(*,*)'Found node II, IPRES==',II,IPRES
                GO TO 10
              ENDIF
            ENDDO
          ENDIF
C  10     PA1=EX(NOR(II))
C         PB1=EY(NOR(II))
   10     PA1=EX(II)
          PB1=EY(II)
          DO J=1,NBP
            JJ=IBC(J)
C           PA2=EX(NOR(JJ))
C           PB2=EY(NOR(JJ))
            PA2=EX(JJ)
            PB2=EY(JJ)
            IF((ABS(PA2-PA1-PDIST).LT.EPS).AND.
     :         (ABS(PB1-PB2)      .LT.EPS))THEN
              INDXF=INDXF+1
              IFBC1(INDXF)=-II
              IFBC2(INDXF)=I
              IFEQV(INDXF)=INDXF+1
              INDXF=INDXF+1
              IFBC1(INDXF)=JJ
              IFBC2(INDXF)=J
              IFEQV(INDXF)=INDXF-1
            ENDIF
          ENDDO                  ! loop on J
        ENDDO                    ! loop on I
        WRITE(LSC,10100)'X',PDIST,INDXF-NFPF3
      ELSE IF(KPBC.EQ.2)THEN      ! periodicity in the Y-direction
        DO I=1,NBP
          II=IBC(I)
          IPRES=0
          IF(NFPF3.GT.0)THEN
            DO K=1,NFPF3
              IF((IFBC1(K).EQ.II).OR.(IFBC1(K).EQ.-II))THEN
                IPRES=IFBC1(K)
C               WRITE(*,*)'Found node II, IPRES==',II,IPRES
                GO TO 20
              ENDIF
            ENDDO
          ENDIF
C  20     PA1=EY(NOR(II))
C         PB1=EX(NOR(II))
   20     PA1=EY(II)
          PB1=EX(II)
          DO J=1,NBP
            JJ=IBC(J)
C           PA2=EY(NOR(JJ))
C           PB2=EX(NOR(JJ))
            PA2=EY(JJ)
            PB2=EX(JJ)
            IF((ABS(PA2-PA1-PDIST).LT.EPS).AND.
     :         (ABS(PB1-PB2)      .LT.EPS))THEN
              INDXF=INDXF+1
              IFBC1(INDXF)=-II
              IFBC2(INDXF)=I
              IFEQV(INDXF)=INDXF+1
              INDXF=INDXF+1
              IFBC1(INDXF)=JJ
              IFBC2(INDXF)=J
              IFEQV(INDXF)=INDXF-1
            ENDIF
          ENDDO
        ENDDO
        WRITE(LSC,10100)'Y',PDIST,INDXF-NFPF3
      ENDIF                    ! X-direction or Y-direction periodicity (KPBC)
C     WRITE(*,*)'IFEQV array follows NFP=',NFP
C     CALL IMATPP(IFEQV,NFP)
C     CALL MITPRT(IFEQV,NFP,1,NFP)
C
C   To complete JFBC1, JFBC2 find triangles with a midpoint in the
C    IFBC1 values defined above
C
      ICOUNT=0
      DO JEL=1,NE
        DO J5=4,6
          K5=LEM(J5,JEL)
          DO II=NFPF3+1,NFP
            IF(ABS(IFBC1(II)).EQ.K5)THEN
              JFBC1(II)=JEL
              JFBC2(II)=J5
              ICOUNT=ICOUNT+1
            END IF
          ENDDO
        ENDDO
      ENDDO
C     WRITE(*,*)'JFBC1 array follows NFP=',NFP
C     CALL IMATPP(JFBC1,NFP)
C     WRITE(*,*)'JFBC2 array follows NFP=',NFP
C     CALL IMATPP(JFBC2,NFP)
C
C    Finally, set IBCTYP = 11 as used for a locked fault (not applied
C    for remesh, where IBCTYP is picked up from retained point attributes)
C
      IF(IREM.EQ.0)THEN
        DO K=NFPF3+1,NFP
          NIBC=IFBC2(K)
          IBCTYP(NIBC)=11
          IBCTYP(NIBC+NBP)=11
        ENDDO
      ENDIF
C     IF(IFLT.EQ.0)IFLT=1
      RETURN
C
10099 FORMAT('A PB command has been issued without necessary ',
     :        'allocation of memory - no action')
10100 FORMAT('Periodic boundaries in ',A1,'-direction now implemen',
     :'ted, lambda = ',F8.4,' with ',I5,' nodes.')
      END
C
      SUBROUTINE PBWIDTH(EX,EY,IFBC1,IFEQV,NOR,NUP,NFP,NFPF3,KPBC,PDIST)
C
C   This routine called only if remeshing required on a problem with
C   periodic boundaries.  The current periodicity is obtained prior
C   to remeshing
C
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION IFBC1(NFP),IFEQV(NFP)
C
      PDSUM=0.0
      DO K=NFPF3+1,NFP
        NP1=ABS(IFBC1(K))
        NP2=ABS(IFBC1(IFEQV(K)))
        IF(KPBC.EQ.1)THEN
          PD1=ABS(EX(NP1)-EX(NP2))
        ELSE IF(KPBC.EQ.2)THEN
          PD1=ABS(EY(NP1)-EY(NP2))
        ENDIF
        PDSUM=PDSUM+PD1
      ENDDO
      PDIST=PDSUM/FLOAT(NFP-NFPF3)
      RETURN
      END
C
      SUBROUTINE FLTINIT(EX,EY,LEM,NOR,IBC,IFBC1,IFBC2,IFEQV,
     :                   JFBC1,JFBC2,NUP,NE,NFP,NFPF3,NBP,
     :                   NSEG,LUW,LSC,IERR)
C
C   this routine is called by the MESH command if IMSH=3 and IFLT > 0 to 
C   establish arrays that are used by SEMFLT to set continuity conditions 
C   on fault-type discontinuities.  FLTINIT replaces part of the function of 
C   FLTBC and, subject to value of IFLT, should be called whenever a fault 
C   discontinuity is present, before boundary conditions are set.
C
C   we first check the entire mesh to find midpoints that are present only in
C   one triangle - they must be either on external boundary or on a fault.
C
C   NFP = the number of fault nodes (2 x the number duplicated)
C         actually = NFPF3 obtained in trimesh.
C   IFBC1(1,NFP) = the node numbers that fall on a fault
C   IFBC2(1,NFP) = the position in the array IBC of the node given in IFBC1
C   IFEQV(1,NFP) = the position in the array IFBC1 of the equivalent node 
C                  on the other side of the fault e.g. if J=IFEQV(I),
C                  IFBC1(I) and IFBC1(J) face each other across the fault
C   JFBC1 = the element number that contains the node in IFBC
C          if the node is a midpoint node, 0 if it is a vertex node.
C   JFBC2 = the position (i.e. 4-6) in the element JFBC1 of the
C          node in IFBC. If JFBC1=0 then JFBC2=0
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP),NTMP(NUP)
      DIMENSION IBC(NBP)
      DIMENSION IFBC1(NFP),IFBC2(NFP),IFEQV(NFP)
      DIMENSION JFBC1(NFP),JFBC2(NFP)
      EPS=1.E-6
C
C     WRITE(*,*)'FLTINIT: NFP = ',NFP,' NFPF3 =',NFPF3
C
C   the arrays IFBC1 and IFBC2 have been established by trimesh if
C   this routine is called. To associate pairs we search the IBC 
C   array for values that have the same location
C
      DO I=1,NFPF3
        N1=IFBC1(I)
        EX1=EX(IABS(N1))
        EY1=EY(IABS(N1))
        DO J=1,NFPF3
          N2=IFBC1(J)
          IF(N2.NE.N1)THEN
            EX2=EX(IABS(N2))
            EY2=EY(IABS(N2))
            DX=EX2-EX1
            DY=EY2-EY1
            DS=SQRT(DX*DX+DY*DY)
            IF(DS.LT.EPS)THEN
              IFEQV(I)=J
              IFEQV(J)=I
C             WRITE(*,*)'Fault node pairs are I,J,N1,N2 = ',I,J,N1,N2
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
C   to associate element numbers and fault midpoints, a fault midpoint
C   can be found in only one element
C
      DO K=1,NFP
        JFBC1(K)=0
        JFBC2(K)=0
      ENDDO
      IFOUND=0
      DO JEL=1,NE
        DO JM=4,6
          JMPN=LEM(JM,JEL)
          DO K=1,NFP
            IF(ABS(IFBC1(K)).EQ.JMPN)THEN
              JFBC1(K)=JEL
              JFBC2(K)=JM
              IFOUND=IFOUND+1
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      WRITE(*,100)IFOUND,NSEG
  100 FORMAT('Fault segments: ',I5,' Total labelled segments: ',I5,
     ://,'*********************************************************',
     :  '*********************',/)
C     WRITE(*,*)'JFBC1 at end of FLTINIT NFP =',NFP
C     CALL IMATPP(JFBC1,NFP)
C     CALL MITPRT(JFBC1,2,NFP/2,NFP,6)
      RETURN
C    
C    consistency checks only activated for debugging
C
      DO I=1,NFP
        IF(JFBC1(I).NE.0)THEN
          MM=LEM(JFBC2(I),JFBC1(I))
          IV1=MOD(JFBC2(I)-2,3)+1
          IV2=MOD(JFBC2(I)-1,3)+1
          JV1=LEM(IV1,JFBC1(I))
          JV2=LEM(IV2,JFBC1(I))
          WRITE(*,*)'3 nodes for each midpoint ',I,JV1,MM,JV2,IFBC1(I)
        ELSE
          WRITE(*,*)'vertex nodes:',I,IFBC1(I)
        ENDIF
      ENDDO
C
      RETURN
      END
C*--------------------------------------------------------------------
C*    Basil / Sybil:   cgbc.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE SEMFLT(BIG,STK2,STK3,IK2,NK2,K2BW,
     :                  EX,EY,LEM,NOR,IBCTYP,IBNGH,QLOAD,
     :                  QBND,IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :                  VELXO,VELYO,IOFF,IFLT,IREM,
     :                  NUP,NE,NROWS,NBP,NFP,NFPF3,LUW,IERR)
C
C    This routine uses the IFBC1 array to identify nodes / segments
C    on a fault and computes the components of the local unit normal to
C    the fault and loads those elements into the STK2 array
C    It must be recalled whenever the geometrical arrangement of
C    of the elements changes, but is independent of changes to
C    rheological coefficients etc.  It is called by CGRUN as part
C    of the matrix assembly procedure
C
      DOUBLE PRECISION BIG,ANX,ANY,QSR,D1,D2,DST
      DOUBLE PRECISION STK2,STK3
      DOUBLE PRECISION QNX,QNY,QDST1,QDST2,QDSTM
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION QLOAD(NROWS),QBND(NBP*2)
      DIMENSION IBCTYP(NBP*2),IBNGH(NBP*2)
      DIMENSION IFBC1(NFP),IFBC2(NFP),IFEQV(NFP)
      DIMENSION JFBC1(NFP),JFBC2(NFP)
      DIMENSION QNX(NFP),QNY(NFP)
      DIMENSION QDST1(NFP),QDST2(NFP),QDSTM(NFP)
C     DIMENSION IDST(NFP)
      DIMENSION LIK(3),LJK(3),LKK(3)
      DIMENSION STK2(NK2),STK3(NFP),IK2(NK2)
C
C  NFP = the number of nodes on the fault
C  IFBC1 = the node numbers that fall on the fault
C         (if INFLT=1, negative are on the side x=0)
C         (if INFLT=2, negative are on the side closer to the origin)
C  IFBC2 = the position in the arrays IBC, IBCTYP of the node given
C          in IFBC1.  IBC includes boundaries and faults
C  IFEQV = gives the position in the array IFBC1 of the equivalent node 
C          on the other side of the fault (e.g. if J=IFEQV(I),
C          then IFBC1(I) and IFBC1(J) are across the fault from each other)
C  JFBC1 = the element number that contains the node in IFBC1
C          if the node is a midpoint node, 0 if it is a vertex node.
C  JFBC2 = the position (i.e. 4-6) in the element JFBC1 of the
C          node in IFBC1. If JFBC1=0 then JFBC2=0
C
      NFPO2=NFP/2
      IERR=0
C     WRITE(*,*)'SEMFLT: NBP, NFP, NFPF3, NK2 =',NBP,NFP,NFPF3,NK2
C
      DO I=1,NK2
        STK2(I)=0.D0
        IK2(I)=0
      ENDDO
      DO I=1,NFP
        STK3(I)=0.D0
      ENDDO
C
C   temporary arrays locally allocated
      DO I=1,NFP
        QNX(I)=0.0
        QNY(I)=0.0
        QDST1(I)=0.0
        QDST2(I)=0.0
        QDSTM(I)=0.0
C       IDST(I)=0
      ENDDO
C
C  for each node on the fault(s), identify if it is a midpint node
C
      DO 10 I=1,NFP
        IP1=IFBC1(I)
        JEL=JFBC1(I)
C
C  If JEL!=0, then IP1 is a midpoint node of the element JEL.
C  LIK contains the absolute node numbers,
C  LJK contains the local node number and 
C  LKK contains the node positions in the IFBC1 array.
C  boundary segment is LIK(1)-LIK(3)-LIK(2), counter-clockwise
C  although all fault nodes are listed in IFBC1,2, the continuity conditions
C  require that the tractions are the same, so we only need to process nodes
C  on one side of the fault. 
C
        IF((JEL.NE.0).AND.(IP1.GT.0)) THEN
          LJK(3)=JFBC2(I)
          LJK(2)=LJK(3)-3
          LJK(1)=MOD(LJK(2)+1,3)+1
          LIK(1)=LEM(LJK(1),JEL)    ! node no. of vertex preceding
          LIK(2)=LEM(LJK(2),JEL)    ! node no. of vertex following
          LIK(3)=IP1                ! midpoint node no.
          LK1=NOR(LIK(1))
          LK2=NOR(LIK(2))
          D1=EX(LK2)-EX(LK1)
          D2=EY(LK2)-EY(LK1)
          DST=DSQRT(D1*D1+D2*D2)
C
C   find the addresses to store D1 and D2; store in QNX and QNY
C   In general we expect 1 (midpoints) or 2 (vertex) contributions
C   to QNX or QNY.  Endpoints of a terminating fault also expect 1.
C
          LKK(1)=0
          LKK(2)=0
          LKK(3)=I
          DO 15 II=1,NFP
            IEX=IFBC1(II)
            IF(LIK(1).EQ.IEX)THEN
              LKK(1)=II
              QDST1(II)=DST
            ENDIF
            IF(LIK(2).EQ.IEX)THEN
              LKK(2)=II
              QDST2(II)=DST
            ENDIF
   15     CONTINUE
          DO 20 JJ=1,3
            IF(LKK(JJ).NE.0) THEN
              QNX(LKK(JJ))=QNX(LKK(JJ))+D2
              QNY(LKK(JJ))=QNY(LKK(JJ))+D1
C             IDST(LKK(JJ))=IDST(LKK(JJ))+1     ! counter for debugging
            END IF
   20     CONTINUE
        END IF
   10 CONTINUE
C     IF(IREM.NE.0)THEN
C       DO K=1,NFP
C         WRITE(*,*)'K,QNX(K),QNY(K)=',K,QNX(K),QNY(K)
C       ENDDO
C     ENDIF
C
C  Load the matrix elements in STK2, IK2 and STK3 if needed
C
      NSI1=0
      DO 100 I=1,NFP
        IP1=IFBC1(I)
C
C   the following section uses only one side of the fault, but stores
C   also the opposite normal vector for nodes on the opposite side
C   of the fault.  It does not matter which side we use.
C
        IF(IP1.GT.0) THEN
          IP2=ABS(IFBC1(IFEQV(I)))
C
C  Evaluate the cosine and sine elements stored in QNX and QNY
C  (where there are contributions from 2 elements we effectively average)
C
          ANX=QNX(I)
          ANY=QNY(I)
          QSR=DSQRT(ANX*ANX+ANY*ANY)
          IF(JFBC1(I).NE.0)QDSTM(I)=QSR
          ANX=ANX/QSR
          ANY=-ANY/QSR
C
C  Now store normal vector components in STK2 Row i1,i2,j1,j2  Col si1,si2
C
          NRI1=IP1
          NRJ1=IP1+NUP
          NRI2=IP2
          NRJ2=IP2+NUP
          NSI1=NSI1+1
          NSI2=NFPO2+NSI1
C
C -----------------------------------------------------------
C    (anx,any) pairs go into 4 rows in the STK2 matrix
C    row i1  col si1  +nx  col si2  -ny  node ip1 x component
C    row i2  col si1  -nx  col si2  +ny  node ip2 x component
C    row j1  col si2  +ny  col si1  +nx  node ip1 y component
C    row j2  col si2  -ny  col si1  -nx  node ip2 y component
C -----------------------------------------------------------
C For the transpose:
C    Row si1 - continuity condition on normal velocity
C    Row si2 - continuity condition on tangential velocity
C -----------------------------------------------------------
C
          CALL MATPN2(NRI1,NSI1, ANX,NSI2,-ANY,STK2,IK2,NK2)
          CALL MATPN2(NRI2,NSI1,-ANX,NSI2, ANY,STK2,IK2,NK2)
          CALL MATPN2(NRJ1,NSI1, ANY,NSI2, ANX,STK2,IK2,NK2)
          CALL MATPN2(NRJ2,NSI1,-ANY,NSI2,-ANX,STK2,IK2,NK2)
          IF(JFBC1(I).EQ.0)THEN
            AAA=(QDST1(I)+QDST2(I))/6.0
          ELSE
            AAA=(2.0*QDSTM(I))/3.0
          ENDIF
C
C   amendments to QLOAD now follow for fault-type conditions
C   first case: friction condition on normal or tangential fault slip
C   QLOAD (and STK3) for faults comprises NFP/2 entries for normal 
C   components of traction followed by NFP/2 entries for tangential component.
C   IBCTYP = 11: fault locked, IBCTYP = 21: friction coefficient applied
C
          IF(IBCTYP(IFBC2(I)).EQ.21)THEN
            IF(QBND(IFBC2(I)).NE.0.0)THEN
              STK3(NSI1)=-1.0/(AAA*QBND(IFBC2(I)))
            ELSE
              STK3(NSI1)=BIG
            ENDIF
          ENDIF
          IF(IBCTYP(IFBC2(I)+NBP).EQ.21)THEN
            IF(QBND(IFBC2(I)+NBP).NE.0.0)THEN
              STK3(NSI2)=-1.0/(AAA*QBND(IFBC2(I)+NBP))
            ELSE
              STK3(NSI2)=BIG
            ENDIF
          ENDIF
C
C  For Periodic Boundary conditions: with optional velocity difference
C  intended for use e.g. in extensional environment.  If no velocity
C  difference these QLOAD elements remain zero, unless changed by next condition
C
          IF(I.GT.NFPF3)THEN             ! only PB nodes
            IF(IOFF.EQ.1) THEN
              VNO=VELXO*ANX+VELYO*ANY    ! needs checking
              VTO=VELYO*ANX-VELXO*ANY
              QLOAD(NSI1)=VNO
              QLOAD(NSI2)=VTO
            END IF
          ENDIF
C
C  if the fault is unlocked, set the tangential stress to YLDSTR
C  with the appropriate sign   Setting of QLOAD to move to B.C.s ?
C        *** under construction  ***
          ISGN=-1
C         IF(IBCTYP(IFBC2(I)).EQ.13.OR.IBCTYP(IFBC2(I)).EQ.14) THEN
C           TSTR=ISGN*Y****R/AAA
C           STK3(NSI2)=BIG
C           QLOAD(NSI2+NUP2)=BIG*TSTR
C         ENDIF
        ENDIF       ! if (IP1 > 0)
  100 CONTINUE
C     IF(IREM.NE.0)THEN
C       WRITE(*,*)'IK2 array:NK2 =',NK2
C       CALL MITPP(IK2,NK2)
C       CALL MITPRT(IK2,2,NK2/2,NFP)
C       WRITE(*,*)'STK2 array:NK2 =',NK2
C       CALL DMATPP(STK2,NK2,LUW)
C       CALL DMATPRT(STK2,2,NK2/2,NFP,6)
C       WRITE(*,*)'STK3 array:NFP =',NFP
C       CALL DMATPP(STK3,NFP,LUW)
C       CALL DMATPRT(STK3,NFP/2,2,NFP)
C     ENDIF
      RETURN
      END
C
      SUBROUTINE MATPN2(IROW,ICOL1,VALU1,ICOL2,VALU2,STK2,IK2,NK2)
C
C     Replaces MATPL2 to put VALU1 and VALU2 into columns
C     ICOL1 and ICOL2 of the STK2 matrix
C
      DOUBLE PRECISION STK2,VALU1,VALU2
      DIMENSION STK2(NK2),IK2(NK2)
C
      KC1=(IROW-1)*2+1
      KC2=KC1+1
      IK2(KC1)=ICOL1
      STK2(KC1)=VALU1
      IK2(KC2)=ICOL2
      STK2(KC2)=VALU2
      RETURN
      END

C
      SUBROUTINE MATPL2(IROW,ICOL,VALU,STK2,IK2,NK2,K2BW,LUW,IERR)
C
C Places the VALU into the STK2 matrix at row IROW and column ICOL
C
      DOUBLE PRECISION STK2,VALU
      DIMENSION STK2(NK2),IK2(NK2)

      KC0=(IROW-1)*K2BW
      KC=KC0
  100 KC=KC+1
      KBS=KC-KC0
C
C    Check that the available length of the row is not exceeded
C
      IF(KBS.GT.K2BW)GO TO 999
      IF(IK2(KC).LE.0)THEN
        IK2(KC)=ICOL
        STK2(KC)=VALU
      ELSE
        GO TO 100
      END IF
      RETURN
C
C    Error trap
C
  999 WRITE(*,10999)K2BW,IROW,ICOL,VALU
      IERR=3
  900 RETURN
10999 FORMAT(' ALLOWED SPACE IN STK2 EXCEEDED : K2BW=',I5,
     1',  ROW=',I5,',  COL=',I5,' VAL=',F9.3,';  INCREASE K2BW !!',//)
      END
