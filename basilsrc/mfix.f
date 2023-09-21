C*--------------------------------------------------------------------
C*    Basil / Sybil:   mfix.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE MFIX(JNONUM,EX,EY,VC,VHB,LEM,NOR,
     :                SSQ,FROT,IELFIX,IBC,
C    :                EX1,EY1,NOMOV,SSQ1,FROT1,IELF1,
     :                IVIS,ITHDI,THDINT,THDISS,
     :                NUP,NE,NBP,PNI,PLI,LSC,IERR)
C
C    routine to move nodes in order to fix a finite
C    element mesh that has become too distorted for 
C    further use.  MFIX should be called at every timestep.
C    If no elements can be moved satisfactorily, then none will
C    be moved.  If the result is not satisfactory, try adjusting
C    the total no. of moves (NMOVED) or the size of the worst
C    element list (currently 10).  The size of JNOALT is
C     2*NMOVED*7
C
      DIMENSION EX(NUP),EY(NUP),SSQ(NUP),FROT(NUP)
      DIMENSION IELFIX(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION IBC(NBP)
      DIMENSION VHB(8,NE)
      DIMENSION PNI(42),PLI(21)
      DIMENSION THDINT(7,NE),THDISS(7,NE)
C       local arrays
      DIMENSION NOMOV(NUP),JNOALT(700),IELF1(NUP)
      DIMENSION EX1(NUP),EY1(NUP),SSQ1(NUP),FROT1(NUP)
      DIMENSION SIDE2(3,10),SD2(3),NODET(3),DEFWST(10),IDEF(10)
C
      JIND=0
      NMOVED=50
      JNORSZ=NMOVED*2*7
      DO JI=1,JNORSZ
        JNOALT(JI)=0
      ENDDO
      DO JJ=1,NUP
        EX1(JJ)=EX(JJ)
        EY1(JJ)=EY(JJ)
      ENDDO
C
C    Establish a list of nodes that should not be moved
C      NOMOV(K)=1 means node K should not be moved
C
      DO K=1,NUP
        NOMOV(K)=0
      ENDDO
C
C    Node will not be moved if it is located on a boundary
C
      DO K=1,NBP
        NOMOV(IBC(K))=1
      ENDDO
C
C    If node belongs to element with anomalous viscosity, do not move
C
      DO N=1,NE
        DO K=1,7
          VIS=VC
          IF(IVIS.NE.0)VIS=VHB(K,N)
          IF(ABS(VIS-VC).GT.1.E-4)THEN
            DO J=1,3
              NOMOV(LEM(J,N))=1
            ENDDO
          END IF
        ENDDO
      ENDDO
C
C    Repeat the entire adjustment procedure NMOVED times
C
      DO 100 JK=1,NMOVED
C
C    define a deformation index for an element and determine
C    a list of the 10 worst elements, i.e. most collapsed
C    The deformation index has a maximum of 1.5 for an equilateral
C    triangle, and approaches 1.0 as the triangle approaches
C    (0.,0.,180.) or (0.,90.,90.)
C
      DO K1=1,10
        IDEF(K1)=0
        DEFWST(K1)=1.6
      ENDDO
      DO 25 N=1,NE
C
C    Calculate lengths of sides (squared) and cosines of
C    angles for the triangle elements
C
        DO K1=1,3
          K2=MOD(K1,3)+1
          K3=MOD(K1+1,3)+1
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          BY=EY(LK2)-EY(LK3)
          CX=EX(LK3)-EX(LK2)
          SD2(K1)=BY*BY + CX*CX
        ENDDO
        DEFIND=0.0
        DO K1=1,3
          K2=MOD(K1,3)+1
          K3=MOD(K1+1,3)+1
          SSS=2.0*SQRT(SD2(K2)*SD2(K3))
          CSA=(SD2(K2) + SD2(K3) - SD2(K1))/SSS
          DEFIND=DEFIND + CSA
        ENDDO
        IF(DEFIND.GE.DEFWST(10))GO TO 25
        DO K1=1,10
          IF(DEFIND.LT.DEFWST(K1))THEN
            KINS=K1
            GO TO 130
          END IF
        ENDDO
  130   IF(KINS.EQ.10)GO TO 150
        KS1=KINS+1
        DO K1=10,KS1,-1
          DEFWST(K1)=DEFWST(K1-1)
          IDEF(K1)=IDEF(K1-1)
          DO K2=1,3
            SIDE2(K2,K1)=SIDE2(K2,K1-1)
          ENDDO
        ENDDO
  150   DEFWST(KINS)=DEFIND
        IDEF(KINS)=N
        DO K1=1,3
          SIDE2(K1,KINS)=SD2(K1)
        ENDDO
   25 CONTINUE
C      WRITE(LSC,10001)IDEF,DEFWST
C10001 FORMAT('Element no.',10I5,'  has deformation index =',/,10f7.4)
C
C    Work on each bad element in order.  Try to move 2 nodes
C    and then revise the worst element list.  If problem has a 
C    symmetry axis, symmetry should be preserved, by fixing 2 elements
C    at a time
C
      IPMOV=0
      IWST=0
   29 IWST=IWST+1
C
C    If none of the elements on the list are amenable
C    to adjustment, then redo IDEF or exit
C
      IF(IWST.GT.10)THEN
        IF(IPMOV.EQ.0)GO TO 101
        IF(IPMOV.NE.0)GO TO 100
      END IF
C
C    Check whether any of the nodes of this element are 
C    on the prohibited list.  If so they don't move
C
      DO J=1,3
        NODET(J)=LEM(J,IDEF(IWST))
        IF(NOMOV(NODET(J)).NE.0)NODET(J)=0
      ENDDO
C
C     Check which side of the triangle is shortest
C
   50 S2MIN=99.9
      K1M=0
      DO J=1,3
        IF((NODET(J).NE.0).AND.(SIDE2(J,IWST).LE.S2MIN))THEN
          K1M=J
          S2MIN=SIDE2(J,IWST)
        END IF
      ENDDO
      IF(K1M.EQ.0)GO TO 29
      NOD1=LEM(K1M,IDEF(IWST))
C
C    Node NOD1 is the node to be moved. Check that thre is no
C    element that contains both NOD1 and IPMOV.  Node mustn't
C    belong to any element of which another node (IPMOV) has been 
C    readjusted since the last update of IDEF
C
      IF(IPMOV.NE.0)THEN
        DO J=1,NE
          ITT=2
          DO K=1,3
            LKJ=LEM(K,J)
            IF(LKJ.EQ.NOD1)ITT=ITT-1
            IF(LKJ.EQ.IPMOV)ITT=ITT-1
          ENDDO
          IF(ITT.EQ.0)THEN
            NODET(K1M)=0
            GO TO 50
            END IF
        ENDDO
      END IF
C
C    Evaluate node location for which triangle would be equilateral (XE,YE)
C    and move node K1 half way towards it
C
  210 K2M=MOD(K1M,3)+1
      K3M=MOD(K1M+1,3)+1
      K4M=K1M+3
      K5M=MOD(K1M,3)+4
      K1NOD=NOR(NOD1)
      K2NOD=NOR(LEM(K2M,IDEF(IWST)))
      K3NOD=NOR(LEM(K3M,IDEF(IWST)))
      X1=EX(K1NOD)
      Y1=EY(K1NOD)
      SSS=SQRT(S2MIN)
      XM=0.5*(EX(K2NOD) + EX(K3NOD))
      YM=0.5*(EY(K2NOD) + EY(K3NOD))
      CSTHET=(EY(K2NOD) - EY(K3NOD))/SSS
      SNTHET=(EX(K3NOD) - EX(K2NOD))/SSS
      R3O2=0.5*SQRT(3.0)
      XE=XM + R3O2*SSS*CSTHET
      YE=YM + R3O2*SSS*SNTHET
      EX(K1NOD)=XE
      EY(K1NOD)=YE
C
C     Check that none of the 6 elements containing this node
C     (a) have been turned inside out, or
C     (b) have a lesser deformation index than DEFWST(IWST)
C     If more than 3 attempts are made to move a node then
C     it is judged that there is no useful move to be made
C     and an attempt is made to move another node in this element
C     before trying to fix the next worst element
C
      ICHCK=0
  220 EX(K1NOD)=0.5*(EX(K1NOD) + X1)
      EY(K1NOD)=0.5*(EY(K1NOD) + Y1)
      ICHCK=ICHCK+1
      IF(ICHCK.GT.3)THEN
        EX(K1NOD)=X1
        EY(K1NOD)=Y1
        NODET(K1M)=0
        GO TO 50
      END IF
      NCNT=0
      DO 90 N=1,NE
        DO K=1,3
          IF(LEM(K,N).EQ.NOD1)THEN
C
C   (a) If triangle inverted, TRIA < 0
C
            TRIA=0.0
            DO K1=1,3
              K2=MOD(K1,3)+1
              K3=MOD(K1+1,3)+1
              LK2=NOR(LEM(K2,N))
              LK3=NOR(LEM(K3,N))
              BY=EY(LK2)-EY(LK3)
              CX=EX(LK3)-EX(LK2)
              SD2(K1)=BY*BY + CX*CX
              TRIA=TRIA + (EX(LK2)*EY(LK3) - EX(LK3)*EY(LK2))
            ENDDO    ! on K1
            IF(TRIA.LT.0)GO TO 220
C
C   (b) If deformation index is < DEFWST(IWST)
C
            DEFIND=0.0
            DO K1=1,3
              K2=MOD(K1,3)+1
              K3=MOD(K1+1,3)+1
              SSS=2.0*SQRT(SD2(K2)*SD2(K3))
              CSA=(SD2(K2) + SD2(K3) - SD2(K1))/SSS
              DEFIND=DEFIND + CSA
            ENDDO     ! on K1
            IF(DEFIND.LE.DEFWST(IWST))GO TO 220
C
C    if 6 elements have been tested, that's all
C
            NCNT=NCNT+1
            IF(NCNT.EQ.6)GO TO 91
          END IF
        ENDDO    ! on K
   90 CONTINUE
   91 JIND=JIND+1
      JNOALT(JIND)=NOD1
C      WRITE(LSC,10004)IDEF(IWST),NOD1,X1,Y1,EX(K1NOD),EY(K1NOD)
C10004 FORMAT(' Element/node =',2I5,'  old (x,y) =',2F7.4,
C     1 '  new (x,y) =',2F7.4)
      IF(IPMOV.EQ.0)THEN
        IPMOV=NOD1
        GO TO 29
      END IF
  100 CONTINUE
  101 CONTINUE
      JNONUM=JIND
C
C    Check which midpoint nodes have been moved, calculate
C    new location.  On each triangle, record node to right of 
C    moved vertex
C
      JIND=JNONUM
      DO 560 N=1,NE
        DO KK=1,3
          NKK=LEM(KK,N)
          DO I=1,JNONUM
            IF(NKK.EQ.JNOALT(I))THEN
              KK3=MOD(KK+1,3)+1
              KK4=KK+3
              NKK3=LEM(KK3,N)
              NKK4=IABS(LEM(KK4,N))
              NJ1=NOR(NKK)
              NJ3=NOR(NKK3)
              NJ4=NOR(NKK4)
              EX(NJ4)=(EX(NJ3)+EX(NJ1))*0.5
              EY(NJ4)=(EY(NJ3)+EY(NJ1))*0.5
              IF(JIND.GE.JNORSZ)THEN
                WRITE(LSC,*)'JNOALT in MFIX is overflowing, ',
     :                      'something wrong'
                IERR=1
                GO TO 600
              END IF
              JIND=JIND+1
              JNOALT(JIND)=NKK4
              GO TO 559
            END IF
          ENDDO   !  on I
  559     CONTINUE
        ENDDO     !  on KK
  560 CONTINUE
  600 JNONUM=JIND
      WRITE(LSC,10005)JNONUM
10005 FORMAT(I6,' nodes have been moved by routine MFIX')
      IF(JNONUM.NE.0.AND.IERR.EQ.0.AND.ICR.GT.0)THEN
       IF(ITHDI.GT.0)
     :    CALL REFIX(JNONUM,EX,EY,LEM,NOR,EX1,EY1,JNOALT,
     :               THDINT,THDISS,NUP,NE,PNI,PLI,LSC,IERR)
        IF(ICR.GT.0)
     :    CALL REFIT(JNONUM,EX,EY,LEM,NOR,SSQ,FROT,IELFIX,
     :                            EX1,EY1,JNOALT,
     :                            SSQ1,FROT1,IELF1,
     :                            NUP,NE,LSC,IERR)
        END IF
      RETURN
      END
      SUBROUTINE REFIX(JNONUM,EX,EY,LEM,NOR,EX1,EY1,JNOALT,
     :                 THDINT,THDIN1,NUP,NE,PNI,PLI,LSC,IERR)
C
C    to be called following MFIX to re-interpolate values defined
C    on the interpolation points, for elements in which any node has
C    been moved.  Values on nodes are handled by REFIT.
C
      PARAMETER (NPLIM=400)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION EX1(NUP),EY1(NUP)
      DIMENSION PNI(42),PLI(21)
C
C   THDIN1 uses same array as THDISS for two distinct purposes
C   in each case usage is temporary and two usages should not overlap
C   thereby avoiding declaration of yet another array.
C
      DIMENSION THDINT(7,NE),THDIN1(7,NE)
C      local arrays
      DIMENSION A0(3),BY(3),CX(3),COOL(3),XPT(3),YPT(3)
      DIMENSION QOOL(6)
      DIMENSION JNOALT(700)
      DIMENSION NELMOV(NPLIM),ELEC(3,NPLIM)
      DIMENSION ARR(7,6),VINT(6)
C
C   the following matrix is the 7x6 matrix that projects from the
C   7 interpolation points to a best fit set of 6 nodal values
C
      SAVE ARR
      DATA ARR/1.9743924,0.1434449,0.1434449,-0.4126757,-0.4126757,
     : 0.2563768,-0.6923076,0.1434449,1.9743924,0.1434449,0.2563768,
     :-0.4126757,-0.4126757,-0.6923076,0.1434449,0.1434449,1.9743924,
     :-0.4126757,0.2563768,-0.4126757,-0.6923076,0.09812581,0.1562660,
     :0.09812581,1.0345624,-0.2822974,-0.2822974,0.1775148,0.09812581,
     :0.09812581,0.1562660,-0.2822974,1.0345624,-0.2822974,0.1775148,
     :0.1562660,0.09812581,0.09812581,-0.2822974,-0.2822974,1.0345624,
     :0.1775148/
C
      EPS=1.E-4
      DO JJ=1,NE
        DO KK=1,7
          THDIN1(KK,JJ)=THDINT(KK,JJ)
        ENDDO
      ENDDO
      RATMAX=0.0
C
C      WRITE(LSC,*)'JNOALT matrix follows, JNONUM =',JNONUM
C     NROW=JNONUM/10
C     IF(NROW*10.LT.JNONUM)NROW=NROW+1
C      IF(JNONUM.NE.0)CALL MITPRT(JNOALT,10,NROW,NROW*10,LSC)
C
C    Make a list of elements affected by node move operations
C    and determine the X, Y limits of the new element location
C
      NP=0
      DO NPP=1,NPLIM
        NELMOV(NPP)=0
      ENDDO
      DO 60 N=1,NE
      DO 15 KK=1,3
        NKK=LEM(KK,N)
        DO I=1,JNONUM
          IF(NKK.EQ.JNOALT(I))THEN
C
C   element N is affected by remeshing.  save its centroid and size
C   (radius of circumscribed circle) in ELEC
C
            NP=NP+1
            IF(NP.GT.NPLIM)THEN
              WRITE(LSC,*)'parameter NPLIM in REFIT should be increased'
              STOP
            END IF
            NORA=NOR(LEM(1,N))
            XCN=EX(NORA)
            YCN=EY(NORA)
            DO KKK=2,3
              NORA=NOR(LEM(KKK,N))
              XCN=XCN+EX(NORA)
              YCN=YCN+EY(NORA)
            ENDDO
            XCN=XCN/3.0
            YCN=YCN/3.0
            NORA=NOR(LEM(1,N))
            XDIF=EX(NORA)-XCN
            YDIF=EY(NORA)-YCN
            REL2M=XDIF*XDIF+YDIF*YDIF
            DO KKK=2,3
              NORA=NOR(LEM(KKK,N))
              XDIF=EX(NORA)-XELCN
              YDIF=EY(NORA)-YELCN
              REL2=XDIF*XDIF+YDIF*YDIF
              IF(REL2.GT.REL2M)REL2M=REL2
            ENDDO
            NELMOV(NP)=N
            ELEC(1,NP)=XCN
            ELEC(2,NP)=YCN
            ELEC(3,NP)=SQRT(REL2M)
            GO TO 60
          END IF
        ENDDO       ! on I
   15 CONTINUE
   60 CONTINUE
      WRITE(LSC,*)NP,' elements are affected by node move operations'
C
C   work through the list of elements again, to check which of the
C   (new) nodes or interpolation points is potentially inside the (old)
C    element (using centroid coordinates and radius of circumscribed circle)
C
      DO 160 N=1,NE
        NORA=NOR(LEM(1,N))
        XELCN=EX1(NORA)
        YELCN=EX1(NORA)
        DO KKK=2,3
          NORA=NOR(LEM(KKK,N))
          XELCN=XELCN+EX1(NORA)
          YELCN=YELCN+EY1(NORA)
        ENDDO
        XELCN=XELCN/3.0
        YELCN=YELCN/3.0
        NORA=NOR(LEM(1,N))
        XDIF=EX1(NORA)-XELCN
        YDIF=EY1(NORA)-YELCN
        REL2M=XDIF*XDIF+YDIF*YDIF
        DO KKK=2,3
          NORA=NOR(LEM(KKK,N))
          XDIF=EX1(NORA)-XELCN
          YDIF=EY1(NORA)-YELCN
          REL2=XDIF*XDIF+YDIF*YDIF
          IF(REL2.GT.REL2M)REL2M=REL2
        ENDDO
C
C   check for possible overlap with list of moved elements
C   aim is to locate points associated with element NNEW (new)
C    inside element N (old)
C
      DO 140 K=1,NP
        NNEW=NELMOV(K)
        DISTX=XELCN-ELEC(1,K)
        DISTY=YELCN-ELEC(2,K)
        DIST=SQRT(DISTX*DISTX+DISTY*DISTY)
C
C    if no overlap of domain go to next element in list
C
        IF(DIST.GT.(SQRT(REL2M)+ELEC(3,K)))GO TO 140
C
C    if possible overlap, get (new) vertex node positions
C
        DO KX=1,3
          NORA=NOR(LEM(KX,NNEW))
          XPT(KX)=EX(NORA)
          YPT(KX)=EY(NORA)
        ENDDO
C
C    then get interpolation coefficients based on (old) coordinates
C
        DO 120 KV=1,3
          K2=MOD(KV,3)+1
          K3=MOD(KV+1,3)+1
          LKV=NOR(LEM(KV,N))
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          X1=EX1(LKV)
          X2=EX1(LK2)
          X3=EX1(LK3)
          Y1=EY1(LKV)
          Y2=EY1(LK2)
          Y3=EY1(LK3)
          A0(KV)=X2*Y3 - X3*Y2
          BY(KV)=Y2    - Y3
          CX(KV)=X3    - X2
  120   CONTINUE
C
C    TRI is twice the area of the triangle element
C
        TRI=X2*Y3 - X3*Y2 + X3*Y1 - X1*Y3 + X1*Y2 - X2*Y1
C
C     For (new) interpolation points = 1,7
C
        DO 110 KP1=1,7
          KP2=KP1+7
          KP3=KP2+7
          XP=XPT(1)*PLI(KP1)+XPT(2)*PLI(KP2)+XPT(3)*PLI(KP3)
          YP=YPT(1)*PLI(KP1)+YPT(2)*PLI(KP2)+YPT(3)*PLI(KP3)
C
C     get the (old) local ccordinates
C
          NODIN=1
          DO 109 KPX=1,3
            COOL(KPX)=(A0(KPX) + XP*BY(KPX) + YP*CX(KPX))/TRI
            IF((COOL(KPX).GT.1.0+EPS).OR.(COOL(KPX).LT.-EPS))NODIN=0
  109     CONTINUE
C
C     If (new) interpolation point is in (old) element, interpolate
C     to get (new) value using (old) interpolation function
C
          IF(NODIN.EQ.1)THEN
C
C     compute the effective nodal values of the function to be
C     interpolated (using old element values)
C
            DO KE=1,6
              VINT(KE)=0.0
              DO KPP=1,7
                VINT(KE)=VINT(KE)+ARR(KPP,KE)*THDIN1(KPP,N)
              ENDDO
            ENDDO
C
C     compute the (old) quadratic interpolation function values
C
            DO KE=1,3
              QOOL(KE)=COOL(KE)*(2.0*COOL(KE)-1.0)
            ENDDO
            DO KE=4,6
              KPF=KE-3
              KPB=MOD(KE+1,3)+1
              QOOL(KE)=4.0*COOL(KPF)*COOL(KPB)
            ENDDO
C
C   the interpolated fuction value then is inner product of VINT and QOOL
C
            SUM=0.0
            DO KE=1,6
              SUM=SUM+QOOL(KE)*VINT(KE)
            ENDDO
            THDINT(KP1,NNEW)=SUM
C
      RAT=ALOG(THDINT(KP1,NNEW)/THDIN1(KP1,NNEW))
      IF(ABS(RAT).GT.RATMAX)RATMAX=ABS(RAT)
      IF(ABS(RAT).GT.0.4)THEN
      WRITE(LSC,10020)NNEW,KP1,XP,YP,THDINT(KP1,NNEW),THDIN1(KP1,NNEW)
10020 FORMAT(' Node/element ',2I5,' (x,y) =',2F9.5,/,
     1' old/new THDINT',2G12.5)
      END IF
C
          END IF
  110   CONTINUE
C
  140 CONTINUE
  160 CONTINUE
      WRITE(LSC,*)'Greatest change to THDINT is factor ',EXP(RATMAX)
      RETURN
      END
      SUBROUTINE REFIT(JNONUM,EX,EY,LEM,NOR,SSQ,FROT,IELFIX,
     :                 EX1,EY1,JNOALT,SSQ1,FROT1,IELF1,
     :                 NUP,NE,LSC,IERR)
C
C    to be called following MFIX to re-interpolate crustal thickness
C    values etc., for any nodes that have been moved.  Values on
C    interpolation points are handled by REFIX
C
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION EX1(NUP),EY1(NUP)
      DIMENSION SSQ(NUP),FROT(NUP)
      DIMENSION IELFIX(NUP)
      DIMENSION JNOALT(700),IELF1(NUP)
C      local arrays
      DIMENSION SSQ1(NUP),FROT1(NUP)
      DIMENSION A0(3),BY(3),CX(3),COOL(3)
C
      EPS=1.E-4
      IF (ISSQACTIVE.EQ.1) THEN
        DO JJ=1,NUP
          IELF1(JJ)=IELFIX(JJ)
          SSQ1(JJ)=SSQ(JJ)
        ENDDO
      END IF
      IF (IROTACTIVE.EQ.1) THEN
        DO JJ=1,NUP
          FROT1(JJ)=FROT(JJ)
        ENDDO
      END IF
C      WRITE(LSC,*)'JNOALT matrix follows, JNONUM =',JNONUM
C     NROW=JNONUM/10
C     IF(NROW*10.LT.JNONUM)NROW=NROW+1
C      IF(JNONUM.NE.0)CALL MITPRT(JNOALT,10,NROW,NROW*10,LSC)
C
C    For each element, check whether any of its nodes has been moved
C
      DO 60 N=1,NE
        DO KK=1,6
          NKK=IABS(LEM(KK,N))
          DO I=1,JNONUM
            IF(NKK.EQ.JNOALT(I))GO TO 16
          ENDDO      ! on I
        ENDDO        ! on KK
        GO TO 60
C
C     If one node in this element has been moved, check all
C     moved nodes, to see whether they are now within this element
C
   16 CONTINUE
      DO 59 KJN=1,JNONUM
      NKJ=JNOALT(KJN)
      IF(NKJ.NE.0)THEN
C
C     First calculate the geometrical coefficients for (old) N
C
        DO K=1,3
          K2=MOD(K,3)+1
          K3=MOD(K+1,3)+1
          LK=NOR(LEM(K,N))
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          X1=EX1(LK)
          X2=EX1(LK2)
          X3=EX1(LK3)
          Y1=EY1(LK)
          Y2=EY1(LK2)
          Y3=EY1(LK3)
          A0(K)=X2*Y3 - X3*Y2
          BY(K)=Y2    - Y3
          CX(K)=X3    - X2
        ENDDO        ! on K
C
C    TRI is twice the area of the triangle element
C
        TRI=X2*Y3 - X3*Y2 + X3*Y1 - X1*Y3 + X1*Y2 - X2*Y1
C
C     Now calculate the natural coordinates of (new) NKJ
C
        LKK=NOR(NKJ)
        XP=EX(LKK)
        YP=EY(LKK)
        DO K=1,3
          CNL=(A0(K) + XP*BY(K) + YP*CX(K))/TRI
          IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 59
          COOL(K)=CNL
        ENDDO          ! on K
C
C     If we reach here the point is within the triangle, so
C     get the interpolated values
C     For IELFIX, use the interpolation function and then round
C     up to 1 or down to 0 with 0.5 as the criterion
C
        SI=0.0
        FI=0.0
        DM=0.0
        SK=DFLTSSQ
        FK=0.0
        IELF=0
        DO K=1,3
          LK=LEM(K,N)
          NLK=NOR(LK)
          IF (ISSQACTIVE.EQ.1) THEN
            SK=SSQ1(LK)
            IELF=IELF1(LK)
          END IF
          IF (IROTACTIVE.EQ.1) FK=FROT1(LK)
          K2=MOD(K,3)+1
          K3=MOD(K+1,3)+1
          CNK=COOL(K)*(COOL(K) - COOL(K2) - COOL(K3))
          DM=DM + FLOAT(IELF)*CNK
          SI=SI + SK*CNK
          FI=FI + FK*CNK
        ENDDO          ! on K
        DO K=4,6
          LK=IABS(LEM(K,N))
          IF (ISSQACTIVE.EQ.1) THEN
            SK=SSQ1(LK)
            IELF=IELF1(LK)
          END IF
          IF (IROTACTIVE.EQ.1) FK=FROT1(LK)
          K2=K-3
          K3=MOD(K+1,3)+1
          CNK=4.0*COOL(K2)*COOL(K3)
          DM=DM + FLOAT(IELF)*CNK
          SI=SI + SK*CNK
          FI=FI + FK*CNK
        ENDDO          ! on K
        IF (ISSQACTIVE.EQ.1) THEN
          SSQ(NKJ)=SI
          IF(DM.LT.0.0)DM=0.0
          IF(DM.GT.1.0)DM=1.0
          IELFIX(NKJ)=INT(DM+0.5)
        END IF
        IF (IROTACTIVE.EQ.1) FROT(NKJ)=FI
C
C     WRITE(LSC,10020)NKJ,N,XP,YP,SSQ1(NKJ),SI,IELF1(NKJ),DM,
C    1 IELFIX(NKJ)
C10020 FORMAT(' Node/element ',2I5,' (x,y) =',2F9.5,/,
C    1' old/new S',2G12.5,'  old/new ELF',I10,F8.5,I10)
C
C    wipe the node no. in JNOALT after interpolating, there may
C    be more than one occurrence.  Exit if JNOALT empty.
C
        IST=0
        DO I=1,JNONUM
          IF(NKJ.EQ.JNOALT(I))JNOALT(I)=0
          IF(JNOALT(I).NE.0)IST=1
        ENDDO
      IF(IST.EQ.0)RETURN
C
C   go back to check other nodes in this element and other
C    elements
C
        ENDIF          !  (NKJ != 0)
   59 CONTINUE
   60 CONTINUE
C
C   the above strategy will get most of the moved nodes in
C   minimal time, but may miss a node if it has been moved
C   right outside original domain of one of its elements,
C   following multiple moves.  For the time being, a check 
C   is put in to ensure this hasn't happened.  Code may need
C   to be extended later.
C
      IST=0
      DO 70 I=1,JNONUM
      IF(JNOALT(I).NE.0)THEN
      WRITE(LSC,10010)JNOALT(I)
      IST=1
10010 FORMAT(' Node no',I5,' has not been reinterpolated')
      END IF
   70 CONTINUE
      IF(IST.NE.0)IERR=1
      RETURN
      END
