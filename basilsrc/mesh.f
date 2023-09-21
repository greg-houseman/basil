C*--------------------------------------------------------------------
C*    Basil / Sybil:   mesh.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE MESH(NX,NY,XLEN,YLEN,EX,EY,LEM,NUP,NE,NN,
     :                 IFLT,IMSHP,LUW,LSC,IERR)
C
C    A regular mesh of triangular elements is defined in a
C   box which is initially rectangular (XLEN*YLEN).  EX and EY
C   are the X and Y coordinates of the nodes and LEM gives the
C   nodes associated with each element (in anticlockwise order
C   from the bottom left : 1, 2 & 3).  The nodes are numbered
C   from 1 in the bottom left hand corner, in columns. ITP = 0
C   applies to original mesh geometry, ITP = 1 refers to rows swapped.
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
C
      ITP=IMSHP+1
      DY=YLEN/FLOAT(NY)
      DX=XLEN/FLOAT(NX)
      NX1=NX+1
      NY1=NY+1
C
C     JEL is an index for the elements
C
      JEL=1
C
C     First the nodes and elements on the left hand margin
C
      DO 10 J=1,NY
      EX(J)=0.
      EY(J)=FLOAT(J-1)*DY
      JD=MOD(J+ITP+1,2)
      LEM(1,JEL)=J
      LEM(2,JEL)=J+NY1+JD
      LEM(3,JEL)=J+1
      JEL=JEL+1
   10 CONTINUE
      EX(NY1)=0.
      EY(NY1)=YLEN
C
C    Now work up the columns, two elements at a time
C
      NLC=NX1*NY1
      HALF=0.5*FLOAT(MOD(ITP+1,2))
      DO 50 I=1,NX
      DO 30 J=1,NY
C
C    NL is the number of the next node in a zig-zag path,
C    up the column
C
      NL=I*NY1+J
      JD=MOD(J+ITP,2)
      EX(NL)=DX*(FLOAT(I)-0.5*FLOAT(JD))
      EY(NL)=FLOAT(J-1)*DY
      JELN=JEL+NY
C
C    For the original row ordering: first a pair of elements on the
C    odd-number row, then a pair on the even number row
C
      IF(ITP.EQ.0)THEN
        IF(JD.EQ.1)THEN
          LEM(1,JEL)=NL
          LEM(2,JEL)=NL+1
          LEM(3,JEL)=NL+1-NY1
          LEM(1,JELN)=NL
          IF(I.EQ.NX)THEN
            NLC=NLC+1
            LEM(2,JELN)=NLC
            EX(NLC)=XLEN
            EY(NLC)=FLOAT(J-1)*DY
          ELSE
            LEM(2,JELN)=NL+NY1
          END IF
          LEM(3,JELN)=NL+1
C
C    For the alternative orientation of the pair
C
        ELSE IF(JD.EQ.0)THEN
          LEM(1,JEL)=NL-NY1
          LEM(2,JEL)=NL
          LEM(3,JEL)=NL+1
          LEM(1,JELN)=NL
          IF(I.EQ.NX)THEN
            LEM(2,JELN)=NLC+1
            IF(J.EQ.NY)THEN
              NLC=NLC+1
              EX(NLC)=XLEN
              EY(NLC)=YLEN
            END IF
          ELSE
            LEM(2,JELN)=NL+1+NY1
          END IF
          LEM(3,JELN)=NL+1
        END IF
C
C    Now if the row order is swapped: again a pair of elements on the
C    odd number row, then a pair on the even number row
C
      ELSE IF(ITP.EQ.1)THEN
        IF(JD.EQ.0)THEN
          LEM(1,JEL)=NL-NY1
          LEM(2,JEL)=NL
          LEM(3,JEL)=NL+1
          LEM(1,JELN)=NL
          IF(I.EQ.NX)THEN
            NLC=NLC+1
            LEM(2,JELN)=NLC
          ELSE
            LEM(2,JELN)=NL+NY1+1
          END IF
          LEM(3,JELN)=NL+1
C
C    For the alternative orientation of the pair
C
          ELSE IF(JD.EQ.1)THEN
            LEM(1,JEL)=NL
            LEM(2,JEL)=NL+1
            LEM(3,JEL)=NL-NY1+1
            LEM(1,JELN)=NL
            IF(I.EQ.NX)THEN
              LEM(2,JELN)=NLC
            EX(NLC)=XLEN
            EY(NLC)=FLOAT(J-1)*DY
            ELSE
              LEM(2,JELN)=NL+NY1
            END IF
            LEM(3,JELN)=NL+1
          END IF
        END IF
C
      JEL=JEL+1
   30 CONTINUE
C
      NT=NL+1
      EX(NT)=DX*(FLOAT(I)-HALF)
      EY(NT)=YLEN
      JEL=JELN+1
   50 CONTINUE
      JEL=JEL-1
C
C    If there is an internal fault, duplicate other side of fault
C
      IF(IFLT.EQ.2) THEN
          DO 100 II=(NN/2+1),NN
                IJ=II-NN/2
                EX(II)=EX(IJ)+EX(NN/2)
                EY(II)=EY(IJ)
  100     CONTINUE
          DO 110 JJ=(NE/2+1),NE
                JI=JJ-NE/2
                LEM(1,JJ)=LEM(1,JI)+NN/2
                LEM(2,JJ)=LEM(2,JI)+NN/2
                LEM(3,JJ)=LEM(3,JI)+NN/2
  110     CONTINUE
          JEL=JEL*2
          NLC=NLC*2
      END IF
C
C    Check the number of elements and nodes
C    note that NN differs by one depending on whether ITP = 0 or 1
C
      IF((JEL.NE.NE).OR.(NLC.NE.NN))THEN
      WRITE(LUW,10050)JEL,NE,NLC,NN
      WRITE(LSC,10050)JEL,NE,NLC,NN
10050 FORMAT(' Error in Subroutine MESH: JEL, NE, NLC, NN =',4I6)
      IERR = 1
      END IF
C
      RETURN
      END
      SUBROUTINE MPNODE(NX,NY,EX,EY,LEM,NUP,NE,NN,IFLT,IMSHP,LUW,LSC)
C
C   To assign node numbers to the midpoint nodes used
C   in the velocity equations. The values are stored in
C   LEM(4-6,NE).  ITP = 0 refers to the original mesh geometry
C   ITP = 1 refers to rows swapped.
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
C
C    Commence with left hand column
C    The nodes are numbered sequentially in vertical columns
C
C     WRITE(*,*)'in MPNODE, NN, NMP, NE =',NN,NMP,NE
      ITP=IMSHP+1
      NMP=NUP-NN
      NM4=NN
      NM5=NM4+NY+1-ITP
      NM6=NM5+NY/2
      DO 10 J=1,NY
      JEL=J
      JD=MOD(J+ITP,2)
      IM5=6-JD
      IM6=5+JD
      NM4=NM4+1
      IF(JD.EQ.0)NM5=NM5+1
      NM6=NM6+1
      LEM(4,JEL)=NM4
      LEM(IM5,JEL)=NM5
      LEM(IM6,JEL)=NM6
   10 CONTINUE
C
C    working up the columns and moving from left to right, two
C    columns of elements are processed at a time.  Within an element
C    local node number convention is anticlockwise ; vertices 1-3
C    then midpoints 4-6, wiith 4 adjacent to 3 in the positive direction.
C    Global node numbers are at this stage organised with all vertex numbers
C    listed first, then vertical columns of midpoint nodes of length
C    for ITP=0: NY, (NY/2+1), NY, (NY/2), NY, (NY/2+1), ...
C    for ITP=1: NY, (NY/2), NY, (NY/2+1), NY, (NY/2), ....
C
      DO 50 I=1,NX
        NM4=NN+NY*3/2+1-ITP+(I-1)*(3*NY+1)
        NM5=NM4+NY+ITP
        NM6=NM5+NY/2
        NM4N=NM6
        NM5N=NM4N+NY+1-ITP
        NM6N=NM5N+NY/2
        JEL=(2*I-1)*NY
        JELN=JEL+NY
C
        DO J=1,NY
          JD=MOD(J+1+ITP,2)
          IM5=6-JD
          IM6=5+JD
          NM4=NM4+1
          IF(JD.EQ.0)NM5=NM5+1
          NM6=NM6+1
          JEL=JEL+1
          LEM(4,JEL)=NM4
          LEM(IM5,JEL)=NM5
          LEM(IM6,JEL)=NM6
          JELN=JELN+1
          NM4N=NM4N+1
          IF(JD.EQ.1)NM5N=NM5N+1
          NM6N=NM6N+1
          LEM(4,JELN)=NM4N
          LEM(IM6,JELN)=NM5N
          LEM(IM5,JELN)=NM6N
        ENDDO     ! on J
   50 CONTINUE    ! on I
      NT=NM6N
C
C    If there is an internal fault, duplicate other side of fault
C
      IF(IFLT.EQ.2) THEN
           DO 110 JJ=(NE/2+1),NE
             JI=JJ-NE/2
             LEM(4,JJ)=LEM(4,JI)+NMP/2
                 LEM(5,JJ)=LEM(5,JI)+NMP/2
             LEM(6,JJ)=LEM(6,JI)+NMP/2
  110      CONTINUE
           NT=NT*2
      END IF
C
C    Check the number of mid-point nodes
C
      IF(NT.NE.NMP.AND.IMSHP.NE.0)WRITE(LUW,10002)NT,NMP
10002 FORMAT(' NO. OF MIDPOINT NODES INCORRECT :',2I6)
C
C       Calculate position of midpoint nodes
C
      DO JEL=1,NE
        NJ1=LEM(1,JEL)
        NJ2=LEM(2,JEL)
        NJ3=LEM(3,JEL)
        NJ4=LEM(4,JEL)
        NJ5=LEM(5,JEL)
        NJ6=LEM(6,JEL)
        EX(NJ4)=(EX(NJ3)+EX(NJ1))*0.5
        EX(NJ5)=(EX(NJ1)+EX(NJ2))*0.5
        EX(NJ6)=(EX(NJ2)+EX(NJ3))*0.5
        EY(NJ4)=(EY(NJ3)+EY(NJ1))*0.5
        EY(NJ5)=(EY(NJ1)+EY(NJ2))*0.5
        EY(NJ6)=(EY(NJ2)+EY(NJ3))*0.5
      ENDDO
      RETURN
      END
      SUBROUTINE NORDER(EX,EY,LEM,NOR,NUP,NE,NN,LSC,IERR)
C
C     Routine to order the nodes for insertion of the matrix
C     elements in such a way as to minimise the bandwidth
C     The nodes are ordered sequentially in columns
C     This version orders sequentially in rows by choosing
C     minimum X-coordinate from all minimum Y-coordinates.
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
C       local arrays
      DIMENSION NOR2(NUP)
C
C      WRITE(LSC,10101)
C10101 FORMAT('Reordering Nodes Begun')
      EPS=1.E-4
C
C          Zero NOR and NOR2 arrays
C
      DO 64 INN=1,NUP
      NOR(INN)=0
      NOR2(INN)=0
   64 CONTINUE
C
C       Construct the NOR and NOR2 arrays
C
      DO 70 IN=1,NUP
        XCC=EX(IN)
        YCC=EY(IN)
        IMIN=1
        DO 65 I=1,NUP
          XC=EX(I)
          YC=EY(I)
          IF(YCC.GT.YC+EPS) THEN
              IMIN=IMIN+1
          ELSE IF(YCC.GT.(YC-EPS).AND.XCC.GT.(XC+EPS))THEN
              IMIN=IMIN+1
          END IF
   65   CONTINUE
   67   IF(NOR(IMIN).EQ.0) GO TO 68
          IMIN=IMIN+1
        GO TO 67
   68   IF(IMIN.GT.NUP) GO TO 80
        NOR(IMIN)=IN
        NOR2(IN)=IMIN
        IMIN=1
   70 CONTINUE
C
C     Reorder the values in LEM
C      CALL MITPRT(LEM,6,NE,6*NE,LSC)
C
      DO I=1,NE
        DO J=1,6
          NI=LEM(J,I)
          K=NOR2(NI)
          LEM(J,I)=K
        ENDDO     ! on J
      ENDDO       ! on I
C      CALL MATPRT(EX,NUP,1,NUP,LSC)
C      CALL MATPRT(EY,NUP,1,NUP,LSC)
C     CALL MITPRT(NOR,NUP,1,NUP,LSC)
C      CALL MITPRT(NOR2,NUP,1,NUP,LSC)
C      CALL MITPRT(LEM,6,NE,6*NE,LSC)
      RETURN
   80 WRITE(LSC,10100)IMIN,IN,NUP
10100 FORMAT(' NORDER IS WRONG, IMIN =',I5,' IN =',I5,'  NUP =',I5,//)
      IERR = 1
      RETURN
      END
C
      SUBROUTINE MESH2(NX,NY,XLEN,YLEN,EX,EY,LEM,NUP,NE,NN,
     :                 IMESH,IFLT,LUW,LSC,IERR)
C
C    A regular mesh of triangular elements is defined in a
C   box which is initially square (or rectangular).  EX and EY
C   are the X and Y coordinates of the nodes and LEM gives the
C   nodes associated with each element (in anticlockwise order
C   from the bottom left : 1, 2 & 3).  The nodes are numbered
C   from 1 in the bottom left hand corner, in columns.
C   IMESH=1 |/|/|
C   IMESH=2 |\|/|
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
C
      DX=XLEN/FLOAT(NX)
      DY=YLEN/FLOAT(NY)
      NX1=NX+1
      NY1=NY+1
      NMP=NUP-NN
C
C     JEL is an index for the elements
C
      JEL=0
C
C     First the nodes
C
      DO 10 I=1,NX1
      DO 20 J=1,NY1
         NL=(I-1)*NY1+J
         EX(NL)=FLOAT(I-1)*DX
         EY(NL)=FLOAT(J-1)*DY
   20 CONTINUE
   10 CONTINUE
C
C    Now the elements
C
      DO 30 I=1,NX
      DO 40 J=1,NY
C
         JEL=JEL+1
C
C vertex nodes
C
         NL1=(I-1)*NY1+J
         NL2=NL1+1
         NL3=NL1+NY1
         NL4=NL3+1
C
C midpoint nodes
C
         NM1=NN+(I-1)*(3*NY+1)+J
         NM2=NM1+NY+J-1
         NM3=NM2+1
         NM4=NM2+2
         NM5=NM1+3*NY+1
C
         IF (IMESH.EQ.1.OR.MOD(J,2).EQ.1) THEN
C
C    the lower side
C
           LEM(1,JEL)=NL1
           LEM(2,JEL)=NL3
           LEM(3,JEL)=NL2
           LEM(4,JEL)=NM1
           LEM(5,JEL)=NM2
           LEM(6,JEL)=NM3
C
C      the upper side
C
           JELN=JEL+NY
           LEM(1,JELN)=NL3
           LEM(2,JELN)=NL4
           LEM(3,JELN)=NL2
           LEM(4,JELN)=NM3
           LEM(5,JELN)=NM5
           LEM(6,JELN)=NM4
         ELSE
C
C    the lower side
C
           LEM(1,JEL)=NL1
           LEM(2,JEL)=NL3
           LEM(3,JEL)=NL4
           LEM(4,JEL)=NM3
           LEM(5,JEL)=NM2
           LEM(6,JEL)=NM5
C
C      the upper side
C
           JELN=JEL+NY
           LEM(1,JELN)=NL4
           LEM(2,JELN)=NL2
           LEM(3,JELN)=NL1
           LEM(4,JELN)=NM3
           LEM(5,JELN)=NM4
           LEM(6,JELN)=NM1
         ENDIF
   40 CONTINUE
         JEL=JEL+NY
   30 CONTINUE
      NMT=NM5-NN
C
C    If there is an internal fault, duplicate other side of fault
C
      IF(IFLT.EQ.2) THEN
           DO 50 II=1,NN/2
              IJ=II+NN/2
C             EX(IJ)=2.0*EX(NN/2)-EX(II)
              EX(IJ)=EX(NN/2)+EX(II)
              EY(IJ)=EY(II)
   50      CONTINUE
           DO 60 II=1,NE/2
              IJ=II+NE/2
              LEM(1,IJ)=LEM(1,II)+NN/2
              LEM(2,IJ)=LEM(2,II)+NN/2
              LEM(3,IJ)=LEM(3,II)+NN/2
              LEM(4,IJ)=LEM(4,II)+NMP/2
              LEM(5,IJ)=LEM(5,II)+NMP/2
              LEM(6,IJ)=LEM(6,II)+NMP/2
   60      CONTINUE
           JEL=JEL*2
           NL=NL*2
           NMT=NMT*2
      END IF
C
C  Find coordinates for the midpoint nodes
C
      DO 120 JJ=1,NE
        NJ1=LEM(1,JJ)
        NJ2=LEM(2,JJ)
        NJ3=LEM(3,JJ)
        NJ4=LEM(4,JJ)
        NJ5=LEM(5,JJ)
        NJ6=LEM(6,JJ)  
        EX(NJ4)=(EX(NJ3)+EX(NJ1))/2
        EX(NJ5)=(EX(NJ1)+EX(NJ2))/2
        EX(NJ6)=(EX(NJ2)+EX(NJ3))/2
        EY(NJ4)=(EY(NJ3)+EY(NJ1))/2
        EY(NJ5)=(EY(NJ1)+EY(NJ2))/2
        EY(NJ6)=(EY(NJ2)+EY(NJ3))/2
  120 CONTINUE
C
C    Check the number of elements and nodes
C
      IF(JEL.NE.NE) THEN
        WRITE(LUW,10050)JEL,NE
        IERR = 1
      END IF
      IF(NL.NE.NN) THEN
        WRITE(LUW,10051)NL,NN
        IERR = 1
      END IF
      IF(NMT.NE.NMP) THEN
        WRITE(LUW,10002)NMT,NMP
        IERR = 1
      END IF
C      WRITE(LSC,10003)NUP,NN,NMP
C      DO 130 JJ=1,NUP
C         WRITE(LSC,10001)JJ,EX(JJ),EY(JJ)
C  130 CONTINUE
C      CALL MITPRT(LEM,6,NE,6*NE,LSC)
C
10050 FORMAT(' NO. OF ELEMENTS INCORRECT :',2I6)
10051 FORMAT(' NO. OF NODES INCORRECT :',2I6)
10001 FORMAT('I=',I6,' EX=',F8.4,' EY=',F8.4)
10002 FORMAT(' NO. OF MIDPOINT NODES INCORRECT :',2I6)
10003 FORMAT('NUP=',I6,' NN=',I6,' NMP=',I6)
C
      RETURN
      END
