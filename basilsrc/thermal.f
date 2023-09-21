      SUBROUTINE THSOLVE
C
C     this routine sets up and solves the equation that gives the
C     rates of change of temperature at each node
C
      CALL THSEMBL(BIG,STK,IK,NK1,KBW,NCOMP,
     :             EX,EY,LEM,NOR,IBC,IBCTYP,
     :             QBND,DIFF,NUP,NE,NBP)
      CALL MATVCT(TEMPT,RHS,STK(NUP*KBW+1),IKB,NK1)
      CALL CONJGR
C
C     save the rates of change
C
      RETURN
      END
      SUBROUTINE THSEMBL(BIG,STK,IK,NK1,KBW,NCOMP,
     :                 EX,EY,LEM,NOR,IBC,IBCTYP,
     :                 QBND,DIFF,NUP,NE,NBP)
C
C     This routine calculates and assembles the entries for the
C     matrices used in the calculation of rate of change of temperature.
C     The entries are stored in STK(NK1).  Compressed format is used
C     so the entries are indexed by the values in IK.
C     Note that this version of
C     SEMBL stores rows for use by Subroutine CONJGR
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION IBC(NBP),IBCTYP(NBP*2)
      DIMENSION QBND(NBP*2)
      SAVE W
      DIMENSION W(7)
      DIMENSION STK(NK1)
      DIMENSION IK(NK1)
      DIMENSION BY(3),CX(3),DNDP(84)
      DOUBLE PRECISION TRIA,DPEX2,DPEX3,DPEY2,DPEY3
C
      DATA W/0.12593918,0.12593918,0.12593918,0.13239415,
     10.13239415,0.13239415,0.225/
C
C     each matrix element is calculted by numerical integration
C    using a 7 point approximation from p421 in Huebner.  For
C    a quadratic integrand (as here) the integral is exact.
C    the seven points are given the preceeding (in W) weights
C     Clear the address arrays for the matrices
C
      DO 20 I=1,NK1
   20 IK(I)=0
C
C    Initialise the diagonal elements
C
      NUP2=NUP*2
      DO 15 I=1,NUP2
      KC=(I-1)*KBW+1
      IK(KC)=I
   15 STK(KC)=0.0
C
C   Cycle through the elements
C
      DO 70 N=1,NE
C
C    Set up the geometrical coefficients for the natural
C    coordinates used in the triangle elements
C
      TRIA=0.0
      DO 10 K1=1,3
      K2=MOD(K1,3)+1
      K3=MOD(K1+1,3)+1
      LK2=NOR(LEM(K2,N))
      LK3=NOR(LEM(K3,N))
      BY(K1)=EY(LK2)-EY(LK3)
      CX(K1)=EX(LK3)-EX(LK2)
      DPEX2 = EX(LK2)
      DPEX3 = EX(LK3)
      DPEY2 = EY(LK2)
      DPEY3 = EY(LK3)
   10 TRIA=TRIA + (DPEX2*DPEY3-DPEX3*DPEY2)
      TRIA=2.0*TRIA
C
C   TRIA is 4 times the area of the triangle element
C
      CALL DNCOM(0,BY,CX,DNDP)
C
C    For each element calculate the 21 independent entries to the
C    thermal matrices M and K.
C
      DO 70 I=1,6
      K1=LEM(I,N)
      IF(K1.LT.0)K1=-K1
      DO 70 J=1,6
      K2=LEM(J,N)
      IF(K2.LT.0)K2=-K2
C
C    The value of the integrand at each of the seven points
C    is calculated, and summed with the appropriate weight
C    each gradient term DNDP is pre-multiplied by a factor of TRIA/2
C    In integrating, the sum is therefore multiplied by TRIA/4 and 
C    divided by TRIA*TRIA/4.
C
      IF(K2.GE.K1)THEN
        SUM1=0.0
        DO 30 K=1,7
          KIN=(I-1)*14 + (K-1)*2 + 1
          KJN=(J-1)*14 + (K-1)*2 + 1
          DNDX2=DNDP(KIN)*DNDP(KJN)
          DNDY2=DNDP(KIN+1)*DNDP(KJN+1)
          AIJ = DNDX2 + DNDY2
          SUM1=SUM1 + W(K)*AIJ
   30   CONTINUE
        SUM1=DIFF*SUM1/TRIA
C
C     Now add in the values of the M matrix
C
        IF(K1.LE.3)THEN
          IF(K2.LE.3)THEN
            IF(K1.EQ.K2)THEN
C   diagonal element, vertex nodes
              TRM=6.0          
            ELSE
C   vertex to vertex, off-diagonal
              TRM=-1.0
            END IF
          ELSE                         !i.e., if (K2.GT.3)
            IF((MOD(K2,3)+1).EQ.K1)THEN
C   vertex to opposite midpoint
              TRM=-4.0
            ELSE
C   vertex to adjacent midpoint node
              TRM=0.0
            END IF
          END IF                       !end of (K2.LE. or .GT.3)
        ELSE                           !i.e., if (K1.GT.3)
          IF(K1.EQ.K2)THEN
C  diagonal element, midpoint nodes
            TRM=32.0
          ELSE
C  midpoint to midpoint, off-diagonal
            TRM=16.0
          END IF
        END IF                         !end of (K1.LE. or .GT.3)
C
C    scale and store the matrix entry
C
        TRM=TRM*TRIA/720.0
C
C     The matrices are assembled and stored in compressed mode,
C            ++++++++++  ROW WISE ++++++++
C     SUM1 goes into row K1, col K2, and into row K2, col K1.
C     TRM is stored in the first NUP rows, SUM1 in the second NUP rows.
C     Only the upper diagonal half is calculated and stored.
C
        CALL MATPAD(K1,K2,TRM,STK,IK,NK1,KBW)
        CALL MATPAD(K1,K2,SUM1,STK(NUP*KBW+1),IK,NK1,KBW)
      END IF                           !end of (K2.GE.K1)
C
   70 CONTINUE
C
C    set the diagonal element for nodes with fixed temperature
C
      DO 80 II=1,NBP
      NR=IBC(II)
      NSX=(NR-1)*KBW+1
      NSY=(NUP+NR-1)*KBW+1
      IY=II+NBP
      IF(IBCTYP(II).EQ.0.OR.IBCTYP(II).EQ.10)STK(NSX)=BIG
      IF(IBCTYP(IY).EQ.0.OR.IBCTYP(IY).EQ.10)STK(NSY)=BIG
      IF(IBCTYP(II).EQ.2)STK(NSX)=STK(NSX)+QBND(II)
      IF(IBCTYP(IY).EQ.2)STK(NSY)=STK(NSY)+QBND(IY)
   80 CONTINUE
C
C     WRITE(LUW,*)'IK matrix indices follow'
C     CALL MITPRT(IK,KBW,NUP2P,KBW*NUP2P,LUW)
C     WRITE(LUW,*)'STK matrix follows'
C     CALL MATPRT(STK,KBW,NUP2P,KBW*NUP2P,LUW)
      RETURN
      END
