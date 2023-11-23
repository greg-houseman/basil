C*--------------------------------------------------------------------
C*    Basil / Sybil:   crust.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE CRUST(DT,SSQTS,FROTTS,SSQ,FROT,UVP,EX,EY,
     :                 PNI,PLI,LEM,NOR,NUP,NE,NROWS,NCOMP)
C
C    Routine to calculate the rate of change in crustal thickness
C    on each mode and update the crustal thickness. The vorticity
C    is also calculated and integrated to give the finite rotation
C
      DOUBLE PRECISION XWRK
      DOUBLE PRECISION PNI,PLI,BY,CX,DNDP,YY,XX
      DOUBLE PRECISION TRIA,DNDPX,DNDPY,UI,VI,UIK,VIK,XXK,DIV,OMEG
      DOUBLE PRECISION YCOTK,ANGL,SI,FI,DLSDX,DLSDY,DFRDX,DFRDY
      DOUBLE PRECISION DB2,DA2,SANGL,CSANGL,PI,A2DT,UPX,UPY
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      DIMENSION SSQ(NUP),SSQTS(NUP),FROT(NUP),FROTTS(NUP)
      DIMENSION UVP(NROWS)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION PNI(7,6),PLI(7,3)
      DIMENSION LEM(6,NE)
      DIMENSION XWRK(NUP*3)
      DIMENSION BY(3),CX(3),DNDP(84)
      DIMENSION YY(6),XX(6)
      SAVE PI
      DATA PI/3.14159265358979/
C
      DO I=1,NUP*3
        XWRK(I)=0.0
      ENDDO
C
C    Set up the geometrical coefficients for the natural
C    coordinates used in the triangle elements
C
      DO 60 N=1,NE
        DO 10 KK=1,3
          K2=MOD(KK,3)+1
          K3=MOD(KK+1,3)+1
          LK1=NOR(LEM(KK,N))
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          XX(KK)=EX(LK1)
          YY(KK)=EY(LK1)
          BY(KK)=EY(LK2)-EY(LK3)
          CX(KK)=EX(LK3)-EX(LK2)
   10   CONTINUE
C
C   TRIA is here 2 times the area of the triangle element
C
        TRIA=CX(3)*BY(2)-CX(2)*BY(3)
        DO 11 KK=4,6
          KK1=MOD(KK+1,3)+1
          KK2=MOD(KK+2,3)+1
          XX(KK)=0.5D0*(XX(KK1)+XX(KK2))
          YY(KK)=0.5D0*(YY(KK1)+YY(KK2))
   11   CONTINUE
C
C    Use argument 5 to get gradients on the boundaries
C
        CALL DNCOM(5,BY,CX,DNDP)
        DO 62 K=1,6
          LK=LEM(K,N)
          UIK=UVP(LK)
          VIK=UVP(LK+NUP)
          XXK=XX(K)
          YCOTK=1.D0/DTAN(YY(K))
          ANGL=0.5D0
C
C    ANGL is the internal angle made by the boundaries of the
C      triangle at this node
C
          K1=MOD(K-1,3)+1
          K2=MOD(K,3)+1
          K3=MOD(K+1,3)+1
          IF(K.LE.3)THEN
            DB2= BY(K2)*BY(K2) + CX(K2)*CX(K2)
            DA2= BY(K3)*BY(K3) + CX(K3)*CX(K3)
            SANGL=TRIA/DSQRT(DB2*DA2)
            CSANGL=1.D0-SANGL*SANGL
            IF(CSANGL.LT.0.D0)CSANGL=0.D0
            CSANGL=DSQRT(CSANGL)
            ANGL=0.5D0*DATAN2(SANGL,CSANGL)/PI
          END IF
          DIV=0.0
          OMEG=0.0
          DLSDX=0.0
          DLSDY=0.0
          DFRDX=0.0
          DFRDY=0.0
          SI=DFLTSSQ
          FI=0.0
C
C    Calculate the interpolated value of the divergence at node K
C     and the crustal thickness gradients for the midpoints
C
          DO 55 I=1,6
            KIN=(I-1)*14 + (K-1)*2 +1
            DNDPX=DNDP(KIN)
            DNDPY=DNDP(KIN+1)
            NI=LEM(I,N)
            UI=UVP(NI)
            VI=UVP(NI+NUP)
            DIV=DIV + UI*DNDPX + VI*DNDPY
            OMEG=OMEG + VI*DNDPX - UI*DNDPY
            IF(NCOMP.EQ.-1)THEN
              DIV=DIV + YCOTK*XXK*VI*DNDPX
              OMEG=OMEG - YCOTK*XXK*UI*DNDPX
            END IF
C
            IF(K.GE.4)THEN
              IF (ISSQACTIVE.EQ.1) SI=SSQTS(NI)
              DLSDX=DLSDX + SI*DNDPX
              DLSDY=DLSDY + SI*DNDPY
              IF (IROTACTIVE.EQ.1) FI=FROTTS(NI)
              DFRDX=DFRDX + FI*DNDPX
              DFRDY=DFRDY + FI*DNDPY
            END IF
   55     CONTINUE
C
C     Store the divergence in XWRK(1-NUP), the angle in
C     XWRK(NUP+1-2*NUP), and the vorticity in XWRK(2*NUP+1-3*NUP).
C
          IF(NCOMP.EQ.-1)THEN
            DIV=DIV + VIK*YCOTK*TRIA
            OMEG=OMEG - UIK*YCOTK*TRIA
          END IF
          A2DT=ANGL/TRIA
          XWRK(LK)=XWRK(LK) + DIV*A2DT
          LKA=LK+NUP
          XWRK(LKA)=XWRK(LKA) + ANGL
          LKV=LKA+NUP
          XWRK(LKV)=XWRK(LKV) + OMEG*A2DT/2.D0
C
C    add the Eulerian correction terms for the midpoint velocities
C
          IF(K.GE.4)THEN
            LK1=LEM(K1,N)
            LK3=LEM(K3,N)
            UPX=UVP(LK)-0.5*(UVP(LK1)+UVP(LK3))
            UPY=UVP(LKA)-0.5*(UVP(LK1+NUP)+UVP(LK3+NUP))
            XWRK(LK)=XWRK(LK) + (UPX*DLSDX + UPY*DLSDY)*A2DT
            XWRK(LKV)=XWRK(LKV) - (UPX*DFRDX + UPY*DFRDY)*A2DT*PI/180.D0
          END IF
   62   CONTINUE
   60 CONTINUE
C
C     SSQ contains the log. of the crustal thickness, which is
C     now updated. FROT contains the rotation angle (in degrees)
C      from the initial orientation
C
      DO 70 N=1,NUP
        NNU=N+NUP
        NNNU=NNU+NUP
        XWRK(N)=XWRK(N)/XWRK(NNU)
        XWRK(NNNU)=XWRK(NNNU)/XWRK(NNU)
        IF (ISSQACTIVE.EQ.1)SSQ(N)=SSQTS(N) - DT*XWRK(N)
        IF (IROTACTIVE.EQ.1)FROT(N)=FROTTS(N) + DT*180.D0*XWRK(NNNU)/PI
   70 CONTINUE
      RETURN
      END
      SUBROUTINE TSTEPS(T,IDT0,TSAVE,TTSV,MPDEF,DVMX,SIGMA,DT)
C
C    Calculate the timestep using the following criteria:
C
C    b) No element should be stretched or compressed by more
C       than (100/MPDEF)% in any given timestep
C    c) The product of growth rate (sigma) times DT < 1/4
C
      EPS=0.0001
      IT=IDT0
      DEFLIM=DVMX*FLOAT(MPDEF)
C     GRLIM=0.25      ! for growth rate criteia not enabled at present
C
C     successively halve time step until it passes deformation criteria
C
      DO WHILE(FLOAT(IT).LT.DEFLIM)
        IT=IT*2
      ENDDO
C
C     successively double time step if it is already smaller than needed
C
      DO WHILE((IT.GE.2*IDT0).AND.(FLOAT(IT).GT.(2.0*DEFLIM)))
        IT=IT/2
      ENDDO
      DT=1.0/FLOAT(IT)
C
C    ensure that solution is calculated and saved on integer multiple
C    of TSAVE, even if it is stored at other times also because of KSAVE
C
      NLEVL=1+INT(T/TSAVE)
      IF((T.LT.(TSAVE*FLOAT(NLEVL)-EPS)).AND.
     :   ((T+DT).GT.(TSAVE*FLOAT(NLEVL)+EPS)))THEN
        DT=TSAVE*FLOAT(NLEVL)-T
      END IF
      IF(ABS(T+DT-TSAVE*FLOAT(NLEVL)).LT.EPS)TTSV=TSAVE*FLOAT(NLEVL-1)
C
      RETURN
      END
      SUBROUTINE DFSTEP(DT,EXTS,EYTS,EX,EY,UVP,NOR,LEM,NUP,NE,NN,
     :                  NCOMP,IOFF,TBXOFF,TBYOFF,VELXO,VELYO,
     :                  SEGMIN,LSEGMIN,SHAPEMIN,LSHAPEMIN,
     :                  ANGMIN,LANGMIN)
C
C    This routine updates the triangular mesh using the present
C     velocity field and returns information on minimum segment length
C     and minmimum shape parameter
C
      DIMENSION EX(NUP),EY(NUP),NOR(NUP),UVP(NUP,2),LEM(6,NE)
      DIMENSION EXTS(NUP),EYTS(NUP)
      DIMENSION SD2(3)
C
      IF (IOFF.NE.0) THEN
        TBXOFF = TBXOFF + DT*VELXO
        TBYOFF = TBYOFF + DT*VELYO
      ENDIF

      DO N=1,NUP
        NI=NOR(N)
        IF((NI.GE.1).AND.(NI.LE.NN))THEN
          EX(NI)=EXTS(NI) + DT*UVP(N,1)
          EY(NI)=EYTS(NI) + DT*UVP(N,2)
C
C     correction for thin spherical sheet
C
          IF(NCOMP.EQ.-1)THEN
            EX(NI)=EX(NI) + DT*EXTS(NI)*UVP(N,2)/TAN(EYTS(NI))
          ENDIF
        ENDIF
      ENDDO
C
C   Calculate new position of midpoint nodes
C
      DO JEL=1,NE
        NJ1=NOR(LEM(1,JEL))
        NJ2=NOR(LEM(2,JEL))
        NJ3=NOR(LEM(3,JEL))
        NJ4=NOR(LEM(4,JEL))
        NJ5=NOR(LEM(5,JEL))
        NJ6=NOR(LEM(6,JEL))
        EX(NJ4)=(EX(NJ3)+EX(NJ1))*0.5
        EX(NJ5)=(EX(NJ1)+EX(NJ2))*0.5
        EX(NJ6)=(EX(NJ2)+EX(NJ3))*0.5
        EY(NJ4)=(EY(NJ3)+EY(NJ1))*0.5
        EY(NJ5)=(EY(NJ1)+EY(NJ2))*0.5
        EY(NJ6)=(EY(NJ2)+EY(NJ3))*0.5
      ENDDO
C
C    Calculate lengths of sides (squared) and cosines of
C    angles for the triangle elements
C    The sum of the cosines has a maximum of 1.5 for an equilateral
C    triangle, and approaches 1.0 as the triangle approaches
C    (0.,0.,180.) or (0.,90.,90.)
C
      SEGSQMIN=999.999
      SHAPEMIN=2.0
      COSMAX = 0
      DO 25 N=1,NE
        DO 10 K1=1,3
          K2=MOD(K1,3)+1
          K3=MOD(K1+1,3)+1
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          BY=EY(LK2)-EY(LK3)
          CX=EX(LK3)-EX(LK2)
          SD2(K1)=BY*BY + CX*CX
          IF (SD2(K1).LT.SEGSQMIN) THEN
            SEGSQMIN=SD2(K1)
            LSEGMIN=N
          END IF
   10   CONTINUE
        DEFIND=0.0
        DEFSEG=0
        DO 20 K1=1,3
          K2=MOD(K1,3)+1
          K3=MOD(K1+1,3)+1
          SSS=2.0*SQRT(SD2(K2)*SD2(K3))
          CSA=(SD2(K2) + SD2(K3) - SD2(K1))/SSS
          DEFIND=DEFIND + CSA
          IF ((CSA.GT.0).AND.(CSA.GE.COSMAX)) THEN
            COSMAX=CSA
            ANGMIN=ACOS(CSA)
            LANGMIN=N
          END IF
   20   CONTINUE
        IF (DEFIND.LT.SHAPEMIN) THEN
          SHAPEMIN=DEFIND
          LSHAPEMIN=N
        END IF
   25 CONTINUE
      SEGMIN=SQRT(SEGSQMIN)
      RETURN
      END

      SUBROUTINE WORKED(THDISS,THDINT,NP,DTDISS,BETA,
     :                  DMIN,DMAX,LUW,LSC)
C
C    This routine integrates the work done by thermal dissipation 
C    at every integration point
C
      DIMENSION THDISS(NP),THDINT(NP)
C
      DMIN=1.E30
      DMAX=0.0
      DO 100 J=1,NP
        THDINT(J)=THDINT(J)+(THDISS(J)-BETA*THDINT(J))*DTDISS
        IF(THDINT(J).LT.DMIN)DMIN=THDINT(J)
        IF(THDINT(J).GT.DMAX)DMAX=THDINT(J)
  100 CONTINUE
      WRITE(LUW,*)'Min and Max of THDINT =',DMIN,DMAX
      WRITE(LSC,*)'Min and Max of THDINT =',DMIN,DMAX
      RETURN
      END

      SUBROUTINE PLOAD(ISET,EX,EY,LEM,NOR,IELFIX,PNI,PLI,
     :                 SSQ,QLOAD,ARGANP,BRGANP,THRESH,
     :                 NROWS,NUP,NE,NCOMP,LSC)
C
C     Routine to calculate the load vector due to horizontal
C     crustal thickness variation...put into QLOAD.
C
      DOUBLE PRECISION PNI,PLI,BY,CX,DNDP,XX,YY,XLA,YCOT,WGT
      DOUBLE PRECISION YLA,DNIDX,DNIDY,EXL,WEXL,SUM1,SUM2
      DOUBLE PRECISION EXZERO,ARGANE
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION NOR(NUP)
      DIMENSION PNI(7,6),PLI(7,3)
      DIMENSION QLOAD(NROWS)
      DIMENSION SSQ(NUP),IELFIX(NUP)
      DIMENSION BY(3),CX(3),DNDP(84)
      DIMENSION WGT(7)
      DIMENSION XX(3),YY(3),XLA(7),YCOT(7)
C
C     DATA W/0.05,0.05,0.05,0.13333333,0.13333333,0.13333333,
C    10.45/
      SAVE WGT
      DATA WGT/0.12593918,0.12593918,0.12593918,0.13239415,
     1       0.13239415,0.13239415,0.225/
C
C    Clear the load array
C
      EXZERO=EXP(2.0*DFLTSSQ)
      DO I=1,NROWS
        QLOAD(I)=0.0
      ENDDO
C
C    Set IELFIX where crustal thickness has exceeded THRESH
C    (set once and leave set for any particular node)
C
      ISET=0
      DO J=1,NUP
        IF(SSQ(J).GT.THRESH)IELFIX(J)=1
        ISET=ISET+IELFIX(J)
      ENDDO
C     WRITE(LSC,10001)ISET
C10001 FORMAT(I5,' Node points have so far been set in PLOAD')
C
C    Cycle over the elements
C
      DO 60 N=1,NE
C
C    Set up the geometrical coefficients for the natural
C    coordinates used in the triangle elements
C
        DO 10 K1=1,3
          K2=MOD(K1,3)+1
          K3=MOD(K1+1,3)+1
          LK1=NOR(LEM(K1,N))
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          BY(K1)=EY(LK2)-EY(LK3)
          CX(K1)=EX(LK3)-EX(LK2)
          IF(NCOMP.EQ.-1)THEN
            XX(K1)=EX(LK1)
            YY(K1)=EY(LK1)
          END IF
   10   CONTINUE
C
C   TRIA (in PLOAD) is twice the area of the triangle element
C
        TRIA=CX(3)*BY(2)-CX(2)*BY(3)
C
C   Evaluate the cot(latitude) term if required at each
C    integration point
C
        IF(NCOMP.EQ.-1)THEN
          DO 26 K7=1,7
            XLA(K7)=XX(1)*PLI(K7,1)+XX(2)*PLI(K7,2)+XX(3)*PLI(K7,3)
            YLA=YY(1)*PLI(K7,1)+YY(2)*PLI(K7,2)+YY(3)*PLI(K7,3)
            YCOT(K7)=1.D0/DTAN(YLA)
   26     CONTINUE
        END IF
C
C    Calculate the values and gradients of the interpolation
C    function at the seven integration points
C
        CALL DNCOM(0,BY,CX,DNDP)
C
C     for node number I
C
        DO I=1,6
          LI=LEM(I,N)
          SUM1=0.D0
          SUM2=0.D0
C
C     Integrate over the element
C
          DO 55 K=1,7
C
C     Calculate EXL=ArS**2  at each integration point(K)
C
            ARGANE=0.D0
            EXL=0.D0
            DO 35 J=1,6
              JN=LEM(J,N)
              ARGANE=ARGANE+(ARGANP+FLOAT(IELFIX(JN))*BRGANP)*PNI(K,J)
              EXL=EXL + 2.D0*SSQ(JN)*PNI(K,J)
   35       CONTINUE
            EXL=ARGANE*(EXP(EXL)-EXZERO)
C
C    perform integration
C
            KIN=(I-1)*14 + (K-1)*2 +1
            DNIDX=DNDP(KIN)
            DNIDY=DNDP(KIN+1)
            IF(NCOMP.EQ.-1)THEN
              DNIDY=DNIDY+YCOT(K)*(PNI(K,I)*TRIA+XLA(K)*DNIDX)
            END IF
            WEXL=WGT(K)*EXL
            SUM1=SUM1 + WEXL*DNIDX
            SUM2=SUM2 + WEXL*DNIDY
   55     CONTINUE
C
          QLOAD(LI)=QLOAD(LI)+0.5*SUM1
          QLOAD(LI+NUP)=QLOAD(LI+NUP)+0.5*SUM2
        ENDDO          !  on node number I
   60 CONTINUE
C
      RETURN
      END
