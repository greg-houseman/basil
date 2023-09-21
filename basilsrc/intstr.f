C*--------------------------------------------------------------------
C*    Basil / Sybil:   intstr.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE INTSTR(EX,EY,NOR,NUP,LEM,NE,SE,VHB,UVP,IVV,
     :                  NROWS,NFP,IFBC,IFEQV,JFBC1,JFBC2)
C
C   This subroutine calculates the stress
C
C     DOUBLE PRECISION UVP
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
      DOUBLE PRECISION TRIA2,DPEX1,DPEX2,DPEY1,DPEY2
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION UVP(NROWS)
      DIMENSION NOR(NUP)
      DIMENSION VHB(8,NE)
      DIMENSION IFBC(NFP),IFEQV(NFP)
      DIMENSION JFBC1(NFP),JFBC2(NFP)
      DIMENSION BY(3),CX(3),DNDP(84)
      DIMENSION AUX(6),AUY(6),AUP(3)
      DIMENSION AFXX(6),AFXY(6),AFYY(6),NODE(6)
      EPS=1.E-4
C
C  Default exponent and coefficient for the viscosity
C   calculation. They will only be modified if IVIS > 0.
C
      SEXP=SE
      SCOF=1.0
C
C Cycle over the fault nodes
C
      NUP2 = NUP*2
      DO 10 JB=1,NFP
         KMD=JFBC2(JB)
        IF(KMD.GT.3) THEN
         MID=IABS(IFBC(JB))
         JEL=JFBC1(JB)
         JEL2=JFBC1(IFEQV(JB))
C
C Integrate the stress
C
           KM2=KMD-3
           KM1=MOD(KM2+1,3)+1
           MI1=LEM(KM1,JEL)
           MI2=LEM(KM2,JEL)
C
C    ADST is the length of the side of the present element,
C    and ANX, ANY are the x and y components of the normal.
C
           LK1=NOR(MI1)
           LK2=NOR(MI2)
           D1=EX(LK2)-EX(LK1)
           D2=EY(LK2)-EY(LK1)
           ADST=SQRT(D1*D1 + D2*D2)
           APHI=ATAN2(D1,D2)
           ANX=COS(APHI)
           ANY=-SIN(APHI)
C
C IF IVIS > 0, modify SEXP and SCOF for each element
C
C           IF(IVIS.GT.0) THEN
C             SEXP=VISEXP(JEL)
C             SCOF=VISCOF(JEL)
C           END IF
C
C Calculate stress components
C
           DO 30 JI=1,6
              NODE(JI)=IABS(LEM(JI,JEL))
              AUX(JI)=UVP(NODE(JI))
              AUY(JI)=UVP(NODE(JI)+NUP)
              IF (JI.LT.4) AUP(JI)=UVP(NOR(NODE(JI))+NUP2)
C for the fault case
C              IF (JI.LT.4) AUP(JI)=UVP(NOR(NODE(JI))+NUP2+NFP)
   30      CONTINUE
C   Set up the arrays for DNCOM
C   Note: TRIA2 is 2 times the area of the triangle element
C
           TRIA2=0.0
           DO 40 K1=1,3
              K2=MOD(K1,3)+1
              K3=MOD(K1+1,3)+1
              LK2=NOR(NODE(K2))
              LK3=NOR(NODE(K3))
              BY(K1)=EY(LK2)-EY(LK3)
              CX(K1)=EX(LK3)-EX(LK2)
              DPEX1 = EX(LK2)
              DPEX2 = EX(LK3)
              DPEY1 = EY(LK2)
              DPEY2 = EY(LK3)
              TRIA2=TRIA2 + (DPEX1*DPEY2-DPEX2*DPEY1)
   40      CONTINUE
           CALL DNCOM(5,BY,CX,DNDP)
           DO 70 JL=1,3
C
C  Only calculate the stresses for the vertex nodes
C
              DUDX=0.0
              DUDY=0.0
              DVDX=0.0
              DVDY=0.0
              DO 50 JI=1,6
                  IA=2*((JI-1)*7+JL-1)
                  DUDX=DUDX+AUX(JI)*DNDP(IA+1)
                  DUDY=DUDY+AUX(JI)*DNDP(IA+2)
                  DVDX=DVDX+AUY(JI)*DNDP(IA+1)
                  DVDY=DVDY+AUY(JI)*DNDP(IA+2)
   50         CONTINUE
              DUDX=DUDX/TRIA2
              DUDY=DUDY/TRIA2
              DVDX=DVDX/TRIA2
              DVDY=DVDY/TRIA2
              DWDZ=-DUDX-DVDY
              EDXY=0.5*(DUDY+DVDX)
              EDOT=SQRT(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ
     1            + 2.0*EDXY*EDXY)
C
C evaluate the local effective viscosity
C
              SELOC = SEXP
              IF (IVV.GE.1) THEN
                SCOF=VHB(JL,JEL)
              ENDIF
              IF (IVV.GE.1) THEN
                IF (IVV.GT.2) SELOC=VHB(8,JEL)
                AVIS = SCOF * EDOT**(1.0/SELOC -1.0)
              ELSE
                AVIS = 1.0
              ENDIF
              IF(AVIS.GT.1.0E+5)AVIS=1.0E+5
              AFXX(JL)=2.0*AVIS*DUDX
              AFXY(JL)=2.0*AVIS*EDXY
              AFYY(JL)=2.0*AVIS*DVDY
   70      CONTINUE
           SUMXX=(AFXX(KM1)+AFXX(KM2))*0.5*ADST
           SUMXY=(AFXY(KM1)+AFXY(KM2))*0.5*ADST
           SUMYY=(AFYY(KM1)+AFYY(KM2))*0.5*ADST
           SUMP=(AUP(KM1)+AUP(KM2))*0.5*ADST
           QA=0.5*(SUMXX+SUMYY)
           QB=0.5*(SUMXX-SUMYY)
           THETA=0.5*ATAN2(-AFXY(JL),QB)
           THPHI=THETA-APHI
           CS2L=COS(2.0*THPHI)
           SN2L=SIN(2.0*THPHI)
           CS2T=COS(2.0*THETA)
           SNOR=QA+QB*CS2L/CS2T
           STAN=-QB*SN2L/CS2T
C endif kmd > 3
       ENDIF
   10 CONTINUE
      END
