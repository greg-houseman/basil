
C*    Basil / Sybil:   cginit.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE CNVERG(ICONV,XNEW,XOLD,DIF,RMSU,RMSV,RMSP,
     :                  KFIN,IT,NCOMP,EX,EY,WFIT,NROWS,
     :                  NUP,NN,NFP,IFLT,KSTEP,ITSTOP,LSC,LUW)
C
C    Tests convergence by comparing the old and new solutions
C    at each iteration.  The maximum change is compared against
C    the RMS value of the velocity, separately for each
C    component of velocity, and pressure if NCOMP.ne.0.
C    NCOMP=0 for the thin viscous sheet case.
C
      DOUBLE PRECISION XNEW,XOLD
      DOUBLE PRECISION STRD,STRN,STRS,STRNEW,STROLD
      DOUBLE PRECISION WFCO,WFITD,RMS,RMSSTR
      DOUBLE PRECISION UN,DUD,UD,UOLD,UNEW,US,RMSUD
      DOUBLE PRECISION VN,DVD,VD,VOLD,VNEW,VS,RMSVD
      DOUBLE PRECISION PN,DPD,PD,POLD,PNEW,PS,RMSPD
      DOUBLE PRECISION DIFT,DIFU,DIFV,DIFX,DIFXP
      INTEGER*4 IDAY(3)
      DIMENSION XNEW(NROWS),XOLD(NROWS)
      DIMENSION EX(NUP),EY(NUP)
      SAVE DIFXP,NLINC
C
      NUP2=NUP*2
      ICONV=0
      DIFT=DIF
      WFITD=WFIT
      IF(IT.EQ.1)THEN
        DIFXP=1.D0
        NLINC=0
      END IF
C
      UD=0.0
      VD=0.0
      US=0.0
      VS=0.0
C
C    Find the maximum difference between old and new solutions
C
      DO 10 I=1,NUP
        NI=I+NUP
        UN=XNEW(I)
        VN=XNEW(NI)
        DUD=UN-XOLD(I)
        DUD=DUD*DUD
        IF(DUD.GT.UD)THEN
          IMX=I
          UD=DUD
          UOLD=XOLD(I)
          UNEW=UN
        END IF
        DVD=VN-XOLD(NI)
        DVD=DVD*DVD
        IF(DVD.GT.VD)THEN
          IMY=I
          VD=DVD
          VOLD=XOLD(NI)
          VNEW=VN
        END IF
C
C     Calculate the RMS value of the two fields
C
        US=US + UN*UN
        VS=VS + VN*VN
   10 CONTINUE
C
C    evaluate RMS velocities
C
      RMSUD=DSQRT(US/DFLOAT(NUP))
      RMSVD=DSQRT(VS/DFLOAT(NUP))
      RMS=DSQRT(RMSUD*RMSUD + RMSVD*RMSVD)
      DIFU=DSQRT(UD)/RMS
      DIFV=DSQRT(VD)/RMS
C
C    reduce step size if convergence stalled
C
      DIFX=DSQRT(DIFU*DIFU+DIFV*DIFV)
      IF((IT.GT.1).AND.(DIFX.GT.DIFXP))THEN
        WFITD=0.5*WFIT
        WRITE(LSC,10020)WFITD
10020   FORMAT('non-monotonic convergence, check AC, WFIT = ',F5.3)
      END IF
      DIFXP=DIFX
      WFCO=1.D0-WFITD
C
C     Update the solution vector in XNEW
C
      DO I=1,NUP2
        XNEW(I)= WFITD*XNEW(I) + WFCO*XOLD(I)
      ENDDO
C
C     Repeat the above calculation for the stress integrals on the fault
C
      STRD=0.0
      STRS=0.0
      STROLD=0.0
      STRNEW=0.0
      RMSSTR=0
      IMSTR=0
Cc    DIFSTR=0
      IF (NFP.GE.1) THEN
        DO 20 I=1,NFP
          NI=NUP2+I
          STRN=XNEW(NI)
          DPD=STRN-XOLD(NI)
          DPD=DPD*DPD
          IF(DPD.GT.STRD)THEN
            IMSTR=I
            STRD=DPD
            STROLD=XOLD(NI)
            STRNEW=STRN
          END IF
C
C     Calculate the RMS value of the stress integrals
C
          STRS=STRS + STRN*STRN
C
C     Calculate the new solution vector
C
          XNEW(NI)= WFITD*STRN + WFCO*XOLD(NI)
   20   CONTINUE
        RMSSTR=DSQRT(STRS/DFLOAT(NFP))
      END IF
C
C     Repeat the above calculation for the pressure field
C
      IMP=0
      IF(NCOMP.GT.0)THEN
        PD=0.0
        PS=0.0
        DO 30 I=1,NN
          NI=NUP2+NFP+I
          PN=XNEW(NI)
          DPD=PN-XOLD(NI)
          DPD=DPD*DPD
          IF(DPD.GT.PD)THEN
            IMP=I
            PD=DPD
            POLD=XOLD(NI)
            PNEW=PN
          END IF
C
C     Calculate the RMS value of the pressure field
C
          PS=PS + PN*PN
C
C     Calculate the new solution vector
C
          XNEW(NI)= WFITD*PN + WFCO*XOLD(NI)
   30   CONTINUE
        RMSPD=DSQRT(PS/DFLOAT(NN))
        RMSP=RMSPD
      END IF
C
C     Watch for unchanging solution but failure to converge
C
      IF(KFIN.LT.5)THEN
        NLINC=NLINC+1
      ELSE
        NLINC=0
      END IF
C
C    Apply the convergence test (based on u and v only)
C
      RMSU=RMSUD
      RMSV=RMSVD
      IF(((DIFV.LE.DIFT).AND.(DIFU.LE.DIFT)).OR.(NLINC.GT.4))ICONV=1
      IF(NLINC.GT.4)WRITE(LSC,10021)DIF
10021 FORMAT('beware non-linear convergence not achieved but solution ',
     :       ' unchanging, ACNL = ',G12.5,' consider rescaling.')
      CALL DATIME(IDAY,ISEC,1)
C
C    ouput convergence data
C
      CALL REPORTCV(LSC,IFLT,NCOMP,IT,KFIN,UOLD,UNEW,RMSU,VOLD,VNEW,
     :               RMSV,STROLD,STRNEW,RMSSTR,POLD,PNEW,RMSP)
      IF((KSTEP.EQ.0).OR.(ICONV.EQ.1).OR.(IT.EQ.ITSTOP))
     :CALL REPORTCV(LUW,IFLT,NCOMP,IT,KFIN,UOLD,UNEW,RMSU,VOLD,VNEW,
     :               RMSV,STROLD,STRNEW,RMSSTR,POLD,PNEW,RMSP)
C     IF(IT.EQ.1)WRITE(LSC,10115)
C     WRITE(LSC,10118) KFIN,ISEC
C         CALL MATPRT(XNEW(NUP2+NFP+1),NN,1,NN,LSC)
C         DO 30999 IJ=1,NN
C           WRITE(LSC,31999)IJ,EX(IJ),EY(IJ),XNEW(IJ)
C31999      FORMAT(I4,1x,3G12.4)
C30999    CONTINUE
      RETURN
      END
C
      SUBROUTINE REPORTCV(LPP,IFLT,NCOMP,IT,KFIN,UOLD,UNEW,RMSU,VOLD,
     :                    VNEW,RMSV,STROLD,STRNEW,RMSSTR,POLD,PNEW,RMSP)
C
      DOUBLE PRECISION UOLD,UNEW,VOLD,VNEW,POLD,PNEW,STRNEW,STROLD
      DOUBLE PRECISION RMSSTR,SFAC
      SFAC=1.D3
C
C     subroutine called from CNVERG to report output of non-linear iterations
C
      IF ((IFLT.EQ.0).AND.(NCOMP.LE.0)) THEN
        IF (IT.EQ.1) WRITE(LPP,10111)
        WRITE(LPP,10116)IT,KFIN,UOLD-UNEW,RMSU,VOLD-VNEW,RMSV
      ELSE IF ((IFLT.NE.0).AND.(NCOMP.LE.0)) THEN
        IF (IT.EQ.1) WRITE(LPP,10112)
        WRITE(LPP,10117)IT,KFIN,UOLD-UNEW,RMSU,VOLD-VNEW,RMSV,
     :                  (STROLD-STRNEW)*SFAC,RMSSTR
      ELSE IF ((IFLT.EQ.0).AND.(NCOMP.GT.0)) THEN
        IF (IT.EQ.1) WRITE(LPP,10113)
        WRITE(LPP,10117)IT,KFIN,UOLD-UNEW,RMSU,VOLD-VNEW,RMSV,
     :                  POLD-PNEW,RMSP
      ELSE IF ((IFLT.NE.0).AND.(NCOMP.GT.0)) THEN
        IF (IT.EQ.1) WRITE(LPP,10114)
        WRITE(LPP,10118)IT,KFIN,UOLD-UNEW,RMSU,VOLD-VNEW,RMSV,
     :                  (STROLD-STRNEW)*SFAC,RMSSTR,POLD-PNEW,RMSP
      END IF
      RETURN
C
10111 FORMAT('  it     KF        dUo        Urms         dVo   ',
     :       '     Vrms  ')
10112 FORMAT('  it     KF        dUo        Urms         dVo   ',
     :       '     Vrms        dSTRo      STRrms ')
10113 FORMAT('  it     KF        dUo        Urms         dVo   ',
     :       '     Vrms      Po        Prms    ')
10114 FORMAT('  it     KF        dUo        Urms         dVo   ',
     :       '     Vrms        dSTRo      STRrms ',
     :       '     dPo         Prms   ')

10116 FORMAT(2I6,2(F13.8,F11.6))
10117 FORMAT(2I6,2(F13.8,F11.6),F13.5,F11.4)
10118 FORMAT(2I6,4(F13.8,F11.6))
      END
C
C----------------------------------------------------------------
C
      SUBROUTINE CGRUN(INTV,RLV,UVP,EX,EY,LEM,NOR,VHB,THDISS,
     :                 THDINT,SSQ,QLOAD,QBND,IBC,IBCTYP,IBNGH,
     :                 IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,PNI,PLI,
     :                 STK,STK2,STK3,STK4,IK,IK2,IK4,
     :                 VISCMIN,VISCMAX,ED2IMIN,ED2IMAX,
     :                 THDIMIN,THDIMAX,NE,NBP,NN,NUP,NFP,NROWS,
     :                 KBWP,K2BWP,K4BWP,NKP,NK2P,NK4P,IREAD,
     :                 IT,RMSU,RMSV,NITER,IOFF,NFPF3,IVERB,IERR)
C
C    routine to control the calculation of the velocity field
C    at each timestep
C
      INCLUDE 'indices.parameters'
C     COMMON/STKDIM/KBWP,K2BWP,NKP,NK2P,K4BWP,NK4P
C     SAVE /STKDIM/
      COMMON/AI/LUW,LSC,LBC,LLG
      DOUBLE PRECISION XLOAD,PDP,XWRK,PWTS
      DOUBLE PRECISION STK,STK2,STK3,STK4,VISN
      DOUBLE PRECISION PNI,PLI
      DOUBLE PRECISION BIG,GAMMAD,TREFD
      DIMENSION LEM(6,NE)
      DIMENSION VHB(8,NE)
      DIMENSION THDISS(7,NE),THDINT(7,NE)
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION SSQ(NUP)
      DIMENSION XLOAD(NROWS),PDP(NROWS),XWRK(NROWS)
      DIMENSION QLOAD(NROWS)
      DIMENSION IBC(NBP),IBCTYP(NBP*2)
      DIMENSION IFBC1(NFP),IFBC2(NFP),IFEQV(NFP),JFBC1(NFP),JFBC2(NFP)
      DIMENSION QBND(NBP*2)
      DIMENSION UVP(NROWS)
      DIMENSION PNI(7,6),PLI(7,3)
      DIMENSION INTV(64)
      DIMENSION RLV(64)
      DIMENSION STK(NKP),STK2(NK2P),STK3(NFP),STK4(NK4P)
      DIMENSION IK(NKP),IK2(NK2P),IK4(NK4P)
      DIMENSION VISN(7,NE)
      DIMENSION PWTS(NN)
      DATA IFIRST/0/
      SAVE IFIRST
C
C     Setting and use of viscosity parameters should be consistent with
C     definition of IVV as follows
C
C       * ivv    visc (vhb[i,1-7])   se (vhb[i,8])
C       *  0       constant           constant(1)    No vhb array
C       *  1       variable           constant(1)
C       *  2       constant           constant(>1)   
C       *  3       variable           constant(>1)
C       *  4       variable           variable
C
C     Assemble the STK2 and STK3 matrices which deal with the fault bc's
C
      BIG=RLV(IBIG)
      GAMMAD=RLV(IGAMMA)
      TREFD=RLV(ITREF)
      IF(INTV(IIFLT).NE.0)THEN
        CALL SEMFLT(BIG,STK2,STK3,IK2,NK2P,K2BWP,
     :              EX,EY,LEM,NOR,IBCTYP,IBNGH,QLOAD,
     :              QBND,IFBC1,IFBC2,IFEQV,JFBC1,JFBC2,
     :              RLV(IVELXO),RLV(IVELYO),IOFF,
     :              INTV(IIFLT),NUP,NE,NROWS,
     :              NBP,NFP,NFPF3,LUW,IERR)
        IF (IERR.NE.0) GO TO 300
C       WRITE(*,*)'STK3 after SEMFLT NFP =',NFP
C       CALL DMATPP(STK3,NFP,LUW)
      ENDIF
C
C  Assemble the STK4 (pressure) matrix
C  NOTE: STK2, STK3, STK4 do not use the viscosity, so they don't need
C    Be included in the following iterative loop
C
      IF(INTV(INCOMP).GT.0) THEN
            CALL SEMBL4(STK4,IK4,NK4P,K4BWP,
     :                  EX,EY,LEM,NOR,INTV(INCOMP),PNI,PLI,QLOAD,
     :                  INTV(IIFLT),NUP,NE,
     :                  NROWS,NFP,LUW,IERR)
        IF (IERR.NE.0) GO TO 300
      END IF
C
C    Copy the UVP array (single) into the XWRK array (double)
C
      DO J=1,NROWS
        XWRK(J)=UVP(J)
      ENDDO
C       WRITE(*,*)'VHB before SEMBL NE =',NE
C       CALL SMATPP(VHB,8*NE)
C       WRITE(*,*)'QBND before SEMBL'
C       CALL SMATPP(QBND,2*NBP)
C       WRITE(*,*)'IBCTYP before SEMBL'
C       CALL MITPP(IBCTYP,2*NBP)
C
C    Commence the iteration cycle for non-linear rheology
C
      IT=0
  200 IT=IT+1
      IF(IT.GT.INTV(IITSTOP)) THEN
        IERR = 1
        WRITE(LUW,10100)INTV(IITSTOP)
10100   FORMAT('Allowed non-linear interations exceeded: ITSTOP =',I6)
        GO TO 300
      END IF
C
C   Initialise VISN:  In general handled by VISK, but not
C   if UVP is still empty or SE = 1.0
C
      IF(((IREAD.EQ.0).AND.(IFIRST.EQ.0)).OR.(INTV(IIVV).LE.1))THEN
        CALL VISCLIN(VISN,VHB,RLV(IVC),NE,INTV(IIVIS))
      ELSE
C
C    XWRK is passed to VISK to compute current viscosity field VISN
C
        CALL VISK(BIG,INTV(IJLV),RLV(IHLENSC),RLV(IAREA),
     :                RLV(IVOLUM),RLV(IVISCD),RLV(IDVMX),RLV(ISE),
     :                EX,EY,LEM,NOR,SSQ,XWRK,VISN,VHB,THDISS,THDINT,
     :                GAMMAD,TREFD,PNI,PLI,INTV(INCOMP),INTV(IIVV),
     :                INTV(IIVIS),RLV(IVC),
     :                VISCMIN,VISCMAX,ED2IMIN,ED2IMAX,
     :                THDIMIN,THDIMAX,NUP,NE,NROWS,
     :                LUW,LSC,IERR)
        IF (IERR.NE.0) GO TO 300
C       WRITE(*,*)'VISN after VISK, NROWS =',7*NE
C       CALL DMATPP(VISN,7*NE,LUW)
      ENDIF
C
C     Assemble the stiffness matrices using viscosity VISN from VISK
C     NCOMP=0 or -1 for thin viscous sheet or shell
C     NCOMP=1 or 2 for plane-strain or axisymmetric solution
C
      CALL SEMBL(BIG,STK,IK,NKP,KBWP,INTV(INCOMP),
     :           EX,EY,PWTS,LEM,NOR,IBC,IBCTYP,PNI,PLI,
     :           QBND,VISN,NUP,NN,NE,
     :           NBP,IERR)
      IF (IERR.NE.0) GO TO 300
C
C     put the load vector into XLOAD and add boundary conditions
C
      CALL LOADBC(BIG,XLOAD,QLOAD,QBND,IBC,IBCTYP,
     :            NUP,NBP,NROWS)
C
C     Solve the matrix equation using the conjugate gradient
C     method, solving U and V simultaneously.  The estimate of
C     the solution obtained from the previous iteration is in XWRK
C     XLOAD initially contains the load vector, but after CONJGR is 
C     called XLOAD contains the new solution.
C
      DO J=1,NROWS
        PDP(J)=XWRK(J)
      ENDDO
      CALL CONJGR(XLOAD,PDP,KFIN,RLV(IAC),INTV(INCOMP),
     :            STK,STK2,STK3,STK4,PWTS,
     :            IK,IK2,IK4,NKP,NK2P,NK4P,KBWP,K2BWP,K4BWP,
     :            NUP,NN,NROWS,NFP,
     :            INTV(IIFLT),NITER,LUW,LSC,IVERB,IERR)
      IF(IERR.NE.0)GO TO 300
C     WRITE(*,*)'XLOAD after CONJGR, NUP =',NUP
C     CALL DMATPP(XLOAD,NUP,LUW)
C
C     Test convergence if non-linear. Copy solution from XLOAD to XWRK
C
      CALL CNVERG(JCONV,XLOAD,XWRK,RLV(IACNL),RMSU,RMSV,RMSP,
     :            KFIN,IT,INTV(INCOMP),EX,EY,RLV(IWFIT),
     :            NROWS,NUP,NN,NFP,
     :            INTV(IIFLT),INTV(IKSTEP),INTV(IITSTOP),LSC,LUW)
      DO I=1,NROWS
        XWRK(I)=XLOAD(I)
      ENDDO
      IFIRST=1                     ! first solution obtained
C
C    go back to 200 for the next iteration
C
      IF((JCONV.EQ.0).AND.(INTV(IIVV).GT.1))GO TO 200 
C
C    Call VISK once more to evaluate dissipation integrals using final soln -> exit
C
      CALL VISK(BIG,INTV(IJLV),RLV(IHLENSC),RLV(IAREA),
     :            RLV(IVOLUM),RLV(IVISCD),RLV(IDVMX),RLV(ISE),
     :            EX,EY,LEM,NOR,SSQ,XLOAD,VISN,VHB,THDISS,THDINT,
     :            GAMMAD,TREFD,PNI,PLI,INTV(INCOMP),INTV(IIVV),
     :            INTV(IIVIS),RLV(IVC),
     :            VISCMIN,VISCMAX,ED2IMIN,ED2IMAX,
     :            THDIMIN,THDIMAX,NUP,NE,NROWS,
     :            LUW,LSC,IERR)
      IF (IERR.NE.0)GO TO 300
C
C     CALL INTSTR(EX,EY,NOR,NUP,LEM,NE,RLV(ISE),
C    :            VHB,XWRK,INTV(IIVV),INTV(INROWS),NFP,
C    :            IFBC1,IFEQV,JFBC1,JFBC2)
C  Update pressure bit of UVP for the thin-viscous sheet case
C  ***Only correct for t=0, when pres=-tzz
C  *****must update so pres=szz-tzz  and szz=AR*S^2 for t>0
C
      IF(INTV(INCOMP).LE.0)THEN
        CALL TVSPRES(EX,EY,NOR,LEM,XLOAD,NUP,NE,NROWS,
     :               NN,NFP,RLV(ISE))
      ENDIF
C
C    iterations finished; solution copied back to UVP (sp)
C
      DO J=1,NROWS
        UVP(J)=XWRK(J)
      ENDDO
C
C    print viscosity array if required
C
      IF(INTV(IIPRINT).NE.0)THEN
        WRITE(LUW,10112)
10112   FORMAT(' the viscosity distribution follows ....',/)
        CALL DMATPRT(VISN,7,NE,7*NE,LUW)
      END IF
      RETURN
C
  300 CONTINUE
      WRITE(*,10113)IERR
      WRITE(LUW,10113)IERR
10113 FORMAT('CGRUN is returning IERR = ',I3,' job stopping')
      RETURN
      END
C
      SUBROUTINE SEMBL(BIG,STK,IK,NK1,KBW,NCOMP,
     :                 EX,EY,PWTS,LEM,NOR,IBC,IBCTYP,PNI,PLI,
     :                 QBND,VISN,NUP,NN,NE,NBP,IERR)
C
C     This routine calculates and assembles the entries for the
C     matrix which previously included K1,K5,K7.
C     The entries are stored in STK(NK1).  Compressed format is used
C     so the entries are indexed by the values in IK.
C     The matrix is not necessarily symmetric as in previous cases.
C     Note that this version of
C     SEMBL stores rows for use by Subroutine CONJGR, whereas
C     the old version stored columns for use by DECOMP.
C
      DOUBLE PRECISION ASCT1,ASCT2,TRIA,TRIX,YLA,RKI,RKJ,VIS,BIG
      DOUBLE PRECISION SUM1,SUM5,SUM7,XR,YL,QNI,QNJ,AIJ,BIJ,CIJ
      DOUBLE PRECISION DNIDX,DNIDY,DNJDX,DNJDY,DNDX2,DNDY2
      DOUBLE PRECISION WGT,XX,YY,BY,CX,XRAD,YCOT,PNI,PLI,DNDP
      DOUBLE PRECISION STK,VISN,PWTS
C
      COMMON/AI/LUW,LSC,LBC,LLG
      SAVE /AI/
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION IBC(NBP),IBCTYP(NBP*2)
      DIMENSION QBND(NBP*2)
      DIMENSION VISN(7,NE)
      DIMENSION WGT(7)
      DIMENSION STK(NK1)
      DIMENSION IK(NK1)
      DIMENSION PWTS(NN)
      DIMENSION BY(3),CX(3),DNDP(84)
      DIMENSION PNI(7,6),PLI(7,3)
      DIMENSION XX(3),XRAD(7),YY(3),YCOT(7)
C
      SAVE WGT
      DATA WGT/0.12593918,0.12593918,0.12593918,0.13239415,
     1         0.13239415,0.13239415,0.225/
C
C     each matrix element is calculted by numerical integration
C    using a 7 point approximation from p421 in Huebner.  For
C    a quadratic integrand (as here) the integral is exact.
C    the seven points are given the preceeding (in WGT) weights
C
C     NCOMP=-1 -> thin viscous sheet on sphere
C     NCOMP=0  -> thin viscous sheet; Cartesian
C     NCOMP=1  -> 2D plane strain
C     NCOMP=2  -> 2D cylindrical axisymmetry (axis: x = 0)
C
      ASCT1=4.D0
      ASCT2=2.D0
      IF(NCOMP.GT.0)THEN
         ASCT1=2.D0
         ASCT2=0.D0
      END IF
C
C     Clear the address arrays for the matrices
C
      DO I=1,NK1
        IK(I)=0
      ENDDO
      NUP2=NUP*2
      DO I=1,NUP2
        KC=(I-1)*KBW+1
        IK(KC)=I
        STK(KC)=0.D0
      ENDDO
      DO I=1,NN
        PWTS(I)=0.D0
      ENDDO
C     CALL DMATPRT(PNI,7,6,42,LSC)
C     CALL DMATPRT(PLI,7,3,21,LSC)
C
C   Cycle through the elements
C
      DO 70 N=1,NE
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
          XX(K1)=EX(LK1)
          YY(K1)=EY(LK1)
          BY(K1)=EY(LK2)-EY(LK3)
          CX(K1)=EX(LK3)-EX(LK2)
   10   CONTINUE
C
C   TRIA is 4 times the area of the triangle element
C
        TRIA=(CX(3)*BY(2)-CX(2)*BY(3))*2.D0
        IF(TRIA.LE.0.D0)THEN
          WRITE(6,*)'Problem in SEMBL, element inverted'
          WRITE(6,*)'N, XX(1), YY(1) =',N, XX(1), YY(1)
          IERR=1
          RETURN
        ENDIF
        TRIX=TRIA*TRIA*0.5D0
C
C   The X-coordinate is interpreted as radius if NCOMP=2
C
        IF(NCOMP.EQ.2)THEN
          DO 25 K7=1,7
            XRAD(K7)=XX(1)*PLI(K7,1)+XX(2)*PLI(K7,2)+XX(3)*PLI(K7,3)
   25     CONTINUE
C
C   The Y-coordinate is interpreted as latitude if NCOMP=-1
C    assumed to be in radians at this point.
C
        ELSE IF(NCOMP.EQ.-1)THEN
          DO 26 K7=1,7
            XRAD(K7)=XX(1)*PLI(K7,1)+XX(2)*PLI(K7,2)+XX(3)*PLI(K7,3)
            YLA     =YY(1)*PLI(K7,1)+YY(2)*PLI(K7,2)+YY(3)*PLI(K7,3)
            YCOT(K7)=1.D0/DTAN(YLA)
   26     CONTINUE
        END IF
C
C  DNCOM computes dNi/dy & dNj/dx (times 2*area)
C
        CALL DNCOM(0,BY,CX,DNDP)
C
C    For each element calculate the twenty-one entries to the
C    stiffness matrices K1 and K7.
C
        DO I=1,6
          K1=LEM(I,N)
C         IF(I.LE.3)PWTS(NOR(K1))=PWTS(NOR(K1))+TRIA  ! for use in DEMEANP
          IF(NCOMP.EQ.2)RKI=EX(NOR(K1))
C
          DO J=1,6
            K2=LEM(J,N)
            IF(NCOMP.EQ.2)RKJ=EX(NOR(K2))
C
C    The value of the integrand at each of the seven points
C    is calculated, and summed with the appropriate weight
C
            SUM1=0.D0
            SUM5=0.D0
            SUM7=0.D0
            DO 30 K=1,7
              KIN=(I-1)*14 + (K-1)*2 + 1
              KJN=(J-1)*14 + (K-1)*2 + 1
              DNIDX=DNDP(KIN)
              DNJDX=DNDP(KJN)
              DNIDY=DNDP(KIN+1)
              DNJDY=DNDP(KJN+1)
              DNDX2=DNIDX*DNJDX
C
C   for spherical thin viscous sheet (NCOMP=-1)
C
              IF(NCOMP.EQ.-1)THEN
                XR=XRAD(K)
                YL=YCOT(K)
                QNI=PNI(K,I)*TRIA*0.5D0
                QNJ=PNI(K,J)*TRIA*0.5D0
                DNIDY=DNIDY+YL*(XR*DNIDX-QNI)
                DNJDY=DNJDY+YL*(XR*DNJDX-QNJ)
                AIJ = ASCT1*DNDX2 + DNIDY*DNJDY
                BIJ = DNIDY*DNJDX
                DNIDY=DNIDY+1.5D0*YL*QNI
                DNJDY=DNJDY+1.5D0*YL*QNJ
                CIJ = DNDX2 + ASCT1*DNIDY*DNJDY + 3.D0*YL*YL*QNI*QNJ
                DNJDY=DNJDY+1.5D0*YL*QNJ
                BIJ = BIJ+ASCT2*DNIDX*DNJDY
C
C    AIJ is the u-u coefficient, 
C    BIJ is the u-v coefficient
C    CIJ is the v-v coefficient
C
              ELSE        !   not spherical sheet
                DNDY2=DNIDY*DNJDY
                AIJ = ASCT1*DNDX2 + DNDY2
                CIJ = DNDX2 + ASCT1*DNDY2
                BIJ = DNIDY*DNJDX + ASCT2*DNIDX*DNJDY
C
C  For axisymmetric geometry (around x = 0)
C
                IF(NCOMP.EQ.2)THEN
                  XR=XRAD(K)
                  BIJ = XR*BIJ
                  CIJ = XR*CIJ
                  IF((DABS(RKI).GT.1.D-6).AND.(DABS(RKJ).GT.1.D-6))
     1                 AIJ = XR*AIJ+PNI(K,I)*PNI(K,J)*TRIX/XR
                END IF
              END IF
C
              VIS=VISN(K,N)
              SUM1=SUM1 + WGT(K)*AIJ*VIS
              SUM5=SUM5 + WGT(K)*BIJ*VIS
              SUM7=SUM7 + WGT(K)*CIJ*VIS
   30       CONTINUE
            SUM1=SUM1/TRIA
            SUM5=SUM5/TRIA
            SUM7=SUM7/TRIA
C
C     The matrix is assembled and stored in compressed mode,
C            ++++++++++  ROW WISE ++++++++
C ----------------------------------------------
C     SUM1 goes into row K1, col K2
C     SUM5 goes into row K1, col K2+NUP
C     SUM5 also goes into row K2+NUP, col K1
C     SUM7 goes into row K1+NUP, col K2+NUP
C ----------------------------------------------
C      CALL MATPAD(K1,K2,SUM1)
C      CALL MATPAD(K1,K2+NUP,SUM5)
C      CALL MATPAD(K2+NUP,K1,SUM5)
C      CALL MATPAD(K1+NUP,K2+NUP,SUM7)
C
C  The above 4 lines are replaced by the following 5 in 
C   order to upper diagonalize the STK matrix
C
            CALL MATPAD(K1,K2+NUP,SUM5,STK,IK,NK1,KBW,IERR)
            IF (IERR.NE.0) GO TO 900
            IF(K1.LE.K2) THEN
              CALL MATPAD(K1,K2,SUM1,STK,IK,NK1,KBW,IERR)
              IF (IERR.NE.0) GO TO 900
              CALL MATPAD(K1+NUP,K2+NUP,SUM7,STK,IK,NK1,KBW,IERR)
              IF (IERR.NE.0) GO TO 900
            END IF
          ENDDO             !  on J, node number in element
        ENDDO               !  on I, interpolation function
   70 CONTINUE              !  on N, element nmber
C
C    set the diagonal element for nodes with fixed velocity
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
C     WRITE(*,*)'IK matrix indices follow'
C     CALL MITPRT(IK,KBW,NUP2,KBW*NUP2,6)
C     WRITE(*,*)'STK matrix follows'
C     CALL DMATPRT(STK,KBW,NUP2,KBW*NUP2,6)
  900 RETURN
      END
C
      SUBROUTINE SEMBL4(STK4,IK4,NK4,K4BW,
     :                  EX,EY,LEM,NOR,NCOMP,PNI,PLI,QLOAD,
     :                  IFLT,NUP,NE,NROWS,NFP,LUW,IERR)
C
C     This routine calculates and assembles the entries for the
C     K4 matrix.  The entries are stored in STK4(NK4P)
C     using compressed format, indexed by the values in IK4.
C
      DOUBLE PRECISION TRIA,XR,BIJ,CIJ,PLJ,SUMA,SUMB
      DOUBLE PRECISION WGT,XX,BY,CX,XRAD,PNI,PLI,DNDP
      DOUBLE PRECISION STK4
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION QLOAD(NROWS)
      DIMENSION PNI(7,6)
      DIMENSION PLI(7,3)
      DIMENSION WGT(7)
      DIMENSION STK4(NK4)
      DIMENSION IK4(NK4)
      DIMENSION BY(3),CX(3),DNDP(84)
      DIMENSION XX(3),XRAD(7)
      SAVE WGT
      DATA WGT/0.12593918,0.12593918,0.12593918,0.13239415,
     10.13239415,0.13239415,0.225/
C
C     each matrix element is calculted by numerical integration
C    using a 7 point approximation from p421 in Huebner.  For
C    a quadratic integrand (as here) the integral is exact.
C    the seven points are given the preceeding (in W) weights
C
C     Clear the address arrays for the matrices
C
      DO I=1,NK4
        IK4(I)=0
      ENDDO
C
C   Cycle through the elements
C
      DO 70 N=1,NE
C
C    Set up the geometrical coefficients for the natural
C    coordinates used in the triangle elements
C
        TRIA=0.D0
        DO 20 K1=1,3
          K2=MOD(K1,3)+1
          K3=MOD(K1+1,3)+1
          LK1=NOR(LEM(K1,N))
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          XX(K1)=EX(LK1)
          BY(K1)=EY(LK2)-EY(LK3)
          CX(K1)=EX(LK3)-EX(LK2)
   20   CONTINUE
C
C   TRIA is 4 times the area of the triangle element
C
        IF(NCOMP.EQ.2)THEN
          TRIA=(CX(3)*BY(2)-CX(2)*BY(3))*2.D0
          DO 25 K7=1,7
            XRAD(K7)=XX(1)*PLI(K7,1)+XX(2)*PLI(K7,2)+XX(3)*PLI(K7,3)
   25     CONTINUE
        END IF
C
        CALL DNCOM(0,BY,CX,DNDP)
C
C    For each element calculate the eighteen entries to the
C    stiffness matrices K3 and K4.
C
        DO I=1,6
          K1=LEM(I,N)
          DO J=1,3
            K2=LEM(J,N)
            NRK2=NOR(K2)
C
C    The value of the integrand at each of the seven points
C    is calculated, and summed with the appropriate weight
C
            SUMA=0.D0
            SUMB=0.D0
            KIN=(I-1)*14 + 1
            KJL=(J-1)*7 + 1
C
C    Do the integration, BIJ is the u-p coefficient, CIJ is the v-p coefficient
C
            DO 30 K=1,7
              BIJ=DNDP(KIN)
              CIJ=DNDP(KIN+1)
              PLJ=PLI(K,J)
C
C   For axisymmetric case
C
              IF(NCOMP.EQ.2)THEN
                XR=XRAD(K)
                BIJ = XR*BIJ+PNI(K,I)*TRIA*0.5D0
                CIJ = XR*CIJ
              END IF
C
              SUMA=SUMA + WGT(K)*BIJ*PLJ
              SUMB=SUMB + WGT(K)*CIJ*PLJ
              KIN=KIN+2
              KJL=KJL+1
   30       CONTINUE
C
C     The matrix is assembled and stored in compressed mode.
C            ++++++++++  COLUMN WISE ++++++++
C
            KC0=(NRK2-1)*K4BW
            KC=KC0
   35       KC=KC+1
            K4BS=KC-KC0
C
C    Check that the available length of the row is not exceeded
            IF(K4BS.GT.K4BW-1)GO TO 102
C
C    Initialise a matrix entry if none already exists
C
            IF(IK4(KC).LE.0)THEN
              IK4(KC)=K1
              STK4(KC)=0.D0
              IK4(KC+1)=K1+NUP
              STK4(KC+1)=0.D0
            END IF
C
C    If correct matrix entry is found, add the nodal contribution
C    from the current element.  TRIA cancels out in integration
C
            IF(IK4(KC).EQ.K1)THEN
              STK4(KC) = STK4(KC) + SUMA*0.5D0
              STK4(KC+1) = STK4(KC+1) + SUMB*0.5D0
C
C    If incorrect matrix entry is found, keep looking
C
            ELSE
              GO TO 35
            END IF
          ENDDO                  ! on J, node number in element
        ENDDO                    ! in I, interpolation fucntion
   70 CONTINUE                   ! on N, element number
C      WRITE(LUW,*)'IK4 array from SEMBL4'
C      CALL MITPRT(IK4,K4BW,NNP,K4BW*NNP,LUW)
C      WRITE(LUW,*)'STK4 array from SEMBL4'
C      CALL DMATPRT(STK4,K4BW,NNP,K4BW*NNP,LUW)
      RETURN
  102 WRITE(LUW,10103)K4BW,K1,K2
      IERR=3
  900 RETURN
10103 FORMAT(' ALLOWED SPACE IN K4 EXCEEDED : K4BW =',I5,
     1',  K1 =',I5,',  K2 =',I5,';  INCREASE K4BW !!',//)
      END
C
      SUBROUTINE MATPAD(IROW,ICOL,VALU,STK,IK,NK,KBW,IERR)
C
C Adds the VALU to the STK matrix at row IROW and column ICOL
C
      COMMON/AI/LUW,LSC,LBC,LLG
      DOUBLE PRECISION STK,VALU
      DIMENSION STK(NK),IK(NK)

      KC0=(IROW-1)*KBW
      KC=KC0
  100 KC=KC+1
      KBS=KC-KC0
C
C    Check that the available length of the row is not exceeded
C
      IF(KBS.GT.KBW)GO TO 999
C
C    Initialise a matrix entry if none already exists
C
      IF(IK(KC).LE.0)THEN
        IK(KC)=ICOL
        STK(KC)=0.D0
      END IF
C
C    If correct matrix entry is found, add the value
C
      IF(IK(KC).EQ.ICOL)THEN
        STK(KC) = STK(KC)+VALU
        RETURN
      ELSE
        GO TO 100
      END IF
C
C    Error trap
C
  999 WRITE(LSC,10999)KBW,KBS,IROW,ICOL,VALU
      WRITE(LUW,10999)KBW,KBS,IROW,ICOL,VALU
      IERR=3
      RETURN
10999 FORMAT(' ALLOWED SPACE IN STK EXCEEDED : KBW=',I5,', KBS=',I5,
     1',  ROW=',I5,',  COL=',I5,' VAL=',F9.3,';  INCREASE KBW !!',//)
      END
C
      SUBROUTINE LNCORD(PNI,PLI)
C
C    This routine sets up the local coordinates for the
C    integration points in the triangle.  These are the
C    same for every element, and this routine need only
C    be called once at start of program
C    For each I, 7 values of L(i) and N(i) are stored in
C    consecutive memory slots.
C
      DOUBLE PRECISION PNI,PLI,ALFA,CL
      DIMENSION PNI(42),PLI(21)
      SAVE ALFA
      DIMENSION CL(3)
      DIMENSION ALFA(10)
      DATA ALFA/0.05971587,0.47014206,0.79742699,0.10128651,
     1          0.33333333,0.0,0.5,1.0,0.0,0.33333333/
C
      INP=0
      INL=0
      I1=1
      I2=2
      I3=3
      DO 30 K=1,7
C
C   set up coordinates for K near vertex
C
      IF(K.LE.3)THEN
      CL(1)=ALFA(4)
      CL(2)=ALFA(4)
      CL(3)=ALFA(4)
      CL(K)=ALFA(3)
      GO TO 20
      END IF
C
C    set up coordinates for K near midpoint
C
      IF(K.LE.6)THEN
      IX=MOD(K,3)+1
      CL(1)=ALFA(2)
      CL(2)=ALFA(2)
      CL(3)=ALFA(2)
      CL(IX)=ALFA(1)
      GO TO 20
      END IF
C
C    K is the centroid
C
      CL(1)=ALFA(5)
      CL(2)=ALFA(5)
      CL(3)=ALFA(5)
C
C    PLI and PNI are arranged in rows of 7, each row consists
C    of interpolation function Li, Ni values at 7 integration point
C    Row 1 is for I=1, through to row 3 for L, row 6 for N.
C
   20 INP=INP+1
      INL=INL+1
      PLI(INL)   =CL(I1)
      PLI(INL+7) =CL(I2)
      PLI(INL+14)=CL(I3)
      PNI(INP)   =CL(I1)*CL(I1) - CL(I1)*(CL(I2)+CL(I3))
      PNI(INP+7) =CL(I2)*CL(I2) - CL(I2)*(CL(I3)+CL(I1))
      PNI(INP+14)=CL(I3)*CL(I3) - CL(I3)*(CL(I1)+CL(I2))
      PNI(INP+21)=4.0*CL(I3)*CL(I1)
      PNI(INP+28)=4.0*CL(I1)*CL(I2)
      PNI(INP+35)=4.0*CL(I2)*CL(I3)
   30 CONTINUE
C
C    An alternative coding for the above group of statements
C    follows - yet untested
C
C      DO 30 I1=1,3
C      INL1=INL+7*(L-1)
C      PLI(INL1)=CL(L)
C      I2=MOD(I1,3)+1
C      I3=MOD(I2,3)+1
C      INP1=INP+7*(L-1)
C      INP2=INP+21+7*(L-1)
C      PNI(INP1)=CL(I1)*CL(I1) - CL(I1)*(CL(I2)+CL(I3))
C      PNI(INP2)=4.0*CL(I3)*CL(I1)
C   30 CONTINUE
C      WRITE(LUW,*)'Interpolation functions Li follow'
C      CALL MATPRT(PLI,7,3,21,LUW)
C      WRITE(LUW,*)'Interpolation functions Ni follow'
C      CALL MATPRT(PNI,7,6,42,LUW)
      RETURN
      END
C
      SUBROUTINE DNCOM(JJ,BY,CX,DNDP)
C
      DOUBLE PRECISION BY,CX,DNDP
      DOUBLE PRECISION ALFA,CL1,CL2,CL3
C
      DIMENSION BY(3),CX(3),DNDP(84)
      SAVE ALFA
      DIMENSION ALFA(10)
      DATA ALFA/0.05971587,0.47014206,0.79742699,0.10128651,
     10.33333333,0.0,0.5,1.0,0.0,0.33333333/
C
C     This routine calculates the gradient of the interpolation
C     function N(I) at a location in the element determined by
C     the parameter J.  The calculation for X and Y gradient is
C     the same, though different geometrical elements in B are
C     used.  For numbering the points (I & J), vertices are 1-3
C     midpoints are 4-6, and the centroid (J only) is 7.  The
C     interpolation functions are given in Huebner (p345).
C
C     The gradient values DNDP produced by this routine need to be
C     divided by the factor of 2*(triangle area) to produce actual
C     spatial gradients.  This factor is introduced in matrix assembly.
C
C     The interpolation function is for a vertex point
C
      IND=1
      DO 30 I=1,3
        DO J=1,7
C
C     J  is a vertex point
C
          IF(J.LE.3)THEN
            KN2=4
            IF(J.EQ.I)THEN
              KN1=3
              KN3=4
            ELSE
              KN1=4
              KN3=3
            END IF
            GO TO 20
          END IF
C
C    J is a midpoint
C
          IF(J.LE.6)THEN
            IX=MOD(I+1,3)+4
            KN3=2
            IF(J.EQ.IX)THEN
              KN1=1
              KN2=2
            ELSE
              KN1=2
              KN2=1
            END IF
            GO TO 20
          END IF
C
C     J is the centroid
C
          KN1=5
          KN2=5
          KN3=5
   20     I2=MOD(I,3)+1
          I3=MOD(I+1,3)+1
          CL1=ALFA(KN1+JJ)
          CL2=ALFA(KN2+JJ)
          CL3=ALFA(KN3+JJ)
          DNDP(IND)  =BY(I)*(2.D0*CL1-CL2-CL3)-CL1*(BY(I2)+BY(I3))
          DNDP(IND+1)=CX(I)*(2.D0*CL1-CL2-CL3)-CL1*(CX(I2)+CX(I3))
          IND=IND+2
        ENDDO          ! on J
   30 CONTINUE         ! on I
C
C    The interpolation function is for a midpoint
C
      DO 230 I=4,6
        DO J=1,7
          I1=MOD(I,3)+1
          I2=MOD(I1,3)+1
          I3=MOD(I1+1,3)+1
C
C    J is a vertex point
C
          IF(J.LE.3)THEN
            IF(J.EQ.I2)THEN
              KN2=3
              KN3=4
            ELSE IF(J.EQ.I3)THEN
              KN2=4
              KN3=3
            ELSE
              KN2=4
              KN3=4
            END IF
            GO TO 220
          END IF
C
C    J is a mid-point node
C
          IF(J.LE.6)THEN
            I5 = 4 + MOD(I,3)
            I6 = 4 + MOD(I+1,3)
            IF(J.EQ.I5)THEN
              KN2=1
              KN3=2
            ELSE IF(J.EQ.I6)THEN
              KN2=2
              KN3=1
            ELSE
              KN2=2
              KN3=2
            END IF
            GO TO 220
          END IF
C
C    J is the centroid
C
          KN2=5
          KN3=5
C
  220     CL2=ALFA(KN2+JJ)
          CL3=ALFA(KN3+JJ)
          DNDP(IND)   = 4.D0*(CL3*BY(I2) + CL2*BY(I3))
          DNDP(IND+1) = 4.D0*(CL3*CX(I2) + CL2*CX(I3))
          IND=IND+2
        ENDDO               !  on J
  230 CONTINUE              !  on I
      RETURN
      END
C
      SUBROUTINE LOADBC(BIG,XLOAD,QLOAD,QBND,IBC,IBCTYP,
     :                  NUP,NBP,NROWS)
C
C    Modifies the load vector obtained from PLOAD by fixing the b.c.'s
C
      DOUBLE PRECISION XLOAD,BIG
      DIMENSION XLOAD(NROWS)
      DIMENSION IBC(NBP),IBCTYP(NBP*2)
      DIMENSION QLOAD(NROWS)
      DIMENSION QBND(NBP*2)
C
      DO I=1,NROWS
        XLOAD(I)=QLOAD(I)
      ENDDO
C
C    Fix the velocity and stress boundary conditions
C
C     WRITE(*,*)'I,QBND(I),QBND(I+NBP),ITYPX,ITYPY'
      DO 200 I=1,NBP
        NNX=IBC(I)
        NNY=NNX+NUP
        ITYPX=IBCTYP(I)
        ITYPY=IBCTYP(I+NBP)
C
C    Applied velocity boundary condition
C
        IF(ITYPX.EQ.0.OR.ITYPX.EQ.10)XLOAD(NNX)=BIG*QBND(I)
        IF(ITYPY.EQ.0.OR.ITYPY.EQ.10)XLOAD(NNY)=BIG*QBND(I+NBP)
C
C    Applied traction boundary condition
C
        IF(ITYPX.EQ.1.OR.ITYPX.EQ.12.OR.ITYPX.EQ.14)
     :                       XLOAD(NNX)=XLOAD(NNX)+QBND(I)
        IF(ITYPY.EQ.1.OR.ITYPY.EQ.12.OR.ITYPY.EQ.14)
     :                       XLOAD(NNY)=XLOAD(NNY)+QBND(I+NBP)
C       WRITE(*,*)I,QBND(I),QBND(I+NBP),ITYPX,ITYPY
  200 CONTINUE
C     WRITE(*,*)'XLOAD from LOADBC, BIG =',BIG
C     CALL DMATPP(XLOAD,NROWS,6)
      RETURN
      END
C
      SUBROUTINE VISK(BIG,JLV,HLENSC,AREASP,VOLUMSP,VISCDSP,DVMX,
     :                SE,EX,EY,LEM,NOR,SSQ,UVP,VISN,VHB,THDISS,
     :                THDINT,GAMMAD,TREFD,PNI,PLI,NCOMP,IVV,
     :                IVIS,VC,
     :                VISCMIN,VISCMAX,ED2IMIN,ED2IMAX,
     :                THDIMIN,THDIMAX,
     :                NUP,NE,NROWS,LUW,LSC,IERR)
C
C    Routine to calculate the viscosity at the seven integration
C    point locations in element no. N.  A non-Newtonian fluid
C    with stress exponent SE is assumed
C
      COMMON/SSQVAL/ISSQACTIVE,IROTACTIVE,DFLTSSQ
C
      DOUBLE PRECISION VISN,UVP     !  UVP is XWRK in CGRUN
      DOUBLE PRECISION WGT,XX,YY,BY,CX,YCOTK,PNI,PLI,DNDP
      DOUBLE PRECISION TRIA,XLA,YLA,DNIDX,DNIDY
      DOUBLE PRECISION BIG,VOLUM,VISCD,AREA,DMIN,BIGV,BIGV10,BIGINV
      DOUBLE PRECISION SELOC,EDPOW,HITE,VISDIS,DUDX,DUDY,DVDX,DVDY
      DOUBLE PRECISION UI,VI,UUI,VVI,EXL,PNIKI,ADVDX,EDXY,DWDZ,EDO,EDOT
      DOUBLE PRECISION AVIS,OVIS,DVIS,DISSIP,HEXL,VO,VN
      DOUBLE PRECISION GAMMAD,TREFD,DENOMR,DNUMER,WRKTRM,BLOGD,TEMDEP
      DOUBLE PRECISION BDISS,DAREA,TOOPI
      DIMENSION VISN(7,NE),UVP(NROWS)
      DIMENSION PNI(7,6),PLI(7,3)
      DIMENSION WGT(7)
      DIMENSION BY(3),CX(3),DNDP(84)
      DIMENSION XX(3),YY(3)
C   next arrays not allocated unless GAMMA .ne. 0
      DIMENSION THDISS(7,NE),THDINT(7,NE)
C
      DIMENSION EX(NUP)
      DIMENSION EY(NUP)
      DIMENSION NOR(NUP)
      DIMENSION LEM(6,NE)
      DIMENSION SSQ(NUP)
      DIMENSION VHB(8,NE)
C
      SAVE WGT
      DATA WGT/0.12593918,0.12593918,0.12593918,0.13239415,
     1         0.13239415,0.13239415,0.225/
C
      TOOPI=2.d0*3.1415926535d0
      VOLUM=0.D0
      VISCD=0.D0
      AREA=0.D0
      DMIN=1.D0
      DVMX=0.0     ! single precision
      HEXL=1.D0
      ED2IMIN=9.999e+35
      ED2IMAX=0.0
      VISCMIN=ED2IMIN
      VISCMAX=ED2IMAX
      THDIMIN=ED2IMIN
      THDIMAX=ED2IMAX
C
C    A maximum viscosity of BIG/10000. is permitted
C
C     BIGV=1.0D+20
      BIGV=BIG/1.0D4
      BIGV10=BIGV*0.01D0
      BIGINV=1.D0/BIGV
      SELOC=SE
      EDPOW=(1.D0/SELOC)-1.D0
      IF(ABS(SE-1.0).LE.1.E-5)EDPOW=0.0
C
C    Loop over elements; look up exponent for each element if required
C
      DO 60 N=1,NE
        IF(IVV.GE.4)THEN
          SELOC=VHB(8,N)
          EDPOW=(1.D0/SELOC)-1.D0
        ENDIF
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
          XX(K1)=EX(LK1)
          YY(K1)=EY(LK1)
          BY(K1)=EY(LK2)-EY(LK3)
          CX(K1)=EX(LK3)-EX(LK2)
   10   CONTINUE
        TRIA=CX(3)*BY(2)-CX(2)*BY(3)
C
C   TRIA is 2 times the area of the triangle element
C   Check for inverted elements indicating mesh distortion
C
        IF (TRIA.LT.0.D0) THEN
          IERR = 2
          WRITE(LUW,11000)LK1,EX(LK1),EY(LK1),EX(LK2),EY(LK2),
     :                      EX(LK3),EY(LK3)
          RETURN
        ENDIF
11000   FORMAT('ELEMENT ',I5,' INVERTED: VERTEX POSITIONS '/,
     :         2G14.5/,2G14.5/,2G14.5/)
C
C   DNCOM computes dNi/dy & dNj/dx (times 2*area)
C
        CALL DNCOM(0,BY,CX,DNDP)
C
        HITE=0.D0
        VISDIS=0.D0
        DO 62 K=1,7
          OVIS=VISN(K,N)
          DUDX=0.D0
          DVDX=0.D0
          DUDY=0.D0
          DVDY=0.D0
          UUI=0.D0
          VVI=0.D0
          EXL=DFLTSSQ
          IF(JLV.NE.0)EXL=0.D0
C
C    Calculate the interpolated value of the strain rate
C    and the crustal thickness (logarithm) at point K
C
          DO 55 I=1,6
            NI=LEM(I,N)
            KIN=(I-1)*14 + (K-1)*2 +1
            PNIKI=PNI(K,I)
C
C   JLV.ne.0 means that viscosity is scaled by layer thickness.
C   ISSQACTIVE should be set = 1 if JLV.ne.0 (or other reasons)
C   DFLTSSQ should be set to 1.0/HLENSC
C
            IF (JLV.NE.0) EXL=EXL + SSQ(NI)*PNIKI
            DNIDX=DNDP(KIN)
            DNIDY=DNDP(KIN+1)
            UI=UVP(NI)
            VI=UVP(NI+NUP)
            DUDX=DUDX + UI*DNIDX
            DUDY=DUDY + UI*DNIDY
            DVDX=DVDX + VI*DNIDX
            DVDY=DVDY + VI*DNIDY
            IF(NCOMP.EQ.-1)THEN
              UUI=UUI+UI*PNIKI
              VVI=VVI+VI*PNIKI
            ELSE IF(NCOMP.EQ.2)THEN
              UUI=UUI+UI*PNIKI
            END IF
   55     CONTINUE
C
C    DVMX is used to check timestep
C
          ADVDX = (DABS(DUDX) + DABS(DVDY))/TRIA
          IF(ADVDX.GT.DVMX)DVMX=ADVDX
          EDXY=0.5D0*(DUDY+DVDX)
C
C    Corrections to the strain-rate expressions for the
C     thin spherical sheet (NCOMP=-1) and
C     cylindrical axisymmetry (NCOMP=2)
C
          IF(NCOMP.EQ.-1)THEN
            XLA=XX(1)*PLI(K,1)+XX(2)*PLI(K,2)+XX(3)*PLI(K,3)
            YLA=YY(1)*PLI(K,1)+YY(2)*PLI(K,2)+YY(3)*PLI(K,3)
            YCOTK=1.D0/DTAN(YLA)
            DUDX=DUDX+VVI*YCOTK*TRIA
            DVDY=DVDY+XLA*YCOTK*DVDX
            EDXY=EDXY+0.5D0*YCOTK*(XLA*DUDX-UUI*TRIA)
          ELSE IF(NCOMP.EQ.2)THEN
            XLA=XX(1)*PLI(K,1)+XX(2)*PLI(K,2)+XX(3)*PLI(K,3)
          END IF
          DWDZ=0.D0
          IF(NCOMP.LE.0)DWDZ=-DUDX-DVDY
          IF((NCOMP.EQ.2).AND.(XLA.GT.0.D0))DWDZ=TRIA*UUI/XLA
          EDO=DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ + 2.D0*EDXY*EDXY
C
C     the second invariant of the strain rate tensor is :
C
          EDOT=(DSQRT(EDO))/TRIA
          IF(EDOT.LT.ED2IMIN)ED2IMIN=EDOT
          IF(EDOT.GT.ED2IMAX)ED2IMAX=EDOT
C
C      the viscosity is computed using the local strain-rate:
C      VHB is interpreted as a viscosity coefficient, eta=B/2
C      with EDOT = sqrt(epsij*epsij); tau_ij = 2.*eta*eps_ij
C      eta = VHB*(EDOT**EDPOW)
C
          AVIS=VC          ! default viscosity if VHB not used
          IF(IVV.GE.2)THEN
            IF(EDOT.GT.0.D0) THEN
              AVIS=EDOT**EDPOW
            ELSE
              AVIS=BIGV
            END IF
          ENDIF
C
C     thermal dissipation per unit volume
C
          DISSIP=2.D0*AVIS*EDOT*EDOT
C 
C     a maximum viscosity is permitted but is approached gradually
C
C         IF(AVIS.GT.BIGV10) AVIS=1.D0/(BIGINV+(1.D0/AVIS))
          IF(AVIS.GT.BIGV) AVIS=BIGV
          VISN(K,N)=AVIS
C
C     Viscosity is multiplied by the normalised crustal thickness
C     if required (JLV.NE.0). To keep matrix elements O(1) Argand no. 
C     is divided by HLENSC in main program
C
          IF(JLV.NE.0) THEN
            HEXL=DEXP(EXL)*HLENSC
            DISSIP=DISSIP*HEXL
            VISN(K,N)=AVIS*HEXL
          ENDIF
C
C    use the integrated work (THDINT) to modulate the strength coefficient
C
          IF(GAMMAD.NE.0.D0)THEN
            WRKTRM=GAMMAD*THDINT(K,N)
            IF(TREFD.EQ.0.D0)THEN       ! if GAMMA is constant
              TEMDEP=DEXP(-WRKTRM)
              BDISS=VHB(K,N)*TEMDEP
            ELSE                        ! if GAMMA is temperature dependent
              BLOGD=LOG(VHB(K,N))
              DNUMER=BLOGD+(TREFD*BLOGD-1.D0)*WRKTRM
              DENOMR=1.D0+TREFD*WRKTRM
              BDISS=DEXP(DNUMER/DENOMR)
            ENDIF
C
C    set effective local viscosity and save thermal dissipation rate
C
            VISN(K,N)=VISN(K,N)*BDISS
            DISSIP=DISSIP*BDISS
            THDISS(K,N)=DISSIP
C
C    if no strain weakening
C
          ELSEIF(IVIS.GT.0)THEN
            VISN(K,N)=VISN(K,N)*VHB(K,N)
            DISSIP=DISSIP*VHB(K,N)
          ENDIF
        IF(VISN(K,N).LT.VISCMIN)VISCMIN=VISN(K,N)
        IF(VISN(K,N).GT.VISCMAX)VISCMAX=VISN(K,N)
        IF(DISSIP.LT.THDIMIN)THDIMIN=DISSIP
        IF(DISSIP.GT.THDIMAX)THDIMAX=DISSIP
C
C    Find the point where the viscosity change is maximum
C
        DVIS=AVIS/OVIS
        IF(DVIS.GT.1.D0)DVIS=1.D0/DVIS
        IF(DVIS.LE.DMIN)THEN
          DMIN=DVIS
          KV=K
          NV=N
          VO=OVIS
          VN=AVIS
        END IF
C
C    Integrate the thermal dissipation, volume and area
C
        VISDIS=VISDIS + WGT(K)*DISSIP
        IF(NCOMP.LE.0)THEN
          HITE=HITE + WGT(K)*HEXL
        ELSEIF(NCOMP.EQ.2)THEN
          HITE=HITE + WGT(K)*XLA
        ENDIF
   62 CONTINUE
      DAREA=TRIA*0.5D0
      VISCD = VISCD + VISDIS*DAREA
      VOLUM = VOLUM + HITE*DAREA
      AREA = AREA + DAREA
   60 CONTINUE
      IF(NCOMP.LE.0)THEN
        VISCD = VISCD/HLENSC
        VOLUM = VOLUM/HLENSC
      ELSEIF(NCOMP.EQ.2)THEN
        VISCD=VISCD*TOOPI
        VOLUM=VOLUM*TOOPI
      ENDIF
C
C   return single precision variables
C
      AREASP = AREA
      VISCDSP = VISCD
      VOLUMSP = VOLUM
      IDBUG=1
      IF(IDBUG.NE.1)THEN
        WRITE(LUW,*)'VISN matrix follows'
        CALL DMATPRT(VISN,7,NE,7*NE,LUW)
      END IF
      RETURN
      END
C
      SUBROUTINE VISCLIN(VISN,VHB,VC,NE,IVIS)
C
C    routine to load VHB into VISN if VISK is not used
C    (linear viscosity) used in matrix assembly
C
      DOUBLE PRECISION VISN
      DIMENSION VISN(7,NE),VHB(8,NE)
      IF(IVIS.GT.0)THEN
        DO J=1,NE
          DO K=1,7
            VISN(K,J)=VHB(K,J)
          ENDDO
        ENDDO
      ELSE
        DO J=1,NE
          DO K=1,7
            VISN(K,J)=VC
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END

