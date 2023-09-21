C*--------------------------------------------------------------------
C*    Basil / Sybil:   conjgr.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE CONJGR(XDP,PDP,KFIN,AC,NCOMP,STK,STK2,STK3,STK4,PWTS,
     :                  IK,IK2,IK4,NK,NK2,NK4,KBW,K2BW,K4BW,
C    :                  D1,D2,D3,R,RTIL,P,WORK,
     :                  NUP,NN,NROWS,NFP,IFLT,NITER,LUW,LSC,IVERB,IERR)
C
C    conjugate gradient method used to solve matrix equation
C    Ax=b. XDP initially holds b, but is replaced by successive
C    versions of the solution XDP. PDP is the initial guess for x.
C    STK and STK2 are stored in compressed format, row-wise, with column
C    indices stored in IK and IK2.  Similarly STK4 is stored in a
C    compressed format, column-wise, with row indices stored in IK4.
C    Each term in STK3 is a diagonal element since there are non-zeros
C    only on the diagonal.  NROWS is the number of rows in complete matrix,
C    and MAXROW is the maximum no. of non-zero elements in any 
C    row of matrix. AC is the accuracy criterion of the convergence test.
C
C    The matrix equation Ax=b has the form:
C
C      [ K     K2    K4 ] [ u ] 
C      [                ] [   ]
C      [ K2T   K3    0  ] [ s ] = X
C      [                ] [   ]
C      [ K4T   0     0  ] [ p ]
C
C     (K, K2, K3, and K4 represent STK, STK2, STK3 and STK4)
C     where u = the velocity components, s = the stress integrals, and
C       p = the pressure
C
C    The second version of this routine uses a primitive
C    form of pre-conditioning, in which the aim is to rescale
C    the diagonal elements of the matrix to unity.  The diagonal 
C    elements are found in the first column of STK.
C
      DOUBLE PRECISION PAPDN,ALPHA,BETA,AC2,RRP,PBAR
      DOUBLE PRECISION RESIUV,RESIST,RESIPR,RESID
      DOUBLE PRECISION PDP,R,XDP,RTIL,WORK,D1,D2,D3,PWTS
      DOUBLE PRECISION STK,STK2,STK3,STK4
      DIMENSION R(NROWS),RTIL(NROWS),PDP(NROWS),XDP(NROWS),WORK(NROWS)
      DIMENSION STK(NK),STK2(NK2),STK3(NFP),STK4(NK4)
      DIMENSION IK(NK),IK2(NK2),IK4(NK4)
      DIMENSION D1(NUP*2),D2(NFP),D3(NN)
      DIMENSION PWTS(NN)
C
C    set up initial vectors
C
      NUP2=NUP*2
      AC2=10.0*AC*AC
C      AC2=1.0*AC*AC
C
C    first matrix-vector multiply
C
      CALL MATVEC(PDP,R,NCOMP,STK,STK2,STK3,STK4,
     :            IK,IK2,IK4,NK,NK2,NK4,KBW,K2BW,K4BW,
     :            NUP,NN,NFP,NROWS,IFLT)
      DO 20 I=1,NROWS
        R(I)=XDP(I)-R(I)
        XDP(I)=PDP(I)
        WORK(I)=0.0
   20 CONTINUE
C
C    normalise residual by condition matrix M
C
      CALL PRECON(RTIL,R,NROWS,1,NCOMP,STK,STK2,STK3,STK4,
     :            D1,D2,D3,IK,IK2,IK4,
     :            NK1,NK2,NK4,KBW,K2BW,K4BW,
     :            NUP,NN,NFP,IFLT)
      DO 30 I=1,NROWS
        PDP(I)=RTIL(I)
   30 CONTINUE
C
C     calculate the norm of the residual
C
C     RESID=0.0
C     DO 40 I=1,NROWS
C       RESID=RESID+RTIL(I)*R(I)
C  40 CONTINUE
      RESIUV=0.0
      RESIST=0.0
      RESIPR=0.0
      DO I=1,NUP2
        RESIUV=RESIUV+RTIL(I)*R(I)
      ENDDO
      IF(NFP.GT.0)THEN
        DO I=NUP2+1,NUP2+NFP
          RESIST=RESIST+RTIL(I)*R(I)
        ENDDO
      ENDIF
      IF(NCOMP.GT.0)THEN
        DO I=NUP2+NFP+1,NROWS
          RESIPR=RESIPR+RTIL(I)*R(I)
        ENDDO
      ENDIF
      RESID=RESIUV+RESIST+RESIPR
C
C    commence the iteration cycle
C
      K=0
   50 K=K+1
C
C     calculate alpha = (r.r)/(p.Ap)
C
      CALL MATVEC(PDP,WORK,NCOMP,STK,STK2,STK3,STK4,
     :            IK,IK2,IK4,NK,NK2,NK4,KBW,K2BW,K4BW,
     :            NUP,NN,NFP,NROWS,IFLT)
      PAPDN=0.0
      DO I=1,NROWS
        PAPDN=PAPDN+PDP(I)*WORK(I)
      ENDDO
C
C     update solution vector x and residual vector r
C
      ALPHA=RESID/PAPDN
      DO I=1,NROWS
        XDP(I)=XDP(I)+ALPHA*PDP(I)
        R(I)=R(I)-ALPHA*WORK(I)
      ENDDO
C
C    if one of the variables is not constrained, remove average.
C    the condition for activating this step needs more thought as
C    whether a variable is adequately constrained depends on the
C    applied boundary conditions. But anyway, idea does not work.
C
C     IF(NCOMP.GE.1)CALL DEMEANP(1,XDP,PBAR,NROWS,NUP,NN,PWTS)
C
C    normalise residual by condition matrix M
C
      CALL PRECON(RTIL,R,NROWS,0,NCOMP,STK,STK2,STK3,STK4,
     :            D1,D2,D3,IK,IK2,IK4,
     :            NK1,NK2,NK4,KBW,K2BW,K4BW,
     :            NUP,NN,NFP,IFLT)
C
C     calculate beta = (r(k+1).r(k+1))/(r(k).r(k))
C
      RRP=RESID
C     RESID=0.0
C     DO I=1,NROWS
C       RESID=RESID+RTIL(I)*R(I)
C     ENDDO
      RESIUV=0.0
      RESIST=0.0
      RESIPR=0.0
      DO I=1,NUP2
        RESIUV=RESIUV+RTIL(I)*R(I)
      ENDDO
      IF(NFP.GT.0)THEN
        DO I=NUP2+1,NUP2+NFP
          RESIST=RESIST+RTIL(I)*R(I)
        ENDDO
      ENDIF
      IF(NCOMP.GT.0)THEN
        DO I=NUP2+NFP+1,NROWS
          RESIPR=RESIPR+RTIL(I)*R(I)
        ENDDO
      ENDIF
      RESID=RESIUV+RESIST+RESIPR
      BETA=RESID/RRP
C
C     if satisfactory convergence return solution
C
C     CALL DATIME(IDAY,ISEC,1)
      IF(IVERB.GT.0)THEN
        IF(MOD(K,IVERB).EQ.0)
     :            WRITE(LSC,10002)K,RESIUV,RESIST,RESIPR,BETA
      ENDIF
      IF((RESIUV.LT.AC2).AND.(RESIPR.LT.AC2).AND.
     :                       (ABS(RESIST).LT.AC2))THEN
        KFIN=K
        RETURN
      END IF
C
C     in case convergence is not happening print last iterations
C
      IF(K.GT.NITER-30)THEN
        WRITE(LSC,10002)K,RESIUV,RESIST,RRP,BETA
      ENDIF
10002 FORMAT('CONJGR: k = ',I5,' RESIUV = ',G12.5,'  RESIST = ',G12.5,
     1' RESIPR =',G12.5,' beta = ',G12.5)
C
C    update the direction vector p
C
      DO I=1,NROWS
        PDP(I)=RTIL(I)+BETA*PDP(I)
      ENDDO
C
C     return for the next iteration unless too many 
C
      IF(K.LE.NITER)GO TO 50
      WRITE(LSC,10001)K,BETA
      WRITE(LUW,10001)K,BETA
10001 FORMAT(' Convergence failed in CONJGR. K =',I5,
     :       ' beta (step) =',G12.5)
      IERR=1
      RETURN
      END
      SUBROUTINE MATVEC(VEC,RES,NCOMP,STK,STK2,STK3,STK4,
     :                  IK,IK2,IK4,NK1,NK2,NK4,KBW,K2BW,K4BW,
     :                  NUP,NN,NFP,NROWS,IFLT)
C
C    this routine called by CONJGR, to post-multiply stiffness matrix 
C    by vector VEC.  Matrix is stored in compressed mode, row-wise
C    with column indices stored in IK.  Only the upper diagonal
C    half of matrix is stored. NROWSP=NUP2P+NFPP+NNP is total no. of rows
C    in matrix.  KBW, K2BW and K4BW define the no of non-zero elements
C    in any row of matrix.
C
C     INCLUDE 'cg.parameters'
      DOUBLE PRECISION RES,VEC,STK,STK2,STK3,STK4
      DIMENSION VEC(NROWS),RES(NROWS)
      DIMENSION STK(NK1),STK2(NK2),STK3(NFP),STK4(NK4)
      DIMENSION IK(NK1),IK2(NK2),IK4(NK4)
C
      NUP2=NUP*2
      DO I=1,NROWS
        RES(I)=0.0
      ENDDO
C
C    for each row in the STK array
C
      DO 200 IROW1=1,NUP2
      INDEX=(IROW1-1)*KBW
C
C    for all entries in this row
C
      DO 150 I=1,KBW
      INDEX=INDEX+1
      ICOL1=IK(INDEX)
      IF(ICOL1.LE.0)GO TO 200
C
C     dot product terms derived from upper half of STK
C
      RES(IROW1)=RES(IROW1)+STK(INDEX)*VEC(ICOL1)
C
C     dot product terms derived from lower half of STK
C
      IF(ICOL1.NE.IROW1)THEN
      RES(ICOL1)=RES(ICOL1)+STK(INDEX)*VEC(IROW1)
      END IF
  150 CONTINUE
  200 CONTINUE
C ----------------------------------------
C FOR THE FAULT CASE (INFLT.NE.0)
C ----------------------------------------
      IF(IFLT.NE.0) THEN
C
C     now add in the cross-component contributions from
C      the STK2 matrix and its transpose
C
        DO 400 IROW2=1,NUP2
          INDEX=(IROW2-1)*K2BW
C
C    for all entries in this row
C
          DO 350 I=1,K2BW
            INDEX=INDEX+1
            ICOL2=IK2(INDEX)
            IF(ICOL2.LE.0)GO TO 400
            ITRAN2=ICOL2+NUP2
C
C     dot product terms derived from K2.v
C
            RES(IROW2)=RES(IROW2)+STK2(INDEX)*VEC(ITRAN2)
C
C     dot product terms derived from (transpose of K2).u
C
            RES(ITRAN2)=RES(ITRAN2)+STK2(INDEX)*VEC(IROW2)
  350     CONTINUE
  400   CONTINUE
C
C   now add in the terms derived from the STK3 (diagonal) matrix
C
        DO 500 I=1,NFP
          IROW3=I+NUP2
          RES(IROW3)=RES(IROW3)+STK3(I)*VEC(IROW3)
  500   CONTINUE
      ENDIF               ! IFLT != 0
C
C ----------------------------------------
C FOR THE INCOMPRESSIBLE CASE (NCOMP.GT.0)
C ----------------------------------------
      IF(NCOMP.LE.0)RETURN
C
C     now add in the cross-component contributions from
C      the STK4 matrix and its transpose.
C      for incompressible flow only.  Note that
C      STK4 is stored columnwise, whereas other matrices
C      are stored rowwise.
C
      DO 600 ICOL4=1,NN
      INDEX=(ICOL4-1)*K4BW
C
C    for all entries in this row
C
      DO 550 I=1,K4BW
      INDEX=INDEX+1
      IROW4=IK4(INDEX)
      IF(IROW4.LE.0)GO TO 600
      ITRAN4=ICOL4+NUP2+NFP
C
C     dot product terms derived from K4.p
C
      RES(IROW4)=RES(IROW4)+STK4(INDEX)*VEC(ITRAN4)
C
C     dot product terms derived from K4(t).u
C
      RES(ITRAN4)=RES(ITRAN4)+STK4(INDEX)*VEC(IROW4)
  550 CONTINUE
  600 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PRECON(RTIL,R,NROWS,JTIME,NCOMP,STK,STK2,STK3,
     :                  STK4,D1,D2,D3,IK,IK2,IK4,
     :                  NK1,NK2,NK4,KBW,K2BW,K4BW,
     :                  NUP,NN,NFP,IFLT)
C
C    normalise residual by condition matrix M
C
      DOUBLE PRECISION R,RTIL,D1,D2,D3
      DOUBLE PRECISION STK,STK2,STK3,STK4
      DIMENSION D1(NUP*2),D2(NFP),D3(NN)
      DIMENSION R(NROWS),RTIL(NROWS)
      DIMENSION STK(NK1),STK2(NK2),STK3(NFP),STK4(NK4)
      DIMENSION IK(NK1),IK2(NK2),IK4(NK4)
C
      NUP2=NUP*2
C
C    if JTIME=1 the diagonal matrices are defined
C
      IF(JTIME.EQ.1)THEN
      DO 10 I1=1,NUP2
        KC=KBW*(I1-1)+1
        IF (STK(KC).NE.0.0) THEN
          D1(I1)=1.0/STK(KC)
        ELSE
          D1(I1)=1.0
        END IF
   10 CONTINUE
C ---------------------------
C FOR THE FAULT CASE
C ---------------------------
      IF(IFLT.NE.0) THEN
        DO I=1,NFP
          D2(I)=STK3(I)
        ENDDO
C  the following structure was not working / helping
C     DO I2=1,NUP2
C       INDEX=(I2-1)*K2BW
C       DO I3=1,K2BW
C         INDEX=INDEX+1
C         K2=IK2(INDEX)
C         IF (K2.LE.0) GOTO 30
C         D2(K2)=D2(K2)+STK2(INDEX)*STK2(INDEX)*D1(I2)
C       ENDDO
C  30   CONTINUE
C     ENDDO
        DO I=1,NFP
          IF(D2(I).NE.0.0) THEN
             D2(I)=1.0/D2(I)
          ELSE 
            D2(I)=1.0
          END IF
        ENDDO
      ENDIF
C ---------------------------
C FOR THE INCOMPRESSIBLE CASE
C ---------------------------
      IF(NCOMP.GT.0)THEN
      DO 200 I2=1,NN
         D3(I2)=0.0
         INDEX=(I2-1)*K4BW
      DO 210 I3=1,K4BW
         INDEX=INDEX+1
         K4=IK4(INDEX)
         IF(K4.LE.0)GO TO 200
         D3(I2)=D3(I2) + STK4(INDEX)*STK4(INDEX)*D1(K4)
  210 CONTINUE
  200 CONTINUE
      DO 250 I2=1,NN
         IF(D3(I2).NE.0.0) THEN
            D3(I2)=1.0/D3(I2)
         ELSE
            D3(I2)=1.0
         END IF
  250 CONTINUE
      END IF
C---------------------------
      END IF
C
C    for thin viscous sheet calculation preconditioner:
C    simply divide by diagonal elements
C
      DO I1=1,NUP2
        RTIL(I1)=R(I1)*D1(I1)
      ENDDO
      IF(IFLT.NE.0) THEN
        DO I1=1,NFP
          I2=NUP2+I1
          RTIL(I2)=R(I2)*D2(I1)
        ENDDO
      ENDIF
C
      IF(NCOMP.LE.0)RETURN
      DO I1=1,NN
        I2=NUP2+NFP+I1
        RTIL(I2)=R(I2)*D3(I1)
      ENDDO
      RETURN
      END
C
C****************************************************************************
C
      SUBROUTINE DEMEANP(IUV,X,PBAR,NROWS,NUP,NN,PWTS)
C
C      calculates mean value of an insufficiently constrained variable
C      and subtracts that value from current entries of the solution vector X.  
C      The weights (PWTS) set up in SEMBL mean that each value is weighted 
C      by the total area of all the elements that contain that vertex node.
C      only vertex node values are used.
C      The total (WTSUM) comes out at 12 times the actual area
C      because TRIA is 4*area and there are 3 nodes per triangle.
C      IUV = 0 for x-velocity, 1 for y-velocity, 2 for pressure
C
      DOUBLE PRECISION X,PWTS,PBAR,WTSUM
      DIMENSION X(NROWS),PWTS(NN)
C
      K0=IUV*NUP
      PBAR = 0.D0
      WTSUM = 0.D0
      DO J=1,NN
        PBAR=PBAR+X(K0+J)*PWTS(J)
        WTSUM=WTSUM+PWTS(J)
      ENDDO
      PBAR=PBAR/WTSUM
C
C     now remove PBAR from the X vector
C
      DO J=1,NN
        X(K0+J)=X(K0+J)-PBAR
      ENDDO
C
      RETURN
      END
