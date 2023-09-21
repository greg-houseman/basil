C*--------------------------------------------------------------------
C*    Basil / Sybil:   plmesh.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
C      contains subroutines PLMESH, STMARK, ELLPAR, INTBND, PLSEG
C
      SUBROUTINE PLMESH(LABEL,NE,NN,NUP,NBP,EX,EY,LEM,NOR,IBC,
     1                  IBCTYP,VHB,TBXOFF,TBYOFF,IFLT,MPE,MINCOL,
     2                  IRANGECOL,XMINC,XMAXC,YMINC,YMAXC,
     3                  VMIN,VMAX,ZOOM,NCOMP,XMID,YMID,IVERBOSE)
C
C    This routine draws the finite element mesh, which is
C    in store, in order to check its validity, or to show
C    the deformation.
C
      INCLUDE "sybilps.parameters"
      CHARACTER NUMBER*80
      DIMENSION EX(NUP),EY(NUP),VHB(8,NE)
      DIMENSION LEM(6,NE),NOR(NUP),IBC(NBP),IBCTYP(NBP*2)
      DIMENSION POLYX(3),POLYY(3),MENV(3)
C
C     LABEL = 0  - for no labels
C           = 1  - for element numbers - not implemented
C           = 2  - for node numbers - not implemented
C           = 3  - for boundary only
C           = 4  - double mesh
C           = 5  - for viscosity shading
C           = 6  - for double mesh with viscosity shading
C           = 7  - for double boundary only
C           = 12  - for se shading
C           = 13  - for double mesh with se shading
C
C     Flags:  ILAB = 0,1,2 - no labels, elements, nodes
C             IBND = 0,1 - not just boundary, external only
C                  = 2,3 - external + fault, external + interal
C             ITWO = 0,1 - single, double mesh
C             ISHD = 0,1 - no shading, shading
C
C     IOUTLINE = 0   elements shaded, no outline
C     IOUTLINE = 1   elements shaded and outlined
C     IOUTLINE = 2   elements outlined only
C
C     IPARAM = 50   elements shaded according to viscosity
C     IPARAM = 55   elements according to stress/strain exponent
C
C     Defaults:
        ILAB=0
        IBND=0
        ITWO=0
        ISHD=0
        IF(LABEL.EQ.1) ILAB=1
        IF(LABEL.EQ.2) ILAB=2
        IF(LABEL.EQ.3) IBND=1
        IF(LABEL.EQ.4) ITWO=1
        IF(LABEL.EQ.5.OR.LABEL.EQ.12) ISHD=1
        IF(LABEL.EQ.6.OR.LABEL.EQ.13) THEN
                       ISHD=1
                       ITWO=1
        END IF
        IF(LABEL.EQ.5.OR.LABEL.EQ.6) IPARAM=50
        IF(LABEL.EQ.12.OR.LABEL.EQ.13) IPARAM=55
        IF(LABEL.EQ.7) THEN
                       IBND=3
                       ITWO=1
        END IF
        IF(IFLT.EQ.2) THEN
           IF(LABEL.EQ.3) IBND=2
        END IF
        IF(IFLT.EQ.1) THEN
           IF(LABEL.EQ.3) IBND=3
           IF(LABEL.EQ.7) IBND=1
        END IF
        IOUTLINE = 2
        ICOL = 1
        EPS = 1.0E-4
C
C     If shading viscosity, find viscosity range
C
      IF (ISHD.EQ.1) THEN
        IF (IPARAM.EQ.50) THEN
          IF (VMIN.EQ.0.0.AND.VMAX.EQ.0.0) THEN
            VMIN = VHB(7,1)
            VMAX = VHB(7,1)
            DO 5 I=2,NE
              IF (VHB(7,I).LT.VMIN) VMIN=VHB(7,I)
              IF (VHB(7,I).GT.VMAX) VMAX=VHB(7,I)
   5        CONTINUE
          ENDIF
        ELSE IF (IPARAM.EQ.55) THEN
          IF (VMIN.EQ.0.0.AND.VMAX.EQ.0.0) THEN
            VMIN = VHB(8,1)
            VMAX = VHB(8,1)
            DO 6 I=2,NE
              IF (VHB(8,I).LT.VMIN) VMIN=VHB(8,I)
              IF (VHB(8,I).GT.VMAX) VMAX=VHB(8,I)
   6        CONTINUE
          ENDIF
        ENDIF
C     print *,VMIN,' ',VMAX
        RANGEV = VMAX-VMIN
        IF (RANGEV.EQ.0) RANGEV=1
        IOUTLINE = 0
      END IF
      MAXCOL = MINCOL+IRANGECOL
      RANGECOL = FLOAT(IRANGECOL)
C
C     Draw each element in sequence, and write the element
C     number inside
C
      YMID2=YMIN+YMAX
      IF(IBND.GT.0) THEN
        IF(IVERBOSE.NE.0)WRITE(6,*)'external boundary data follow:'
C
C boundary plot only
C
        KB=1
        DO 10 I=1,NE
          DO 20 K1=1,3
            K5=MOD(K1,3)+4
            K2=MOD(K1,3)+1
            IYES=0
            LM=LEM(K5,I)
C
C check if LM is on boundary
C
            DO 30 KK=KB,NBP
              LBC=IBC(KK)
              IF(LM.EQ.LBC) THEN
                IYES=1
                IF(IBND.EQ.2) THEN
C
C check if LM is on fault only if IBND=2
C
                     IF(IBCTYP(KK).EQ.11
     1                            .OR.IBCTYP(KK).EQ.12) THEN
C
C LM is on unlocked part of fault (don't plot locked portion if IBND=2
C
                       IYES=0
                   END IF
                END IF

C  only plot if IYES=1
C
                IF(IYES.EQ.1) THEN
                  X1=EX(NOR(LEM(K1,I)))
                  Y1=EY(NOR(LEM(K1,I)))
                  X2=EX(NOR(LEM(K2,I)))
                  Y2=EY(NOR(LEM(K2,I)))
                  IF ((ABS(ZOOM-1.0).LT.EPS).OR.
     :               (X1.GE.XMINC-EPS.AND.X2.GE.XMINC-EPS.AND.
     :                X1.LE.XMAXC+EPS.AND.X2.LE.XMAXC+EPS.AND.
     :                Y1.GE.YMINC-EPS.AND.Y2.GE.YMINC-EPS.AND.
     :                Y1.LE.YMAXC+EPS.AND.Y2.LE.YMAXC+EPS)) THEN
                    CALL PLOTU(X1,Y1,PENUP)
                    CALL PLOTU(X2,Y2,PENDN)
C
C    output segment end locations to stdout, projected if spherical
C
                    IF(IVERBOSE.NE.0)THEN
                      IF(NCOMP.LT.0)THEN
                        CALL PROJECTDEG(X1,Y1,XMID,YMID,1,NCOMP,IERROR)
                        CALL PROJECTDEG(X2,Y2,XMID,YMID,1,NCOMP,IERROR)
                      ENDIF
                      WRITE(6,10003)KK,X1,Y1,X2,Y2
                    ENDIF
                  END IF
                  IF(ITWO.EQ.1) THEN
C    doesn't check if segment inside cell
                    CALL PLOTU(X1+TBXOFF,Y1+TBYOFF,PENUP)
                    CALL PLOTU(X2+TBXOFF,Y2+TBYOFF,PENDN)
                  END IF
                  GO TO 20
                END IF
              END IF
   30       CONTINUE
   20     CONTINUE
   10   CONTINUE
10003   FORMAT(I6,2F13.5,3X,2F13.5)
C
C    output a list of boundary node numbers
C
C       IF(IVERBOSE.NE.0)THEN
C         WRITE(6,*)'boundary vertex node numbers and coordinates',
C    1              ' follow'
C         DO 15 J=1,NBP
C           NODNO=NOR(IBC(J))
C           IF(NODNO.LE.NN)WRITE(6,*)NODNO,EX(NODNO),EY(NODNO)
C  15     CONTINUE
C       END IF
      ELSE
C
C If IBND = 0, this is not a boundary plot only
C
      IF(IVERBOSE.NE.0)WRITE(6,*)'triangle node allocation follows'
      DO 100 I=1,NE
        KS=3
C
C     Define vertex to commence drawing
C
        M=NOR(LEM(KS,I))
        X=EX(M)
        Y=EY(M)
        IF ((ABS(ZOOM-1.0).LT.EPS).OR.
     :         (X.GE.XMINC-EPS.AND.X.GE.XMINC-EPS.AND.
     :          X.LE.XMAXC+EPS.AND.X.LE.XMAXC+EPS.AND.
     :          Y.GE.YMINC-EPS.AND.Y.GE.YMINC-EPS.AND.
     :          Y.LE.YMAXC+EPS.AND.Y.LE.YMAXC+EPS)) THEN
          CALL PLOTU(X,Y,PENUP)
          XM=0.0
          YM=0.0
C
C     Draw the 3 sides in sequence
C
          DO 50 J=1,3
            KS=MOD(KS,3)+1
            LM=LEM(KS,I)
            M=NOR(LM)
            MENV(J)=M
            X=EX(M)
            Y=EY(M)
            XM = XM+X
            YM = YM+Y
C
C    Draw the line
C
C           IF(ISHD.EQ.1) THEN
              POLYX(J)=X
              POLYY(J)=Y
              IF(J.EQ.3) THEN
                  IF (IPARAM.EQ.50)
     :              ICOL = MINCOL + (VHB(7,I)-VMIN) * RANGECOL/RANGEV
                  IF (IPARAM.EQ.55)
     :              ICOL = MINCOL + (VHB(8,I)-VMIN) * RANGECOL/RANGEV
                  IF (ICOL.LT.MINCOL) ICOL=MINCOL
                  IF (ICOL.GT.MAXCOL) ICOL=MAXCOL
C    fill the element with icol and outline it
                  CALL FILLPOLY(POLYX,POLYY,J,ICOL,IOUTLINE)
              END IF
C           ELSE
C             CALL PLOTU(X,Y,PENDN)
C           END IF
   50     CONTINUE
C
C    output mesh description to stdout (triangle configuration
C    here; see below for node coordinates)
C
            IF(IVERBOSE.NE.0)WRITE(6,10001)I,(MENV(J),J=1,3)
10001       FORMAT(4I5)
C
C Draw a 2nd mesh which is offset by tbxoff,tbyoff
C
          KS=3
          IF(ITWO.EQ.1) THEN
            M=NOR(LEM(KS,I))
            X=EX(M)+TBXOFF
            Y=EY(M)+TBYOFF
            CALL PLOTU(X,Y,PENUP)
            DO 150 J=1,3
             KS=MOD(KS,3)+1
             LM=LEM(KS,I)
             M=NOR(LM)
             X=EX(M)+TBXOFF
             Y=EY(M)+TBYOFF
C
C    Draw the line
C
C            IF(ISHD.EQ.1) THEN
               POLYX(J)=X
               POLYY(J)=Y
               IF(J.EQ.3) THEN
                  IF (IPARAM.EQ.50)
     :              ICOL = MINCOL + (VHB(7,I)-VMIN) * RANGECOL/RANGEV
                  IF (IPARAM.EQ.55)
     :              ICOL = MINCOL + (VHB(8,I)-VMIN) * RANGECOL/RANGEV
                 CALL FILLPOLY(POLYX,POLYY,J,ICOL,IOUTLINE)
               END IF
C            ELSE
C              CALL PLOTU(X,Y,PENUP)
C            END IF
  150      CONTINUE
         END IF
C
         IF(ILAB.EQ.1.AND.((MOD(I,MPE).EQ.1).OR.(ZOOM.GT.1)))THEN
C
C     Write element number
C
C           XW=0.5*(XL+XM)
C           YW=0.5*(YL+YM)
            XW = XM/3.0
            YW = YM/3.0
            RI=FLOAT(I)
C           IF (I.EQ.1) CALL SETFONT('Helvetica',9,8)
            CALL FORMATNUMBER(NUMBER,RI,50)
            LEN = 0
 200        LEN = LEN+1
            N=RI/(10**LEN)
            IF (N.GT.0) GO TO 200
            CALL DRAWLABEL(XW,YW,0,-1,1,NUMBER,LEN)
          END IF
        END IF
  100 CONTINUE
C
C    output node coordinates to stdout
C
C       IF(IVERBOSE.NE.0)THEN
C         WRITE(6,*)'node locations follow'
C         DO 300 J=1,NN
C           WRITE(6,10002)J,EX(J),EY(J)
C 300     CONTINUE
C       END IF
C10002       FORMAT(I5,2F9.5)
      END IF
      RETURN
      END
      SUBROUTINE INTBND(IMP,LEM,NOR,IMAT,IBC,IBNGH,NMP,NE,NN,NUP,NBP,
     :                  EX,EY,BNDS,DENS,DENSREG,ISVP,IREG,NCOMP,
     :                  XMID,YMID,IVERBOSE)
C
C    This routine locates internal boundaries, as defined by two adjoining
C    elements having different IMAT numbers.  If called with ISVP = 0, the
C    boundary is plotted but not saved.  If ISVP>1 boundary data are saved
C    for use in gravity calculation, but not plotted here. On exit ISVP is
C    set to the number of boundary segments
C
      INCLUDE "sybilps.parameters"
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION LEM(6,NE),NOR(NUP)
      DIMENSION IMP(2,NMP),IMAT(NE)
      DIMENSION IBC(2*NBP),IBNGH(2*NBP)
      DIMENSION BNDS(4,ISVP),DENS(7,NE)
C
C     WRITE(6,*)'In INTBND, ISVP =',ISVP
      DO K=1,NMP
        IMP(1,K)=0
        IMP(2,K)=0
      ENDDO
C
C   First we find the two element numbers that are found on either side
C   of every midpoint node
C
      DO J=1,NE
        DO K=4,6
          NMN=NOR(LEM(K,J))-NN
          IF(IMP(1,NMN).EQ.0)THEN
            IMP(1,NMN)=J
          ELSE IF(IMP(2,NMN).EQ.0)THEN
            IMP(2,NMN)=J
          ELSE
            WRITE(*,*)'problem 1 in INTBND; J, K, NMN =',J,K,NMN
            STOP
          END IF
        ENDDO
      ENDDO
C
C   Next we compare the IMAT numbers of every midpoint node, looking for
C    inequality to define internal boundary segment
C
      IF(IVERBOSE.NE.0)
     1    WRITE(6,*)'Boundary segments for region ',IREG
      KSV=0
      DENSREG=0.0
      DO 300 J=1,NMP
        NEL1=IMP(1,J)
        NEL2=IMP(2,J)
        IMNL1=IMAT(NEL1)
        IMNL2=IMAT(NEL2)
C
C   check first for nodes on the external boundary, find boundary segment for 
C    element with IMAT=IREG.  There may be 1 or 2 external segments
C
        IF((IMNL1.EQ.IREG).OR.(IMNL2.EQ.IREG))THEN
          IF(NEL2.EQ.0)THEN
            IF(ISVP.GT.0)DENSREG=DENS(1,NEL1)
            DO M=4,6
              LEMNL=LEM(M,NEL1)
              NMN=NOR(LEMNL)-NN
            DO K=1,NBP
              IF((IBC(K).EQ.LEMNL).AND.(J.EQ.NMN))THEN
                LADD=NOR(IBNGH(K))
                X1=EX(LADD)
                Y1=EY(LADD)
                LADD=NOR(IBNGH(K+NBP))
                X2=EX(LADD)
                Y2=EY(LADD)
                IF(ISVP.EQ.0)THEN
                  CALL PLOTU(X1,Y1,PENUP)
                  CALL PLOTU(X2,Y2,PENDN)
C
C    output segment end locations to stdout, projected if spherical
C
                  IF(IVERBOSE.NE.0)THEN
                    IF(NCOMP.LT.0)THEN
                      CALL PROJECTDEG(X1,Y1,XMID,YMID,1,NCOMP,IERROR)
                      CALL PROJECTDEG(X2,Y2,XMID,YMID,1,NCOMP,IERROR)
                    ENDIF
                    WRITE(6,10003)K,X1,Y1,X2,Y2
                  ENDIF
C                 IF(IVERBOSE.NE.0)WRITE(6,1001)X1,Y1,X2,Y2
C
C    save these coordinates for other use
C
                ELSE
                  KSV=KSV+1
                  IF(KSV.GT.ISVP)GO TO 555
                  BNDS(1,KSV)=X1
                  BNDS(2,KSV)=Y1
                  BNDS(3,KSV)=X2
                  BNDS(4,KSV)=Y2
                END IF
              END IF
            ENDDO           ! on K
            ENDDO           ! on M
C
C  now check for internal boundary for a region with IMAT = IREG
C
          ELSE 
            IF(IMNL1.NE.IMNL2)THEN
              IF((DENSREG.EQ.0.0).AND.(ISVP.GT.0))DENSREG=DENS(1,NEL)
              DO K1=4,6
                DO K2=4,6
                  IF(LEM(K1,NEL1).EQ.LEM(K2,NEL2))THEN
                    IF(IMNL1.EQ.IREG)THEN
                      KMID=K1
                      NEL=NEL1
                    ELSE
                      KMID=K2
                      NEL=NEL2
                    END IF
                    GO TO 450
                  END IF
                ENDDO       ! on K2
              ENDDO         ! on K1
              WRITE(*,*)'problem 3 in INTBND'
              STOP
  450         CONTINUE
              IF(ISVP.GT.0)DENSREG=DENS(1,NEL)
              KV1=MOD(KMID+1,3)+1
              KV2=MOD(KV1,3)+1
              LADD=NOR(LEM(KV1,NEL))
              X1=EX(LADD)
              Y1=EY(LADD)
              LADD=NOR(LEM(KV2,NEL))
              X2=EX(LADD)
              Y2=EY(LADD)
              KSV=KSV+1
              IF(ISVP.EQ.0)THEN
                CALL PLOTU(X1,Y1,PENUP)
                CALL PLOTU(X2,Y2,PENDN)
C
C    output segment end locations to stdout, projected if spherical
C
                IF(IVERBOSE.NE.0)THEN
                  IF(NCOMP.LT.0)THEN
                    CALL PROJECTDEG(X1,Y1,XMID,YMID,1,NCOMP,IERROR)
                    CALL PROJECTDEG(X2,Y2,XMID,YMID,1,NCOMP,IERROR)
                  ENDIF
                  WRITE(6,10003)NEL,X1,Y1,X2,Y2
                  ENDIF
C             IF(IVERBOSE.NE.0)WRITE(6,1001)X1,Y1,X2,Y2
              ELSE
                IF(KSV.GT.ISVP)GO TO 555
                BNDS(1,KSV)=X1
                BNDS(2,KSV)=Y1
                BNDS(3,KSV)=X2
                BNDS(4,KSV)=Y2
              END IF
            END IF
          END IF
        END IF
  300 CONTINUE
 1001 FORMAT(4F12.6)
C
C     IF(ISVP.GT.0)ISVP=KSV
      ISVP=KSV
      RETURN
  555 WRITE(*,*)'Insufficient space provided in BNDS in INTBND'
10003 FORMAT(I6,2F13.5,3X,2F13.5)
      RETURN
      END

      SUBROUTINE PLSEG(LABEL,NUP,NBP,NSEG,EX,EY,NOR,IBC,
     1                 IBCTYP,ISEG,
     2                 IRANGECOL,XMINC,XMAXC,YMINC,YMAXC,
     3                 NCOMP,XMID,YMID,IVERBOSE)
C
C    This routine is only used if the ISEG array is found in the solution
C    This routine draws the mesh segments which are not marked
C    as external boundary segments
C    i.e. the internal boundaries and faults are drawn.
C
      INCLUDE "sybilps.parameters"
      CHARACTER NUMBER*80
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION NOR(NUP),IBC(NBP),IBCTYP(NBP*2)
      DIMENSION ISEG(3,NSEG)
C
      IERROR=0
      IF(IVERBOSE.NE.0)WRITE(6,*)'internal segments follow:'
      DO 10 I=1,NSEG
        IF (ISEG(3,I).GT.1000 .AND. ISEG(3,I).LE.3000) THEN
          M=NOR(ISEG(1,I))
          X1=EX(M)
          Y1=EY(M)
          IF(X1.GE.XMINC.AND.X1.LE.XMAXC.AND.
     :          Y1.GE.YMINC.AND.Y1.LE.YMAXC) THEN
            M=NOR(ISEG(2,I))
            X2=EX(M)
            Y2=EY(M)
            CALL PLOTU(X1,Y1,PENUP)
            CALL PLOTU(X2,Y2,PENDN)
C
C    output segment end locations to stdout, projected if spherical
C
            IF(IVERBOSE.NE.0)THEN
              IF(NCOMP.LT.0)THEN
                CALL PROJECTDEG(X1,Y1,XMID,YMID,1,NCOMP,IERROR)
                CALL PROJECTDEG(X2,Y2,XMID,YMID,1,NCOMP,IERROR)
              ENDIF
              WRITE(6,10003)I,X1,Y1,X2,Y2
            ENDIF
          ENDIF
        ENDIF
 10   CONTINUE
10003 FORMAT(I6,2F13.5,3X,2F13.5)
C
C    output segment array to stdout
C
C       IF(IVERBOSE.NE.0)THEN
C         WRITE(6,*)'segment data follows'
C         DO 300 J=1,NSEG
C           WRITE(6,10002)J,(ISEG(I,J),I=1,3)
C 300     CONTINUE
C       END IF
C10002   FORMAT(4I6)
      RETURN
      END
      SUBROUTINE STMARK(STELPX,STELPY,NSM,NPM,NRM,R0,IVERBOSE,
     :                   XMINC,XMAXC,YMINC,YMAXC)
      DIMENSION STELPX(NPM,NSM),STELPY(NPM,NSM)

      INCLUDE "sybilps.parameters"
      IF (IVERBOSE.NE.0) WRITE(6,10001)

      DO 10 J=1,NSM
         XX=STELPX(1,J)
         YY=STELPY(1,J)
         CALL PLOTU(XX,YY,PENUP)
         DO 20 I=2,NPM
            XX=STELPX(I,J)
            YY=STELPY(I,J)
            IPENSTATE=PENUP
            IF(XX.GE.XMINC.AND.XX.LE.XMAXC.AND.
     :          YY.GE.YMINC.AND.YY.LE.YMAXC) IPENSTATE=PENDN
            CALL PLOTU(XX,YY,IPENSTATE)
   20    CONTINUE
         IF (IVERBOSE.NE.0) THEN
             CALL ELLPAR(STELPX(1,J),STELPY(1,J),NPM-NRM-1,R0,J)
             ANG=0
             CALL ANGLECALC(STELPX(NPM,J),STELPY(NPM,J),
     :                      STELPX(NPM,J)+1.0,STELPY(NPM,J),
     :                      STELPX(NPM-NRM,J),STELPY(NPM-NRM,J),ANG)
             print *,STELPX(NPM,J),STELPY(NPM,J),
     :               STELPX(NPM-NRM,J),STELPY(NPM-NRM,J),ANG
         ENDIF
   10 CONTINUE
10001 FORMAT(' No. Qual  Centre (xc,yc)      axis a    axis b    axis c',
     1'  orientation')
      RETURN
      END

      SUBROUTINE ELLPAR(X,Y,N,RZERO,NUM)
C
C    this subroutine evaluates the parameters for a best-fit
C    ellipse to match a set of points X(N),Y(N)
C    The parameters are (XC,YC) - centre coordinates
C                       (AA,BB) - semi major and minor axes
C                        ALPHA - orientation of semi-major axis
C
      DIMENSION X(N),Y(N)
      PI=3.14159265359
      IQUAL=0
C
C    1.  compute central coordinates, using spatial average
C
      XC=0.0
      YC=0.0
      DO 100 J=1,N
      XC=XC+X(J)
      YC=YC+Y(J)
  100 CONTINUE
      XCC=XC/FLOAT(N)
      YCC=YC/FLOAT(N)
C
C    1a  modified algorithm: locate 2 opposite points that are closest
C        then use midpoint.
C
      DIAM=2.0*RZERO
      NHALF=N/2
      JMIN=1
      DO 101 J=1,NHALF
      XDIF=X(J+NHALF)-X(J)
      YDIF=Y(J+NHALF)-Y(J)
      DTEST=SQRT(XDIF*XDIF+YDIF*YDIF)
      IF(DTEST.LT.DIAM)THEN
        DIAM=DTEST
        JMIN=J
      END IF
  101 CONTINUE
      XC=0.5*(X(JMIN+NHALF)+X(JMIN))
      YC=0.5*(Y(JMIN+NHALF)+Y(JMIN))
C     WRITE(6,*)'XCC,XC, YCC,YC',XCC,XC,YCC,YC
C
C      initialise variables for integration
C
      SUMSN=0.0
      SUMCS=0.0
      SUMTH=0.0
      RJM1SQ=(X(N)-XC)*(X(N)-XC)+(Y(N)-YC)*(Y(N)-YC)
      RJM1=SQRT(RJM1SQ)
      DXJM1Q=(X(1)-X(N))*(X(1)-X(N))+(Y(1)-Y(N))*(Y(1)-Y(N))
C
C   2.  integrate around the set of points to get alpha
C
      DO 200 J=1,N
      RJSQ=(X(J)-XC)*(X(J)-XC)+(Y(J)-YC)*(Y(J)-YC)
      RJ=SQRT(RJSQ)
      SN2TH=2.0*(Y(J)-YC)*(X(J)-XC)/RJSQ
      CS2TH=((X(J)-XC)*(X(J)-XC)-(Y(J)-YC)*(Y(J)-YC))/RJSQ
      IF(J.EQ.N)THEN
      RJP1SQ=(X(1)-XC)*(X(1)-XC)+(Y(1)-YC)*(Y(1)-YC)
      DXJP1Q=(X(1)-X(N))*(X(1)-X(N))+(Y(1)-Y(N))*(Y(1)-Y(N))
      ELSE
      RJP1SQ=(X(J+1)-XC)*(X(J+1)-XC)+(Y(J+1)-YC)*(Y(J+1)-YC)
      DXJP1Q=(X(J+1)-X(J))*(X(J+1)-X(J))+
     1                           (Y(J+1)-Y(J))*(Y(J+1)-Y(J))
      END IF
      RJP1=SQRT(RJP1SQ)
      DTHJP1=(SQRT(DXJP1Q-(RJP1-RJ)*(RJP1-RJ)))/(RJP1+RJ)
      DTHJM1=(SQRT(DXJM1Q-(RJ-RJM1)*(RJ-RJM1)))/(RJM1+RJ)
      DTH=DTHJP1+DTHJM1
C     WRITE(6,*)J, RJ,DTH,DTHJM1,DTHJP1
C
C    sums for orientation
C
      SUMSN=SUMSN+SN2TH*DTH/RJSQ
      SUMCS=SUMCS+CS2TH*DTH/RJSQ
      SUMTH=SUMTH+DTH
C
C    set up RJM1, DXJM1Q for next pass of this loop
C
      RJM1=RJ
      DXJM1Q=DXJP1Q
  200 CONTINUE
C
C   2 possible values of alpha; ambiguity resolved by values of a,b
C
      ALPHA=0.5*ATAN2(SUMSN,SUMCS)
      CSALP=COS(ALPHA)
      SNALP=SIN(ALPHA)
C     WRITE(6,*)'CSALP,SNALP,SUMTH',CSALP,SNALP,SUMTH
C
C    2nd integration to get best-fit estimates of a and b, assumes alpha
C
      SUMA11=0.0
      SUMA12=0.0
      SUMA22=0.0
      SUMB1=0.0
      SUMB2=0.0
      RJM1SQ=(X(N)-XC)*(X(N)-XC)+(Y(N)-YC)*(Y(N)-YC)
      RJM1=SQRT(RJM1SQ)
      DXJM1Q=(X(1)-X(N))*(X(1)-X(N))+(Y(1)-Y(N))*(Y(1)-Y(N))
      DO 300 J=1,N
      RJSQ=(X(J)-XC)*(X(J)-XC)+(Y(J)-YC)*(Y(J)-YC)
      RJ=SQRT(RJSQ)
      CSTMA=((X(J)-XC)*CSALP+(Y(J)-YC)*SNALP)/RJ
      SNTMA=((Y(J)-YC)*CSALP-(X(J)-XC)*SNALP)/RJ
      SNQTMA=SNTMA*SNTMA
      CSQTMA=CSTMA*CSTMA
      SN2QTM=SNQTMA*CSQTMA
      SNQQTM=SNQTMA*SNQTMA
      CSQQTM=CSQTMA*CSQTMA
      IF(J.EQ.N)THEN
      RJP1SQ=(X(1)-XC)*(X(1)-XC)+(Y(1)-YC)*(Y(1)-YC)
      DXJP1Q=(X(1)-X(N))*(X(1)-X(N))+(Y(1)-Y(N))*(Y(1)-Y(N))
      ELSE
      RJP1SQ=(X(J+1)-XC)*(X(J+1)-XC)+(Y(J+1)-YC)*(Y(J+1)-YC)
      DXJP1Q=(X(J+1)-X(J))*(X(J+1)-X(J))+
     1                           (Y(J+1)-Y(J))*(Y(J+1)-Y(J))
      END IF
      RJP1=SQRT(RJP1SQ)
C    no check for SQRT arg being negative
      DTHJP1=(SQRT(DXJP1Q-(RJP1-RJ)*(RJP1-RJ)))/(RJP1+RJ)
      DTHJM1=(SQRT(DXJM1Q-(RJ-RJM1)*(RJ-RJM1)))/(RJM1+RJ)
      DTH=DTHJP1+DTHJM1
C
C    sums for matrix equation
C
      SUMA11=SUMA11+CSQQTM*DTH
      SUMA12=SUMA12+SN2QTM*DTH
      SUMA22=SUMA22+SNQQTM*DTH
      SUMB1=SUMB1+CSQTMA*DTH/RJSQ
      SUMB2=SUMB2+SNQTMA*DTH/RJSQ
C
C    set up RJM1, DXJM1Q for next pass of this loop
C
      RJM1=RJ
      DXJM1Q=DXJP1Q
  300 CONTINUE
C
C   calculate determinant and solve matrix
C
      DET=SUMA12*SUMA12-SUMA11*SUMA22
      IF(ABS(DET).LT.1.E-10)THEN
        A2INV=SUMB1/(SUMA11+SUMA12)
        B2INV=A2INV
      ELSE
        A2INV=(SUMA12*SUMB2-SUMA22*SUMB1)/DET
        B2INV=(SUMA12*SUMB1-SUMA11*SUMB2)/DET
      END IF
      IF (B2INV.LE.0.OR.A2INV.LE.0) IQUAL = 1
      AA=1.0/SQRT(ABS(A2INV))
      BB=1.0/SQRT(ABS(B2INV))
C
C    resolve ambiguity with alpha
C
      IF(BB.GT.AA)THEN
      ALPHA=ALPHA+0.5*PI
      IF(ALPHA.GT.PI)ALPHA=ALPHA-PI
      CC=AA
      AA=BB
      BB=CC
      END IF
      ALPDEG=ALPHA*180.0/PI
C
C    normalise to initial size of circle
C
      AA=AA/RZERO
      BB=BB/RZERO
      CC=1.0/(AA*BB)
C
C   write best-fit parameters
C
      WRITE(6,10002)NUM,IQUAL,XC,YC,AA,BB,CC,ALPDEG
10002 FORMAT(I3,I5,5F10.5,F11.5)
      RETURN
      END

      SUBROUTINE ANGLECALC(X,Y,X1,Y1,X2,Y2,ANG)

      PI=3.14159265359
      ANG = 0
      B1 = Y1-Y
      A1 = X1-X
      A2 = Y2-Y
      B2 = X2-X
      ANGPHI = ATAN2(A2,B2)
      SINPHI = SIN(ANGPHI)
      COSPHI = COS(ANGPHI)
      YVAL = B1*COSPHI - A1*SINPHI
      XVAL = A1*COSPHI + B1*SINPHI

      IF (XVAL.EQ.0.AND.YVAL.EQ.0) THEN
          print *,"Both zero"
      ELSE
          ANG = ATAN2(YVAL,XVAL)
      ENDIF

      RETURN
      END
