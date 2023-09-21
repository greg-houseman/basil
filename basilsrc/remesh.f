C*--------------------------------------------------------------------
C*    Basil / Sybil:   remesh.f  1.0  1 October 2009
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE REMESH(EX,EY,IBCTYP,IBC,IBNGH,SEGSQLIM,ANGLIM,
     :                  PTLOC,ISEG,LEM,NOR,INTV,IVNEW,IFLT,NFP,
     :                  NFPF3,IFBC1,IFEQV,PDIST,MSNODE,MEASUR,
     :                  LUW,LSC,IERR)
C
      INCLUDE "indices.parameters"
      INCLUDE "limits.parameters"
      INTEGER INTV(64),IVNEW(64)
      DIMENSION MEASUR(MAXMEAS), MSNODE(MAXMEAS), PTLOC(3,MAXMEAS)
      REAL, DIMENSION(:), POINTER :: EX, EY
      INTEGER, DIMENSION(:,:), POINTER :: ISEG,LEM
      INTEGER, DIMENSION(:), POINTER :: NOR,IBC,IFBC1,IFEQV
      INTEGER, DIMENSION(:,:), POINTER :: IBCTYP,IBNGH
      INTEGER, ALLOCATABLE :: IPROTECTED(:)
      LOGICAL, ALLOCATABLE :: FOUND(:)
      INTERFACE
        SUBROUTINE FINDPROTECTED(NUMPROTECTED,FOUND,
     :                        IBCTYP,IBC,IBNGH,ISEG,
     :                        NOR,NBP,NN,NUP,NSEG)
        DIMENSION IBC(NBP)
        DIMENSION IBCTYP(NBP*2)
        DIMENSION IBNGH(NBP*2)
        DIMENSION NOR(NUP), ISEG(3,NSEG)
        LOGICAL FOUND(NBP)
        END SUBROUTINE
        SUBROUTINE SETPROTECTED(NUMPROTECTED,IPROTECTED,FOUND,
     :                        IBC,NOR,NBP,NUP,IERR)
        DIMENSION IBC(NBP)
        DIMENSION NOR(NUP)
        INTEGER IPROTECTED(NUMPROTECTED)
        LOGICAL FOUND(NBP)
        END SUBROUTINE
        SUBROUTINE BEMESH(EX,EY,SEGSQLIM,ANGMIN,
     :                  PTLOC,ISEG,IBC,LEM,NOR,IFBC1,IV,IVNEW,
     :                  IPROTECTED,NUMPROTECTED,NUP,NSEG,NBP,
     :                  IFLT,NFPF3,MSNODE,MEASUR,LUW,LSC,IERR)
        INCLUDE "indices.parameters"
        INCLUDE "limits.parameters"
C       REAL, DIMENSION(:), POINTER :: EX, EY
C       INTEGER, DIMENSION(:,:), POINTER :: ISEG
        INTEGER, DIMENSION(:,:), POINTER :: LEM
        INTEGER, DIMENSION(:), POINTER :: NOR
        DIMENSION IPROTECTED(NUMPROTECTED)
        DIMENSION ISEG(3,NSEG)
        DIMENSION EX(NUP),EY(NUP)
        DIMENSION IBC(NBP)
        DIMENSION IFBC1(NFPF3)
        DIMENSION PTLOC(3,MAXMEAS)
        DIMENSION MEASUR(MAXMEAS), MSNODE(MAXMEAS)
        DIMENSION IV(64), IVNEW(64)
        END SUBROUTINE
      END INTERFACE
C
C    this is an umbrella routine that basically calls FINDPROTECTED,
C    SETPROTECTED, and BEMESH (BEMESH was previously called REMESH)
C
      NUP=INTV(INUP)
      NSEG=INTV(INSEG)
      NBP=INTV(INBP)
      NN=INTV(INN)
C
C   if periodic boundaries apply, confirm current value of PDIST
C
      KPBC=INTV(IKPBC)
      IF((KPBC.GT.0).AND.(NFP.GT.NFPF3))
     :  CALL PBWIDTH(EX,EY,IFBC1,IFEQV,NOR,NUP,NFP,NFPF3,KPBC,PDIST)
      NPBPT=NFP-NFPF3
C
C    identify nodes that should be protected in the REMESH procedure
C
      ALLOCATE(FOUND(NBP),STAT=IERR)
      CALL FINDPROTECTED(NUMPROTECTED,FOUND,
     :                   IBCTYP,IBC,IBNGH,ISEG,NOR,
     :                   NBP,NN,NUP,NSEG)
      IF (NUMPROTECTED.GT.0)
     :  ALLOCATE(IPROTECTED(NUMPROTECTED),STAT=IERR)
      IF(IERR.NE.0)THEN
        WRITE(LSC,10200)
        WRITE(LUW,10200)
        RETURN
      ENDIF
      CALL SETPROTECTED(NUMPROTECTED,IPROTECTED,FOUND,
     :                  IBC,NOR,NBP,NUP,IERR)
      IF(IERR.NE.0)THEN
        WRITE(LSC,10201)
        WRITE(LUW,10201)
        RETURN
      END IF
C
      CALL BEMESH(EX,EY,SEGSQLIM,ANGLIM,
     :            PTLOC,ISEG,IBC,LEM,NOR,IFBC1,INTV,IVNEW,
     :            IPROTECTED,NUMPROTECTED,NUP,NSEG,NBP,
     :            IFLT,NFPF3,MSNODE,MEASUR,LUW,LSC,IERR)
C
      IF (IERR.NE.0)THEN
        WRITE(LSC,10202)
        WRITE(LUW,10202)
        RETURN
      END IF
C
      DEALLOCATE(FOUND,STAT=IERR)
      IF (NUMPROTECTED.GT.0) DEALLOCATE(IPROTECTED,STAT=IERR)
C
10200 FORMAT("Memory allocation problem for protected nodes.")
10201 FORMAT("Problem finding protected nodesin SETPROTECT.")
10202 FORMAT("Problem in remeshing operation (BEMESH).")
      RETURN
      END
C
      SUBROUTINE BEMESH(EX,EY,SEGSQLIM,ANGMIN,
     :                  PTLOC,ISEG,IBC,LEM,NOR,IFBC1,IV,IVNEW,
     :                  IPROTECTED,NUMPROTECTED,NUP,NSEG,NBP,
     :                  IFLT,NFPF3,MSNODE,MEASUR,LUW,LSC,IERR)
      INCLUDE "indices.parameters"
      INCLUDE "limits.parameters"
      INCLUDE "input.parameters"
C     REAL, DIMENSION(:), POINTER :: EX, EY
C     INTEGER, DIMENSION(:,:), POINTER :: ISEG
      INTEGER, DIMENSION(:,:), POINTER :: LEM
      INTEGER, DIMENSION(:), POINTER :: NOR
      DIMENSION IPROTECTED(NUMPROTECTED)
      DIMENSION ISEG(3,NSEG)
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION IBC(NBP)
      DIMENSION IFBC1(NFPF3)
      DIMENSION PTLOC(3,MAXMEAS)
      DIMENSION MEASUR(MAXMEAS), MSNODE(MAXMEAS)
      DIMENSION IV(64), IVNEW(64)
      INTERFACE
        SUBROUTINE BNDSEGS(EX,EY,ISEG,IFBC1,IPROTECTED,ICNTPRO,
     :                     NUP,NSEG,IFLT,NFPF3,IBPSIZE,IERR)
C
C      Remove small segments around the external boundary
C      Assumes each segment in ISEG is vertex-vertex?
C
        INCLUDE "limits.parameters"
        DIMENSION EX(NUP),EY(NUP)
        DIMENSION ISEG(3,NSEG)
        DIMENSION IFBC1(NFPF3)
        INTEGER IPROTECTED(ICNTPRO)
        END SUBROUTINE
      END INTERFACE
C      
C    Convert Element/local ccordinates into (X,Y) so that
C    reconversion can be done after remeshing
C
      MSINDX=IV(IMSINDX)
      NSEG=IV(INSEG)
      CALL RESETSERIES(LEM,NOR,IV(INE),NUP,MEASUR,MSNODE,
     :                 MSINDX,PTLOC,EX,EY)
C
C  Find nodes that appear only once in the segment list to set IBPCNT
C  used for array allocation in BNDSEGS; IBPCNT is obtained exactly
C  in BNDSEGS, but allow for extra nodes for terminating boundaries.
C  multi-segment nodes will cause IBPCNT to be overestimated here
C
      ILONE=0
      DO I1=1,NSEG
        DO K=1,2
          N1=ISEG(K,I1)
          INUM=0
          DO I2=1,NSEG
            IF(I1.NE.I2)THEN
              IF((N1.EQ.ISEG(1,I2)).OR.
     :           (N1.EQ.ISEG(2,I2)))INUM=INUM+1
            ENDIF
          ENDDO
        ENDDO
        IF(INUM.EQ.0)ILONE=ILONE+1
      ENDDO
      IBPCNT=NSEG+ILONE/2
C     WRITE(*,*)'BEMESH: NSEG,IBPCNT,ILONE,NUMPROTECTED=',
C    :           NSEG,IBPCNT,IMULT,NUMPROTECTED
C
C    find boundary segments and nodes to be passed to remeshing routine
C    delete nodes if necessary, and remesh
C
      CALL BNDSEGS(EX,EY,ISEG,IFBC1,IPROTECTED,NUMPROTECTED,
     :             NUP,NSEG,IFLT,NFPF3,IBPCNT,IERR)
      DO I=1,64
        IVNEW(I) = IV(I)
      END DO
   5  RETURN
      END SUBROUTINE BEMESH
C      
      SUBROUTINE BNDSEGS(EX,EY,ISEG,IFBC1,IPROTECTED,ICNTPRO,
     :                   NUP,NSEG,IFLT,NFPF3,IBPSIZE,IERR)
C
C   This routine examines the list of segments used to define the 
C   current mesh build, and amends it to remove excessively short 
C   segments. The modified segment list, ISEGCP(3,NSEGNEW) and corresponding
C   set of vertex nodes, IBP(NEWNBP) is passed to remesh_data to build
C   a new mesh. Each segment in ISEG is vertex-vertex, but sequence of
C   segments is randomised by the build process originally used in
C   triangle, and order of two nodes in ISEG is not systematic.
C
      INCLUDE "limits.parameters"
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION ISEG(3,NSEG)
      DIMENSION IPROTECTED(ICNTPRO)
      DIMENSION IFBC1(NFPF3)
C
C    local arrays
C
      DIMENSION ISEGCP(3,NSEG)
      DIMENSION IBP(IBPSIZE),NBUSE(5,IBPSIZE)
      DIMENSION SEGPRO(2,IBPSIZE)
C
C  IBP is a list of nodes on boundary segments (external + internal)
C  IBPCNT is the number of unique vertex points in these segments
C  for a simple closed boundary, IBPCNT = NSEG. More than 2 segs
C  connected to a node may result in IBPCNT < IBPSIZE - no problem
C  First find the set of unique nodes in ISEG and store them in IBP
C
      IBP(1)=ISEG(1,1)
      IBP(2)=ISEG(2,1)
      IBPCNT=2
      DO I=2,NSEG
        DO J=1,IBPCNT
          IF(IBP(J).EQ.ISEG(1,I))GOTO 10
        ENDDO
        IBPCNT=IBPCNT+1
        IBP(IBPCNT)=ISEG(1,I)
   10   CONTINUE
        DO J=1,IBPCNT
          IF(IBP(J).EQ.ISEG(2,I))GOTO 11
        ENDDO
        IBPCNT=IBPCNT+1
        IBP(IBPCNT)=ISEG(2,I)
   11   CONTINUE
      ENDDO
C
C   IBP can contain a mix of internal and external boundary points
C
      IF(IBPCNT.GT.IBPSIZE)THEN
        WRITE(*,*)'WARNING: IBPCNT > IBPSIZE ',IBPCNT,IBPSIZE,
     :            'in BNDSEGS, check setting'
        ERR=1
        STOP
      ENDIF
      WRITE(*,10101)NSEG,IBPCNT,ICNTPRO
10101 FORMAT('Mesh rebuild uses ',I5,' segments and ',I5' nodes with',
     :       I4,' protected')
C
C    copy the segments into ISEGCP, changing indices for ISEGCP
C
      DO I=1,NSEG
        ISEGCP(3,I)=ISEG(3,I)
      END DO
      DO I=1,NSEG
        DO J=1,IBPCNT
          IF (ISEG(1,I).EQ.IBP(J)) ISEGCP(1,I)=J
          IF (ISEG(2,I).EQ.IBP(J)) ISEGCP(2,I)=J
        END DO
      END DO
C
C    what follows is preparation for eliminating points too close
C    first flag nodes that are endpoints or junctions, eligible 
C    nodes for deletion are found in 2 segments, no more, no less.
C
      DO J=1,IBPCNT
        DO K=1,5
          NBUSE(K,J)=0
        ENDDO
      ENDDO
      DO I=1,NSEG
        NBUSE(1,ISEGCP(1,I))=NBUSE(1,ISEGCP(1,I))+1
        NBUSE(1,ISEGCP(2,I))=NBUSE(1,ISEGCP(2,I))+1
      ENDDO
C
C    Store neighbour node addresses in NBUSE(2,J) and NBUSE(3,J)
C    allowing that segments are in random order and either direction
C    and store the SEG addresses in NBUSE(4,J) and NBUSE(5,J)
C    beware that protected nodes NBUSE values are not set - because
C    this method is unreliable if more than two segs connected to a node
C
      DO I=1,NSEG
        DO J=1,IBPCNT
          IF(NBUSE(1,J).EQ.2)THEN
            IF(ISEGCP(1,I).EQ.J)THEN
              IF(NBUSE(2,J).EQ.0)THEN
                NBUSE(2,J)=ISEGCP(2,I)
                NBUSE(4,J)=I
              ELSE
                NBUSE(3,J)=ISEGCP(2,I)
                NBUSE(5,J)=I
              ENDIF
            ELSE IF(ISEGCP(2,I).EQ.J)THEN
              IF(NBUSE(2,J).EQ.0)THEN
                NBUSE(2,J)=ISEGCP(1,I)
                NBUSE(4,J)=I
              ELSE
                NBUSE(3,J)=ISEGCP(1,I)
                NBUSE(5,J)=I
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
C    at this point NBUSE(1,J) contains 2 for nodes eligible for deletion
C    deletion of any nodes may be avoided by setting NBUSE(1,J)=0
C    First make nodes listed in IPROTECTED ineligible for deletion.
C    These are identified in FINDPROTECTED and mark changes in bc's
C
      IF(ICNTPRO.GT.0)THEN
        DO K=1,ICNTPRO
          IPK=IPROTECTED(K)
          DO J=1,IBPCNT
            IF(IPK.EQ.IBP(J))NBUSE(1,J)=0
          ENDDO
        ENDDO
      ENDIF
C  
C   next flag for removal duplicate nodes created for a fault
C   i.e., if node in position J of IBP is found in IFBC1 (negative)
C   if either node in a seg is in this group the seg also flagged
C
      NDUP=0
      NDREM=0
      IF((IFLT.GT.0).AND.(NFPF3.GT.0))THEN
        DO J=1,IBPCNT
          DO K=1,NFPF3
            IF(IBP(J).EQ.-IFBC1(K))THEN
              NDUP=NDUP+1
              NBUSE(1,J)=-1
            ENDIF
          ENDDO
        ENDDO
        DO J=1,IBPCNT
          IF(NBUSE(1,J).EQ.-1)THEN
            ISEGCP(3,NBUSE(4,J))=-1
            ISEGCP(3,NBUSE(5,J))=-1
          ENDIF
        ENDDO
        DO I=1,NSEG
          IF(ISEGCP(3,I).EQ.-1)NDREM=NDREM+1
        ENDDO
        IF((NDUP.GT.0).OR.(NDREM.GT.0))WRITE(*,10200)NDUP,NDREM
      ENDIF
10200 FORMAT('Duplicated nodes / segments used in fault(s) removed: ',
     :       I5,' / ',I5)
C
C    calculate neighbouring segment lengths and angle subtended
C
      ICOUNT=0
      SEGAVL=0.0
      SEGMIN=99.9E+9
      SEGMAX=0.0
      JCOUNT=0
      JMIN=0
      DO J=1,IBPCNT
        IF(NBUSE(1,J).EQ.2)THEN
C         WRITE(*,*)J,NBUSE(2,J),NBUSE(3,J),EX(IBP(J)),
C    :              EX(IBP(NBUSE(2,J))),EX(IBP(NBUSE(3,J)))
          XDIF1=EX(IBP(NBUSE(2,J)))-EX(IBP(J))
          YDIF1=EY(IBP(NBUSE(2,J)))-EY(IBP(J))
          XDIF2=EX(IBP(J))-EX(IBP(NBUSE(3,J)))
          YDIF2=EY(IBP(J))-EY(IBP(NBUSE(3,J)))
          XDIF3=XDIF1+XDIF2
          YDIF3=YDIF1+YDIF2
          ASQUARE=XDIF1*XDIF1 + YDIF1*YDIF1
          BSQUARE=XDIF2*XDIF2 + YDIF2*YDIF2
          CSQUARE=XDIF3*XDIF3 + YDIF3*YDIF3
          AAA=SQRT(ASQUARE)
          BBB=SQRT(BSQUARE)
          COSSUB=(ASQUARE+BSQUARE-CSQUARE)/(2.0*AAA*BBB)
          SEGPRO(1,J)=AAA
          SEGPRO(2,J)=BBB
          ICOUNT=ICOUNT+1
          AANDB=AAA+BBB
          SEGAVL=SEGAVL+AANDB
          IF(AANDB.GT.SEGMAX)SEGMAX=AAA+BBB
          IF((AANDB).LT.SEGMIN)THEN
            SEGMIN=AANDB
            JMIN=J
          ENDIF
C
C   if angle subtended < 170 degrees, prevent deletion
C
          IF(COSSUB.GT.-0.98480775)THEN
            NBUSE(1,J)=0
            JCOUNT=JCOUNT+1
C           WRITE(*,*)'BNDSEGS: retained nodes J,IBP(J),COSSUB=',
C    :                JCOUNT,J,IBP(J),COSSUB
          ENDIF
        ENDIF
      ENDDO
      SEGAVL=SEGAVL/FLOAT(ICOUNT)
      THRESH=0.4*SEGAVL
      WRITE(*,10102)SEGMIN,SEGMAX,SEGAVL,THRESH
10102 FORMAT('Segments, min: ',F12.5,', max: ',F12.5,', average:',
     :        F12.5,', threshold: ',F12.5)
      IF(SEGMIN.GT.THRESH)GO TO 20
C
C   Start removing nodes for which the seg length a+b is less than
C   THRESH - using 0.4 * average a+b.  First go to minimum (a+b)
C   set flags to remove one node JMIN, and one seg, then remake 
C   connections. Remove right side segment - either way should work
C
   15 CONTINUE
      IF(JMIN.NE.0)THEN
        NBUSE(1,JMIN)=-1              !  removed node flag
        JLEFT=NBUSE(2,JMIN)           !  neighbour nodes
        JRIGHT=NBUSE(3,JMIN)
        NBUSE(1,JLEFT)=0              !  protect neighbour nodes
        NBUSE(1,JRIGHT)=0
        KLEFT=NBUSE(4,JMIN)           !  neighbour segments
        KRIGHT=NBUSE(5,JMIN)
        ISEGCP(3,KRIGHT)=-1           !  removed seg flag
        WRITE(*,10103)JMIN,KRIGHT,KLEFT,JLEFT,JRIGHT,SEGMIN
10103   FORMAT('Removing node ',I5,' and seg ',I5,' Retaining Seg ',I5,
     :         '  (nodes ',I5,' - ',I5,' length ',F10.4,')')
C       WRITE(*,*)'Removing Seg ',KRIGHT,' nodes ',
C    :            ISEGCP(1,KRIGHT),ISEGCP(2,KRIGHT)
C       WRITE(*,*)'Retaining Seg ',KLEFT,' nodes ',
C    :            ISEGCP(1,KLEFT),ISEGCP(2,KLEFT)
C
C   connect retained segment to right node.  Note that left / right
C   nodes may appear in either sequence in ISEGCP
C
        IF(ISEGCP(1,KLEFT).EQ.JLEFT)THEN
          ISEGCP(2,KLEFT)=JRIGHT
        ELSEIF(ISEGCP(2,KLEFT).EQ.JLEFT)THEN
          ISEGCP(1,KLEFT)=JRIGHT
        ENDIF
C
C   re-assign left and right neighbours to JLEFT and JRIGHT
C
        IF(NBUSE(2,JLEFT).EQ.JMIN)THEN
          NBUSE(2,JLEFT)=JRIGHT
        ELSEIF(NBUSE(3,JLEFT).EQ.JMIN)THEN
          NBUSE(3,JLEFT)=JRIGHT
        ENDIF
        IF(NBUSE(2,JRIGHT).EQ.JMIN)THEN
          NBUSE(2,JRIGHT)=JLEFT
        ELSEIF(NBUSE(3,JRIGHT).EQ.JMIN)THEN
          NBUSE(3,JRIGHT)=JLEFT
        ENDIF
        SEGPRO(1,JRIGHT)=SEGPRO(1,JMIN)+SEGPRO(2,JMIN)
C       WRITE(*,*)'JLEFT, neighbours',JLEFT,
C    :             NBUSE(2,JLEFT),NBUSE(3,JLEFT)
C       WRITE(*,*)'JRIGHT, neighbours',JRIGHT,
C    :             NBUSE(2,JRIGHT),NBUSE(3,JRIGHT)
C       WRITE(*,*)'seg lengths',SEGPRO(1,JMIN),SEGPRO(2,JMIN),
C    :              SEGPRO(1,JRIGHT)
C
C   find new minimum seglength for unprotected nodes
C
      SEGMIN=SEGMAX
      DO J=1,IBPCNT
        IF(NBUSE(1,J).EQ.2)THEN
          SEGL=SEGPRO(1,J)+SEGPRO(2,J)
          IF(SEGL.LT.SEGMIN)THEN
            SEGMIN=SEGL
            JMIN=J
          END IF
        ENDIF
      ENDDO
      ENDIF                            ! JMIN .n. 0
C
C     WRITE(*,*)'New shortest seg: ',JMIN,'length ',SEGMIN,
C    :          ' average seg length: ',SEGAVL
C
C   Repeat the node removal process if SEGMIN small enough
C
      IF(SEGMIN.LT.THRESH)GO TO 15
   20 CONTINUE
      MAXN=0
      DO I = 1,NSEG
        IF(ISEGCP(3,I).LT.0)MAXN=MAXN+1
      ENDDO
      MAXS=0
      DO I = 1,IBPCNT
        IF(NBUSE(1,I).LT.0)MAXS=MAXS+1
      ENDDO
      WRITE(*,10113)MAXN,MAXS
10113 FORMAT('For removal: total nodes =',I5,' total segs =',I5)
C
C     WRITE(*,*)'ISEGCP array'
C     DO I=1,NSEG
C       WRITE(*,*)'ISEGCP',ISEGCP(1,I),ISEGCP(2,I),ISEGCP(3,I)
C     ENDDO
C
C     Finally remove deleted nodes from the IBP and deleted SEGS
C     from the ISEGCP arrays
C
C     Warning: if faults intersected boundaries or other faults
C     we need to back out the replaced nodes at the intersection
C     points in ISEGCP array.  This is not taken care of yet.
C
      NODDEL=0
      DO J=1,IBPCNT
        IF(NBUSE(1,J).EQ.-1)THEN
          NODDEL=NODDEL+1
        ELSE
          IBP(J-NODDEL)=IBP(J)
          DO I=1,NSEG
            IF(ISEGCP(1,I).EQ.J)ISEGCP(1,I)=J-NODDEL
            IF(ISEGCP(2,I).EQ.J)ISEGCP(2,I)=J-NODDEL
          ENDDO
        ENDIF
      END DO
      NEWNBP=IBPCNT-NODDEL
      NSEGDEL=0
      DO I=1,NSEG
        IF(ISEGCP(3,I).EQ.-1)THEN
          NSEGDEL=NSEGDEL+1
        ELSE
          DO K=1,3
            ISEGCP(K,I-NSEGDEL)=ISEGCP(K,I)
          ENDDO
        ENDIF
      END DO
      NSEGNEW=NSEG-NSEGDEL
C
C    final check that all nodes in ISEG present in IBP
C
      MAXIM=0
      DO I = 1,NSEGNEW
        IF(ISEGCP(1,I).GT.MAXIM)MAXIM=ISEGCP(1,I)
        IF(ISEGCP(2,I).GT.MAXIM)MAXIM=ISEGCP(2,I)
      ENDDO
C     WRITE(*,*)'End BNDSEGS, NSEGNEW,MAXIM = ',NSEGNEW,MAXIM
      IF(MAXIM.GT.NEWNBP)THEN
         WRITE(*,*)'Something wrong in BNDSEGS'
         IERR=1
         RETURN
      ENDIF
C
C     create poly struct, based now on NEWNBP nodes.  NEWNBP is the
C     number of vertex nodes that are passed to triangle for mesh
C     building (no duplicates)
C
      CALL REMESH_DATA(IBP,ISEGCP,NEWNBP,NSEGNEW,IERR)
      RETURN
      END SUBROUTINE BNDSEGS
C
C**************************************************************************
C
      SUBROUTINE LOCNOD(EX,EY,EXNEW,EYNEW,LEM,NOR,
     :                  LENO,ELAD,NE,NUP,NEWNUP,NPTOUT)
C
C    Routine to identify locations of new nodes in old mesh, prior
C    to interpolation. Locations are identified by 1) element number
C    LENO(NEWNUP) and 2) local coordinates ELAD(3,NEWNUP)
C
      DIMENSION LENO(NEWNUP),ELAD(3,NEWNUP)
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE),NOR(NUP)
      DIMENSION EXNEW(NEWNUP),EYNEW(NEWNUP)
      DOUBLE PRECISION COOL(3),QOOL(6),XCOLD,YCOLD,XOLD(3),YOLD(3)
      DOUBLE PRECISION TRI,X3M1,Y3M1,X1M2,Y1M2
      DOUBLE PRECISION EPS
C
      EPS=1.D-4
      NPTOUT=0
C
C   initialise the new arrays
C
      DO I=1,NEWNUP
        LENO(I)=0
        ELAD(1,I)=0.D0
        ELAD(2,I)=0.D0
        ELAD(3,I)=0.D0
      ENDDO
C
C    For each set of coordinates in the new mesh
C
      DO 50 IP=1,NEWNUP
        XP=EXNEW(IP)
        YP=EYNEW(IP)
        DELMIN=1.E35
        NNEAR=0
        IFOUND=0
C
C    Test proximity to element in old mesh
C
        DO 80 N=1,NE
C
C     Calculate the centroid coordinates
C
          XCOLD=0.D0
          YCOLD=0.D0
          DO K=1,3
            LK=NOR(LEM(K,N))
            XOLD(K)=EX(LK)
            YOLD(K)=EY(LK)
            XCOLD=XCOLD+XOLD(K)
            YCOLD=YCOLD+YOLD(K)
          ENDDO
C
C      find max distance from centroid of vertices
C
          XCOLD=XCOLD/3.D0
          YCOLD=YCOLD/3.D0
          RDOLD=(XOLD(1)-XCOLD)**2 + (YOLD(1)-YCOLD)**2
          DO J=2,3
            RDIST=(XOLD(J)-XCOLD)**2 + (YOLD(J)-YCOLD)**2
            IF(RDIST.GT.RDOLD)RDOLD=RDIST
          ENDDO
          RDOLD=SQRT(RDOLD)
C
C      nearest old centroid to current new point
C
          EDIST=SQRT((XCOLD-XP)**2+(YCOLD-YP)**2)
          IF(EDIST.LT.DELMIN)THEN
            DELMIN=EDIST
            NNEAR=N
          ENDIF
C
C  calculate local coordinates of new point in this triangle
C  first translate to (X1,Y1) as a local origin, for computing local 
C  coordinates and determinants
C
          IF(EDIST.LE.(RDOLD*(1.+EPS)))THEN
            Y3M1=YOLD(3)-YOLD(1)
            X3M1=XOLD(3)-XOLD(1)
            Y1M2=YOLD(1)-YOLD(2)
            X1M2=XOLD(1)-XOLD(2)
C
C    TRI is twice the area of the triangle element
C
            TRI=Y1M2*X3M1-X1M2*Y3M1
C
C     Calculate the natural coordinates of the (new) vertex points
C     within the (old) element
C
            COOL(2)=(XP-XOLD(1))*Y3M1-(YP-YOLD(1))*X3M1
            COOL(3)=(XP-XOLD(1))*Y1M2-(YP-YOLD(1))*X1M2
            COOL(1)=TRI-COOL(2)-COOL(3)
            TPEPS=-TRI*EPS
            IF((COOL(1)*(TRI-COOL(1)).GE.TPEPS).AND.
     :         (COOL(2)*(TRI-COOL(2)).GE.TPEPS).AND.
     :         (COOL(3)*(TRI-COOL(3)).GE.TPEPS))THEN
C
C     stop checking elements if circumscribing element found
C
              IFOUND=1
              GO TO 49
            ENDIF   ! if the point is in the old element
          ENDIF     ! if the point is near the old element
   80   CONTINUE ! loop on elements (N)
C
C   if new point is not found in any element, use average vertex
C   values from nearest element
C
   49   CONTINUE
        IF(IFOUND.EQ.0)THEN
          NPTOUT=NPTOUT+1          ! NPTOUT reported in tri_remesh
          IF(NPTOUT.LE.20)WRITE(6,1011)XP,YP,NNEAR
 1011     FORMAT('Point: ',F12.5,', ',F12.5,'  near element ',I5,
     :         ' outside old mesh')
          IF(NPTOUT.EQ.20)WRITE(6,1012)
 1012     FORMAT('Further output of problem points is suppressed')
          LENO(IP)=NNEAR
          DO JJ=1,3
            COOL(JJ)=1.D0/3.D0
          ENDDO
        ELSE
C
C     store the location parameters
C
          LENO(IP)=NNEAR
          DO JJ=1,3
            ELAD(JJ,IP)=COOL(JJ)/TRI
          END DO
        ENDIF
   50 CONTINUE   ! loop on points (IP)
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE NTRPLT(UVPO,UVPN,BGU,LEM,NOR,LENO,ELAD,NE,NUP,
     :                  NEWNUP,IQ)
C
C    Routine to interpolate the nodal values from the old finite
C     element mesh onto a new mesh using the location parameters
C     obtained by LOCNOD
C
      DIMENSION UVPO(NUP),UVPN(NEWNUP)
      DIMENSION LENO(NEWNUP),ELAD(3,NEWNUP)
      DIMENSION LEM(6,NE),NOR(NUP)
      DOUBLE PRECISION COOL(3),QOOL(6),VPO(6)
      IF((IQ.LE.0).OR.(IQ.GT.2))THEN
        WRITE(*,*)'Invalid value of IQ used in NTRPLT; stopping'
        STOP
      END IF
C
C   initialise the new array
C
      DO I=1,NEWNUP
        UVPN(I)=BGU
      ENDDO
C
C    For each nodal point in the new mesh
C
      DO 50 IP=1,NEWNUP
        NNEAR=LENO(IP)
C
C     get the (old) nodal values for this element
C   IQ=2 for quadratic interpolation, =1 for linear interpolation
C
        DO K=1,3*IQ
          LK=NOR(LEM(K,NNEAR))
          VPO(K)=UVPO(LK)
        ENDDO
C
C     Interpolate using the quadratic interpolation function
C
        IF(IQ.EQ.2)THEN
          DO KE=1,3
            COOL(KE)=ELAD(KE,IP)
            QOOL(KE)=COOL(KE)*(2.D0*COOL(KE)-1.D0)
          ENDDO
          DO KE=4,6
            KPF=KE-3
            KPB=MOD(KE+1,3)+1
            QOOL(KE)=4.D0*COOL(KPF)*COOL(KPB)
          ENDDO
C
C    or linear interpolation function
C
        ELSE
          DO KE=1,3
            QOOL(KE)=ELAD(KE,IP)
          ENDDO
        ENDIF
C
C    interpolate and store required node values
C
        VPSUM=0.D0
        DO KE=1,3*IQ
          VPSUM=VPSUM+VPO(KE)*QOOL(KE)
        ENDDO
        UVPN(IP)=VPSUM
   50 CONTINUE   ! loop on new nodal points (IP)
      RETURN
      END

      SUBROUTINE NTRPLTD(IDEN,ODMESH,RDMESH,BGD,IVV,OVMESH,RVMESH,BGV,
     :                  ITHDI,OTHDINT,THDINT,BGDI,
     :                  EX,EY,EXNEW,EYNEW,LEM,LEMNEW,NOR,IMAT,IMATNEW,
     :                  NE,NUP,NEWNE,NEWNUP,IMREG)
C
C    Routine to interpolate the values from integration points of
C    the finite element mesh onto a new mesh, used by remesh_data
C
      DOUBLE PRECISION XOLD(3),YOLD(3),COOL(3),CNL
      DOUBLE PRECISION VNODD(6),VNODV(6),VNODT(6)
      DOUBLE PRECISION XNEW(3),YNEW(3),XCOLD,YCOLD,XCNEW,YCNEW
      DOUBLE PRECISION TRI,X3M1,Y3M1,X1M2,Y1M2
      DOUBLE PRECISION XPOS(7),YPOS(7),QOOL(6),ARR(7,6),WGT(4)
      DOUBLE PRECISION XVRT,YVRT,XMID,YMID,EPS
      COMMON/AI/LUW,LSC,LBC,LLG
      DIMENSION ODMESH(7,NE),RDMESH(7,NEWNE)
      DIMENSION OTHDINT(7,NE),THDINT(7,NEWNE)
      DIMENSION OVMESH(8,NE),RVMESH(8,NEWNE)
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE)
      DIMENSION EXNEW(NEWNUP),EYNEW(NEWNUP),LEMNEW(6,NEWNE)
      DIMENSION IMAT(NE),IMATNEW(NEWNE)
      DIMENSION NOR(NUP)
C
C   the following matrix is the 7x6 matrix that projects from the
C   7 interpolation points to a best fit set of 6 nodal values
C   WGT is used to determine location of integration points
C
      SAVE ARR,WGT
      DATA ARR/1.9743924,0.1434449,0.1434449,-0.4126757,-0.4126757,
     : 0.2563768,-0.6923076,0.1434449,1.9743924,0.1434449,0.2563768,
     :-0.4126757,-0.4126757,-0.6923076,0.1434449,0.1434449,1.9743924,
     :-0.4126757,0.2563768,-0.4126757,-0.6923076,0.09812581,0.1562660,
     :0.09812581,1.0345624,-0.2822974,-0.2822974,0.1775148,0.09812581,
     :0.09812581,0.1562660,-0.2822974,1.0345624,-0.2822974,0.1775148,
     :0.1562660,0.09812581,0.09812581,-0.2822974,-0.2822974,1.0345624,
     :0.1775148/
      DATA WGT/0.69614048,0.30385953,0.82085238,0.17914761/
C
      EPS=1.D-4
C
C    min and max of old arrays (for debugging)
C
      IF(IDEN.GT.0)THEN
        CALL REPORT(ODMESH,7,NE,1,7,DEMINB,DEMAXB)
      ENDIF
      IF(IVV.GT.0)THEN
        CALL REPORT(OVMESH,8,NE,1,7,VEMINB,VEMAXB)
        CALL REPORT(OVMESH,8,NE,8,8,SEMINB,SEMAXB)
      ENDIF
C
C   initialise the new arrays
C
      DO I=1,NEWNE
        IF(IDEN.GT.0)THEN
          DO J=1,7
            RDMESH(J,I)=BGD
          ENDDO
        ENDIF
        IF(IVV.GT.0)THEN
          DO J=1,7
            RVMESH(J,I)=BGV
          ENDDO
          RVMESH(8,I)=1.0
        ENDIF
        IF(ITHDI.GT.0)THEN
          DO J=1,7
            THDINT(J,I)=BGDI
          ENDDO
        ENDIF
      ENDDO
C
C   remeshing should require that region numbers are set for all parts of 
C   the solution domain.  If region numbers aren't defined it will treat 
C   all parts of the solution domain as belonging to the same region.
C
C   Warning: a dicontinuity in region numbers that is not aligned with a 
C   boundary will cause problems in the interpolation of arrays
C
C    Look at each element in new array
C
      IREGN=1
      DO 80 I=1,NEWNE
        IF (IMREG.NE.0) IREGN=IMATNEW(I)
C
C     Calculate the geometrical coefficients
C
        XCNEW=0.D0
        YCNEW=0.D0
        DO J=1,3
          XNEW(J)=EXNEW(LEMNEW(J,I))
          YNEW(J)=EYNEW(LEMNEW(J,I))
          XCNEW=XCNEW+XNEW(J)
          YCNEW=YCNEW+YNEW(J)
        ENDDO
C
C      centroid, and max distance from centroid
C
        XCNEW=XCNEW/3.D0
        YCNEW=YCNEW/3.D0
        RDNEW=(XNEW(1)-XCNEW)**2 + (YNEW(1)-YCNEW)**2
        DO J=2,3
          RDIST=(XNEW(J)-XCNEW)**2 + (YNEW(J)-YCNEW)**2
          IF(RDIST.GT.RDNEW)RDNEW=RDIST
        ENDDO
        RDNEW=SQRT(RDNEW)
C
C      the integration points on the new element are at:
C      for vertex node,   (X,Y)=0.69614048(XV,YV)+0.30385953(XC,YC)
C      for midpoint node, (X,Y)=0.82085238(XM,YM)+0.17914761(XC,YC)
C      (obtained from manipulations of data on p421 of Huebner)
C
        XVRT=WGT(2)*XCNEW
        YVRT=WGT(2)*YCNEW
        XMID=WGT(4)*XCNEW
        YMID=WGT(4)*YCNEW
        DO K=1,3
          XPOS(K)=EXNEW(LEMNEW(K,I))*WGT(1) + XVRT
          YPOS(K)=EYNEW(LEMNEW(K,I))*WGT(1) + YVRT
        ENDDO
        DO K=4,6
          XPOS(K)=EXNEW(LEMNEW(K,I))*WGT(3) + XMID
          YPOS(K)=EYNEW(LEMNEW(K,I))*WGT(3) + YMID
        ENDDO
        XPOS(7)=XCNEW
        YPOS(7)=YCNEW
C
C    look for old triangles that may overlap current (new) triangle
C
        DELMIN=1.D35
        NNEAR=0
        IREG=1
        DO 50 N=1,NE
C
C     don't interpolate to a different region
C
          IF(IMREG.NE.0) IREG=IMAT(N)
          IF(IREGN.EQ.IREG)THEN
C
C    find the centroid of the old element
C
            XCOLD=0.0
            YCOLD=0.0
            DO J=1,3
              XOLD(J)=EX(LEM(J,N))
              YOLD(J)=EY(LEM(J,N))
              XCOLD=XCOLD+XOLD(J)
              YCOLD=YCOLD+YOLD(J)
            ENDDO
C
C      centroid, and max distance from centroid
C
            XCOLD=XCOLD/3.D0
            YCOLD=YCOLD/3.D0
            RDOLD=(XOLD(1)-XCOLD)**2 + (YOLD(1)-YCOLD)**2
            DO J=2,3
              RDIST=(XOLD(J)-XCOLD)**2 + (YOLD(J)-YCOLD)**2
              IF(RDIST.GT.RDOLD)RDOLD=RDIST
            ENDDO
            RDOLD=SQRT(RDOLD)
C
C      check for possible overlap and identify nearest old element
C      must be on same region number for this to work
C
            EDIST=SQRT((XCNEW-XCOLD)**2+(YCNEW-YCOLD)**2)
            IF(EDIST.LT.DELMIN)THEN
              DELMIN=EDIST
              NNEAR=N
            ENDIF
            IF(EDIST.LE.((RDNEW+RDOLD)*(1.+EPS)))THEN
C
C  aiming to find local coordinates of new int.point in old element
C  first translate to (X1,Y1) as a local origin, for computing local 
C  coordinates and determinants
C
              Y3M1=YOLD(3)-YOLD(1)
              X3M1=XOLD(3)-XOLD(1)
              Y1M2=YOLD(1)-YOLD(2)
              X1M2=XOLD(1)-XOLD(2)
C
C    TRI is twice the area of the triangle element
C
              TRI=Y1M2*X3M1-X1M2*Y3M1
C
C     The best-fit equivalent node values (old triangle) are:
C     Interpolation of VHB is done on logarithm of VHB, folled by EXP
C
              IF(IDEN.GT.0)THEN
                DO KE=1,6
                  VNODD(KE)=0.0
                  DO KPP=1,7
                    VNODD(KE)=VNODD(KE)+ARR(KPP,KE)*ODMESH(KPP,N)
                  ENDDO
                ENDDO
              ENDIF
              IF(IVV.GT.0)THEN
                DO KE=1,6
                  VNODV(KE)=0.0
                  DO KPP=1,7
                    VNODV(KE)=VNODV(KE)+ARR(KPP,KE)*LOG(OVMESH(KPP,N))
                  ENDDO
                ENDDO
              ENDIF
              IF(ITHDI.GT.0)THEN
                DO KE=1,6
                  VNODT(KE)=0.0
                  DO KPP=1,7
                    VNODT(KE)=VNODT(KE)+ARR(KPP,KE)*OTHDINT(KPP,N)
                  ENDDO
                ENDDO
              ENDIF
C
C     Calculate the natural coordinates of the (new) integration points
C     within the (old) element
C
              DO 60 J=1,7
                COOL(2)=(XPOS(J)-XOLD(1))*Y3M1-(YPOS(J)-YOLD(1))*X3M1
                COOL(3)=(XPOS(J)-XOLD(1))*Y1M2-(YPOS(J)-YOLD(1))*X1M2
                COOL(1)=TRI-COOL(2)-COOL(3)
                TPEPS=-TRI*EPS
                IF((COOL(1)*(TRI-COOL(1)).GE.TPEPS).AND.
     :             (COOL(2)*(TRI-COOL(2)).GE.TPEPS).AND.
     :             (COOL(3)*(TRI-COOL(3)).GE.TPEPS))THEN
C
C     If we reach here the point is within the triangle. Compute the 
C     quadratic interpolation fn values at the given local coordinates.
C
                  DO KE=1,3
                    COOL(KE)=COOL(KE)/TRI
                    QOOL(KE)=COOL(KE)*(2.D0*COOL(KE)-1.D0)
                  ENDDO
                  DO KE=4,6
                    KPF=KE-3
                    KPB=MOD(KE+1,3)+1
                    QOOL(KE)=4.D0*COOL(KPF)*COOL(KPB)
                  ENDDO
C
C   the interpolated function value then is inner product of VNOD and QOOL
C
                  IF(IDEN.GT.0)THEN
                    SUMP=0.D0
                    DO KE=1,6
                      SUMP=SUMP+QOOL(KE)*VNODD(KE)
                    ENDDO
                    RDMESH(J,I)=SUMP
                  ENDIF
                  IF(IVV.GT.0)THEN
                    SUMP=0.D0
                    DO KE=1,6
                      SUMP=SUMP+QOOL(KE)*VNODV(KE)
                    ENDDO
                    RVMESH(J,I)=EXP(SUMP)
                  ENDIF
                  IF(ITHDI.GT.0)THEN
                    SUMP=0.D0
                    DO KE=1,6
                      SUMP=SUMP+QOOL(KE)*VNODT(KE)
                    ENDDO
                    THDINT(J,I)=SUMP
                  ENDIF
                ENDIF  ! in old element or not
   60         CONTINUE ! J (interpolation points of the new triangle)
            ENDIF      ! (possible overlap)
          ENDIF        ! (same region? IREG.EQ.IREGN ?)
   50   CONTINUE       ! I (elements of old mesh)
C
C   if overlap not found, set to value of centroid of nearest element
C
        IF(NNEAR.NE.0)THEN
          DO J=1,7
            IF (IDEN.GT.0) THEN
               IF (ABS(RDMESH(J,I)-BGD).LT.EPS)
     :           RDMESH(J,I)=ODMESH(7,NNEAR)
            ENDIF
            IF (IVV.GT.0) THEN
               IF (ABS(RVMESH(J,I)-BGV).LT.EPS)
     :           RVMESH(J,I)=OVMESH(7,NNEAR)
            ENDIF
            IF (ITHDI.GT.0) THEN
               IF (ABS(THDINT(J,I)-BGDI).LT.EPS)
     :           THDINT(J,I)=OTHDINT(7,NNEAR)
            ENDIF
          ENDDO
C
C    also set the strain-rate index
C
          IF(IVV.GT.0)RVMESH(8,I)=OVMESH(8,NNEAR)
        ELSE
          WRITE(6,*)'warning: NTRPLTD: element ',I,' location problem'
        ENDIF
C
   80 CONTINUE         ! end loop N (elements of new mesh)
C
C    min and max of new array (for debugging)
C
      IF(IDEN.GT.0)THEN
        CALL REPORT(ODMESH,7,NE,1,7,DEMINA,DEMAXA)
        WRITE(LSC,1002)'Density:  ',DEMINB,DEMINA,DEMAXB,DEMAXA
      ENDIF
      IF(IVV.GT.0)THEN
        CALL REPORT(RVMESH,8,NEWNE,1,7,VEMINA,VEMAXA)
        WRITE(LSC,1002)'Viscosity:',VEMINB,VEMINA,VEMAXB,VEMAXA
        CALL REPORT(RVMESH,8,NEWNE,8,8,SEMINA,SEMAXA)
        WRITE(LSC,1002)'Exponent: ',SEMINB,SEMINA,SEMAXB,SEMAXA
      ENDIF
 1002 FORMAT(A10,' before and after, minima ',2G14.5,' maxima ',2G14.5)
C
C    find any unset points; alert user
C    unset points may have been caused by a discontinuity in IMAT
C    numbers not coinciding with a protected region boundary.
C    but there may be other causes. Unset points are set to
C    the value BGD/BGV : set this outside the expected range of values.
C
      INOTSET=0
      DO N=1,NEWNE
        IF(IDEN.GT.0)THEN
          DO J=1,7
            IF (ABS(RDMESH(J,N)-BGD).LT.EPS)THEN
              INOTSET=INOTSET+1
              WRITE(LSC,1022)'density  ',INOTSET,J,N,BGD
            ENDIF
          ENDDO
        ENDIF
        IF(IVV.GT.0)THEN
          DO J=1,7
            IF (ABS(RVMESH(J,N)-BGV).LT.EPS)THEN
              INOTSET=INOTSET+1
              WRITE(LSC,1022)'viscosity',INOTSET,J,N,BGD
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF(INOTSET.NE.0)THEN
        WRITE(LSC,1001)INOTSET
        WRITE(LUW,1001)INOTSET
 1001   FORMAT('Warning (NTRPLTD): ',I5,' interpolation points ',
     :           'remain unset')
 1022   FORMAT(A9,' value unset for (J,N)= (',I2,I6') now set to',G12.5)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE REPORT(ARRAY,N1,NE,IV,IW,XMIN,XMAX)
C
C     routine to find and report Min and Max of an array
C
      DIMENSION ARRAY(N1,NE)
      XMIN=ARRAY(IV,1)
      XMAX=ARRAY(IV,1)
      DO I=2,NE
        DO J=IV,IW
          IF(ARRAY(J,I).LT.XMIN)XMIN=ARRAY(J,I)
          IF(ARRAY(J,I).GT.XMAX)XMAX=ARRAY(J,I)
        ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE SETMAT(EX,EY,EXNEW,EYNEW,LEM,LEMNEW,
     :                  NOR,IMAT,IMATNEW,
     :                  NE,NUP,NEWNE,NEWNUP,IERR)
C
C    Routine to set the region IMAT values for the elements
C    onto a new mesh.  Calling of this routine follows call to
C    tri_remesh, within which there is a check that all IMAT
C    values are defined.
C
      DOUBLE PRECISION COOL(3),XOLD(3),YOLD(3),XNEW(3),YNEW(3)
      DOUBLE PRECISION X3M1,Y3M1,X1M2,Y1M2,TRI,CNL,EPS
      DOUBLE PRECISION XCENT,YCENT,XEMIN,XEMAX,YEMIN,YEMAX,XCO,YCO
      COMMON/AI/LUW,LSC,LBC,LLG
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE)
      DIMENSION EXNEW(NEWNUP),EYNEW(NEWNUP),LEMNEW(6,NEWNE)
      DIMENSION IMAT(NE), IMATNEW(NEWNE)
      DIMENSION NOR(NUP)
      DATA PI,PION2,DTOR/3.141592653,1.570796327,1.7453292519943295E-2/
C
      EPS=1.D-12    ! used to determine if point is inside or just outside
      MWARN=0
      IERR=0
C
C   initialise the new mat array
C
      DO 10 I=1,NEWNE
        IMATNEW(I)=0
  10  CONTINUE
C
C    Look at each element in new array
C
      DO 80 N=1,NEWNE
C
C     Find centroid of element
C
        XCENT=0.D0
        YCENT=0.D0
        DO J=1,3
          LKJ=LEMNEW(J,N)
          XCENT=XCENT+EXNEW(LKJ)
          YCENT=YCENT+EYNEW(LKJ)
        ENDDO
        XCENT=XCENT/3.D0
        YCENT=YCENT/3.D0
C
C     For each element in old array, get geometrical coefficients
C
        DMIN=999.D20
        DO 50 I=1,NE
          XCO=0.0
          YCO=0.0
          DO J=1,3
            LKJ=NOR(LEM(J,I))
            XOLD(J)=EX(LKJ)
            YOLD(J)=EY(LKJ)
            XCO=XCO+XOLD(J)
            YCO=YCO+YOLD(J)
          ENDDO
          XCO=XCO/3.D0
          YCO=YCO/3.D0
          DIST2=(XCENT-XCO)**2 + (YCENT-YCO)**2
          IF(DIST2.LT.DMIN)THEN
            DMIN=DIST2
            ICLO=I
            XCLO=XCO
            YCLO=YCO
          ENDIF
C
C    find approx limits on old element location
C
          XEMIN=XOLD(1)
          XEMAX=XOLD(1)
          YEMIN=YOLD(1)
          YEMAX=YOLD(1)
          DO J=2,3
            IF(XOLD(J).GT.XEMAX)XEMAX=XOLD(J)
            IF(XOLD(J).LT.XEMIN)XEMIN=XOLD(J)
            IF(YOLD(J).GT.YEMAX)YEMAX=YOLD(J)
            IF(YOLD(J).LT.YEMIN)YEMIN=YOLD(J)
          ENDDO
C
C   if the centroid is in this box
C
          IF(((XEMIN-XCENT)*(XEMAX-XCENT).LE.EPS).AND.
     :       ((YEMIN-YCENT)*(YEMAX-YCENT).LE.EPS))THEN
C
C  to get determinant first translate to (X1,Y1) as a local origin, for
C     computing local coordinates and determinants
C
            Y3M1=YOLD(3)-YOLD(1)
            X3M1=XOLD(3)-XOLD(1)
            Y1M2=YOLD(1)-YOLD(2)
            X1M2=XOLD(1)-XOLD(2)
C
C    TRI is twice the area of the triangle element
C
            TRI=Y1M2*X3M1-X1M2*Y3M1
C
C     Calculate the natural coordinates (new centroid in old element)
C
            COOL(2)=(XCENT-XOLD(1))*Y3M1-(YCENT-YOLD(1))*X3M1
            COOL(3)=(XCENT-XOLD(1))*Y1M2-(YCENT-YOLD(1))*X1M2
            COOL(1)=TRI-COOL(2)-COOL(3)
            IF(((COOL(1)-TRI)*COOL(1).LE.0.D0).AND.
     :         ((COOL(2)-TRI)*COOL(2).LE.0.D0).AND.
     :         ((COOL(3)-TRI)*COOL(3).LE.0.D0))THEN
C
C     If we reach here the point is within the (old) triangle
C     first check that point is not already set
C
              IF(IMATNEW(N).NE.0)THEN
                MWARN=MWARN+1
                WRITE(LSC,1001)N,IMATNEW(N),IMAT(I)
 1001           FORMAT('Warning: new element ',I5,' material number',
     :                 ' already set to ',I4,' now ',I4)
                IF(MWARN.GT.5)GO TO 150
              ENDIF
              IMATNEW(N)=IMAT(I)
            ENDIF   ! centroid in this element
          ENDIF     ! centroid in surrounding box
C
   50   CONTINUE    ! index I of old element
C
C    If new element is outside old mesh, set to number of nearest old element
C
        IF(IMATNEW(N).EQ.0)THEN
          IMATNEW(N)=IMAT(ICLO)
          WRITE(6,1003)N,XCENT,YCENT
 1003     FORMAT('Element ',I5,' at ',2G12.5,' outside old mesh')
        ENDIF
   80 CONTINUE      ! index N of new element
C
C      check that all elements have material number set
C
      NOTSET=0
      DO I=1,NEWNE
        IF(IMATNEW(I).EQ.0)THEN
          NOTSET=NOTSET+1
        ENDIF
      ENDDO
      IF(NOTSET.NE.0)THEN
        WRITE(LUW,1002)NOTSET
        WRITE(LSC,1002)NOTSET
 1002   FORMAT('Warning: ',I5,' elements have no material ',
     :            'number after SETMAT')
      ENDIF
C
  150 CONTINUE
      IF(MWARN.GT.5)THEN
        WRITE(LUW,1004)
        WRITE(LSC,1004)
 1004   FORMAT('Appears to be significant overlay of parts of the mesh.',
     :         ' Consider changing boundary conditions')
        IERR=1          ! exit the job
      ENDIF
C
      RETURN
      END
      SUBROUTINE FINDPROTECTED(NUMPROTECTED,FOUND,
     :                        IBCTYP,IBC,IBNGH,ISEG,
     :                        NOR,NBP,NN,NUP,NSEG)
C      
C      this routine will only find external vertices which
C       mark a change in boundary conditions
C       intersect with an internal boundary
C     
      INCLUDE "limits.parameters"
      DIMENSION IBC(NBP)
      DIMENSION IBCTYP(NBP*2)
      DIMENSION IBNGH(NBP*2)
      DIMENSION NOR(NUP),ISEG(3,NSEG)
      LOGICAL FOUND(NBP)
C
      DO I=1,NBP
        FOUND(I)=.FALSE.
      END DO
      NUMPROTECTED=0
      DO 25 I=1,NBP
        IF (IBC(I).GT.NN) THEN
          IBNX = IBCTYP(I)
          IBNY = IBCTYP(I+NBP)
          NODENB1= IBNGH(I)
          NODENB2= IBNGH(I+NBP)
          DO J=1,NBP
            IF (IBC(J).EQ.NODENB1) THEN
              IBNX1 = IBCTYP(J)
              IBNY1 = IBCTYP(J+NBP)
              IF ((FOUND(J).EQV..FALSE.).AND.
     :                    (IBNX.NE.IBNX1.OR.IBNY.NE.IBNY1)) THEN
                NUMPROTECTED = NUMPROTECTED+1
                FOUND(J) = .TRUE.
              END IF
              GO TO 25
            ELSE IF (IBC(J).EQ.NODENB2) THEN
              IBNX2 = IBCTYP(J)
              IBNY2 = IBCTYP(J+NBP)
              IF ((FOUND(J).EQV..FALSE.).AND.
     :                    (IBNX.NE.IBNX2.OR.IBNY.NE.IBNY2)) THEN
                NUMPROTECTED = NUMPROTECTED+1
                FOUND(J) = .TRUE.
              END IF
              GO TO 25
            END IF
          END DO
        END IF
  25  CONTINUE
      DO I=1,NSEG
        NNODE1 = ISEG(1,I)
        NNODE2 = ISEG(2,I)
        DO K=1,NBP
          IF (IBC(K).EQ.NNODE1) INDX1=K
          IF (IBC(K).EQ.NNODE2) INDX2=K
        END DO
        DO J=I+1,NSEG
          IF ((FOUND(INDX1).EQV..FALSE.).AND.(ISEG(1,J).EQ.NNODE1.OR.
     :                      ISEG(2,J).EQ.NNODE1)) THEN
            IF ((ISEG(3,J).NE.ISEG(3,I)).AND.
     :          (ISEG(3,J).GE.EXTMIN.AND.
     :              ISEG(3,J).LE.INTLMAX)) THEN
              NUMPROTECTED = NUMPROTECTED+1
C     print *,ISEG(1,I),ISEG(3,I),ISEG(1,J),ISEG(2,J),ISEG(3,J)
              FOUND(INDX1)=.TRUE.
            END IF
          END IF
          IF ((FOUND(INDX2).EQV..FALSE.).AND.(ISEG(1,J).EQ.NNODE2.OR.
     :                      ISEG(2,J).EQ.NNODE2)) THEN
            IF ((ISEG(3,J).NE.ISEG(3,I)).AND.
     :          (ISEG(3,J).GE.EXTMIN.AND.
     :              ISEG(3,J).LE.INTLMAX)) THEN
              NUMPROTECTED = NUMPROTECTED+1
C     print *,ISEG(2,I),ISEG(3,I),ISEG(1,J),ISEG(2,J),ISEG(3,J)
              FOUND(INDX2)=.TRUE.
            END IF
          END IF
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE SETPROTECTED(NUMPROTECTED,IPROTECTED,FOUND,
     :                        IBC,NOR,NBP,NUP,IERR)
C      
C      this routine sets the array IPROTECTED, for use in the
C        remesh operation.
C     
      DIMENSION IBC(NBP)
      DIMENSION NOR(NUP)
      INTEGER IPROTECTED(NUMPROTECTED)
      LOGICAL FOUND(NBP)
C
      ICNT=0
      DO I=1,NBP
        IF (FOUND(I).EQV..TRUE.) THEN
          ICNT = ICNT+1
          IF (ICNT.GT.NUMPROTECTED) THEN
            IERR = 1
            GO TO 5
          END IF
          IPROTECTED(ICNT)=NOR(IBC(I))
        END IF
      END DO
  5   RETURN
      END
C
      REAL FUNCTION ANGLECALC(X,Y)
C      Calculate the angle between segment X(1),Y(1) to X(2),Y(2) and 
C        segment X(1),Y(1) to X(3),Y(3)
      REAL X(3), Y(3), SLEN(3)
      EPS=1E-5
      SLEN(1) = (X(2)-X(3)) * (X(2)-X(3)) + (Y(2)-Y(3)) * (Y(2)-Y(3))
      SLEN(2) = (X(1)-X(3)) * (X(1)-X(3)) + (Y(1)-Y(3)) * (Y(1)-Y(3))
      SLEN(3) = (X(1)-X(2)) * (X(1)-X(2)) + (Y(1)-Y(2)) * (Y(1)-Y(2))
      ANUM = SLEN(2) + SLEN(3) - SLEN(1)
      DENOM = 2*SQRT(SLEN(2)*SLEN(3))
      IF (DENOM.EQ.0) THEN
        ANGLECALC=0
        IF (ANUM.LT.0) ANGLECALC=3.141592654
      ELSE 
        VAL=ANUM/DENOM
        IF (VAL.GT.1.0) VAL=1.0
        IF (VAL.LT.-1.0) VAL=-1.0
        ANGLECALC= ACOS(VAL)
      ENDIF
      RETURN
      END
      REAL FUNCTION POLARANGLE(X,Y)
      INCLUDE "limits.parameters"
      REAL X, Y
      EPS=1E-5
      IF (X.EQ.0.0) THEN
        IF (Y.GE.0.0) THEN
          POLARANGLE = PION2
        ELSE
          POLARANGLE = PI * 1.5
        END IF
      ELSE
        POLARANGLE = ATAN(Y/X)
        IF (X.LT. 0.0) POLARANGLE = POLARANGLE + PI
        IF (POLARANGLE.LT.-EPS) THEN
          POLARANGLE = POLARANGLE + 2.0*PI
        ELSE IF (POLARANGLE.GT.(2.0*PI-EPS)) THEN
          POLARANGLE = POLARANGLE - 2.0*PI
        END IF 
      END IF 
      RETURN
      END
      SUBROUTINE ELEMENTVAR(RV,AREA,QUALITY)
      INCLUDE "indices.parameters"
      DIMENSION RV(64)
C
      AREA = RV(IELAREA)
      QUALITY = RV(IELQUAL)
      RETURN
      END
