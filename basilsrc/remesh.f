C*--------------------------------------------------------------------
C*    Basil / Sybil:   remesh.f  1.0  1 October 2009
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE REMESH(EX,EY,IBCTYP,IBC,QBND,SEGSQLIM,ANGLIM,
     :                  PTLOC,ISEG,LEM,INTV,IVNEW,IFLT,NFP,
     :                  NFPF3,IFBC1,IFEQV,PDIST,MSNODE,MEASUR,
     :                  LUW,LSC,IVERB,IERR)
C
      INCLUDE "indices.parameters"
      INCLUDE "limits.parameters"
      INTEGER INTV(64),IVNEW(64)
      DIMENSION MEASUR(MAXMEAS), MSNODE(MAXMEAS), PTLOC(3,MAXMEAS)
      REAL, DIMENSION(:), POINTER :: EX,EY,QBND
      INTEGER, DIMENSION(:,:), POINTER :: ISEG,LEM
      INTEGER, DIMENSION(:), POINTER :: IBC,IFBC1,IFEQV
      INTEGER, DIMENSION(:), POINTER :: IBCTYP,IBNGH
      INTERFACE
        SUBROUTINE EXTERNAL_SEG_X_(ISEG,EX,EY,NSEG,NUP,IERR)    BIND(C)
          DIMENSION EX(NUP),EY(NUP)
          DIMENSION ISEG(3,NSEG)
          INTEGER NSEG,IERR
        END SUBROUTINE
      END INTERFACE
      INTERFACE
        SUBROUTINE BNDSEGS(EX,EY,QBND,ISEG,IBC,IBCTYP,IFBC1,
     :                     NN,NBP,NUP,NSEG,IFLT,NFPF3,IVERB,IERR)
        INCLUDE "limits.parameters"
C
        INTEGER, DIMENSION(:,:), POINTER :: ISEGCP
        INTEGER, DIMENSION(:,:), POINTER :: IBPNBR
        INTEGER, DIMENSION(:), POINTER :: IBP
        INTEGER, DIMENSION(:), POINTER :: IPROT
        INTEGER, DIMENSION(:), POINTER :: IBPTYP
        REAL, DIMENSION(:), POINTER :: RBND
        REAL, DIMENSION(:,:), POINTER :: SEGPRO
        DIMENSION EX(NUP),EY(NUP),QBND(NBP*2)
        DIMENSION ISEG(3,NSEG)
        DIMENSION IBC(NBP),IBCTYP(NBP*2)
        DIMENSION IFBC1(NFPF3)
        END SUBROUTINE
      END INTERFACE
C
C    this routine calls BNDSEGS after a couple of preliminaries
C
      NUP=INTV(INUP)
      NSEG=INTV(INSEG)
      NBP=INTV(INBP)
      NN=INTV(INN)
      NE=INTV(INE)
      MSINDX=INTV(IMSINDX)
      KPBC=INTV(IKPBC)

C     WRITE(6,*) NSEG
C     DO I=1,10
C       WRITE(6,*) ISEG(1,I), ISEG(2,I), ISEG(3,I)
C     ENDDO
      CALL EXTERNAL_SEG_X_(ISEG,EX,EY,NSEG,NUP,IERR)
      IF (IERR.NE.0)THEN
        WRITE(LSC,10203)
        WRITE(LUW,10203)
        RETURN
      END IF
C
C   if periodic boundaries apply, confirm current value of PDIST
C
      IF((KPBC.GT.0).AND.(NFP.GT.NFPF3))
     :  CALL PBWIDTH(EX,EY,IFBC1,IFEQV,NUP,NFP,NFPF3,KPBC,PDIST)
      NPBPT=NFP-NFPF3
C      
C    Convert Element/local ccordinates into (X,Y) so that
C    tracked values can be re-located after remeshing
C
      CALL RESETSERIES(LEM,NE,NUP,MEASUR,MSNODE,MSINDX,PTLOC,EX,EY)
C
C reconstruct the SEG and node arrays used to define a new mesh
C
      CALL BNDSEGS(EX,EY,QBND,ISEG,IBC,IBCTYP,IFBC1,
     :             NN,NBP,NUP,NSEG,IFLT,NFPF3,IVERB,IERR)
      DO I=1,64
        IVNEW(I) = INTV(I)
      END DO
C
      IF (IERR.NE.0)THEN
        WRITE(LSC,10202)
        WRITE(LUW,10202)
        RETURN
      END IF
C
10202 FORMAT("Problem in remeshing operation (BNDSEGS).")
10203 FORMAT('External boundary segment crossing, not remeshing.'
     :         ' Consider changing boundary conditions')
      RETURN
      END
C
      SUBROUTINE BNDSEGS(EX,EY,QBND,ISEG,IBC,IBCTYP,IFBC1,
     :                   NN,NBP,NUP,NSEG,IFLT,NFPF3,IVERB,IERR)
C
C   This routine examines the list of segments used to define the 
C   current mesh build, and amends it to remove excessively short 
C   segments. The modified segment list, ISEGCP(3,NSEGNEW) and corresponding
C   set of vertex nodes, IBP(NEWNBP) is passed to remesh_data to build
C   a new mesh. Each segment in ISEG is vertex-vertex. On input, sequence
C   of segments is randomised by the build process originally used in
C   triangle, and order of two nodes in ISEG is not systematic.
C
      INCLUDE "limits.parameters"
C
C  ISEGCP, IBP allocated for use in REMESH_DATA, SEGPRO local
C
      INTEGER, DIMENSION(:,:), POINTER :: ISEGCP
      INTEGER, DIMENSION(:,:), POINTER :: IBPNBR
      INTEGER, DIMENSION(:), POINTER :: IBP
      INTEGER, DIMENSION(:), POINTER :: IPROT
      INTEGER, DIMENSION(:), POINTER :: IBPTYP
      REAL, DIMENSION(:), POINTER :: RBND
      REAL, DIMENSION(:,:), POINTER :: SEGPRO
C
C  arrays received from calling routine BEMESH
C
      DIMENSION EX(NUP),EY(NUP),QBND(NBP*2)
      DIMENSION ISEG(3,NSEG)
      DIMENSION IBC(NBP),IBCTYP(NBP*2)
      DIMENSION IFBC1(NFPF3)
      INTERFACE
        SUBROUTINE REMESH_DATA_(IBP,ISEGCP,IBPTYP,RBND,NEWNBP,
     :                          NSEGNEW,IBPSIZE,IVERB,IERR)    BIND(C)
          INTEGER, DIMENSION(:,:) :: ISEGCP
          INTEGER, DIMENSION(:) :: IBP,IBPTYP
          REAL, DIMENSION(:) :: RBND
          INTEGER NEWNBP,NSEGNEW,IBPSIZE,IVERB,IERR
        END SUBROUTINE
      END INTERFACE

C
C  IBP is a list of nodes on boundary segments (external + internal)
C  IBPCNT is the number of unique vertex points in these segments
C  for a simple closed boundary, IBPCNT = NSEG. More than 2 segs
C  connected to a node may result in IBPCNT < IBPSIZE - no problem
C  First find the set of unique nodes in ISEG and store them in IBP
C
      ALLOCATE(ISEGCP(3,NSEG),STAT=IERR)
C
C  first reorganise the segments in the same kind of sequence 
C  used in the initial setup with label numbers in ascending order
C
      CALL REORDR(ISEG,ISEGCP,NSEG,NBEXT,INBTOT)
C
C   the following call removes segments that may form a vertical 
C   needle - needed to avoid entrainment problems
C
      NSEGT=NSEG        ! NSEGT may be reset by TNDRL
      CALL TNDRL(ISEGCP,NSEG,NSEGT,EX,EY,NUP,NBEXT)
C
C   allocate space for IBP, noting that requirement is same as NSEG if
C   all boundary segments form closed loop, but possibly has one more
C   node for each unattached fault included, or one less for a doubly
c   attached fault.  Typically a few points at end of IBP will be unused.
C
      IBPSIZE=NSEGT+INBTOT
      ALLOCATE(IBP(IBPSIZE),STAT=IERR)
      ALLOCATE(IBPNBR(2,IBPSIZE),STAT=IERR)
      ALLOCATE(IBPTYP(IBPSIZE*2),STAT=IERR)
      ALLOCATE(RBND(IBPSIZE*2),STAT=IERR)
      ALLOCATE(IPROT(IBPSIZE),STAT=IERR)
      ALLOCATE(SEGPRO(2,IBPSIZE),STAT=IERR)
      DO J=1,IBPSIZE
        IBP(J)=0
        IBPNBR(1,J)=0
        IBPNBR(2,J)=0
        IPROT(J)=0
        IBPTYP(J)=-1
        IBPTYP(J+IBPSIZE)=-1
        RBND(J)=0.0
        RBND(J+IBPSIZE)=0.0
      ENDDO
C
C  next set up the IBP array, external boundary: one node per segment
C  but flag points where label changes as protected
C
      IBP(1)=ISEGCP(1,1)
      IBPNBR(1,1)=ISEGCP(1,NBEXT)       !closing the loop
      IBPNBR(2,1)=ISEGCP(2,1)
      IPROT(1)=1
      DO K=2,NBEXT
        JBP=K
        IBP(K)=ISEGCP(1,K)
        IBPNBR(1,K)=ISEGCP(1,K-1)
        IBPNBR(2,K)=ISEGCP(2,K)
        IF(ISEGCP(3,K-1).NE.(ISEGCP(3,K)))IPROT(K)=1
      ENDDO
C     DO K=1,NBEXT
C       WRITE(*,*)'IBPS ',K,IBPNBR(1,K),IBP(K),IBPNBR(2,K),ISEGCP(3,K)
C     ENDDO
C
C   for internal segments, check for prior occurrances of first and last
C
      NXTBEG=NBEXT+1            ! seg number of first seg in next boundary
      JBP=NBEXT
      DO WHILE(NXTBEG.LT.NSEGT)
        JBPS=JBP
        NXTEND=NSEGT             ! applies if last boundary
        K=NXTBEG
        IFOUND=0
C
C    boundary sections are treated in sequence delimited by [NXTBEG,NXTEND]
C
        DO WHILE((IFOUND.EQ.0).AND.(K.LT.NSEGT))
          IF(ISEGCP(3,K+1).NE.ISEGCP(3,K))THEN
            NXTEND=K             ! last seg number in current boundary
            IFOUND=1
          ENDIF
          K=K+1
        ENDDO
C       WRITE(*,*)'Boundary ',ISEGCP(3,K),' NXTBEG =',NXTBEG,
C    :            ' NXTEND =',NXTEND
C
C  set up IBPNBR arrays to identify neighbour points and protect points
C  on which the boundary label changes
C
        IFND1=0
        IFND2=0
        DO L=1,JBPS
          IF(ISEGCP(1,NXTBEG).EQ.IBP(L))THEN
            IFND1=L
            IPROT(L)=1
            IBPNBR(1,L+1)=IBP(L)
C         WRITE(*,*)'Pro1ecting node IBP(L),X,Y',
C    :               IBP(L),EX(IBP(L)),EY(IBP(L))
          ENDIF
          IF(ISEGCP(2,NXTEND).EQ.IBP(L))THEN
            IFND2=L
            IPROT(L)=1
            IBPNBR(2,L)=IBP(L)
C         WRITE(*,*)'Pro2ecting node IBP(L),X,Y',
C    :               IBP(L),EX(IBP(L)),EY(IBP(L))
          ENDIF
        ENDDO
        IF(IFND1.EQ.0)THEN       ! save unconnected first point
          JBP=JBP+1
          IBP(JBP)=ISEGCP(1,NXTBEG)
          IBPNBR(2,JBP)=ISEGCP(2,NXTBEG)
          IPROT(JBP)=1
C         WRITE(*,*)'Pro3ecting node JBP,IBP(JBP),X,Y',
C    :               JBP,IBP(JBP),EX(IBP(JBP)),EY(IBP(JBP))
        ENDIF
        DO M=NXTBEG,NXTEND-1     ! save intermediate points
          JBP=JBP+1
          IBP(JBP)=ISEGCP(2,M)
          IBPNBR(1,JBP)=ISEGCP(1,M)
          IBPNBR(2,JBP)=ISEGCP(2,M+1)
        ENDDO
        IF(IFND2.EQ.0)THEN       ! save unconnected last point
          JBP=JBP+1
          IBP(JBP)=ISEGCP(2,NXTEND) ! save unconnected last point
          IBPNBR(1,JBP)=ISEGCP(1,NXTEND)
          IPROT(JBP)=1
C         WRITE(*,*)'Pro3ecting node JBP,IBP(JBP),X,Y',
C    :               JBP,IBP(JBP),EX(IBP(JBP)),EY(IBP(JBP))
        ENDIF
        NXTBEG=NXTEND+1
      ENDDO    !   DO WHILE(NXTBEG.LT.NSEGT)
      IBPCNT=JBP           ! total number of points used in mesh spec.
      IF(IBPCNT.GT.IBPSIZE)THEN
        WRITE(*,*)'WARNING: IBPCNT > IBPSIZE ',IBPCNT,IBPSIZE,
     :            'in BNDSEGS, check setting'
        ERR=1
        STOP
      ENDIF
C
C    also protect any nodes where the boundary condition type changes
C    first associate current boundary conditions with IBP arrays
C
      NSET=0
      DO LX=1,IBPCNT
        IB=IBP(LX)
        DO K=1,NBP            ! find entries for LX in IBC arrays
          NK=IBC(K)
          IF(NK.EQ.IB)THEN
            LY=LX+IBPSIZE
            IBPTYP(LX)=IBCTYP(K)
            IBPTYP(LY)=IBCTYP(K+NBP)
            RBND(LX)=QBND(K)
            RBND(LY)=QBND(K+NBP)
            NSET=NSET+1
          ENDIF
        ENDDO
      ENDDO
      WRITE(*,10100)NSET,IBPCNT
10100 FORMAT('Boundary conditions saved for ',I5,
     :       ' nodes in IBP array of size ',I5)
C
C     then find adjacent nodes and compare IBPTYP values
C     depending on BC's this test may find some points twice or once
C
      DO LX=1,IBPCNT
        IB=IBP(LX)
        IBM1=IBPNBR(1,LX)
        DO LM=1,IBPCNT          ! find entries for LX in IBC arrays
          IF(IBP(LM).EQ.IBM1)LXM=LM
        ENDDO
        LY=LX+IBPSIZE
        LYM=LXM+IBPSIZE
        IF(((IBPTYP(LX).EQ.0).AND.(IBPTYP(LXM).NE.0)).OR.
     :     ((IBPTYP(LY).EQ.0).AND.(IBPTYP(LYM).NE.0)))THEN
          IPROT(LX)=1
C         WRITE(*,*)'Protecting node LX,IBP(LX),X,Y',
C    :               LX,IBP(LX),EX(IBP(LX)),EY(IBP(LX))
        ENDIF
        IF(((IBPTYP(LX).NE.0).AND.(IBPTYP(LXM).EQ.0)).OR.
     :     ((IBPTYP(LY).NE.0).AND.(IBPTYP(LYM).EQ.0)))THEN
          IPROT(LXM)=1
C         WRITE(*,*)'Protecting node LXM,IBP(LX),X,Y',
C    :               LXM,IBP(LXM),EX(IBP(LXM)),EY(IBP(LXM))
        ENDIF
      ENDDO
C
C    count how many nodes are protected
C
      JCOUNT=0                               ! no. of nodes to be removed
      ICNTPRO=0
      DO J=1,IBPCNT
        ICNTPRO=ICNTPRO+IPROT(J)
      ENDDO
      WRITE(*,10101)NSEGT,IBPCNT,ICNTPRO
10101 FORMAT('Mesh rebuild examines ',I5,' segments and ',I5,
     :       ' nodes with',I4,' protected')
C
C    to decide which segments/nodes can be deleted, calculate neighbouring 
C    segment lengths and angle subtended.  This section avoids start
C    and end points of structures because such points are already protected
C
      ICOUNT=0                    ! no. of segs to include in average
      SEGAVL=0.0
      SEGMIN=99.9E+9
      SEGMAX=0.0
      JMIN=0
      DO J=2,IBPCNT-1
        IB=IBP(J)
        IBM1=IBPNBR(1,J)
        IBP1=IBPNBR(2,J)
        IF(IABS(IPROT(J)).NE.1)THEN     ! avoid nodes protected or for deletion
C         WRITE(*,*)'IBM1,IB,IBP1=',IBM1,IB,IBP1
          XDIF1=EX(IBP1)-EX(IB)
          YDIF1=EY(IBP1)-EY(IB)
          XDIF2=EX(IB)-EX(IBM1)
          YDIF2=EY(IB)-EY(IBM1)
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
C   if angle subtended < 170 degrees, node is protected
C
          IF(COSSUB.GT.-0.98480775)IPROT(J)=1
        ENDIF       ! if IPROT(J) .ne. 1
      ENDDO         ! J=2,IBPCNT-1
      SEGAVL=SEGAVL/FLOAT(ICOUNT)
      THRESH=0.4*SEGAVL
      WRITE(*,10102)SEGMIN,SEGMAX,SEGAVL,THRESH
10102 FORMAT('Segments, min: ',F12.5,', max: ',F12.5,', average:',
     :        F12.5,', threshold: ',F12.5)
C
C   Start removing nodes for which the seg length a+b is less than
C   THRESH - using 0.4 * average a+b.  First go to minimum (a+b)
C   set flags to remove one node JMIN, and one seg, then remake 
C   connections. Remove right side segment - either way should work
C
      DO WHILE(SEGMIN.LT.THRESH)      ! repeated removals
        IF(JMIN.EQ.0)THEN
          WRITE(*,*)'There appears to be a problem in BNDSEGS: JMIN=0'
          STOP
        ENDIF
        JCOUNT=JCOUNT+1               ! removing another node (JMIN)
        IPROT(JMIN)=-1                ! flag for removal of IBP(JMIN)
        IPROT(JMIN-1)=1               ! protect neighbour nodes
        IPROT(JMIN+1)=1
        SEGPRO(2,JMIN-1)=SEGMIN       ! update seg lengths either side
        SEGPRO(1,JMIN+1)=SEGMIN
C
C   find new minimum seglength for unprotected nodes
C
        SEGMIN=SEGMAX
        DO J=2,IBPCNT-1
          IF(IABS(IPROT(J)).NE.1)THEN    ! not protected or removed already
            SEGL=SEGPRO(1,J)+SEGPRO(2,J)
            IF(SEGL.LT.SEGMIN)THEN
              SEGMIN=SEGL
              JMIN=J
            END IF
          ENDIF
        ENDDO
C       WRITE(*,*)'New shortest seg: ',JMIN,' IBP(JMIN) ',IBP(JMIN),
C    :            'length ',SEGMIN
C
C   Repeat node removal process if SEGMIN small enough
C
      ENDDO         ! DO WHILE(SEGMIN.LT.THRESH)
C
C    In each case where a node is removed identify segment to be removed 
C    (where IPROT(J)=-1 and IBP(JN)=ISEGCP(2,JRM)).
C
      DO JN=1,IBPCNT
        IF(IPROT(JN).EQ.-1)THEN
          IFOUND=0
          JRM=0
          DO WHILE ((IFOUND.EQ.0).AND.(JRM.LE.NSEGT))
            JRM=JRM+1
            IF(ISEGCP(2,JRM).EQ.IBP(JN))IFOUND=1
          ENDDO
          IF(JRM.NE.0)THEN
            ISEGCP(3,JRM)=-1               ! flag to remove seg
            ISEGCP(1,JRM+1)=ISEGCP(1,JRM)  ! reconnect kept segment
          ELSE
            WRITE(*,*)'BNDSEGS problem: seg with deleted node not found'
            STOP
          ENDIF
          WRITE(*,10103)IBP(JN),JRM,JRM+1,
     :                ISEGCP(1,JRM+1),ISEGCP(2,JRM+1)
10103     FORMAT('Removing node ',I5,'  and segment ',I5,
     :          '   Retaining Seg ',I5,'  (nodes ',I5,' - ',I5,')')
        ENDIF                  ! IPROT(JN) = -1
      ENDDO                    ! JN = 1,IBPCNT
C
C   Check the node and segment deletions add up
C
      MAXN=0
      DO I = 1,NSEGT
        IF(ISEGCP(3,I).LT.0)MAXN=MAXN+1
      ENDDO
      MAXS=0
      DO I = 1,IBPCNT
        IF(IPROT(I).LT.0)MAXS=MAXS+1
      ENDDO
      IF((JCOUNT.NE.MAXN).OR.(JCOUNT.NE.MAXS))THEN
        WRITE(*,10114)MAXN,MAXS,JCOUNT
10114   FORMAT('Problem with node removal: MAXN =',I4,' MAXS =',I4,
     :         '  JCOUNT - ',I4)
      ELSE
        WRITE(*,10113)MAXN,MAXS
10113   FORMAT('A total of ',I4,' nodes and ',I4,' segments will be ',
     :         'removed from the present SEG array')
      ENDIF
C
C     Remove deleted nodes from the IBP and deleted SEGS
C     from the ISEGCP arrays
C
      INEW=0
      DO I=1,NSEGT
        IF(ISEGCP(3,I).NE.-1)THEN
          INEW=INEW+1
          DO K=1,3
            ISEGCP(K,INEW)=ISEGCP(K,I)
          ENDDO
        ENDIF
      END DO
      NSEGNEW=INEW
      JNEW=0
      DO J=1,IBPCNT
        IF(IPROT(J).NE.-1)THEN
          JNEW=JNEW+1
          IBP(JNEW)=IBP(J)
          IPROT(JNEW)=IPROT(J)
          IBPTYP(JNEW)=IBPTYP(J)
          IBPTYP(JNEW+IBPSIZE)=IBPTYP(J+IBPSIZE)
          RBND(JNEW)=RBND(J)
          RBND(JNEW+IBPSIZE)=RBND(J+IBPSIZE)
        ENDIF
      END DO
      NEWNBP=JNEW
C     WRITE(*,*)'NSEG,NSEGNEW,NEWNBP',NSEGT,NSEGNEW,NEWNBP
C
C   finally swap the node numbers in ISEG for the index positions in IBP
C   sign-swapping to avoid more than one change in each ISEG entry
C
      DO J=1,NEWNBP
        NODN=IBP(J)
        DO I=1,NSEGNEW
          IF(ISEGCP(1,I).EQ.NODN)ISEGCP(1,I)=-J
          IF(ISEGCP(2,I).EQ.NODN)ISEGCP(2,I)=-J
        ENDDO
      ENDDO
      DO I=1,NSEGNEW
        ISEGCP(1,I)=-ISEGCP(1,I)
        ISEGCP(2,I)=-ISEGCP(2,I)
        IF((ISEGCP(1,I).LT.0).OR.(ISEGCP(2,I).LT.0))THEN
          WRITE(*,*)'Swapping nodes for indices in BNDSEGS appears ',
     :          'to have failed somehow, ISEG =',(ISEGCP(J,I),J=1,3)
          IERR=1
          STOP
        ENDIF
      ENDDO
C
C     write the arrays for checking if verbose mode set
C
      IF(IVERB.NE.0)THEN
        OPEN(12,FILE='SEGr')
        WRITE(12,*)'IBP array: K IBP(K) EX EY IBPTYP(x,y),IPROT'
        DO K=1,NEWNBP
          NKB=IBP(K)
          WRITE(12,'(I4,I6,2F10.5,2I4,2F10.5,I4)')K,NKB,
     :      EX(NKB),EY(NKB),IBPTYP(K),IBPTYP(K+IBPSIZE),
     :      RBND(K),RBND(K+IBPSIZE),IPROT(K)
        ENDDO
        WRITE(12,*)'SEG array prior to new mesh generation'
        DO K=1,NSEGNEW
          WRITE(12,'(4I6)')K,ISEGCP(1,K),ISEGCP(2,K),ISEGCP(3,K)
        ENDDO
        CLOSE(12)
      ENDIF
      IF(IERR.NE.0)RETURN
C
C     create poly struct, based now on NEWNBP nodes.  NEWNBP is the
C     number of vertex nodes that are passed to triangle for mesh rebuild
C
      CALL REMESH_DATA(IBP,ISEGCP,IBPTYP,RBND,NEWNBP,NSEGNEW,
     :                 IBPSIZE,IVERB,IERR)
C
C    seems to be a question whether the DEALLOCATES are necessary
C
      DEALLOCATE(IBP)
      DEALLOCATE(IBPTYP)
      DEALLOCATE(RBND)
      DEALLOCATE(IPROT)
      DEALLOCATE(SEGPRO)
      DEALLOCATE(ISEGCP)
      DEALLOCATE(IBPNBR)
      RETURN
      END SUBROUTINE BNDSEGS
C
C**************************************************************************
C
      SUBROUTINE REORDR(ISEG,ISEGCP,NSEG,NBEXT,INBTOT)
C
C    The purpose of this subroutine is to reorganise the SEG array so
C    that we can check continuity etc.
C
      DIMENSION ISEG(3,NSEG),ISEGCP(3,NSEG)
C
C   Count the number of SEGS with A. Label < 1000 (external), or
C    B. in next label number, until all SEGS accounted for
C
      IEMAX=1000          ! initially for external boundary
      LABMAX=3000
      IEMIN=0
      NLABMIN=IEMAX
      NLABMAX=0
      NKL=0
      ISP=0
      IEXT=1
      INBTOT=0
      WRITE(*,10001)NSEG
10001 FORMAT('A total of ',I5,' segments are transcribed in sequence')
      DO WHILE((ISP.LT.NSEG).AND.(NLABMIN.LT.LABMAX))
C
C  external boundary: find min and max label numbers and count segments.
C
        NXT=0
        IF(IEXT.EQ.1)THEN
          DO K=1,NSEG
            NL=ISEG(3,K)
            IF(NL.LT.IEMAX)THEN
              NXT=NXT+1            ! no. of SEGS in current boundary section
              IF(NL.GT.NLABMAX)NLABMAX=NL
              IF(NL.LT.NLABMIN)NLABMIN=NL
            ENDIF
          ENDDO
C
C   internal boundary: find next smallest label number and count segments
C
        ELSE
          NPLAB=NLABMAX+1
          IF(NPLAB.LT.IEMAX)NPLAB=IEMAX
          NLABMIN=LABMAX
          DO K=1,NSEG
            NL=ISEG(3,K)
            IF((NL.GE.NPLAB).AND.(NL.LT.NLABMIN))NLABMIN=NL
          ENDDO
          NLABMAX=NLABMIN
          DO K=1,NSEG
            NL=ISEG(3,K)
            IF(NL.EQ.NLABMIN)NXT=NXT+1
          ENDDO
        ENDIF
C       WRITE(*,*)'for NLABMIN / MAX =',NLABMIN,NLABMAX,' NXT =',NXT
C
C   find the smallest node number that is unmatched in NLABMIN
C   so presumably marks the start of a boundary section
C
        NODMIN=0
        DO K=1,NSEG
          IF(ISEG(3,K).EQ.NLABMIN)THEN
            DO L=1,2
              IPAIR=0
              DO J=1,NSEG
                IF((ISEG(3,J).EQ.NLABMIN).AND.(J.NE.K))THEN
                  IF((ISEG(L,K).EQ.ISEG(1,J)).OR.
     :               (ISEG(L,K).EQ.ISEG(2,J)))IPAIR=1
                ENDIF     !   is a matching node number found
              ENDDO
              IF(IPAIR.EQ.0)THEN
C               WRITE(*,*)'unpaired: L =',L,' node =',ISEG(L,K),
C    :          ' label =',ISEG(3,K)
                IF((ISEG(L,K).LT.NODMIN).OR.(NODMIN.EQ.0))THEN
                  NODMIN=ISEG(L,K)
                  NODSEG=K
                ENDIF
              ENDIF       !  no pair found
            ENDDO
          ENDIF           !  on section LABMIN
        ENDDO
C
C   if no unmatched node is found, assume the boundary is a closed loop
C   and just choose the smallest node number to start with
C
        IF(NODMIN.EQ.0)THEN
          DO K=1,NSEG
            IF((ISEG(3,K).GE.NLABMIN).AND.(ISEG(3,K).LE.NLABMAX))THEN
              IF(ISEG(1,K).LT.NODMIN)THEN
                NODMIN=ISEG(1,K)
                NODSEG=K
              ENDIF
              IF(ISEG(2,K).LT.NODMIN)THEN
                NODMIN=ISEG(2,K)
                NODSEG=K
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
C  we begin to transcribe the segments, starting with the segment 
C  SEGMIN and node NODMIN, just determined
C
        J=ISP+1
        NKL=NKL+1
        ISEGCP(1,J)=NODMIN
        ISEGCP(3,J)=NLABMIN
        IF(ISEG(1,NODSEG).EQ.NODMIN)THEN
          ISEGCP(2,J)=ISEG(2,NODSEG)
        ELSE
          ISEGCP(2,J)=ISEG(1,NODSEG)
        ENDIF
        ISEG(3,NODSEG)=-ISEG(3,NODSEG)    ! to flag it as done
C       WRITE(*,*)'entry J',J,ISEGCP(1,J),ISEGCP(2,J),ISEGCP(3,J)
C
C  progressively transcribe the remaining entries for SEGS on the
C  current boundary section
C
        DO WHILE(J.LT.(ISP+NXT))  ! while we still have boundary segments to fill
          IFOUND=0
          JP1=J+1
          K=0
          DO WHILE((IFOUND.EQ.0).AND.(K.LE.NSEG))
            K=K+1
            NL=ISEG(3,K)
            IF((NL.GE.0).AND.(NL.LE.NLABMAX))THEN
              IF(ISEG(1,K).EQ.ISEGCP(2,J))THEN
                ISEGCP(1,JP1)=ISEG(1,K)
                ISEGCP(2,JP1)=ISEG(2,K)
                ISEGCP(3,JP1)=ISEG(3,K)
                ISEG(3,K)=-ISEG(3,K)
                IFOUND=1
                NKL=NKL+1
              ELSE IF(ISEG(2,K).EQ.ISEGCP(2,J))THEN
                ISEGCP(1,JP1)=ISEG(2,K)
                ISEGCP(2,JP1)=ISEG(1,K)
                ISEGCP(3,JP1)=ISEG(3,K)
                ISEG(3,K)=-ISEG(3,K)
                IFOUND=1
                NKL=NKL+1
              END IF      !  if next segment found
            ENDIF         !  if label for external segment
          ENDDO
          J=JP1
C         IF(IEXT.EQ.1)
C      :  WRITE(*,*)'entry J',J,ISEGCP(1,J),ISEGCP(2,J),ISEGCP(3,J)
        ENDDO
C
C   check if this boundary section is closed or open
C
        IF(IEXT.EQ.1)THEN     ! if external boundary
          NBEXT=NXT
          IF(ISEGCP(2,J).EQ.ISEGCP(1,ISP+1))THEN
            WRITE(*,10002)ISEGCP(1,ISP+1),NBEXT
10002       FORMAT('The external boundary is closed at node ',I6,
     :             ' with ',I4,' segments')
          ELSE
            WRITE(*,10003)
10003       FORMAT('REORDR: Warning, external boundary appears to',
     :                 'not close')
          ENDIF
        ELSE                 ! if internal boundary section
          INBTOT=INBTOT+1
          IF(ISEGCP(2,J).EQ.ISEGCP(1,ISP+1))THEN
            WRITE(*,10004)NLABMIN,ISEGCP(1,ISP+1),NXT
10004       FORMAT('Boundary ',I5,' is closed at node ',I7,' with ',I4,
     :              'segments')
          ELSE
            WRITE(*,10005)NLABMIN,ISEGCP(1,ISP+1),ISEGCP(2,J),NXT
10005       FORMAT('Boundary ',I5,' starts and ends at nodes ',I6,
     :             ' and ',I6,' with ',I4,' segments')
           ENDIF
        ENDIF
        ISP=ISP+NXT         ! running total of resequenced segments
        IEXT=0              ! further boundary segments are internal
      ENDDO                 ! end of the main loop on sections
C
      IF(ISP.NE.NSEG)THEN
        WRITE(*,10006)ISP
10006   FORMAT('Segment re-assignment in REORDR appears to have ',
     :         'a problem, ISP = ',I6)
      ENDIF
      RETURN
      END
C
C**************************************************************************
C
      SUBROUTINE LOCNOD(EX,EY,EXNEW,EYNEW,LEM,
     :                  LENO,ELAD,NE,NUP,NEWNUP,NPTOUT)
C
C    Routine to identify locations of new nodes in old mesh, prior
C    to interpolation. Locations are identified by 1) element number
C    LENO(NEWNUP) and 2) local coordinates ELAD(3,NEWNUP)
C
      DIMENSION LENO(NEWNUP),ELAD(3,NEWNUP)
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE)
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
            LK=LEM(K,N)
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
      SUBROUTINE NTRPLT(UVPO,UVPN,BGU,LEM,
     :                  LENO,ELAD,NE,NUP,NEWNUP,IQ)
C
C    Routine to interpolate the nodal values from the old finite
C     element mesh onto a new mesh using the location parameters
C     obtained by LOCNOD
C
      DIMENSION UVPO(NUP),UVPN(NEWNUP)
      DIMENSION LENO(NEWNUP),ELAD(3,NEWNUP)
      DIMENSION LEM(6,NE)
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
          LK=LEM(K,NNEAR)
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
     :                  EX,EY,EXNEW,EYNEW,LEM,LEMNEW,
     :                  IMAT,IMATNEW,NE,NUP,NEWNE,NEWNUP,IMREG)
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
        CALL REPORT(RDMESH,7,NEWNE,1,7,DEMINA,DEMAXA)
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
     :                  IMAT,IMATNEW,NE,NUP,NEWNE,NEWNUP,IERR)
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
            LKJ=LEM(J,I)
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

      SUBROUTINE TNDRL(ISEGCP,NSEG,NSEGT,EX,EY,NUP,NBEXT)
C
C   This routine is designed to remove segments from an internal
C   boundary that has been folded to form a vertical "needle".  The object
C   is to remove the needle without badly affecting the adjoining boundary
C   Beware that this action can transfers minor material from one domain 
C   to another, which could cause other problems, but it is in place here
C   because any mesh created with a needle-like tendril will introduce
C   efficiency and accuracy problems, generally causing the code to fail
C   in convergence. If causing problems, comment out the CALL to TNDRL.
C
      INTEGER, DIMENSION(:), POINTER :: IREM
      DIMENSION ISEGCP(3,NSEG)
      DIMENSION EX(NUP),EY(NUP)
      ALLOCATE(IREM(NSEG),STAT=IERR)
C
C    For each internal segment (excluding first and last) find the
C    gradient.  If 'too large' flag that segment for deletion
C
      GMAX=10.0
      NODCON=0
      DO J=NBEXT+1,NSEG
        IREM(J)=0
      ENDDO
      JREM=0
      EXMIN=9999.99
      EXMAX=-9999.99
      NODCON=0
      DO J=NBEXT+2,NSEG-1
        IF(((ISEGCP(3,J)-ISEGCP(3,J+1)).EQ.0).AND.
     :     ((ISEGCP(3,J)-ISEGCP(3,J-1)).EQ.0))THEN
          J1=ISEGCP(1,J)
          J2=ISEGCP(2,J)
          GL=(EY(J2)-EY(J1))/(EX(J2)-EX(J1))
          IF(ABS(GL).GT.GMAX)THEN
            IREM(J)=1
            JREM=JREM+1             ! counter of segs removed
            IF(EX(J2).GT.EXMAX)EXMAX=EX(J2)
            IF(EX(J2).LT.EXMIN)EXMIN=EX(J2)
            IF((IREM(J-1).NE.1).AND.(NODCON.EQ.0))THEN
              NODCON=ISEGCP(2,J-1)
              EXL=EX(NODCON)
              EYL=EY(NODCON)
C             WRITE(*,*)J,ISEGCP(1,J-1),ISEGCP(2,J-1),ISEGCP(3,J-1)
            ENDIF
C
C    the node that replaces the removed section is moved to average
C    coordinates of the two proximal retained nodes
C
          ELSE
            IF((NODCON.NE.0).AND.(IREM(J-1).EQ.1))THEN
              EXR=EX(ISEGCP(1,J))
              EYR=EY(ISEGCP(1,J))
              ISEGCP(1,J)=NODCON
              EX(NODCON)=0.5*(EXL+EXR)
              EY(NODCON)=0.5*(EYL+EYR)
C             WRITE(*,*)J,ISEGCP(1,J),ISEGCP(2,J),ISEGCP(3,J)
              NODCON=0
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF(JREM.GT.0)WRITE(*,10099)JREM,NSEG,EXMIN,EXMAX,GMAX
10099 FORMAT('TNDRL Removed ',I5,' nodes of ',I5,' between X =',F9.5,
     : ' and X =',F9.5,' for slope > ',F6.3)
C
C    condense the segment array, ignoring removed segments
C
      JN=NBEXT
      DO J=NBEXT+1,NSEG
        IF(IREM(J).EQ.0)THEN
          JN=JN+1
          DO K=1,3
            ISEGCP(K,JN)=ISEGCP(K,J)
          ENDDO
        ENDIF
      ENDDO
      NSEGT=JN             ! the new segment array size
C     DO J=NBEXT+1,NSEGT
C       WRITE(*,*)J,ISEGCP(1,J),ISEGCP(2,J),ISEGCP(3,J),
C    :  EX(ISEGCP(1,J)),EX(ISEGCP(2,J))
C     ENDDO
      DEALLOCATE(IREM)
      RETURN
      END
