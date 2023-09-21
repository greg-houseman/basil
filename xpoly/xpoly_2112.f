C
C   Routine 'xpoly' is used to define a poly file as used in
C   the triangle routine used to construct a 2D finite element
C   mesh comprised of triangles.  xpoly produces an out file 
C   called x.poly.  It looks for 3 input files as follows:

C   "boundary"
C   a set of coordinate pairs, assumed to form a closed polygon 
C   provided in sequence.  Anticlockwise sequence around the
C   boundary is preferred. The final point will be joined to 
C   the first point to close the polygon. Do not repeat the first
C   point or duplicate any points. Points that are too close to
C   each other will only cause problems in the mesh generation.
C   Parts of the boundary can be labelled with boundary numbers
C   by including a line with '#' in column 1 followed by an 
C   integer label.  The label is applied to the boundary segment 
C   that follows each point.  These labels may be used to set
C   boudary conditions in a basil run.
C
C   "internal"
C   allows for multiple internal stuctures that may be open or 
C   closed, but are assumed to be non-intersecting and either
C   entriely within or connecting to the external boundary. 
C   Sets of coordinate pairs for distinct structures are separated 
C   by a single line with '#' in column 1.  'C' in column 2 of 
C   that line means the structure will be closed.  Do not repeat 
C   the first point for a closed structure. An  integer that follows 
C   on the same line can provide a label for that structure.
C
C   "connections"
C   allows for internal boundary segments that are connected to 
C   points that are already defined in "boundary" or "internal".  
C   Each connecting segment starts and ends with a point
C   that is already present - as determined by having a location
C   sufficiently close to an existing point in (x,y) space.
C   Points are listed in consecutive pairs in file "connections"
C   Distinct connecting segments separated by # in column 1 
C   of the preceding line, with a integer label following.
C
C   The segment labels may be used in a basil program to in
C   setting boundary conditions or fault conditions
C
C   revised by G. Houseman, December 2021.
C   
      PARAMETER (NBND=500, MSEG=20)
      DIMENSION XBND(2*NBND),YBND(2*NBND),XCON(NBND),YCON(NBND)
      DIMENSION ICB(MSEG),ICB1(MSEG),ICE(MSEG),ICE1(MSEG)
      DIMENSION IBN(MSEG),IIN(MSEG),ICN(MSEG),ICLOSE(MSEG)
      DIMENSION INTB1(MSEG),INTBL(MSEG),INTEB(MSEG),INTEL(MSEG)
      CHARACTER(80) INSTRING,HEADER
      CHARACTER(6) STR
      CHARACTER(10) CLOSED
C
C    MSEG is the number of separately labelled segments
C    IBN(MSEG), IIN(MSEG) are the label numbers that
C    get attached to each group of segments, for segments that
C    belong to external boundary, internal boundary, and connections
C
      DO J=1,MSEG
        INTEB(J)=0
        INTEL(J)=0
      ENDDO
C    read the sets of points on the external boundary
C
      OPEN(8,FILE='boundary',STATUS='OLD',ERR=900)
      NPT=0
      NOBND=0
      INTEB(1)=1
   10 CONTINUE
      READ(8,'(A80)',END=20)INSTRING(1:80)
C
C    if this condition entered, a new boundary segment is initialised
C
      IF(INSTRING(1:1).EQ.'#')THEN
        NOBND=NOBND+1
        IF((NOBND).GT.MSEG)THEN
          WRITE(6,*)'Insufficient memory in xpoly, increase MSEG'
          STOP
        ENDIF
        READ(INSTRING(3:20),*,END=20)IBN(NOBND)
        INTEB(NOBND)=NPT+1         ! node index of first point in segment
      ELSE                         ! read coordinate pairs
        READ(INSTRING(1:80),*,END=20)XB2,YB2
        ISEG=ISEG+1
        NPT=NPT+1
        IF(NPT.GT.NBND)THEN
          WRITE(6,*)'Insufficient memory allocated, increase NBND'
          STOP
        ENDIF
        XBND(NPT)=XB2
        YBND(NPT)=YB2
        INTEL(NOBND)=NPT       ! node index of last point in segment
      ENDIF
      GO TO 10
   20 CONTINUE                 ! end of boundary file
      CLOSE(8)
      NBT=NPT
      NTN=NPT                                ! running total: no. of points
      IF(INTEB(NOBND).EQ.0)NOBND=NOBND-1
C
C     join last point of boundary to first point without repeat
C
      INTEL(NOBND)=INTEB(1)
      IF(NPT.LT.3)THEN
        WRITE(6,*)'Insufficient points are available'
        STOP
      ENDIF
C
C    report external boundary points
C
      WRITE(6,10021)NPT,NOBND
      DO I=1,NOBND
        KF=INTEL(I)
        IF(KF.EQ.1)KF=NPT
        NIE=KF-INTEB(I)+1                    ! points added by this segment
        WRITE(6,10019)'Ex',I,IBN(I),NIE
        KO=INTEB(I)
        DO WHILE (KO.LE.KF)
          KE=KO+3
          IF(KE.GT.KF)KE=KF
          WRITE(6,1050)(XBND(J),YBND(J),J=KO,KE)
          KO=KO+4
        ENDDO
      ENDDO
10019 FORMAT(A2,'ternal structure ',I5,' label ',I5,' has ',I3,
     :           ' points',A10)
10021 FORMAT('External boundary loop has ',13X,I4,' points in',
     :       I3,' segments')
C
C    if there is a file called internal, then read this also.
C    distinct internal structures are separated by single character #
C    in column 1, preceding pairs of coordinates on subsequent lines.
C
C    A structure may be open or closed.  If closed, first and last points
C    are connected to each other, otherwise they are treated as end points
C    of the structure. A closed structure is indicated by C following #.
C    In absence of C, the structure is treated as open.
C
C    NOINT = number of distinct structures; INTB1(NOINT) and INTBL(NOINT)
C    are addresses of first and last points of structure NOINT.
C    Probably best to avoid interesting situations like intersecting
C    structures.
C
      NOINT=0
      NINTRNL=0
      OPEN(8,FILE='internal',STATUS='OLD',ERR=50)
      INTB1(1)=NPT+1
   30 CONTINUE
      READ(8,'(A80)',END=40)INSTRING(1:80)
C
C    close last loop, start next loop
C
      IF(INSTRING(1:1).EQ.'#')THEN
        HEADER=INSTRING
        IF(NOINT.GT.MSEG)THEN
          WRITE(6,*)'Insufficient memory in xpoly, increase MSEG'
          STOP
        ENDIF
        NOINT=NOINT+1                             ! structure index
        ICLOSE(NOINT)=0
        IF(INSTRING(2:2).EQ.'C')ICLOSE(NOINT)=1   ! switch for closed loop
        INTB1(NOINT)=NPT+1                        ! identify first point
      ELSE
        READ(INSTRING,*,ERR=40,END=40)XI2,YI2
        IF(NPT.GE.NBND)THEN
          WRITE(6,*)'Insufficient memory in xpoly, increase NBND'
          STOP
        ENDIF
        NPT=NPT+1
        IF(NPT.EQ.(INTB1(NOINT)+1))READ(HEADER(3:20),*,END=90)IIN(NOINT)
        XBND(NPT)=XI2
        YBND(NPT)=YI2
        INTBL(NOINT)=NPT                        ! last point unless replaced
      ENDIF
C
C   go and read next line
C
      GO TO 30
   40 CONTINUE
      CLOSE(8)
C
C   decrement NOINT if there was a final # without more points
C
      IF(INTB1(NOINT).GT.NPT)NOINT=NOINT-1
C
C    check points in each structure
C
      IF(NOINT.GT.0)THEN
        DO I=1,NOINT
          NIB=INTBL(I)-INTB1(I)+1                    ! points added by this segment
          NTN=NTN+NIB                                ! accummulate total points
          NINTRNL=NINTRNL+NIB+(ICLOSE(I)-1)          ! accummulate total segments
          IF(NIB.LT.(2+ICLOSE(I)))THEN 
            WRITE(10029)I
10029       FORMAT('Too few points in structure ',I3,
     :           ' Check file: internal')
            STOP
          ENDIF
C
C    report information for internal structures
C
          CLOSED='          '
          IF(ICLOSE(I).EQ.1)CLOSED=', (closed)'
          WRITE(6,10019)'In',I,IIN(I),NIB,CLOSED
          KO=INTB1(I)
          DO WHILE (KO.LE.INTBL(I))
            KE=KO+3
            IF(KE.GT.INTBL(I))KE=INTBL(I)
            WRITE(6,1050)(XBND(J),YBND(J),J=KO,KE)
 1050       FORMAT(8F9.3,' &')
            KO=KO+4
          ENDDO
        ENDDO
      ENDIF
C
C     jump to here if no internal file
C
   50 CONTINUE
C
C     if there is a file called connections, then read this also
C     this list of points assumed to begin and end on points already
C     defined in previous files.  Multiple segments may be included
C     sequentially provided each segment begins and ends on a pre-
C     defined point.  NTN is running total of distinct nodes
C
      NCSEG=0
      NSEG=0
      NCON=0
      OPEN(8,FILE='connections',STATUS='OLD',ERR=80)
   60 CONTINUE
      READ(8,'(A80)',END=80)INSTRING(1:80)
      IF((INSTRING(1:1).EQ.'#').OR.(NCON.EQ.0))THEN
      HEADER=INSTRING
C
C    if this condition entered, then operation uses previous set of points read
C
        IF((NSEG+1).GT.MSEG)THEN
          WRITE(6,*)'Insufficient memory in xpoly, increase MSEG'
          STOP
        ENDIF
C
C    sort out details of segment just finished
C
        IF(NCON.NE.0)THEN
          ICB1(NSEG)=NTN+1         ! node index of 1st internal point
          ICE1(NSEG)=NTN+KE-2      ! node index of last internal point
C
C    segments in the connections file should begin and end with points
C    already defined.  confirm first point already present
C
          CALL PREZENT(XCON(1),YCON(1),XBND,YBND,NPT,IPR)
          IF(IPR.EQ.0)THEN
            WRITE(*,10020)'first',XCON(1),YCON(1)
10020       FORMAT(A5,' point of connections (',G12.5,
     :             ',',G12.5,') not previously defined',/,
     :             'xpoly ignores the rest of the connections file')
            STOP
          ELSE
            ICB(NSEG)=IPR     ! node index of first point
          END IF
C
C    confirm last point already present
C
          CALL PREZENT(XCON(KE),YCON(KE),XBND,YBND,NPT,IPR)
          IF(IPR.EQ.0)THEN
            WRITE(*,10020)'last ',XCON(KE),YCON(KE)
            STOP
          ELSE
            ICE(NSEG)=IPR     ! node index of last point
          ENDIF
C
C   points of current segment are in XCON(1:KE), YCON(1:KE)
C   add new points to XBND and YBND
C
   75     CONTINUE
          ICE(NSEG)=IPR      ! node index of last point
          WRITE(*,2020)NSEG,ICN(NSEG),KE-2
 2020     FORMAT('Connection ',I4,' Label ',I5,' adds ',7X,I4,' points')
          IF(KE.GT.2)THEN
             DO K=1,KE-2
                XBND(NTN+K)=XCON(K+1)
                YBND(NTN+K)=YCON(K+1)
             ENDDO
             NPT=NPT+KE-2       ! update NPT used by PREZENT
C
C    write output for a REG statement as used by basil
C
             WRITE(6,1050)XBND(ICB(NSEG)),YBND(ICB(NSEG)),
     :                   (XBND(J),YBND(J),J=NTN+1,NTN+KE-2),
     :                    XBND(ICE(NSEG)),YBND(ICE(NSEG))
             WRITE(6,2021)ICB(NSEG),XBND(ICB(NSEG)),YBND(ICB(NSEG)),
     :                    ICE(NSEG),XBND(ICE(NSEG)),YBND(ICE(NSEG))
 2021        FORMAT(' Connected to points: ',I4,' (',F9.3,',',F9.3,') ',
     :              ' and ',I4,' (',F9.3,',',F9.3,')')
            NTN=NTN+KE-2
            NCSEG=NCSEG+KE-1
          ELSE IF(KE.LT.2)THEN
            ICB1(NSEG)=0            ! key for no internal nodes on segment
            WRITE(6,1011)NSEG
 1011       FORMAT('Insufficient points in segment ',I2,' of ',
     :             'connections')
            STOP
          ELSE                       ! KE = 2
            NCSEG=NCSEG+1
          ENDIF
        ENDIF        ! if NCON.ne.0
C
C    initiate a new segment
C
        KE=0
        NSEG=NSEG+1
        NCON=NCON+1
        ICB(NSEG)=0
C
C    points for a new segment are first written to XCON,YCON(NCT)
C
      ELSE
        KE=KE+1
        IF(KE.GT.1)READ(HEADER(3:20),*,END=90)ICN(NSEG)
        IF(KE.GT.NBND)THEN
          WRITE(6,*)'Insufficient memory in xpoly, increase NBND'
          STOP
        ENDIF
        READ(INSTRING,*,END=80)XCON(KE),YCON(KE)
      ENDIF
C
C     look for more lines to read
C
      GO TO 60
      CLOSE(8)
C
C     skip to here if no connections file; begin writing poly file
C
   80 CONTINUE
      IF(ICB(NSEG).EQ.0)NSEG=NSEG-1
      STR=' 2 0 0'
      OPEN(10,FILE='x.poly')
      WRITE(10,1010)NTN,STR
 1010 FORMAT(I5,A6)
C
C   write set of all saved points
C
      DO K=1,NTN
        WRITE(10,1020)K,XBND(K),YBND(K)
      ENDDO
 1020 FORMAT(I5,X,2F12.5)
C
C   write total segment number
C
      WRITE(10,1030)NBT+NINTRNL+NCSEG
 1030 FORMAT(I5,' 1')
C
C   write external boundary segments, closing back to first point
C
      DO K=1,NBT-1
        IBLAB=1
        DO J=1,NOBND
          KF=INTEL(J)
          IF(KF.EQ.1)KF=NBT
          IF((K.GE.INTEB(J)).AND.(K.LE.KF))IBLAB=IBN(J)
        ENDDO
        WRITE(10,1040)K,K,K+1,IBLAB
      ENDDO
      WRITE(10,1040)NBT,INTEB(NOBND),1,IBN(NOBND)   ! external last segment
 1040 FORMAT(3I4,I3)
C
C   write internal structure segments
C 
C     WRITE(*,*)'writing internal structure, NOINT=',NOINT
      NS=NBT                                ! segment number for output file
      IF(NOINT.GT.0)THEN
        DO J=1,NOINT
          DO K=INTB1(J),INTBL(J)-1
            NS=NS+1
            WRITE(10,1041)NS,K,K+1,IIN(J)
          ENDDO
          IF(ICLOSE(J).EQ.1)THEN
            NLST=INTBL(J)
            NS=NS+1
            WRITE(10,1041)NS,NLST,INTB1(J),IIN(J)
 1041       FORMAT(3I4,I5)
          ENDIF
        ENDDO
      ENDIF
C
C   write additional connecting segments as needed
C
C     WRITE(*,*)'writing connections, NSEG =',NSEG 
      IF(NSEG.GT.0)THEN
        DO J=1,NSEG
          NS=NS+1
          IF(ICB1(J).NE.0)THEN
            WRITE(10,1041)NS,ICB(J),ICB1(J),ICN(J)    ! 1st segment
            KA=ICB1(J)
            DO WHILE (KA.LT.ICE1(J))
              KB=KA+1
              NS=NS+1
              WRITE(10,1041)NS,KA,KB,ICN(J)           ! internal segments
              KA=KA+1
            ENDDO
            NS=NS+1
            WRITE(10,1041)NS,ICE1(J),ICE(J),ICN(J)    ! last segment
          ELSE
            WRITE(10,1041)NS,ICB(J),ICE(J),ICN(J)     ! one segment only
          ENDIF
        ENDDO
      END IF
C
C    finish and close file
C
      WRITE(10,*)'0'
      WRITE(10,*)'0'
      CLOSE(10)
      WRITE(*,*)'A new poly file has been written to x.poly'
C
      STOP
   90 WRITE(*,*)'please refer specification of segment labels ',
     :   ' needed in files internal, connections, following #_'
      STOP
  900 WRITE(*,*)'xpoly requires at least a file called boundary'
      STOP
      END
      SUBROUTINE PREZENT(XX,YY,XBND,YBND,NPT,IPR)
C
C   this subroutine checks if the point (XX,YY) is already present
C   in the list of points in (XBND(NPT), YBND(NPT))
C
      DIMENSION XBND(NPT),YBND(NPT)
      IPR=0
      DO J=1,NPT
        IF(SQRT((XX-XBND(J))*(XX-XBND(J)) + 
     :            (YY-YBND(J))*(YY-YBND(J))).LE.0.0001)THEN
          IPR=J
        ENDIF
      ENDDO
      RETURN
      END
