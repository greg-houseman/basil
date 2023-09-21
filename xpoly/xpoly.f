C
C   Routine 'xpoly' is used to define a poly file as used in
C   the triangle routine used to construct a 2D finite element
C   mesh comprised of triangles.  xpoly produces an out file 
C   called x.poly.  It looks for 2 input files as follows:

C   "boundary"
C   a set of coordinate pairs, assumed to form a closed polygon 
C   provided in sequence.  Anticlockwise sequence around the
C   boundary is preferred. The final point should be the same
C   location as the first point.  The external boundary can be
C   divided into successive segments labelled with boundary numbers.  
C   To do so begin each segment with # in column 1 followed by the 
C   label (an integer).  The first point in a new segment should be
C   the same as the last point in the preceding segment. These
C   labels may be used to set boundary conditions in a basil run.
C
C   "internal"
C   allows for multiple internal stuctures that may be open or 
C   closed, but are assumed to be non-intersecting and either
C   entriely within or connecting to the external boundary. 
C   Sets of coordinate pairs for distinct structures are separated 
C   by a single line with '#' in column 1.  For a closed segment
C   the last point in the segment should be the same as the first.
C   that line means the structure will be closed.  As in the
C   boundar file each segment may have a label, provided as an
C   integer followung the # on the segment header.
C
C   The segment labels may be used in a basil execution to
C   set boundary conditions or fault conditions.  To be activated
C   fault segments should have an integer label > 2000
C
C   Warning: the polyfile should not contain duplicate points.  The 
C   begining and end of each segment are checked for duplication
C   in order to allow closing a loop or reconnection with a 
C   previously defined structure.  You may get some complicated
C   problems if there is other duplicated locations present.
C
C   revised by G. Houseman, March 2022.
C   
      PARAMETER (NBND=1000, MSEG=100, NSPT=500)
      DIMENSION XBND(2*NBND),YBND(2*NBND)
      DIMENSION NPSEG(MSEG),NBSEG(MSEG)
      DIMENSION INODE(MSEG,NSPT),IBODE(MSEG,NSPT)
      DIMENSION IBN(MSEG),IIN(MSEG),ICLOSE(MSEG)
      CHARACTER(80) INSTRING,HEADER
      CHARACTER(6) STR
      CHARACTER(10) CLOSED
C
C    MSEG is the number of separately labelled segments
C    IBN(MSEG), IIN(MSEG) are the label numbers that
C    get attached to each group of segments, for segments that
C    belong to external boundary, or internal boundary segments
C
      DO J=1,MSEG
        DO K=1,NSPT
          INODE(J,K)=0
          IBODE(J,K)=0
        ENDDO
        NBSEG(J)=0
      ENDDO
C
C    read the sets of points on the external boundary
C
      OPEN(8,FILE='boundary',STATUS='OLD',ERR=900)
      NPT=0
      NOBND=0
C
C    read the header on the next section of the boundary
C
   10 CONTINUE
      READ(8,'(A80)',END=20)INSTRING(1:80)
      IF(INSTRING(1:1).EQ.'#')THEN
        NOBND=NOBND+1                                      ! segment number
        IF(NOBND.GT.MSEG)THEN
          WRITE(6,*)'Insufficient memory in xpoly, increase MSEG'
          STOP
        ENDIF
        READ(INSTRING(2:20),*,ERR=44,END=20)IBN(NOBND)      ! segment label
      ELSE
        READ(INSTRING(1:80),*,ERR=44,END=20)XB2,YB2
        NBSEG(NOBND)=NBSEG(NOBND)+1                         ! number of points in segment
      ENDIF
      GO TO 10
   20 CONTINUE
      IF(NBSEG(NOBND).EQ.0)NOBND=NOBND-1
      REWIND 8
      WRITE(6,10131)'boundary',NOBND
C
C    read the point coordinate pairs on this section of the boundary
C
      DO J=1,NOBND
        READ(8,'(A80)',END=40)INSTRING(1:80)
        IF(NBSEG(J).GT.NSPT)THEN
          WRITE(6,10123)'Ex',J
          STOP
        ENDIF
        IF(NBSEG(J).GE.2)THEN
          DO M=1,NBSEG(J)
            READ(8,*,ERR=40,END=40)XB2,YB2
            IPR=0
C
C     check if point already defined (only check segment ends)
C
            IF(((M.EQ.1).OR.(M.EQ.NBSEG(J))).AND.(NPT.GE.1))
     :        CALL PREZENT(XB2,YB2,XBND,YBND,NPT,IPR)
            IF(IPR.NE.0)THEN
              WRITE(6,10022)M,J,IPR,XB2,YB2
10022         FORMAT('Point ',I4,' of segment ',I4,
     :             ' is node ',I4,' at ',2G12.5)
              IBODE(J,M)=IPR
            ELSE
              IF(NPT.GE.NBND)THEN
                WRITE(6,*)'Insufficient memory in xpoly, increase NBND'
                STOP
              ENDIF
              NPT=NPT+1
              XBND(NPT)=XB2
              YBND(NPT)=YB2
              IBODE(J,M)=NPT
            ENDIF
          ENDDO
        ELSE
          WRITE(10029)IBN(J),'boundary'
10029     FORMAT('Too few points in sub-structure ',I3,
     :           ' Check file: ',A8)
          STOP
        ENDIF
C
C    report information for external structure
C
        CLOSED='          '
        WRITE(6,10019)'Ex',J,IBN(J),NBSEG(J)-1,CLOSED
        WRITE(6,10050)(XBND(IBODE(J,M)),YBND(IBODE(J,M)),M=1,NBSEG(J))
      ENDDO
10050 FORMAT(8F9.3,' &')
10019 FORMAT(A2,'ternal segment ',I5,' label ',I5,' has ',I3,
     :           ' sub-segments',A10)
      CLOSE(8)
      NBT=NPT                          ! total: external points
      WRITE(6,10130)'Ex',NOBND,NBT
10130 FORMAT(A2,'ternal boundaries defined for',I5' sections total ',
     :        I5,' segments = points',/)
      IF(NPT.LT.3)THEN
        WRITE(6,*)'Insufficient points are available'
        STOP
      ENDIF
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
C    NOINT = number of distinct internal structures; 
C    NPSEG(J) is the number of nodes in structure J.
C    Avoid having node locations too close together, or the mesh will
C    contain numerous small triangles.
C
      NOINT=0
      DO J=1,MSEG
        NPSEG(J)=0
        IIN(J)=0
        ICLOSE(J)=0
      ENDDO
      OPEN(8,FILE='internal',STATUS='OLD',ERR=50)
C
C    first read this file to find number of segments and length of segments
C
   30 CONTINUE
      READ(8,'(A80)',END=40)INSTRING(1:80)
      IF(INSTRING(1:1).EQ.'#')THEN
        NOINT=NOINT+1                                      ! segment number
        IF(NOINT.GT.MSEG)THEN
          WRITE(6,*)'Insufficient memory in xpoly, increase MSEG'
          STOP
        ENDIF
        READ(INSTRING(2:20),*,ERR=45,END=40)IIN(NOINT)      ! segment label
      ELSE
        READ(INSTRING(1:80),*,ERR=45,END=40)XI1,YI1
        NPSEG(NOINT)=NPSEG(NOINT)+1                         ! number of points in segment
      ENDIF
      GO TO 30
   40 CONTINUE
      IF(NPSEG(NOINT).EQ.0)NOINT=NOINT-1
      REWIND 8
      WRITE(6,10131)'internal',NOINT
10131 FORMAT(A8,' file read with ',I5,' sections')
C
C    read each section of the internal file
C
      DO J=1,NOINT
        CLOSED='          '
        READ(8,'(A80)',END=40)INSTRING(1:80)
        IF(NPSEG(J).GT.NSPT)THEN
          WRITE(6,10123)'In',J
10123     FORMAT('Too many points in ',A2,'ternal segment ',I5,
     :           ' increase NSPT and recompile xpoly')
          STOP
        ENDIF
        IF(NPSEG(J).GE.2)THEN
          DO M=1,NPSEG(J)
            READ(8,*,ERR=40,END=40)XI2,YI2
            IPR=0
C
C     check if point already defined (only check segment ends)
C
            CALL PREZENT(XI2,YI2,XBND,YBND,NPT,IPR)
            IF(IPR.NE.0)THEN
              IF((M.EQ.1).OR.(M.EQ.NPSEG(J)))THEN
                WRITE(*,10022)M,J,IPR,XI2,YI2
                INODE(J,M)=IPR
              ELSE
                WRITE(*,*)'Connections to existing nodes are only ',
     :          'permitted at ends of segment'
                WRITE(*,10022)M,J,IPR,XI2,YI2
                WRITE(*,*)'Fix internal entries and rerun xpoly'
                STOP
              ENDIF
            ELSE
              IF(NPT.GE.NBND)THEN
                WRITE(6,*)'Insufficient memory in xpoly, increase NBND'
                STOP
              ENDIF
              NPT=NPT+1
              XBND(NPT)=XI2
              YBND(NPT)=YI2
              INODE(J,M)=NPT                         ! last point unless replaced
            ENDIF
            IF(INODE(J,1).EQ.INODE(J,NPSEG(J)))CLOSED=', (closed)'
          ENDDO
        ELSE
          WRITE(10029)IIN(J),'internal'
          STOP
        ENDIF
C
C    report information for internal structure
C
        WRITE(6,10019)'In',J,IIN(J),NPSEG(J),CLOSED
        WRITE(6,10050)(XBND(INODE(J,M)),YBND(INODE(J,M)),M=1,NPSEG(J))
      ENDDO
      CLOSE(8)
      WRITE(6,10129)'In',NOINT,NPT-NBT
10129 FORMAT(A2,'ternal boundaries defined for',I5' sections adding',
     :        I5,' points',/)
C
C     jump to here if no internal file; start writing the poly file
C
   50 CONTINUE
      STR=' 2 0 0'
      OPEN(10,FILE='x.poly')
      WRITE(10,1010)NPT,STR
 1010 FORMAT(I5,A6)
C
C   write set of all saved points
C
      DO K=1,NPT
        WRITE(10,1020)K,XBND(K),YBND(K)
      ENDDO
 1020 FORMAT(I5,X,2F12.5)
C
C   number of segments of closed structure = no. of nodes
C   for open structure it is one less than no. of nodes
C
      NINTRNL=0
      DO J=1,NOINT
        NINTRNL=NINTRNL+NPSEG(J)
      ENDDO
      NINTRNL=NINTRNL-NOINT
C     WRITE(*,*)'writing segments, NBT, NINTRNL =',NBT,NINTRNL
      WRITE(10,1030)NBT+NINTRNL
 1030 FORMAT(I5,' 1')
C
C   write external boundary segments, closing back to first point
C   the label is applied to the segment that follows each point
C
      NS=0
      DO J=1,NOBND
        DO M=2,NBSEG(J)
          NS=NS+1
          WRITE(10,1041)NS,IBODE(J,M-1),IBODE(J,M),IBN(J)
        ENDDO
      ENDDO
C
C   write internal structure segments
C 
C     NS=NBT           ! causing error      ! segment number for output file
      IF(NOINT.GT.0)THEN
        DO J=1,NOINT
          DO M=2,NPSEG(J)
            NS=NS+1
            WRITE(10,1041)NS,INODE(J,M-1),INODE(J,M),IIN(J)
          ENDDO
 1041     FORMAT(3I4,I5)
        ENDDO
      ENDIF
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
   44 WRITE(*,*)'Problem reading boundary file, check contents'
      STOP
   45 WRITE(*,*)'Problem reading internal file, check contents'
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
