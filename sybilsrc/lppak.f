C*--------------------------------------------------------------------
C*    Basil / Sybil:   lppak.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C     Some of these subroutines, originally written by Greg Houseman
C     and Brian Kennett, were modified for postscript output
C     by Kerry Gallagher
C---------------------------------------------------------------------
C    
C     PLOTPAK   - low level plotting routines   
C   
C     Version for MAC II to Laser Printer - Mar. 1990  
C-----------------------------------------------------------------bk
      SUBROUTINE ASPECT(RASP)   
C   
C     sets width to height ratio for characters to RASP 
C   
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      ASP = RASP
      return
      end   
C-----------------------------------------------------------------  
      SUBROUTINE CIRCLE(RADIUS,NSIDES)  
C   
C    draws circle centred at current pen location, with 
C    circumference divided into NSIDES straight segments
C      (NSIDES=0 gives smooth curve) See p.3-12 
C   
      CHARACTER IWORD*40
      SAVE PSCA,TOOPI
      DATA PSCA,TOOPI/28.346457,6.283185307/
C
C    make circle using arc command (convert cm to points)
C
      IF(NSIDES.EQ.0)THEN
      RA=RADIUS*PSCA
      WRITE(IWORD,101)RA
  101 FORMAT('currentpoint ',F6.2,' 0 360 arc ')
      CALL PSEND(IWORD,30,0)
C
C    make circle using straight line segments
C
      ELSE IF(NSIDES.NE.0)THEN
      RA=RADIUS
      DTH=TOOPI/FLOAT(NSIDES)
C
C   RA is nominal radius, RD is distance to vertex from centre
C
      RD=2.0*RA/(1.0+COS(DTH*0.5))
      CALL PLOT(0.0,RD,-3)
      ANGL=TOOPI-DTH*0.5
      XL=2.0*RD*SIN(DTH*0.5)
      DO 60 J=1,NSIDES
      DX=XL*COS(ANGL)
      DY=XL*SIN(ANGL)
      CALL PLOT(DX,DY,-2)
          ANGL=ANGL-DTH
   60 CONTINUE
      CALL PLOT(0.0,-RD,-3)
      END IF
      RETURN
      END   
C---------------------------------------------------------------------- 
      SUBROUTINE CSYMBL(X,Y,IP,SIZE,INT)
C   
C      writes a centered symbol at location (X,Y). The symbol is
C      is selected from the list below by INT for 1<INT<10  
C      and is circle,triangle,square,pentagon,hexagon,heptagon, 
C      octagon for 11<INT<17
C   
      CHARACTER ISYM*1  
      CHARACTER IWORD*40,ICOPY*40
      DIMENSION ISYM(10)
      DIMENSION ICIR(10)
      SAVE PSCA,ISYM,ICIR
      DATA PSCA/28.346457/
      DATA ISYM/'O','X','*','+','#','$','@','8','H','Z'/
      DATA ICIR/20,3,4,5,6,7,8,9,10,11/ 
      IF(INT.GE.11)GO TO 20 
C
C     select character size
C
   40 SZ=ABS(SIZE)*PSCA*1.5
      WRITE(ICOPY,50)SZ
   50 FORMAT(F6.2,' hV')
      CALL PSEND(ICOPY,9,1)
C   
C    turn symbol mode on
C   
C      WRITE(IWORD,55)ISYM(INT) 
C   55 FORMAT('SM',A1,';')   
C      CALL PSEND(IWORD,4,1) 
C   
C     move pen to symbol location, symbol is written after move 
C   
   20 CALL PLOTU(X,Y,IP)
C      IF(INT.LE.10)GO TO 30 
      CALL CIRCLE(SIZE*0.75,ICIR(INT-10))   
C   
C     turn symbol mode off  
C   
C   30 IWORD(1:3)='SM;'  
C      CALL PSEND(IWORD,3,0) 
      RETURN
      END  
C---------------------------------------------------------------bk- 
      SUBROUTINE DASHLN(LDASH,LPAT) 
      CHARACTER IWORD*40
      DIMENSION IDA(4)
C   
C     defines  the style for line drawing   
C      ldash < 0   -  reset to solid line   
C      ldash = 8   -  reset to solid line   
C   
C      ldash = 0   -  dots at calling points
C            = 1   -  dots  
C            = 2   -  half dash 
C            = 3   -  long dash 
C            = 4   -  chain dotted  
C            = 5   -  long and short
C            = 6   -  long and two short
C   
C       lpat - percentage of diagonal of paper used for 
C              a pattern (approx.eq.1000 pts)
C   
C     IF(LPAT.EQ.0 )LPAT=2
C     WRITE(IWORD,100)
C 100 FORMAT('stroke newpath ')
C     CALL PSEND(IWORD,15,1)
C
      IF(LDASH.LE.0.OR.LDASH.GE.8)THEN
        WRITE(IWORD,101)
  101   FORMAT(' stroke [] 0 setdash ')  
        CALL PSEND(IWORD,21,1)
C
      ELSEIF(LDASH.EQ.1)THEN
        IDA(1)=LPAT*12
        IDA(2)=LPAT*4
        WRITE(IWORD,102)(IDA(J),J=1,2)
        CALL PSEND(IWORD,28,1)
C
      ELSEIF(LDASH.EQ.2)THEN
        IDA(1)=LPAT*5
        IDA(2)=IDA(1)
        WRITE(IWORD,102)(IDA(J),J=1,2)
        CALL PSEND(IWORD,28,1)
C   
      ELSEIF(LDASH.EQ.3)THEN
        IDA(1)=LPAT*8
        IDA(2)=LPAT*2
        WRITE(IWORD,102)(IDA(J),J=1,2)
        CALL PSEND(IWORD,28,1)
C
      ELSEIF(LDASH.EQ.7)THEN
        IDA(1)=LPAT*2
        IDA(2)=IDA(1)
        WRITE(IWORD,102)(IDA(J),J=1,2)
  102   FORMAT(' stroke [',I3,' ',I3,'] 0 setdash ')
        CALL PSEND(IWORD,28,1)
C
      ELSEIF(LDASH.GE.4)THEN
        IDA(1)=LPAT*5
        IDA(2)=LPAT*2
        IDA(3)=LPAT*1
        IDA(3)=LPAT*2
        WRITE(IWORD,103)(IDA(J),J=1,4)
  103   FORMAT(' stroke [',I3,' ',I3,' ',I3,' ',I3,'] 0 setdash ')
        CALL PSEND(IWORD,36,1)
C
      END IF
      RETURN
      END
C---------------------------------------------------------------bk- 
      SUBROUTINE ITALIC(THETA)  
C   
C     defines angle of slant for labels in degrees  
C     useful range -60 < theta <60  
C   
      character*40 iword
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
c   
      thet = theta  
      xt = tan(theta/57.2957795)
      WRITE(iword,55) xt  
 55   format('SL',f4.2,';') 
      call psend(iword,7,0) 
C   
      return
      end   
C------------------------------------------------------------------ 
      SUBROUTINE NUMBER(X,Y,SIZE,RN,ANGL,NSF)   
C   
C     writes a number on the plot: if NSF=klm, format is Fkl.m  
C      if NSF=-lm, RN is fixed to an integer and format is Ilm. 
C   
      CHARACTER IFORM*8
      CHARACTER IWORD*80
C   
C    create the format expression in IFORM  
C   
      IF(NSF.LT.0)GO TO 20  
      ITOT=NSF/10   
      IDPL=MOD(NSF,10)  
      WRITE(IFORM,55)ITOT,IDPL 
   55 FORMAT('(F',I2,'.',I1,')')
      WRITE(IWORD,IFORM)RN 
      GO TO 30  
C   
C    for integer format 
C   
   20 ITOT=-NSF 
      WRITE(IFORM,65)ITOT  
   65 FORMAT('(I',I2,')   ')
      IR=IFIX(RN)   
      WRITE(IWORD,IFORM)IR 
C   
C     encode number and send to plotter 
C   
   30 CALL SYMBOL(X,Y,SIZE,IWORD,ANGL,ITOT) 
      RETURN
      END   
C-------------------------------------------------------------le-   
      SUBROUTINE CONVERTUTOXPTS(VAL,PTS)
      INCLUDE "sybilps.parameters"
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET

      VAL_CM = VAL * A
      PTS = VAL_CM*PTSPERCM
      RETURN
      END
C-------------------------------------------------------------le-   
      SUBROUTINE CONVERTUTOYPTS(VAL,PTS)
      INCLUDE "sybilps.parameters"
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET

      VAL_CM = VAL * C
      PTS = VAL_CM*PTSPERCM
      RETURN
      END
C-------------------------------------------------------------bk-   
      SUBROUTINE PEN(IPEN,IVEL)
      CHARACTER IWORD*80
C
C     This routine is not called by Sybil
C     SETPENCOLOR and SETLINEWIDTH perform these functions
C
C     IPEN sets line type and colour for laser printer
C     IPEN=1 is a thin black line
C     IPEN=2 is a thick black line
C     IPEN=3,4,5,6 is a thin line of gray shade (IPEN-2)/5
C
C     IVEL is a switch for the HP plotter pen velocity
C     it has no effect within this library
C
      PIN=FLOAT(IPEN-2)/5.0
      IF(PIN.LT.0.0)PIN=0.0
      IF(PIN.GT.1.0)PIN=1.0
C     WRITE(IWORD,101)PIN
  101 FORMAT('stroke newpath ',F4.2,' setgray ')
C     CALL PSEND(IWORD,28,1)
      IWIDE=1
C     IF(IPEN.GE.2)IWIDE=2
      IF(IPEN.GE.1)IWIDE=IPEN
      WRITE(IWORD,102)IWIDE
  102 FORMAT(I2,' setlinewidth ')
      CALL PSEND(IWORD,16,1)
      RETURN
      END
C--------------------------------------------------------------le-
      SUBROUTINE SETPENCOLOR(IPEN)
      CHARACTER IWORD*80
C   
C     IPEN sets the current colour for laser printer
C
      INCLUDE "sybilps.parameters"
      WRITE(IWORD,101)IPEN
  101 FORMAT(' stroke ',I4,' setcolor ')
      CALL PSEND(IWORD,24,1)
      RETURN
      END   
C--------------------------------------------------------------le-
      SUBROUTINE SETRGBCOLOR(RED,GREEN,BLUE)
      CHARACTER IWORD*80
C   
C   this routine is suspect because of stroke command at beginning
C   of output command seems to conflict with stroke in other contexts
C   use SETRGB instead.  *** this routine to be deleted ***
C
      INCLUDE "sybilps.parameters"
      IF(RED.GT.1.0)RED=1.0
      IF(GREEN.GT.1.0)GREEN=1.0
      IF(BLUE.GT.1.0)BLUE=1.0
      IF(RED.LT.0.0)RED=0.0
      IF(GREEN.LT.0.0)GREEN=0.0
      IF(BLUE.LT.0.0)BLUE=0.0
      WRITE(IWORD,101)RED,GREEN,BLUE
      WRITE(*,101)RED,GREEN,BLUE
  101 FORMAT(' stroke ',F5.3,' ',F5.3,' ',F5.3,' setrgbcolor ')
      CALL PSEND(IWORD,38,1)
      RETURN
      END   
C--------------------------------------------------------------le-
      SUBROUTINE SETLINEWIDTH(WDTH)
      CHARACTER IWORD*80
C   
C     WDTH sets line width (in pts) for laser printer
C
      INCLUDE "sybilps.parameters"
      IF (WDTH.EQ.0.0) THEN
        WIDE = 1.0
      ELSE IF (WDTH.LT.0.0) THEN
        WIDE = (0.0-WDTH)*(PTSPERCM/10.0)
      ELSE
        WIDE = WDTH
      END IF
      WRITE(IWORD,102)WIDE
  102 FORMAT(F7.3,' setlinewidth ')
      CALL PSEND(IWORD,21,1)
      RETURN
      END   
C----------------------------------------------------------------   
      SUBROUTINE PLOT(X,Y,I)
C   
C     Raises (I=3) or lowers (I=2) pen and moves to coordinates 
C      (X,Y) if I>0 or to current position plus (X,Y) if I<0
C   
      INCLUDE "sybilps.parameters"
      CHARACTER IWORD*40
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
C   
C     rotate plot by 90 degrees if necessary
C   
   20 CONTINUE  
      XP=X  
      YP=Y  
C     IF(IROT.NE.0)GO TO 30 
C     YP=X  
C     XP=-Y 
C     IF(I.GT.0)XP=21.0-Y   
C
C    convert cm to points
C
   30 XP=PTSPERCM*XP
      YP=PTSPERCM*YP
C
C     plot
C
      IF(XP.GT.9999.99)XP=9999.99
      IF(YP.GT.9999.99)YP=9999.99
      IF(XP.LT.-999.99)XP=-999.99
      IF(YP.LT.-999.99)YP=-999.99
      IF(I.EQ. 2)WRITE(IWORD,60)XP,YP
      IF(I.EQ. 3)WRITE(IWORD,61)XP,YP
      IF(I.EQ.-2)WRITE(IWORD,62)XP,YP
      IF(I.EQ.-3)WRITE(IWORD,63)XP,YP
      IF(I.EQ. 4)WRITE(IWORD,64)XP,YP
      IF(I.EQ.-4)WRITE(IWORD,65)XP,YP
   60 FORMAT(F7.2,' ',F7.2,' pL ')
   61 FORMAT(F7.2,' ',F7.2,' pM ')
   62 FORMAT(F7.2,' ',F7.2,' rL ')
   63 FORMAT(F7.2,' ',F7.2,' rM ')
   64 FORMAT(F7.2,' ',F7.2,' pF ')
   65 FORMAT(F7.2,' ',F7.2,' rF ')
      CALL PSEND(IWORD,19,1)
      RETURN
      END
C-------------------------------------------------------------le
      SUBROUTINE CLIPREGIONU(X,Y,HGT,WDTH)

      CHARACTER BUF*30
      WRITE(BUF,30)
   30 FORMAT('gsave ')
      CALL PSEND(BUF,6,1)
      WRITE(BUF,40)
   40 FORMAT('newpath ')
      CALL PSEND(BUF,8,1)
      CALL PLOTU(X,Y,3)
      CALL PLOTU(X+WDTH,Y,2)
      CALL PLOTU(X+WDTH,Y+HGT,2)
      CALL PLOTU(X,Y+HGT,2)
      WRITE(BUF,50)
   50 FORMAT('closepath clip newpath ')
      CALL PSEND(BUF,25,1)
C     WRITE(BUF,60)
C  60 FORMAT('grestore ')
C     CALL PSEND(BUF,9,1)
      RETURN
      END
C-------------------------------------------------------------le
      SUBROUTINE UNSETCLIPREGIONU()

      CHARACTER BUF*30
      WRITE(BUF,60)
   60 FORMAT('grestore ')
      CALL PSEND(BUF,9,1)
      RETURN
      END

C-------------------------------------------------------------le
      SUBROUTINE FILEOPEN(NAME,LEN,ID)
      INCLUDE "sybilps.parameters"
      CHARACTER NAME*255
      CHARACTER FNAME*255

      ID = 3
      FNAME(1:ISYB_FILENAME_MAX) = ' '
      FNAME(1:LEN) = NAME(1:LEN)
      OPEN(ID,FILE=FNAME,STATUS='old',IOSTAT=IOS)
      IF(IOS.EQ.0)CLOSE(ID,STATUS='delete')
      OPEN(ID,FILE=FNAME,STATUS='new',IOSTAT=IOS)
      IF(IOS.NE.0) ID = 0

      RETURN
      END
C-----------------------------------------------------------gh--bk  
      SUBROUTINE HPLOTS(ION,IRO,LPL,ILS)
C   
C    Initialises plotter (ION.eq.1):
C                 if ILS .eq. 0  -  mapped from A4 paper
C                 if ILS .eq. 1  -  mapped from A3 paper
C    Terminates plot file (ION.ne.1) 
C   
      INCLUDE "sybilps.parameters"
      CHARACTER IWORD*80
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/PRS001/IPROJ,TM2(2,2),TC2(2)
      SAVE
C   
      IF(ION.NE.1)GO TO 80  
      A=1.0 
      B=0.0 
      C=1.0 
      D=0.0 
      IPROJ=0   
      TM2(1,1)=1.0  
      TM2(2,1)=0.0  
      TM2(1,2)=0.0  
      TM2(2,2)=1.0  
      TC2(1)=0.0
      TC2(2)=0.0
      IL34 = ils
      LPLOT=LPL 
      IROT=IRO  
C
C     Postcript initialisation
C
      WRITE(IWORD,99)
   99 FORMAT('%!PS-Adobe-3.0')
      CALL PSEND(IWORD,14,2)
C     WRITE(IWORD,100)
C 100 FORMAT('initmatrix ')
C     CALL PSEND(IWORD,11,1)
      WRITE(IWORD,101)
  101 FORMAT('/pM {stroke moveto} def /rM {rmoveto} def ')
      CALL PSEND(IWORD,42,2)
      WRITE(IWORD,102)
  102 FORMAT('/pL {lineto} def /rL {rlineto} def 2 setlinejoin ')
      CALL PSEND(IWORD,49,2)
      WRITE(IWORD,103)
  103 FORMAT('/hV {/Helvetica findfont exch scalefont setfont} def ')
      CALL PSEND(IWORD,53,2)
      WRITE(IWORD,104)
  104 FORMAT('/sF {findfont exch scalefont setfont} def ')
      CALL PSEND(IWORD,43,2)
      WRITE(IWORD,105)
  105 FORMAT('/right {stringwidth pop -1 mul } def ')
      CALL PSEND(IWORD,37,2)
      WRITE(IWORD,106)
  106 FORMAT('/xcentre {stringwidth pop -2 div } def ')
      CALL PSEND(IWORD,39,2)
      WRITE(IWORD,201)
  201 FORMAT('/pF {lineto gsave fill grestore} def ')
      CALL PSEND(IWORD,37,2)
      WRITE(IWORD,202)
  202 FORMAT('/rF {rlineto gsave fill grestore} def ')
      CALL PSEND(IWORD,38,2)
      TMP=DEFAULTWIDTH/10.0*PTSPERCM
      WRITE(IWORD,203)TMP
  203 FORMAT('/pgw ',F7.2,' def ')
      CALL PSEND(IWORD,17,2)
      TMP=DEFAULTHEIGHT/10.0*PTSPERCM
      WRITE(IWORD,204)TMP
  204 FORMAT('/pgh ',F7.2,' def ')
      CALL PSEND(IWORD,17,2)
      WRITE(IWORD,109)
  109 FORMAT('newpath ')
      CALL PSEND(IWORD,8,1)
C    allow for 0.25in margin around A4 page
      WRITE(IWORD,110)72/4,72/4
  110 FORMAT(I3,I3,' translate ')
      CALL PSEND(IWORD,17,2)
      IF (IROT.NE.0) THEN
        WRITE(IWORD,111) 
  111   FORMAT('559 0 translate ')
        CALL PSEND(IWORD,16,2)
        WRITE(IWORD,112)
  112   FORMAT('90 rotate')
        CALL PSEND(IWORD,11,2)
      END IF
C
C    select the default pen
C
C     CALL PEN(1,0)
      RETURN
C
C     Clear buffer, finish plot and close plotfile
C
   80 CONTINUE
      IWORD(1:16)='stroke showpage '
      CALL PSEND(IWORD,16,2)  
      IWORD(1:9)='%%Trailer'
      CALL PSEND(IWORD,9,2)  
      IWORD(1:5)='%%EOF'
      CALL PSEND(IWORD,5,2)  
      CLOSE(LPLOT)
      RETURN  
      END
C------------------------------------------------------------------ 
      SUBROUTINE PLOTU(X,Y,II)  
C   
C     scales user coordinates to plotter coordinates
C   
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/PRS001/IPROJ,TM(2,2),TC(2)
      save
C   
C   Next section is only relevant to 3-D projections. See   
C    subroutines PERSPC, PRPLAN 
C   
      XPRO=X
      YPRO=Y
      IF(IPROJ.NE.0)THEN
      XPRO=TM(1,1)*X + TM(2,1)*Y + TC(1)
      YPRO=TM(1,2)*X + TM(2,2)*Y + TC(2)
      END IF
C   
C    Scale coordinates using transformation defined in SCALE
C   
      XP=A*XPRO 
      YP=C*YPRO 
      IF(II.LT.0)GO TO 10   
      XP=XP+B   
      YP=YP+D   
   10 CALL PLOT(XP,YP,II)   
      RETURN
      END   
C------------------------------------------------------------------ 
      SUBROUTINE PLOT3U(X,Y,Z,II)   
C   
C    scales user coordinates to plotter coordinates 
C      projects 3-D coordinates onto 2-D using view point of PERSPC 
C   
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/PRS000/TM(3,2) 
      COMMON/PRS001/IPROJ,TM2(2,2),TC2(2)
      save
C   
      XPRO=X
      YPRO=Y
      IF(IPROJ.NE.0)THEN
      XPRO=TM(1,1)*X + TM(2,1)*Y + TM(3,1)*Z
      YPRO=TM(1,2)*X + TM(2,2)*Y + TM(3,2)*Z
      END IF
C   
C    Scale coordinates using transformation defined in SCALE
C   
      XP=A*XPRO 
      YP=C*YPRO 
      IF(II.LT.0)GO TO 10   
      XP=XP+B   
      YP=YP+D   
   10 CALL PLOT(XP,YP,II)   
      RETURN
      END   
C------------------------------------------------------------------ 
      SUBROUTINE PSEND(ISEND,NC,IFIN)   
C   
C     accumulates a buffer full of plot commands, and when full 
C     (or if IFIN.eq.2) sends them to the plotfile
C     IFIN=0 : the string is added to the buffer as is
C     IFIN>=1: string is added only if room remains on present
C              line, else string starts on next new line
C   
      CHARACTER ISEND*1,IBLANK*1
      CHARACTER IWORD*480
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/BUF1/IWORD
      DIMENSION ISEND(NC)
      save
      DATA INDEX/1/ 
      DATA IBLANK/' '/  
C   
C     INDEX is the pointer to the next empty buffer position
C     copy the received word into the buffer
C   
      IF(IFIN.LT.1)GO TO 50 
C   
C     force a new line if necessary (e.g. for labelling)
C   
      NCR=((INDEX-1)/80+1)*80 - INDEX + 1   
      IF(NCR.GE.NC)GO TO 50 
      DO 17 I=1,NCR 
        IWORD(INDEX:INDEX)=IBLANK   
        INDEX=INDEX+1 
   17 CONTINUE
C   
   50 CONTINUE  
      IF(NC.LE.0)GO TO 60   
      DO 16 I=1,NC  
        IWORD(INDEX:INDEX)=ISEND(I) 
        INDEX=INDEX+1 
   16 CONTINUE
C     WRITE(3,65)(ISEND(I),I=1,NC)  
C  65 FORMAT(' ',40A1)  
C   
C     if there is still room in the buffer return   
C   
   60 IF((INDEX.LE.400).AND.(IFIN.LT.2))RETURN  
C   
C     write out five full lines (more or less if IFIN.eq.2) 
C   
      JN=(INDEX-2)/80 + 1   
      IF(JN.LE.0)RETURN 
      IF((IFIN.LT.2).AND.(JN.GT.5))JN=5 
      DO 18 I=1,JN  
      I1=(I-1)*80+1 
      I2=I*80   
      IF(I2.GE.INDEX)I2=INDEX-1 
      WRITE(LPLOT,55)(IWORD(K:K),K=I1,I2) 
   55 FORMAT(80A1)  
   18 CONTINUE  
C   
C    restore overflow to beginning of buffer
C   
      INDEX=INDEX-JN*80 
      IF(INDEX.LT.1)INDEX=1 
      IN1=INDEX-1   
      IF(IN1.LE.0)RETURN
      DO 20 I=1,IN1 
        IWORD(I:I)=IWORD(400+I:400+I) 
        IWORD(400+I:400+I)=IBLANK   
   20 CONTINUE
      RETURN
      END   
C---------------------------------------------------------gh--bk-   
      SUBROUTINE SCALE(XMIN,XMAX,PX1,PX2,YMIN,YMAX,PY1,PY2) 
C   
C      sets up scale factors used in PLOTU and other routines   
C   
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      common/p00001/a1,a2,b1,b2,c1,c2,d1,d2
      save
      A=(PX2-PX1)/(XMAX-XMIN)   
      B=PX1-A*XMIN  
      C=(PY2-PY1)/(YMAX-YMIN)   
      D=PY1-C*YMIN  
c   
      a1 = xmin 
      a2 = xmax 
      b1 = px1  
      b2 = px2  
      c1 = ymin 
      c2 = ymax 
      d1 = py1  
      d2 = py2  
c   
      RETURN
      END   
C---------------------------------------------------------le-
      SUBROUTINE ORIGIN(XMIN,PX1,YMIN,PY1) 
C   
C      sets up origin used in PLOT
C   
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET

      B = PX1 - A*XMIN
      D = PY1 - C*YMIN

      RETURN
      END   
C-----------------------------------------------------------------  
      SUBROUTINE SYMBOL(X,Y,SIZE,IWORD,ANGL,NCHAR)  
C   
C     writes a Hollerith string on the plot 
C   
      INCLUDE "sybilps.parameters"
      CHARACTER IWORD*80,ICOPY*80
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      NCH=NCHAR
      IF(NCHAR.GT.73)NCH=73
C
C    check for chars needing "escape" (unbalanced parentheses,
C    and backslash)
C
      CALL CHECKSTRING(IWORD,NCH)
C
C     select character size
C
   40 SZ=ABS(SIZE)*PTSPERCM*1.5
      WRITE(ICOPY,50)SZ
   50 FORMAT(F6.2,' hV ')
      CALL PSEND(ICOPY,10,1)
C   
C      move pen to symbol location  
C   
      IP=3  
      IF(SIZE.LT.0.0)IP=-3  
      CALL PLOT(X,Y,IP) 
C
C     select character orientation
C
      AB=0.0
      IF(IROT.EQ.0)AB=90.0
      ANG=ANGL+AB
      WRITE(ICOPY,45)ANG
   45 FORMAT(F7.2,' rotate ')
      CALL PSEND(ICOPY,15,1)
C
C     write character string
C
      ICOPY(1:1)='('
      ICOPY(2:NCH+1)=IWORD(1:NCH)
      ICOPY(NCH+2:NCH+8)=') show '
      CALL PSEND(ICOPY,NCH+8,1)
C
C     reset character orientation if necessary
C
      IF(ANG.NE.0.0) THEN
        ANG=-ANG 
        WRITE(ICOPY,45)ANG
        CALL PSEND(ICOPY,15,1)
      END IF
      RETURN
      END
C-----------------------------------------------------------le-
      SUBROUTINE SETFONT(NAME,NCHAR,ISIZE)  
C   
C     sets the font to the requested font and fontsize
C   
      CHARACTER NAME*80,ICOPY*80

      WRITE(ICOPY,50)ISIZE
   50 FORMAT(I4,' /')
      I = 7
      ICOPY(I:I+NCHAR) = NAME(1:NCHAR)
      I = I+NCHAR+1
      ICOPY(I:I+4)=' sF '
      CALL PSEND(ICOPY,I+4,1)
      WRITE(ICOPY,60)ISIZE
   60 FORMAT('/hgt',I4,' def ')
      CALL PSEND(ICOPY,13,1)

      RETURN
      END
C------------------------------------------------------------------ 
      SUBROUTINE FORMATNUMBER(NUMBER,VAL,NSF)   

      CHARACTER IFORM*8
      CHARACTER NUMBER*12
C   
C    create the format expression in IFORM  
C   
      ITOT=NSF/10   
      IDPL=MOD(NSF,10)  
      WRITE(IFORM,55)ITOT,IDPL 
   55 FORMAT('(F',I2,'.',I1,')')
      WRITE(NUMBER,IFORM)VAL 

      RETURN
      END
C-----------------------------------------------------------le-
      SUBROUTINE DRAWAUTOLABELS(IWORD,LEN,VALUE,IFRMT,PWINDO,IPEN)  

      INCLUDE "sybilps.parameters"
      CHARACTER IWORD*80
      CHARACTER BUF*80
      CHARACTER NUMBER*80
      DIMENSION PWINDO(PWINDO_ENTRIES)
      INTEGER OFFX, OFFY
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET

      IOLDPEN = IPEN
C
C    select the default pen
C
      CALL SETPENCOLOR(1)
      ISIZE = DFLTFNT
      BUF(1:9) = 'Helvetica'
      CALL SETFONT(BUF,9,ISIZE)
      OFFX = J_LEFT
      OFFY = J_TOP
      MODE = 1
      CALL DRAWLABEL(PWINDO(XCMIN),PWINDO(YCMIN),OFFX,OFFY,MODE,
     1                          IWORD,LEN)
      IF (IFRMT.NE.0) THEN
        CALL FORMATNUMBER(NUMBER,VALUE,IFRMT)
        OFFX = J_RIGHT
        LEN = IFRMT/10
        CALL DRAWLABEL(PWINDO(XCMAX),PWINDO(YCMIN),OFFX,OFFY,MODE,
     1                          NUMBER,LEN)
      END IF
      CALL SETPENCOLOR(IOLDPEN)
      RETURN
      END
C-----------------------------------------------------------le-
      SUBROUTINE CHECKSTRING(IWORD,NCH)
C
C    check for chars needing "escape" (unbalanced parentheses,
C    and backslash)
C
      CHARACTER IWORD*80,TMPWORD*80

      ICNT = 0
      DO 10 I=1,NCH
        IF (IWORD(I:I).EQ.'(') ICNTLEFT = ICNTLEFT+1
        IF (IWORD(I:I).EQ.')') ICNTRIGHT = ICNTRIGHT+1
        IF (IWORD(I:I).EQ.'\\') ICNTBACK = ICNTBACK+1
  10  CONTINUE
      IF (ICNTBACK.GT.0.OR.ICNTLEFT.NE.ICNTRIGHT) THEN
        TMPWORD(1:NCH) = IWORD(1:NCH)
        ICNT = 1
        DO 20 I=1,NCH
          IF ((TMPWORD(I:I).EQ.'(').OR.(TMPWORD(I:I).EQ.')').OR.
     :             (TMPWORD(I:I).EQ.'\\')) THEN  
            IF (ICNT.LT.80) IWORD(ICNT:ICNT) = '\\'
            ICNT = ICNT+1
            IF (ICNT.LE.80) IWORD(ICNT:ICNT) = TMPWORD(I:I)
          ELSE
            IF (ICNT.LE.80) IWORD(ICNT:ICNT) = TMPWORD(I:I)
          END IF
          ICNT = ICNT+1
  20    CONTINUE
        NCH = ICNT-1
      END IF
      RETURN
      END
C-----------------------------------------------------------le-
      SUBROUTINE DRAWLABEL(X,Y,OFFX,OFFY,MODE,IWORD,NCHAR)  
C   
C     writes a Hollerith string on the plot 
C   
C     X and Y in user units
C
      INCLUDE "sybilps.parameters"
      CHARACTER IWORD*80,ICOPY*80
      CHARACTER BUF*20
      INTEGER OFFX, OFFY
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      NCH=NCHAR
      IF(NCHAR.GT.73)NCH=73
      DO 10 I=NCHAR,1,-1
        IF (IWORD(I:I).NE.' ') GOTO 15
  10  CONTINUE
  15  NCH=I
C
C    check for chars needing "escape" (unbalanced parentheses,
C    and backslash)
C
      CALL CHECKSTRING(IWORD,NCH)
C   
C      move pen to symbol location  
C   
      IP=PENUP  
      CALL PLOTU(X,Y,IP) 
C
C     write character string
C
      ICOPY(1:1)='('
      ICOPY(2:NCH+1)=IWORD(1:NCH)
      ICOPY(NCH+2:NCH+3)=') '
      INDX=NCH+3
      CALL PSEND(ICOPY,INDX,1)
      INDX=1
      IF (OFFX.EQ.J_RIGHT) THEN
        ICOPY(INDX:INDX+9)='dup right '
        INDX=INDX+10
      ELSE IF (OFFX.EQ.J_CENTRE) THEN
        ICOPY(INDX:INDX+11)='dup xcentre '
        INDX=INDX+12
      ELSE 
        ICOPY(INDX:INDX+2)='0 '
        INDX=INDX+2
      END IF
      IF (OFFY.EQ.J_CENTRE) THEN 
        WRITE(BUF,20)-0.4
      ELSE IF (OFFY.EQ.J_TOP) THEN
        WRITE(BUF,20)-1.0
      ELSE IF (OFFY.EQ.J_BASE) THEN
        WRITE(BUF,20)0.25
      ELSE
        WRITE(BUF,30)0.0
      END IF
   20 FORMAT('hgt ',F4.1,' mul ')
   30 FORMAT(F4.1,'         ')
      ICOPY(INDX:INDX+12) = BUF(1:13)
      INDX=INDX+13
      WRITE(BUF,50)
   50 FORMAT(' rM ')
      ICOPY(INDX:INDX+3) = BUF(1:4)
      INDX = INDX+4
      ICOPY(INDX:INDX+4)='show '
      CALL PSEND(ICOPY,INDX+4,1)
      RETURN
      END
C-----------------------------------------------------------le-
      SUBROUTINE DRAWLABELCM(X,Y,OFFX,OFFY,MODE,IWORD,NCHAR)  
C   
C     writes a Hollerith string on the plot 
C   
C     X and Y in cms
C
      INCLUDE "sybilps.parameters"
      CHARACTER IWORD*80,ICOPY*80
      CHARACTER BUF*20
      INTEGER OFFX, OFFY
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      NCH=NCHAR
      IF(NCHAR.GT.73)NCH=73
C
C    check for chars needing "escape" (unbalanced parentheses,
C    and backslash)
C
      CALL CHECKSTRING(IWORD,NCH)
C   
C      move pen to symbol location  
C   
      IP=PENUP  
      CALL PLOT(X,Y,IP) 
C
C     write character string
C
      ICOPY(1:1)='('
      ICOPY(2:NCH+1)=IWORD(1:NCH)
      ICOPY(NCH+2:NCH+3)=') '
      INDX=NCH+3
      CALL PSEND(ICOPY,INDX,1)
      INDX=1
      IF (OFFX.EQ.J_RIGHT) THEN
        ICOPY(INDX:INDX+9)='dup right '
        INDX=INDX+10
      ELSE IF (OFFX.EQ.J_CENTRE) THEN
        ICOPY(INDX:INDX+11)='dup xcentre '
        INDX=INDX+12
      ELSE 
        ICOPY(INDX:INDX+1)='0 '
        INDX=INDX+2
      END IF
      IF (OFFY.EQ.J_CENTRE) THEN 
        WRITE(BUF,20)-0.5
      ELSE IF (OFFY.EQ.J_TOP) THEN
        WRITE(BUF,20)-1.0
      ELSE
        WRITE(BUF,30)0.0
      END IF
   20 FORMAT('hgt ',F4.1,' mul ')
   30 FORMAT(F4.1,'         ')
      ICOPY(INDX:INDX+12) = BUF(1:13)
      INDX=INDX+13
      WRITE(BUF,50)
   50 FORMAT(' rM ')
      ICOPY(INDX:INDX+3) = BUF(1:4)
      INDX = INDX+4
      ICOPY(INDX:INDX+4)='show '
      CALL PSEND(ICOPY,INDX+4,1)
      RETURN
      END
c-----------------------------------------------------------------le
      SUBROUTINE FILLPOLY(XPTS,YPTS,NUMPTS,ICOLOUR,OUTLINE)   
C   
C     draws the polygon and fills it with ICOLOUR
C     OUTLINE=0  fill only
C     OUTLINE=1  fill + outline polygon
C     OUTLINE=2  outline only
C   
      CHARACTER BUF*30
      INTEGER OUTLINE
      DIMENSION XPTS(NUMPTS),YPTS(NUMPTS)
      
      DO 10 J = 1,NUMPTS
        IF (J.EQ.1) THEN
          CALL PLOTU(XPTS(J),YPTS(J),3)
        ELSE IF (J.EQ.NUMPTS) THEN
          IF (OUTLINE.EQ.0) THEN
            CALL SETCOLOUR(ICOLOUR)
            CALL PLOTU(XPTS(J),YPTS(J),4)
          ELSE
            CALL PLOTU(XPTS(J),YPTS(J),2)
            WRITE(BUF,50)
   50       FORMAT('closepath ')
            CALL PSEND(BUF,10,1)
            IF (OUTLINE.EQ.1) THEN
              WRITE(BUF,60)
   60         FORMAT('gsave ')
              CALL PSEND(BUF,6,1)
              CALL SETCOLOUR(ICOLOUR)
              WRITE(BUF,70)
   70         FORMAT('fill grestore stroke')
              CALL PSEND(BUF,21,1)
            END IF
          END IF
        ELSE
          CALL PLOTU(XPTS(J),YPTS(J),2)
        ENDIF
  10  CONTINUE
      RETURN
      END
c-----------------------------------------------------------------le
      SUBROUTINE DRAWPOLY(XPTS,YPTS,NUMPTS)   
C   
C     draws the polygon and fills it with ICOLOUR
C   
      CHARACTER BUF*30
      DIMENSION XPTS(NUMPTS),YPTS(NUMPTS)
      
      DO 10 J = 1,NUMPTS
        IF (J.EQ.1) THEN
          CALL PLOTU(XPTS(J),YPTS(J),3)
        ELSE IF (J.EQ.NUMPTS) THEN
          CALL PLOTU(XPTS(J),YPTS(J),2)
          WRITE(BUF,50)
   50     FORMAT('closepath stroke')
          CALL PSEND(BUF,16,1)
        ELSE
          CALL PLOTU(XPTS(J),YPTS(J),2)
        ENDIF
  10  CONTINUE
      RETURN
      END
c-----------------------------------------------------------------bk
      SUBROUTINE CREATECOLORMAP(IRED,IGREEN,IBLUE,IRGBV,
     1                           ILIMITS,ISTEPS,ISIZE,NCOLORS)   
C
C     Calculates rgb values and writes them to the ps file
C     rgb values are written as 3 components in the range 0->255
C     The values are written in hex - 2 characters each
C
      DIMENSION ILIMITS(ISTEPS)
      DIMENSION IRED(ISIZE),IBLUE(ISIZE),IGREEN(ISIZE)
      DIMENSION IRGBV(ISTEPS*3)
      CHARACTER BUF*1024
      CHARACTER HEXNUM*7

      MAX = ISIZE + NCOLORS
C        max valid index in the range 0->val
      WRITE(BUF,1000)MAX-1
 1000 FORMAT('[/Indexed /DeviceRGB ',I4,' <')
      CALL PSEND(BUF,27,1)
      INDX = 1
      DO 5 I=1,ISIZE
        WRITE(HEXNUM,1001)IRED(I),IGREEN(I),IBLUE(I)
        DO 25 J=1,6
          IF (HEXNUM(J:J).EQ.' ') HEXNUM(J:J)='0'
   25   CONTINUE
        BUF(INDX:INDX+6) = HEXNUM(1:7)
        INDX = INDX+7
  5   CONTINUE
      CALL PSEND(BUF,INDX-1,1)
      INDX=1
      INTERVALS = ISTEPS-1
      STRIDE = NCOLORS/INTERVALS
      ILIMITS(1) = ISIZE
      DO 10 I=2,ISTEPS
        ILIMITS(I) = (I-1)*STRIDE+ISIZE
 10   CONTINUE
      ILIMITS(ISTEPS) = MAX
      DO 15 I=1,INTERVALS
        NC0 = ILIMITS(I)
        NC1 = ILIMITS(I+1)
        IF (NC0.LT.0) NC0 = 0
        IF (NC1.GT.MAX) NC1 = MAX
        IF (I.EQ.INTERVALS.AND.NC1.NE.MAX) NC1 = MAX
        ICOLR = (I-1)*3+1
        ICOLG = (I-1)*3+2
        ICOLB = (I-1)*3+3
        INCR = (IRGBV(ICOLR+3) - IRGBV(ICOLR))/(NC1-NC0)
        INCG = (IRGBV(ICOLG+3) - IRGBV(ICOLG))/(NC1-NC0)
        INCB = (IRGBV(ICOLB+3) - IRGBV(ICOLB))/(NC1-NC0)
        NXTR = IRGBV(ICOLR)
        NXTG = IRGBV(ICOLG)
        NXTB = IRGBV(ICOLB)
        DO 20 K = NC0,NC1-1
          WRITE(HEXNUM,1001)NXTR,NXTG,NXTB
          DO 30 J=1,6
            IF (HEXNUM(J:J).EQ.' ') HEXNUM(J:J)='0'
   30     CONTINUE
          BUF(INDX:INDX+6) = HEXNUM(1:7)
          INDX = INDX + 7
          IF (INDX.GE.70) THEN
            CALL PSEND(BUF,INDX-1,1)
            INDX=1
          END IF
          NXTR = NXTR+INCR
          NXTG = NXTG+INCG
          NXTB = NXTB+INCB
 20     CONTINUE
 15   CONTINUE
      CALL PSEND(BUF,INDX-1,1)
      WRITE(BUF,1002)
      CALL PSEND(BUF,19,1)
 1001 FORMAT(3Z2,' ')
 1002 FORMAT(' > ] setcolorspace ')
      RETURN
      END
c-----------------------------------------------------------------bk
      SUBROUTINE DRAWTICMARKS(SIZE,FINISH,VAL,START,RINTRVL,
     1                            IJUSTIFY,IHORIZ)   
C
C     draws ticmarks at intervals of RINTRVL, beginning with
C     START. RJUSTIFY controls whether tics are inside the axis
C
      INCLUDE "sybilps.parameters"   ! for J_CENTRE etc.
C
C    amended to allow FINISH < START; previously not working in that case
C
      ROFFSET = 0
      EPS=1.0E-5
      PINT=ABS(RINTRVL)
      TMPSIZE = SIZE   !  don't change the value passed as SIZE
C
      IF (IJUSTIFY.EQ.J_CENTRE) ROFFSET = -SIZE/2
      IF (IJUSTIFY.EQ.J_TOP.OR.IJUSTIFY.EQ.J_RIGHT) TMPSIZE = -SIZE
      NUM = INT((ABS((FINISH-START))*(1.0+EPS)/PINT)) + 1
C
C    draw marks on X axis
C
      IF (IHORIZ.NE.0) THEN
        XVAL = START
        IF(FINISH.LT.START)XVAL=FINISH
        DO 10 I=1,NUM
          YVAL = VAL + ROFFSET
          CALL PLOTU(XVAL,YVAL,3)
          YVAL = VAL + TMPSIZE
          CALL PLOTU(XVAL,YVAL,2)
          XVAL = XVAL + PINT
  10    CONTINUE
C
C    draw marks on Y axis
C
      ELSE
        YVAL = START
        IF(FINISH.LT.START)YVAL=FINISH
        DO 20 I=1,NUM
          XVAL = VAL + ROFFSET
          CALL PLOTU(XVAL,YVAL,3)
          XVAL = VAL + TMPSIZE
          CALL PLOTU(XVAL,YVAL,2)
          YVAL = YVAL + PINT
  20    CONTINUE
      END IF
      RETURN
      END

c-----------------------------------------------------------------le
      SUBROUTINE DRAWCIRCLE(CXMIN,CYMIN,DIAMETER,ICOL,OUTLINE)
      INCLUDE "sybilps.parameters"
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET

      CHARACTER IWORD*80
      RAD=DIAMETER/2
      X=CXMIN+RAD
      Y=CYMIN-RAD
      CALL PLOTU(X,Y,3)
      CALL CONVERTUTOXPTS(RAD,RADCM)
      WRITE(IWORD,60)RADCM,ICOL
   60 FORMAT('currentpoint ',F7.3,' ',I3,' circle ')
      CALL PSEND(IWORD,34,1)
      END

c-----------------------------------------------------------------bk
      SUBROUTINE DRAWRECTANGLEU(XU1,YU1,XU2,YU2)   
C
C     draws a rectangle, given the bottom left and top right
C     co-ordinates in user units
C
      CHARACTER BUF*30
      CALL PLOTU(XU1,YU1,3)
      CALL PLOTU(XU2,YU1,2)
      CALL PLOTU(XU2,YU2,2)
      CALL PLOTU(XU1,YU2,2)
      WRITE(BUF,50)
   50 FORMAT('closepath stroke ')
      CALL PSEND(BUF,17,1)
      RETURN
      END
c-----------------------------------------------------------------bk
      SUBROUTINE FILLTYP(it,spac,ian)   
C   
C     specifies fill type for shading   
C   
C     it    - 3,4  parallel, cross-hatch
C     spac  - line spacing in cm
C     ian   - 0,1,2,3 horizontal, 45 deg, vertical, 135 deg 
C   
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      CHARACTER IWORD*40
C   
      spec = spac*10.0  
      igle = ian*45 
      WRITE(iword,50) it,spec,igle 
   50 format('FT',I2,',',F4.2,',',I3,';')
      call psend(iword,13,0)
C   
      return
      end   
C---------------------------------------------------------------bk- 
      SUBROUTINE LABELBG(X0,Y0,ILABEL,LEN)  
C   
C     Shades rectangle defined by bottom-left corner at (X0,Y0) and
C     size determined by coordinate increments xinc,yinc   
C   
      CHARACTER BUF*80,IWORD*80,ILABEL*80
C
      IWORD(1:1) = '('
      IWORD(2:LEN+1) = ILABEL(1:LEN)
      IWORD(LEN+2:LEN+2) = ')'
      WRITE(BUF,15)
   15 FORMAT('/wdth ')  
      NUMCHAR = 6
      BUF(NUMCHAR+1:NUMCHAR+1+LEN+2) = IWORD(1:LEN+2)
      NUMCHAR = NUMCHAR+LEN+2
      BUF(NUMCHAR+1:NUMCHAR+20) = ' stringwidth pop def'
      NUMCHAR = NUMCHAR+20
      CALL PSEND(BUF,NUMCHAR,1)  
      CALL PLOT(X0,Y0,3)
      WRITE(BUF,35)
   35 FORMAT(' wdth 0 rL')  
      NUMCHAR = 10
      BUF(NUMCHAR+1:NUMCHAR+9) = ' 0 hgt rL'
      NUMCHAR = NUMCHAR+9
      BUF(NUMCHAR+1:NUMCHAR+17) = ' wdth -1 mul 0 rL'
      NUMCHAR = NUMCHAR+17
      BUF(NUMCHAR+1:NUMCHAR+16) = ' closepath fill'
      NUMCHAR = NUMCHAR+16
      CALL PSEND(BUF,NUMCHAR,1)  
C   
      RETURN
      END
C---------------------------------------------------------------bk- 
      SUBROUTINE SHADRT(X0,Y0,XI,YI)  
C   
C     Shades rectangle defined by bottom-left corner at (X0,Y0) and
C     size determined by coordinate increments xinc,yinc   
C   
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET  
      CHARACTER IWORD*40
C
      CALL PLOT(X0,Y0,3)
      CALL PLOT(X0+XI,Y0,2)
      CALL PLOT(X0+XI,Y0+YI,2)
      CALL PLOT(X0,Y0+YI,2)
      WRITE(IWORD,55)
   55 FORMAT('closepath fill ')  
      CALL PSEND(IWORD,15,1)  
C   
      RETURN
      END
C---------------------------------------------------------------bk- 
      SUBROUTINE EDGERT(XI,YI)  
C   
C     Edges rectangle defined by coordinate increments xinc,yinc
C   
        COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET  
        character*40 iword  
C   
        xin = 10.0*a*xi 
        yin = 10.0*c*yi 
        WRITE(iword,55) xin,yin   
   55   format('ER',F7.2,',',F7.2,';')  
        call psend(iword,18,0)  
C   
        return  
        end 
      SUBROUTINE PERSPC(X,Y,Z,DELT) 
C   
C    Sets up parameters for projection transformation : 
C     3-D to 2-D perspective.   
C     (X,Y,Z) are the coordinates of the view point.
C     Only the angles relative to the 3 axes are relevant,  
C     Parallel light projection is assumed.  The resulting  
C     2-D plot is shown with the Z axis oriented at angle   
C     DELT (in degrees anticlockwise) from the vertical.
C     The orientation of the X and Y axes is determined by  
C     the view point coordinates.   
C   
      COMMON/PRS000/TM(3,2) 
      COMMON/PRS001/IPROJ,TM2(2,2),TC2(2)
      save
      DATA PI/3.141592654/  
C   
      IPROJ=1   
      ELL=SQRT(X*X + Y*Y + Z*Z) 
      DEL=DELT*PI/180.0 
C   
C    Direction cosines of viewpoint 
C   
      COSALP=X/ELL  
      COSBET=Y/ELL  
      COSGAM=Z/ELL  
      ALP=ACOS(COSALP)  
      BET=ACOS(COSBET)  
      GAM=ACOS(COSGAM)  
C   
C    Theta = angle -Y^O^X, Phi = angle Z^O^-Y, Psi = angle X^O^-Z   
C   
C     SECTHE=TAN(ALP)*TAN(BET)  
      SECPHI=TAN(BET)*TAN(GAM)  
      SECPSI=TAN(ALP)*TAN(GAM)  
C     THETA=ACOS(1.0/SECTHE)
      PHI=ACOS(1.0/SECPHI)  
      PSI=ACOS(1.0/SECPSI)  
C   
C    Calculate the elements of the 3-D transformation matrix
C   
      TM(1,1)=-SIN(ALP)*SIN(PSI-DEL)
      TM(2,1)= SIN(BET)*SIN(PHI+DEL)
      TM(3,1)=-SIN(GAM)*SIN(DEL)
      TM(1,2)=-SIN(ALP)*COS(PSI-DEL)
      TM(2,2)=-SIN(BET)*COS(PHI+DEL)
      TM(3,2)= SIN(GAM)*COS(DEL)
C   
      RETURN
      END   
      SUBROUTINE PRPLAN(IXYZ,XYZCD) 
C   
C    Specifies a projection plane for subsequent PLOTU commands.
C    IXYZ is either 'X', 'Y' or 'Z' and XYZCD is the coordinate 
C    on the appropriate axis. Coordinates from PLOTU are then   
C    interpreted as lying in the plane IXYZ=XYZCD and projected 
C    using the perspective previously defined in PERSPC 
C   
      CHARACTER IXYZ*1
      COMMON/PRS000/TM(3,2) 
      COMMON/PRS001/IPROJ,TM2(2,2),TC2(2)
      save
      IF(IXYZ.EQ.'X')THEN   
      I1=2  
      I2=1  
      I3=1  
      ELSE IF(IXYZ.EQ.'Y')THEN  
      I1=1  
      I2=2  
      I3=2  
      ELSE IF(IXYZ.EQ.'Z')THEN  
      I1=1  
      I2=1  
      I3=3  
      ELSE  
      RETURN
      END IF
      DO 18 J=1,2   
        DO 16 I=1,2   
          K=I1 + (I-1)*I2   
          TM2(I,J)=TM(K,J)  
   16   CONTINUE
        TC2(J)=XYZCD*TM(I3,J) 
   18 CONTINUE
      RETURN
      END   
C-----------------------------------------------------------tb-
      SUBROUTINE SETRGB(RED,GREEN,BLUE)  
C   
C     sets the RGB color for the folowing plot commands
C     colour components should all be between 0.0 and 1.0
C   
      CHARACTER IWORD*80
      IF(RED.GT.1.0)RED=1.0
      IF(GREEN.GT.1.0)GREEN=1.0
      IF(BLUE.GT.1.0)BLUE=1.0
      IF(RED.LT.0.0)RED=0.0
      IF(GREEN.LT.0.0)GREEN=0.0
      IF(BLUE.LT.0.0)BLUE=0.0
      WRITE(IWORD,101)RED,GREEN,BLUE
  101 FORMAT(' stroke ',F5.2,' ',F5.2,' ',F5.2,' setrgbcolor ')
      CALL PSEND(IWORD,39,1)
C 101 FORMAT(' ',F5.3,' ',F5.3,' ',F5.3,' setrgbcolor ')
C     CALL PSEND(IWORD,31,1)
C
      RETURN
      END
C-----------------------------------------------------------tb-
      SUBROUTINE SETHSB(X,Y,Z)  
C   
C     sets the HSB color for the folowing plot commands
C       X,Y,Z should all be between 0.0 and 1.0
C   
      CHARACTER IWORD*80
      WRITE(IWORD,101)X,Y,Z
  101 FORMAT('stroke ',F5.3,' ',F5.3,' ',F5.3,' sethsbcolor')
      CALL PSEND(IWORD,36,2)
C
      RETURN
      END
C-----------------------------------------------------------tb-
      SUBROUTINE SETCOLOUR(INDX)
C   
C     sets the colour to the values at the INDX position in the
C     array 
C   
      CHARACTER IWORD*80
      WRITE(IWORD,101) INDX
  101 FORMAT(I4,' setcolor ')
      CALL PSEND(IWORD,13,2)
C
      RETURN
      END
C-----------------------------------------------------------tb-
      SUBROUTINE SETGRY(X)  
C   
C     sets the gray scale for the folowing plot commands
C       X should be between 0.0 and 1.0
C       X=1.0 is white and X=0.0 is black
C   
      CHARACTER IWORD*80
      WRITE(IWORD,101)X
  101 FORMAT(F5.3,' setgray')
      CALL PSEND(IWORD,13,2)
C
      RETURN
      END
C-----------------------------------------------------------tb-
      SUBROUTINE SETBLK
C   
C     sets the gray scale to solid black for the folowing plot commands
C   
      CHARACTER IWORD*80
      WRITE(IWORD,101)
  101 FORMAT('stroke 0.0 setgray')
      CALL PSEND(IWORD,18,2)
C
      RETURN
      END
C-----------------------------------------------------------tb-
      SUBROUTINE SETWHT
C   
C     sets the gray scale to solid white for the folowing plot commands
C   
      CHARACTER IWORD*80
      WRITE(IWORD,101)
  101 FORMAT('stroke 1.0 setgray')
      CALL PSEND(IWORD,18,2)
C
      RETURN
      END
C-----------------------------------------------------------tb-

C    dummy routines for consistency with xpak.c and menus.c
C
      SUBROUTINE SYBFLUSH
      END

      SUBROUTINE COLMAP
      END

C     SUBROUTINE C3CODE(AMAT,IWORK,N2,M2,NV,NH,IC0,IC2,DELX,DELY,
C    1 TMIN,TMAX,ZZSCL,XLEV,YLEV0,YHIT,
C    2 XCMIN,XCMAX,YCMIN,YCMAX,BXWD)
C     END

C     SUBROUTINE LABELCOLOURBAR(CMIN,CMAX,SCALE,BARHGT,XLEV,YLEV)
C     END

C     SUBROUTINE MARKCOLOURBAR(ZCMIN,ZCMAX,NUMCNTRS,STEP,
C    1                                         XLEV,YLEV0,YHIT)
C     END
