      SUBROUTINE CNTR(ARRAY,MDIM,NDIM,IDIMA,LOCNTR,HICNTR,                      
     X  NCNTR,WKAREA,EXSC,EXLF)
C                                                                               
C                     MODIFIED FROM                                             
C                                                                               
C ********  THE SLAC UNIFIED GRAPHICS SYSTEM (VAX-11 VERSION)  ********         
C *                      CONTOUR PLOT SUBROUTINE                      *         
C *                                                                   *         
C *                                                                   *         
C *  THE CALLING SEQUENCE IS:                                         *         
C*    CALL ARRAY,MDIM,NDIM,LOCCNTR,HICNTR,WKAREA,EXSC,EXLF)                
C *  ARRAY IS A RECTANGULAR FLOATING POINT ARRAY OF DIMENSION MDIM    *         
C *  AND NDIM; MDIM AND NDIM ARE FIXED POINT; LOCNTR AND HICNTR ARE   *         
C *  FLOATING POINT; NCTR IS FIXED POINT; WKAREA IS A FIXED POINT     *         
C *  ARRAY OF APPROXIMATE DIMENSION MDIM*NDIM/7;         *                      
C *                                                                   *         
C *                          ROBERT C. BEACH                          *         
C *                    COMPUTATION RESEARCH GROUP                     *         
C *                STANFORD LINEAR ACCELERATOR CENTER                 *         
C *                                                                   *         
C *********************************************************************         
C                                                                               
      INTEGER       MDIM,NDIM                                                   
      REAL          ARRAY(IDIMA,1)                                              
      REAL          LOCNTR,HICNTR                                               
      INTEGER       NCNTR                                                       
      INTEGER*4     WKAREA(1)                                                   
      REAL*4        EXXR,EXXL,EXXV,EXYV                                         
      INTEGER*2     EXSC,EXLF                                                   
      INTEGER       NCNT,ICNT,MCNT,IROW,ICOL,ISID                               
      REAL          ZCNT                                                        
      REAL          XPNT,YPNT,TPNT,XLBL,YLBL                                    
      INTEGER       JROW,JCOL,JSID                                              
      INTEGER       KROW,KCOL,KSID                                              
      INTEGER       BBIT,BDFG,INFG,MKFG                                         
      REAL          IVR1,IVR2,DVR1,DVR2                                         
      REAL          XSAV,YSAV                                                   
      INTEGER*4     PWRS(31)                                                    
      INTEGER       XMNB,XMNW,XMNO,XMWD,XMLP,XMUP,XMP0,XMP1                     
      INTEGER       IRET,JRET,KRET                                              
C                                                                               
      INTEGER       INT1,INT2                                                   
      SAVE
C                                                                               
      DATA          PWRS/                                                       
     X            1,          2,          4,          8,                        
     X           16,         32,         64,        128,                        
     X          256,        512,       1024,       2048,                        
     X         4096,       8192,      16384,      32768,                        
     X        65536,     131072,     262144,     524288,                        
     X      1048576,    2097152,    4194304,    8388608,                        
     X     16777216,   33554432,   67108864,  134217728,                        
     X    268435456,  536870912, 1073741824/                                    
      EXXR=-0.015                                                               
      EXXL=-0.120                                                               
      EXXV=-0.060                                                               
      EXYV=0.020                                                                
      IF((NDIM.LT.3).OR.(MDIM.LT.3))THEN                                   
      WRITE(*,302) NDIM,MDIM                                                    
 302  FORMAT(1X,'NO CONTOURING: NDIM = ',I3,' MDIM = ',I3)  
      RETURN
      END IF
C  COMPUTE THE CONTOUR COUNT.                                                   
      KDIM = ((MDIM-2)*(NDIM-1)+NDIM+13)/15                                     
      INT1=MAX0(NCNTR,2)                                                        
      NCNT=INT1+(INT1-1)*EXSC                                                   
C                                                                               
C  LOOP FOR EACH CONTOUR.                                                       
C     
      ICNT=1
  112 ZCNT=LOCNTR+FLOAT(ICNT-1)*(HICNTR-LOCNTR)/FLOAT(NCNT-1)
C     WRITE(*,11111)ZCNT
11111 FORMAT(' Contour level=',G12.5)
      MCNT=MOD(ICNT-1,EXSC+1)                                                 
C  GENERATE LABEL IF NECESSARY.                                                 
C  CLEAR SEGMENT BIT MAP.                                                       
        DO 101 INT1=1,KDIM                                                      
          WKAREA(INT1)=0                                                        
  101   CONTINUE                                                                
C  SET FLAG TO INDICATE BOUNDARYS BEING PROCESSED.                              
        BDFG=1                                                                  
C
C  PROCESS LOWER AND UPPER BOUNDARY.                                            
C
          ICOL=3
  104     IROW=3                                                                
          ISID=1                                                                
          IRET=1                                                                
          GO TO 401                                                             
  102     IROW=MDIM                                                             
          ISID=3                                                                
          IRET=2                                                                
          GO TO 401                                                             
  103     CONTINUE                                                              
          ICOL=ICOL+1
	  IF(ICOL.LE.NDIM)GO TO 104
C
C  PROCESS LEFT AND RIGHT BOUNDARY.                                             
C
          IROW=3
  107     ICOL=3                                                                
          ISID=0                                                                
          IRET=3                                                                
          GO TO 401                                                             
  105     ICOL=NDIM                                                             
          ISID=2                                                                
          IRET=4                                                                
          GO TO 401                                                             
  106     CONTINUE                                                              
          IROW=IROW+1
	  IF(IROW.LE.MDIM)GO TO 107
C
C  SET FLAG TO INDICATE BOUNDARYS NOT BEING PROCESSED.                          
        BDFG=0                                                                  
C  PROCESS INTERIOR SIDES OF SURFACE PATCHES.                                   
C
          IF (NDIM.LT.4) GO TO 111                                                
          ICOL=4
  110     CONTINUE
            IROW=3
  109       CONTINUE
            ISID=0                                                              
            IRET=5                                                              
            GO TO 401                                                           
  108       CONTINUE                                                            
            IROW=IROW+1
	    IF(IROW.LE.MDIM)GO TO 109                   
          ICOL=ICOL+1
	  IF(ICOL.LE.NDIM)GO TO 110
  111 CONTINUE                                                                  
      ICNT=ICNT+1
      IF((NCNTR.EQ.1).OR.(ICNT.GT.NCNT))RETURN
      GO TO 112
C                                                                               
C  INTERNAL ROUTINE TO PROCESS THE ISID-TH SIDE OF THE                          
C  (IROW,ICOL)-TH SURFACE PATCH.  IF THE SIDE HAS NOT BEEN                      
C  CHECKED BEFORE, THEN THE CONTOUR IS EXAMINED TO SEE IF IT                    
C  CROSSES THE SIDE.  IF IT DOES NOT, THE SIDE IS MARKED AS                     
C  HAVING BEEN CHECKED.  IF THE CONTOUR CROSSES THE SIDE, THE                   
C  CONTOUR IS FOLLOWED UNTIL IT IS COMPLETE AND ALL AFFECTED                    
C  SIDES ARE MARKED AS HAVING BEEN CHECKED.                                     
  401 JROW=IROW                                                                 
      JCOL=ICOL                                                                 
      JSID=ISID                                                                 
C  DO ANY CONTOURS BEGIN AT THIS SIDE?                                          
      KRET=1                                                                    
      GO TO 851                                                                 
  402 IF (MKFG.EQ.1) GO TO 471                                                  
      JRET=1                                                                    
      GO TO 501                                                                 
  403 IF (INFG.EQ.0) GO TO 463                                                  
C  START DRAWING THE CONTOUR CURVE.                                             
      JRET=1                                                                    
      GO TO 601                                                                 
  404 BBIT=3                                                                    
      JRET=1                                                                    
      GO TO 701                                                                 
  405 BBIT=2                                                                    
C  FIND THE OTHER SIDE OF THE PATCH.                                            
  411 JSID=MOD(JSID+2,4)                                                        
      KRET=2                                                                    
      GO TO 851                                                                 
  412 IF (MKFG.EQ.1) GO TO 414                                                  
      JRET=2                                                                    
      GO TO 501                                                                 
  413 IF (INFG.NE.0) GO TO 421                                                  
      KRET=3                                                                    
      GO TO 801                                                                 
  414 JSID=MOD(JSID+1,4)                                                        
      KRET=4                                                                    
      GO TO 851                                                                 
  415 IF (MKFG.EQ.1) GO TO 417                                                  
      JRET=3                                                                    
      GO TO 501                                                                 
  416 IF (INFG.NE.0) GO TO 421                                                  
      KRET=5                                                                    
      GO TO 801                                                                 
  417 JSID=MOD(JSID+2,4)                                                        
      KRET=6                                                                    
      GO TO 851                                                                 
 418  JRET=4                                                                    
      GO TO 501                                                                 
 419  CONTINUE                                                                  
C  DRAW CURRENT PART OF THE CONTOUR.                                            
  421 JRET=2                                                                    
      GO TO 701                                                                 
C  MARK THE LINE PROCESSED.                                                     
  431 KRET=7                                                                    
      GO TO 801                                                                 
C  FIND THE ADJACENT SURFACE PATCH.                                             
  441 IF (JSID.EQ.3) GO TO 444                                                  
      IF (JSID.EQ.2) GO TO 443                                                  
      IF (JSID.EQ.1) GO TO 442                                                  
      IF (JCOL.LE.3) GO TO 461                                                  
      JCOL=JCOL-1                                                               
      JSID=2                                                                    
      GO TO 451                                                                 
  442 IF (JROW.LE.3) GO TO 461                                                  
      JROW=JROW-1                                                               
      JSID=3                                                                    
      GO TO 451                                                                 
  443 IF (JCOL.GE.NDIM) GO TO 461                                               
      JCOL=JCOL+1                                                               
      JSID=0                                                                    
      GO TO 451                                                                 
  444 IF (JROW.GE.MDIM) GO TO 461                                               
      JROW=JROW+1                                                               
      JSID=1                                                                    
C  CHECK FOR CLOSURE OF THE CONTOUR LINE.                                       
  451 IF ((JROW.EQ.IROW).AND.(JCOL.EQ.ICOL).AND.(JSID.EQ.ISID))                 
     X  GO TO 471                                                               
      GO TO 411                                                                 
C  FINISH AN OPEN CURVE.                                                        
  461 JRET=2                                                                    
      GO TO 601                                                                 
  462 JROW=IROW                                                                 
      JCOL=ICOL                                                                 
      JSID=ISID                                                                 
  463 KRET=8                                                                    
      GO TO 801                                                                 
C  THE CONTOUR LINE IS NOW COMPLETE.                                            
  471 GO TO (102,103,105,106,108),IRET                                          
C                                                                               
C  INTERNAL ROUTINE TO DETERMINE IF THE (JROW,JCOL,JSID)-TH                     
C  SIDE IS INTERSECTED BY THE CURRENT CONTOUR LINE.  IF THE                     
C  ANSWER IS YES, INFG IS SET TO ONE AND THE COORDINATES                        
C  OF THE INTERSECTION ARE SAVED IN (XPNT,YPNT).                                
C  FIRST, OBTAIN THE INDEPENDENT AND DEPENDENT VARIABLES.                       
 501  CONTINUE                                                                  
      IF (MKFG.NE.0) GO TO 510                                                  
      IF (JSID.EQ.3) GO TO 505                                                  
      IF (JSID.EQ.2) GO TO 503                                                  
      IF (JSID.EQ.1) GO TO 502                                                  
      DVR1=ARRAY(JROW-1,JCOL-1)                                                 
      DVR2=ARRAY(JROW,JCOL-1)                                                   
      XPNT=ARRAY(1,JCOL-1)                                                      
      GO TO 504                                                                 
  502 DVR1=ARRAY(JROW-1,JCOL-1)                                                 
      DVR2=ARRAY(JROW-1,JCOL)                                                   
      YPNT=ARRAY(JROW-1,1)                                                      
      GO TO 506                                                                 
  503 DVR1=ARRAY(JROW-1,JCOL)                                                   
      DVR2=ARRAY(JROW,JCOL)                                                     
      XPNT=ARRAY(1,JCOL)                                                        
  504 IVR1=ARRAY(JROW-1,1)                                                      
      IVR2=ARRAY(JROW,1)                                                        
      GO TO 507                                                                 
  505 DVR1=ARRAY(JROW,JCOL-1)                                                   
      DVR2=ARRAY(JROW,JCOL)                                                     
      YPNT=ARRAY(JROW,1)                                                        
  506 IVR1=ARRAY(1,JCOL-1)                                                      
      IVR2=ARRAY(1,JCOL)                                                        
C  CHECK FOR AN INTERSECTION.                                                   
 507  CONTINUE                                                                  
      IF (((DVR1.LT.ZCNT).AND.(DVR2.LT.ZCNT)).OR.                               
     X  ((DVR1.GE.ZCNT).AND.(DVR2.GE.ZCNT))) GO TO 510                          
C  COMPUTE THE OTHER COORDINATE.                                                
      TPNT=IVR1+(ZCNT-DVR1)*(IVR2-IVR1)/(DVR2-DVR1)                             
      IF (MOD(JSID,2).EQ.0) GO TO 508                                           
      XPNT=TPNT                                                                 
      GO TO 509                                                                 
  508 YPNT=TPNT                                                                 
C  RETURN WITH INTERSECTION.                                                    
  509 INFG=1                                                                    
      GO TO 511                                                                 
C  RETURN WITHOUT INTERSECTION.                                                 
  510 INFG=0                                                                    
 511  CONTINUE                                                                  
      GO TO (403,413,416,419),JRET                                              
C                                                                               
C  INTERNAL ROUTINE TO LABEL A CONTOUR IF A LABEL IS                            
C  REQUIRED.  THE LABEL IS POSITIONED AT THE PROPER OFFSET                      
C  FROM (XPNT,YPNT).                                                            
  601 IF ((BDFG.EQ.0).OR.(EXLF.NE.0).OR.(MCNT.NE.0))                            
     X  GO TO 606                                                               
      XLBL=XPNT                                                                 
      YLBL=YPNT                                                                 
C  GET THE COORDINATES OF THE LABEL.                                            
      IF (JSID.EQ.3) GO TO 604                                                  
      IF (JSID.EQ.2) GO TO 603                                                  
      IF (JSID.EQ.1) GO TO 602                                                  
      XLBL=XLBL+EXXL                                                            
      GO TO 605                                                                 
  602 XLBL=XLBL+EXXV                                                            
      YLBL=YLBL-EXYV                                                            
      GO TO 605                                                                 
  603 XLBL=XLBL+EXXR                                                            
      GO TO 605                                                                 
  604 XLBL=XLBL+EXXV                                                            
      YLBL=YLBL+EXYV                                                            
C  PRODUCE THE LABEL.                                                           
 605  SIZE = 1.0                                                                
      ROT = 0.0                                                                 
C*****      STILL TO BE WRITTEN                                                 
  606 GO TO (404,462),JRET                                                      
C                                                                               
C  INTERNAL PROCEDURE TO DRAW A LINE TO THE POINT (XPNT,YPNT).                  
  701 IF (MCNT.NE.0) GO TO 703                                                  
      CALL PLOTU(XPNT,YPNT,BBIT)                                                
      GO TO 704                                                                 
 703  CALL PLOTU(XPNT,YPNT,4)                                                   
 704  XSAV=XPNT                                                                 
      YSAV=YPNT                                                                 
      GO TO (405,431),JRET                                                      
C                                                                               
C  801 IS AN INTERNAL ROUTINE TO MARK THE (JROW,JCOL,JSID)-TH                   
C  LINE AS HAVING BEEN PROCESSED.  851 IS AN INTERNAL ROUTINE                   
C  TO TEST IF THE (JROW,JCOL,JSID)-TH LINE HAS BEEN MARKED.  IF                 
C  IT HAS, MKFG IS SET TO ONE.                                                  
  801 MKFG=0                                                                    
      GO TO 852                                                                 
  851 MKFG=1                                                                    
  852 KROW=JROW                                                                 
      KCOL=JCOL                                                                 
      KSID=JSID                                                                 
      IF (KSID.NE.2) GO TO 853                                                  
      KSID=0                                                                    
      KCOL=KCOL+1                                                               
      GO TO 854                                                                 
  853 IF (KSID.NE.3) GO TO 854                                                  
      KSID=1                                                                    
      KROW=KROW+1                                                               
  854 XMNB=2*((KROW-3)*(NDIM-1)+(KCOL-3))+KSID                                  
      XMNW=1+XMNB/30                                                            
      XMNO=1+MOD(XMNB,30)                                                       
      XMP0=PWRS(XMNO)                                                           
      XMWD=WKAREA(XMNW)                                                         
      IF (MKFG.EQ.0) GO TO 855                                                  
      XMUP=XMWD/XMP0                                                            
      IF (MOD(XMUP,2).EQ.0) MKFG=0                                              
      GO TO 856                                                                 
  855 XMP1=PWRS(XMNO+1)                                                         
      XMLP=MOD(XMWD,XMP0)                                                       
      XMUP=XMWD/XMP1                                                            
      WKAREA(XMNW)=XMUP*XMP1+XMP0+XMLP                                          
  856 GO TO (402,412,414,415,417,418,441,471),KRET                              
C                                                                               
      END                                                                       
