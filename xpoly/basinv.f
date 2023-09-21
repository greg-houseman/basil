      PROGRAM BASINV
C
C     program to find the best fit set of parameters in a basil
C     calculation, using the criterion of minimum RMS misfit 
C     relative to a set of GPS vectors that are held in selvect.out
C
C    written by Greg Houseman, thanks to Richard Rigby for advice on
C    how to interface fortran and shell script.
C
C     the program uses routines amoeba and amotry provided in Numerical 
C     Recipes library.  In principle this method determines an unknown
C     vector P of NDIM parameters.  The maximum value of NDIM (NMAX) is 
C     primarily restricted by increasingly slow execution. NMAX is here 
C     set to 10, but experimentation with greater values may be justified.
C
C     amoeba calls the routine FUNK which calls the BASRUN shell script
C     in order to do the forward calculation
C
C     to compile: gfortran -o basinv basinv.f
C
C     input parameters
C     PAX(MP,NP) initial set of trial vectors in vicinity of solution
C     Y(MP) vector of misfits for MP trial vectors
C     NDIM - dimension of the vector to be determined
C     NMAX - (maximum permissible NDIM - used in DIMENSION statements)
C     MP=NDIM+1 - number of trial vectors at each iteration
C     FTOL - numerical precision (fractional decrease in misfit
C            at last iteration)
C     FUNK - external subroutine called to evaluate function
C
C     to run: basinv NDIM PAR1 PAR2 ... PAR-NDIM
C
C     on return:
C     PAX(MP) contains the MP best estimates of the solution vector
C     Y(MP) contains the corresponding set of misfits
C     ITER - holds iteration count on exit
C
      PARAMETER (NMAX=10)        !  if you change this, change also in FUNK
      COMMON/XP/KEY,PLIM
      CHARACTER(20) PARK
      CHARACTER(120) OUTLINE
      DIMENSION Y(NMAX+1),PAX(NMAX+1,NMAX),X(NMAX),PAV(NMAX)
      DIMENSION KEY(NMAX),PLIM(2,NMAX)
      PI=3.14159265359
      FTOL=0.0005
      DO K=1,NMAX
        KEY(K)=0
        PLIM(1,K)=0.
        PLIM(2,K)=0.
      ENDDO
C
C    at the end of the calculation inv_record contains a list of all
C    the forward calculations done during the inversion.  This file
C    is initialised here, as is .itnumber, which keeps track of the
C    number of iterations
C
      CALL SYSTEM('/bin/rm inv_record .itnumber')
C
C     read initial parameter estimates from command line.  Up to
C     NMAX parameters will be processed.  In general the search is
C     in the linear space (unrestricted); if parameter is preceded
C     by "L" search in in the logarithmic space.  If parameter is
C     preceded by "T" search is constrained in the tangential space
C     between the values of the next two entries - starting with
C     average of two entries.
C
      NDIM=0
      KP=0
      DO K=1,NMAX
        KP=KP+1
        CALL GETARG(KP,PARK)
        IF(PARK.NE."")THEN
C
C     for search in the logarithm space
C
          IF(PARK.EQ.'L')THEN
            KP=KP+1
            CALL GETARG(KP,PARK)
            NDIM=NDIM+1
            READ(PARK,*)PAX(1,NDIM)
            KEY(NDIM)=1
            IF(PAX(1,NDIM).GT.0.0)THEN
              PAX(1,NDIM)=LOG(PAX(1,NDIM))
            ELSE
              WRITE(6,99)NDIM,PAX(1,NDIM)
              STOP
            ENDIF
C
C     for search in the tangent space 
C
          ELSEIF(PARK.EQ.'T')THEN
            KP=KP+1
            CALL GETARG(KP,PARK)
            READ(PARK,*)PX1
            KP=KP+1
            CALL GETARG(KP,PARK)
            READ(PARK,*)PX2
            PAT=0.5*(PX1+PX2)
            NDIM=NDIM+1
            KEY(NDIM)=2
            IF(ABS(PX2-PX1).GT.1.E-5*PAT)THEN
              PLIM(1,NDIM)=PI/(PX2-PX1)
              PLIM(2,NDIM)=0.5*PI-PLIM(1,NDIM)*PX2
              PAX(1,NDIM)=TAN(PLIM(1,NDIM)*PAT+PLIM(2,NDIM))
            ELSE 
              WRITE(6,98)NDIM,PX1,PX2
              STOP
            ENDIF
C
C     for search in the linear space
C
          ELSE
            NDIM=NDIM+1
            READ(PARK,*)PAX(1,NDIM)
          ENDIF
        ENDIF
      ENDDO
      MP=NDIM+1
   98 FORMAT('parameter ',I2,' : limits ',2G12.5,' should differ')
   99 FORMAT('parameter ',I2,' = ',G12.5,' should be > 0')
      WRITE(6,100)NDIM
  100 FORMAT('basinv inversion proceeds with ',I3,' parameters')
C
C
C     define NDIM+1 points in NDIM space, perturbed in different directions
C     perturbation = 5% of initial value, or 0.05 if initial value is zero
C
      DO I=1,MP
        DO J=1,NDIM
          PAX(I,J)=PAX(1,J)
        ENDDO
      ENDDO
      DELT=1.05
      DO J=2,MP
        PAX(J,J-1)=PAX(J,J-1)*DELT
        IF(PAX(J,J-1).EQ.0.0)PAX(J,J-1)=0.05
      ENDDO
C
C    to begin, we evaluate Y for each of the MP vectors
C
      WRITE(*,10001)MP
10001 FORMAT('Initial choice of ',I2,' vectors:')
      DO K=1,MP
        DO J=1,NDIM
          X(J)=PAX(K,J)
        ENDDO
        Y(K)=FUNK(X,NDIM)
      ENDDO
      WRITE(*,10006)
10006 FORMAT('Now commencing iterations')
C
C    carry out the inversion using the Numerical Recipes algorithm
C    described under section 10.4 of 2nd edition
C
      CALL AMOEBA(PAX,Y,NMAX+1,NMAX,NDIM,FTOL,FUNK,ITER)
C
C    print result, current MP estimates, then average of MP estimates
C
      DO I=1,NDIM
        PAV(I)=0.0
      ENDDO
      YAV=0.0
      N1=NDIM*10+10
      N2=N1+20
      WRITE(*,10003)ITER+MP
10003 FORMAT('After ',I5,' iterations, best fit vectors are:')
10004 FORMAT('For P = ',10F10.5)
10008 FORMAT(' Misfit:',F12.6)
      DO K=1,MP
        DO J=1,3
          IF(KEY(J).EQ.1)THEN
            PAX(K,J)=EXP(PAX(K,J))
          ELSEIF(KEY(J).EQ.2)THEN
            PAX(K,J)=(ATAN(PAX(K,J))-PLIM(2,J))/PLIM(1,J)
          ENDIF
        ENDDO
        WRITE(OUTLINE,10004)(PAX(K,J),J=1,NDIM)
        WRITE(OUTLINE(N1:120),10008)Y(K)
        WRITE(*,*)OUTLINE(1:N2)
        DO I=1,NDIM
          PAV(I)=PAV(I)+PAX(K,I)
        ENDDO
        YAV=YAV+Y(K)
      ENDDO
      DO I=1,NDIM
        PAV(I)=PAV(I)/FLOAT(MP)
      ENDDO
      YAV=YAV/FLOAT(MP)
      WRITE(OUTLINE,10004)(PAV(J),J=1,NDIM)
      WRITE(OUTLINE(N1:120),10008)YAV
      WRITE(*,10005)
      WRITE(*,*)OUTLINE(1:N2)
10005 FORMAT('Centroid best-fit solution is:')
      CALL SYSTEM('/bin/rm .itnumber')
C
      STOP
      END
C
      FUNCTION FUNK(X,NDIM)
C
C     function that runs basil and returns Misfit relative to GPS data
C     called by amoeba and amotry
C
      COMMON/XP/KEY,PLIM
      PARAMETER (IBUF=10,NMAX=10)
      DIMENSION X(NDIM),XAVE(NDIM),XDEV(NDIM),XF(NDIM)
      DIMENSION XSAVE(NMAX,IBUF),KEY(NMAX),PLIM(2,NMAX)
      CHARACTER(128) SCRIPTIN
      DATA ICOUNT/0/
      SAVE ICOUNT,XSAVE
C
C    recent sets of X parameters are saved in XSAVE, which
C    works on a circular buffer of length 10.  Purpose is
C    to avoid situation where iterations loop continuously
C    with no significant changes to parameters
C
      ICN=MOD(ICOUNT,IBUF)+1
      ICOUNT=ICOUNT+1
      DO K=1,NDIM
        IF(KEY(K).EQ.1)THEN
          XF(K)=EXP(X(K))
        ELSEIF(KEY(K).EQ.2)THEN
          XF(K)=(ATAN(X(K))-PLIM(2,K))/PLIM(1,K)
        ELSE
          XF(K)=X(K)
        ENDIF
C       WRITE(*,*)K,KEY(K),X(K),XF(K)
        XSAVE(K,ICN)=XF(K)
      ENDDO
C
C     for each parameter find max deviation from average of
C     IBUF entries
C
      IF(ICOUNT.GT.IBUF)THEN
        DO K=1,NDIM
          XAVE(K)=0.0
          DO J=1,IBUF
            XAVE(K)=XAVE(K)+XSAVE(K,J)
          ENDDO
          XAVE(K)=XAVE(K)/FLOAT(IBUF)
        ENDDO
        DO K=1,NDIM
          XDEV(K)=0.0
          DO J=1,IBUF
            XDEV1=ABS(XSAVE(K,J)-XAVE(K))
            IF(XDEV1.GT.XDEV(K))XDEV(K)=XDEV1
          ENDDO
          IF(XAVE(K).NE.0.)XDEV(K)=XDEV(K)/ABS(XAVE(K))
        ENDDO
C
C     if none of the parameters are changing significantly, STOP
C
        CMAX=0.0
        DO K=1,NDIM
          IF(XDEV(K).GT.CMAX)CMAX=XDEV(K)
        ENDDO
        IF(CMAX.LT.1.E-4)THEN
          WRITE(*,*)'No significant parameter variation: basinv is ',
     :              'stopping with ICOUNT =',ICOUNT
          STOP
        ENDIF
      ENDIF
C
C    execute the basil / sybil / mdcomp routines
C
C      WRITE(SCRIPTIN,10001)(X(J),J=1,3)
C10001 FORMAT('BASRUN ',3F13.5)
C     WRITE(*,*)'starting iteration ',ICOUNT
      WRITE(SCRIPTIN,*)'./BASRUN ',(XF(J),J=1,NDIM)
      ISTAT=SYSTEM(SCRIPTIN)
      IF(ISTAT.NE.0)THEN
        WRITE(*,*)'BASRUN has returned non-zero code: ',ISTAT,
     :            ' after ',ICOUNT,' iterations'
        STOP
      ENDIF
C
C    read the misfit saved in currentMF by BASRUN
C
      OPEN(10, FILE='currentMFSF', STATUS='OLD', ACTION='READ')
      READ(10,*,ERR=50)FUNK,SFAC
      CLOSE(10)
      IF((SFAC.LT.0).OR.(SFAC.GT.1.e5))THEN
        WRITE(*,*)'Scale-factor ',SFAC,' is outside acceptable range, ',
     :            'stopping.'
        STOP
      ENDIF
      RETURN
C
   50 CONTINUE
      WRITE(*,*)'Problem in reading Misfit or Scale-factor, stopping'
      STOP
      END
C
C   the following FORTRAN routines (amoeba and amotry) are obtained from
C   The NUMERICAL RECIPES library described by Press et al. They provide 
C   the implementation of the Downhill Simplex method described in 
C   Section 10.4 of Numerical Recipes (1992 edition).
C   routines amended to pass NDIM to FUNK, and remove PAUSE statement
C
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      REAL ftol,p(mp,np),y(mp),funk
      PARAMETER (NMAX=20,ITMAX=1000)
      EXTERNAL funk
CU    USES amotry,funk
      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      iter=0
1     do 12 n=1,ndim
        sum=0.
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
C
C    amendment to remove PAUSE command
C
      IF (ITER.ge.ITMAX)THEN
        WRITE(*,*)'ITMAX exceeded in amoeba'
        STOP
      ENDIF
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum,ndim)
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END
C
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.
C
      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      INTEGER ihi,mp,ndim,np,NMAX
      REAL amotry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk
CU    USES funk
      INTEGER j
      REAL fac1,fac2,ytry,ptry(NMAX)
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry,ndim)
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.
