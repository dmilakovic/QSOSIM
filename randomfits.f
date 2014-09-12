c-------------------------------------------------------------------------------
      PROGRAM writefits
c     PURPOSE: create a fits file with parameters needed for qsosim fits version
c     OUTPUT:  sin.fits file to replace sin.dat
c-------------------------------------------------------------------------------
c GENERAL DECLARATIONS
      IMPLICIT NONE
      EXTERNAL ran3
c Declare variables
      INTEGER :: unit,i,j,status,blocksize,bitpix,naxis,naxes
      INTEGER :: nrows, tfields, varidat, colnum, idum, inoise
      INTEGER, PARAMETER :: nobj=25
      REAL*4 :: ramin, ramax, dra, decmin, decmax, ddec, zmin, zmax, dz
      REAL*4 :: alphamin, alphamax, dalpha 
      REAL*4 :: wstart, wend, dw, nc, nuplim
      REAL :: ran3
      REAL,DIMENSION(nobj) :: ra, dec, zqso, alpha, vmag, dvavoid, 
     +                        sigblur, s2n
      INTEGER,DIMENSION(nobj) :: numlls
      REAL,DIMENSION(nobj) :: nhills1,blls1,zlls1,nhills2,blls2,zlls2
      CHARACTER*20 :: ttype1(7), tunit1(7), tform1(7)
      CHARACTER*20 :: ttype2(14), tunit2(14), tform2(14)
      CHARACTER*20 :: filename, extname 
      CHARACTER*30 :: errtext
      LOGICAL :: simple, extend
c Define parameters
      filename='sin.fits'
      blocksize=1
      status=0
      simple=.true.
      bitpix=16
      naxis=0
      naxes=0
      extend=.true.
c-------------------------------------------------------------------------------
c GENERATE DATA 
c data: ra,dec,zqso,alpha,vmag,wstart,wend,dw,nc,nuplim,sigblur,s2n,inoise,
c       numlls,dvavoid,nhills(i),blls(i),zlls(i)
      ramin=0.0
      ramax=360.0
      dra=ramax-ramin
     
      decmin=-90.0
      decmax=90.0
      ddec=decmax-decmin
     
      zmin=2.5
      zmax=6.0
      dz=zmax-zmin
      
      alphamin=-1
      alphamax=1
      dalpha=alphamax-alphamin

      wstart=3000
      wend=5200
      dw=0.03
      nc=3.0
      nuplim=1.00
      sigblur=3.0
      s2n=100
      inoise=1

      numlls=0
      dvavoid=100
      nhills1=1.0e21
      blls1=10.0
      zlls1=2.86618097
      nhills2=1.0e21
      blls2=10.0
      zlls2=2.86618097
      
c Generate random coordinates      
      idum=time()
      DO i=1,nobj
          ra(i)=ramin+dra*ran3(idum)
          dec(i)=decmin+ddec*ran3(idum)
          zqso(i)=zmin+dz*ran3(idum)
          alpha(i)=alphamin+dalpha*ran3(idum)
      END DO
c-------------------------------------------------------------------------------
      unit=1
c Check if fits file exists and delete if does
      call DELETEFILE(filename,status)
      call ftgerr(status,errtext)
      print *,status,' ',errtext
c Create a fits file
      call FTINIT(UNIT,filename,blocksize,status)
      call ftgerr(status,errtext)
      print *,status,' ',errtext
c Define primary array parameters
      call FTPHPR(UNIT,simple,bitpix,naxis,naxes,0,1,extend,status)
      call ftgerr(status,errtext)
      print *,status,' ',errtext
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c Define data to be inputted into fits file
      DATA ttype1/'wstart','wend','dw','nc','nuplim','inoise','dvavoid'/
      DATA ttype2/'RA','DEC','zqso','alpha','vmag','sigblur',
     &            's2n','numlls','nhills1',
     &            'blls1','zlls1','nhills2','blls2','zlls2'/
      
      DATA tform1/'E','E','E','E','E','J','E'/
      DATA tform2/'E','E','E','E','E','E','E','J','E','E','E','E','E',
     +            'E'/
      DATA tunit1/'','','','e12','e16','',''/
      DATA tunit2/'','','','','','','','','','','','','',''/
c Create the first binary table HDU
      nrows=1
      tfields=7
      varidat=0
      extname='GENERAL'
      call FTIBIN(unit,nrows,tfields,ttype1,tform1,tunit1,
     &            extname,varidat,status)
c Rename header column names and assign values
      call FTMNAM(unit,'ttype1','wstart',status)
      call FTMKYF(unit,'wstart',wstart,2,'Start wavelength',status)
      call FTMNAM(unit,'ttype2','wend',status)
      call FTMKYF(unit,'wend',wend,2,'End wavelength',status)
      call FTMNAM(unit,'ttype3','dw',status)
      call FTMKYF(unit,'dw',dw,2,'Pixel size',status)
      call FTMNAM(unit,'ttype4','nc',status)
      call FTMKYF(unit,'nc',nc,2,'N(HI) lower cut-off',status)
      call FTMNAM(unit,'ttype5','nuplim',status)
      call FTMKYF(unit,'nuplim',nuplim,2,'N(HI) upper cut-off',status)
      call FTMNAM(unit,'ttype6','inoise',status)
      call FTMKYJ(unit,'inoise',inoise,'inoise parameter',status)
      call FTMNAM(unit,'ttype7','dvavoid',status)
      call FTMKYF(unit,'dvavoid',dvavoid,3,'Avoidance zone',status)
      call ftgerr(status,errtext)
      print *,status,' ',errtext
c Create the second binary table HDU with data pertaining each QSO
      nrows=nobj
      tfields=14
      extname='QSO'
      call FTIBIN(unit,nrows,tfields,ttype2,tform2,tunit2,
     &            extname,varidat,status)
      call FTPCLE(unit,1,1,1,nrows,ra,status)
      call FTPCLE(unit,2,1,1,nrows,dec,status)
      call FTPCLE(unit,3,1,1,nrows,zqso,status)
      call FTPCLE(unit,4,1,1,nrows,alpha,status)
      call FTPCLE(unit,5,1,1,nrows,vmag,status)
      call FTPCLE(unit,6,1,1,nrows,sigblur,status)
      call FTPCLE(unit,7,1,1,nrows,s2n,status)
      call FTPCLJ(unit,8,1,1,nrows,numlls,status)
      call FTPCLE(unit,9,1,1,nrows,nhills1,status)
      call FTPCLE(unit,10,1,1,nrows,blls1,status)
      call FTPCLE(unit,11,1,1,nrows,zlls1,status)
      call FTPCLE(unit,12,1,1,nrows,nhills2,status)
      call FTPCLE(unit,13,1,1,nrows,blls2,status)
      call FTPCLE(unit,14,1,1,nrows,zlls2,status)
c Close fits file
      call FTCLOS(unit,status)
      call ftgerr(status,errtext)
      if (status.eq.0)then 
         print *,status,' File closed'
      else 
         print *,status,' ',errtext
      end if
c End program
      END PROGRAM writefits

c-------------------------------------------------------------------------------
c SUBROUTINES
      subroutine deletefile(filename,status)
C  A simple little routine to delete a FITS file
      integer status,unit,blocksize
      character*20 filename
      character*30 errtext
C  Simply return if status is greater than zero
      if (status .gt. 0)return
C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)
C  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)
      if (status .eq. 0)then
C         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if
C  Free the unit number for later reuse
      call ftfiou(unit, status)
      end subroutine deletefile
c-------------------------------------------------------------------------------
c Numercial Recipes subroutine
      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END FUNCTION ran3
