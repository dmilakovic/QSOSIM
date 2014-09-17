c-------------------------------------------------------------------------------
      SUBROUTINE writefits(outfile,ra,dec,zqso,alpha,npts,lambda,flux,
     &                     flerr,nnflux)
c     PURPOSE: create a fits file with parameters returned by qsosim9
c     OUTPUT:  spectrum.fits file with lambda,flux,sigma,no-noise-flux
c-------------------------------------------------------------------------------
c GENERAL DECLARATIONS
      IMPLICIT NONE
c Declare variables
      INTEGER :: unit,i,j,status,blocksize,bitpix,naxis,naxes,npts
      INTEGER :: nrows, tfields, varidat, colnum, idum, inoise
      REAL*8 :: lambda(npts),flux(npts), flerr(npts), nnflux(npts)
      REAL*8 :: ra,dec,zqso,alpha
      CHARACTER*20 :: ttype1(4), tunit1(4), tform1(4)
      CHARACTER*20 :: ttype2(4), tunit2(4), tform2(4)
      CHARACTER*20 :: outfile, extname 
      CHARACTER*30 :: errtext
      LOGICAL :: simple, extend
c      ALLOCATE(lambda(npts),flux(npts),flerr(npts),nnflux(npts))
c Define parameters
      blocksize=1
      status=0
      simple=.true.
      bitpix=16
      naxis=0
      naxes=0
      extend=.true.
c-------------------------------------------------------------------------------
      unit=1
c Check if fits file exists and delete if does
      call DELETEFILE(outfile,status)
c      call ftgerr(status,errtext)
c      print *,status,' ',errtext
c Create a fits file
      call FTINIT(UNIT,outfile,blocksize,status)
      if (status.eq.0)then 
         print *,status,' Output file initialized'
      else 
         print *,status,' ',errtext
      end if
c Define primary array parameters
      call FTPHPR(UNIT,simple,bitpix,naxis,naxes,0,1,extend,status)
c      call ftgerr(status,errtext)
c     print *,status,' ',errtext
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c Define data to be inputted into fits file
      DATA ttype1/'RA','DEC','Z_QSO','ALPHA'/
      DATA ttype2/'LAMBDA','FLUX','IVAR','NNFLUX'/
      
      DATA tform1/'d','d','d','d'/
      DATA tform2/'d','d','d','d'/
      DATA tunit1/'','','',''/
      DATA tunit2/'','','',''/
c Create the first binary table HDU
      nrows=1
      tfields=4
      varidat=0
      extname='GENERAL'
      call FTIBIN(unit,nrows,tfields,ttype1,tform1,tunit1,
     &            extname,varidat,status)
c Rename header column names and assign values
      call FTMNAM(unit,'ttype1','RA',status)
      call FTMKYD(unit,'RA',ra,5,'',status)
      call FTMNAM(unit,'ttype2','DEC',status)
      call FTMKYD(unit,'DEC',dec,5,'',status)
      call FTMNAM(unit,'ttype3','Z_QSO',status)
      call FTMKYD(unit,'ZQSO',zqso,5,'',status)
      call FTMNAM(unit,'ttype4','ALPHA',status)
      call FTMKYD(unit,'ALPHA',alpha,5,'',status)
c      call ftgerr(status,errtext)
c      print *,status,' ',errtext
c Create the second binary table HDU with data pertaining each QSO
      nrows=npts
 250  format(f10.5,2x,f10.5,2x,f10.5,2x,f10.5)
c      write (*,250)flerr(1:10)
      tfields=4
      extname='QSO'
      call FTIBIN(unit,nrows,tfields,ttype2,tform2,tunit2,
     &            extname,varidat,status)
      call FTPCLD(unit,1,1,1,nrows,lambda,status)
      call FTPCLD(unit,2,1,1,nrows,flux,status)
      call FTPCLD(unit,3,1,1,nrows,flerr,status)
      call FTPCLD(unit,4,1,1,nrows,nnflux,status)
c Close fits file
      call FTCLOS(unit,status)
      call ftgerr(status,errtext)
      if (status.eq.0)then 
         print *,status,' File closed'
      else 
         print *,status,' ',errtext
      end if
c End program
      RETURN
      END SUBROUTINE writefits

c-------------------------------------------------------------------------------
c SUBROUTINES
      subroutine deletefile(outfile,status)
C  A simple little routine to delete a FITS file
      integer status,unit,blocksize
      character*20 outfile
      character*30 errtext
C  Simply return if status is greater than zero
      if (status .gt. 0)return
C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)
C  Try to open the file, to see if it exists
      call ftopen(unit,outfile,1,blocksize,status)
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
