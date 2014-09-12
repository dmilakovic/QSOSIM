c-------------------------------------------------------------------------------
      PROGRAM dat2fits
c     PURPOSE: take spec.dat file created by QSOSIM9 and convert it to spec.fits
c     INPUT: spec.dat (lambda,flux,sigma,no-noise-flux)
c            sin.dat (redshift)
c-------------------------------------------------------------------------------
      IMPLICIT NONE
c Declare variables
      INTEGER :: unit,i,j,ios,nr,status,blocksize,bitpix,naxis,naxes
      INTEGER :: nrowsll, tfields, nrows_z, tfields_z, varidat, colnum
      INTEGER, PARAMETER :: maxrecs=134000
      REAL :: z
      REAL,DIMENSION(*),ALLOCATABLE ::  wl(:), fl(:), flerr(:), nnfl(:)
      CHARACTER :: junk(1), extname(10)
      CHARACTER*20 :: ttype_z(1), tunit_z(1), tform_z(1)
      CHARACTER*20 :: ttype(4), tunit(4), tform(4)
      CHARACTER*20 :: filename
      LOGICAL :: simple, extend
c Define parameters
      nr=0
      filename='spec.fits'
      blocksize=1
      status=0
      simple=.true.
      bitpix=16
      naxis=0
      naxes=0
      extend=.true.
c-------------------------------------------------------------------------------
c Read sin.dat
      open(UNIT=1,FILE='sin.dat')
      read(1,*)z
      print *,z
c Read spec.dat
      OPEN(UNIT=2,FILE='spec.dat') 
      DO j=1,maxrecs 
         READ(2,*,IOSTAT=ios) junk 
         IF (ios.ne.0) EXIT
         IF (j.eq.maxrecs) THEN
            write(*,*) 'Error: Maximum number of records exceeded...'
            write(*,*) 'Exiting program now...'
            STOP 
         ENDIF 
         nr=nr+1
      END DO 
      REWIND(2) 
   !Allocate data variables 
      ALLOCATE(wl(nr)) 
      ALLOCATE(fl(nr))
      ALLOCATE(flerr(nr))
      ALLOCATE(nnfl(nr))
   !Read data 
      DO j=1,nr 
         READ(2,*) wl(j),fl(j),flerr(j),nnfl(j)
      ENDDO 
      CLOSE(2)
c-------------------------------------------------------------------------------
      unit=3
c Check if fits file exists and delete if does
      call DELETEFILE(filename,status)
c Create a fits file
      call FTINIT(UNIT,filename,blocksize,status)
c Define primary array parameters
      call FTPHPR(UNIT,simple,bitpix,naxis,naxes,0,1,extend,status)
c-------------------------------------------------------------------------------
c Define data to be inputted into fits file
      nrowsll=nr
      nrows_z=1
      tfields=4
      tfields_z=1
      varidat=0
      extname='binary'
      DATA ttype_z/'z'/
      DATA tunit_z/' '/
      DATA tform_z/'D'/
      DATA ttype/'lambda','flux','flux_error','no_noise_flux'/
      DATA tunit/'Angst','J','J','J'/
      DATA tform/'E','E','E','E'/
c Create a binary table HDU with redshift
      call FTIBIN(UNIT,nrows_z,tfields_z,ttype_z,tform_z,tunit_z,
     &            extname,varidat,status)
      call FTPCLE(unit,1,1,1,1,z,status)
      call FTMNAM(unit,'ttype1','Z',status)
      call FTMKYF(unit,'Z',2.345,3,'redshift',status)
c Create a binary table HDU with lambda, flux, flux error, no-noise-flux
      call FTIBIN(UNIT,nrowsll,tfields,ttype,tform,tunit,
     &            extname,varidat,status)
c Input data into fits file
      call FTPCLE(unit,1,1,1,nrowsll,wl,status)
      call FTPCLE(unit,2,1,1,nrowsll,fl,status)
      call FTPCLE(unit,3,1,1,nrowsll,flerr,status)
      call FTPCLE(unit,4,1,1,nrowsll,nnfl,status)
c Close fits file
      call FTCLOS(unit,status)
c End program
      END PROGRAM dat2fits

c-------------------------------------------------------------------------------
c SUBROUTINES
      subroutine deletefile(filename,status)
C  A simple little routine to delete a FITS file
      integer status,unit,blocksize
      character*20 filename
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
