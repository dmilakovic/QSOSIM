c-------------------------------------------------------------------------------
      SUBROUTINE readfits(filename, wstart,wend,dw,nc,nuplim,inoise,
     +                   dvavoid,ra,dec,zqso,alpha,vmag,sigblur,s2n,
     +                   numlls)
c     PURPOSE: read fits file containing data for qsosim9
c     INPUT: filename  name of the fits file
c     OUTPUT: 
*            wstart    starting wavelength
*            wend      ending wavelength
*            dw        pixel size
*            nc        N(HI) lower cut-off
*            nuplim    N(HI) upper cut-off
*            inoise    inoise parameter
*            dvavoid   avoidance zone around each additional system in km/h
*            ra        RA coordinates (array)
*            dec       DEC coordinates (array)
*            zqso      Z (array)
*            alpha     spectral index (array)
*            vmag      V magnitude (array)
*            sigblur   spectral resolution (array)
*            numlls    number of user specified additional absorption systems
c------------------------------------------------------------------------------
c GENERAL DECLARATIONS
      IMPLICIT NONE
c Declare variables
      INTEGER :: status,unit,readwrite,blocksize,hdutype,ntable,inoise
      INTEGER,PARAMETER :: nrows=25
      INTEGER :: felem,nelems,nullj,nfound,irow,colnum
      REAL :: nulle,density
      REAL :: wstart,wend,dw,nc,dvavoid,nuplim
      INTEGER,DIMENSION(nrows):: numlls
      REAL,DIMENSION(nrows) :: ra,dec,zqso,alpha,vmag,
     &                                 sigblur,s2n
      character :: filename*20,ttype(20)*10!, nhills,blls,zlls
      logical :: anynull
      character :: errtext*30,card*50, comment*30
c Define parameters
      status=0
      unit=2
      readwrite=0
c------------------------------------------------------------------------------
c Open fits file to read
      call ftopen(unit,filename,readwrite,blocksize,status)
      call ftgerr(status,errtext)
      if (status.eq.0) then 
         write (*,*)status,'File opened'
      end if
c------------------------------------------------------------------------------
c Read contents of 'GENERAL' (ntable=2) binary table
c Read data from columns
      ntable=2
      call ftmahd(unit,ntable,hdutype,status)
      call FTGKYE(unit,'wstart',wstart,comment,status)
      call FTGKYE(unit,'wend',wend,comment,status)
      call FTGKYE(unit,'dw',dw,comment,status)
      call FTGKYE(unit,'nc',nc,comment,status)
      call FTGKYE(unit,'nuplim',nuplim,comment,status)
      call FTGKYJ(unit,'inoise',inoise,comment,status)
      call FTGKYE(unit,'dvavoid',dvavoid,comment,status)
 100  format(2x,a10,f9.3,a4)
 150  format(2x,a10,i3)
      write (*,100)'wstart =',wstart
      write (*,100)'wend =',wend
      write (*,100)'dw =',dw
      write (*,100)'nc =',nc,'xE12'
      write (*,100)'nuplim =',nuplim,'xE16'
      write (*,150)'inoise =',inoise
      write (*,100)'dvavoid =',dvavoid
c------------------------------------------------------------------------------
c Read contents of 'QSO' (ntable=3) binary table
      ntable=3
      call ftmahd(unit,ntable,hdutype,status)
c      ALLOCATE(ra(nrows),dec(nrows),zqso(nrows),alpha(nrows))
c      ALLOCATE(vmag(nrows),sigblur(nrows))
c      ALLOCATE(numlls(nrows),s2n(nrows))
c Read column data, one row at a time, and print them out
      felem=1
      nelems=1
      nulle=0.      
      nullj=0
      write (*,300)'RA','DEC','Z','alpha','num LLS'
      do irow=1,nrows
            call FTGCVE(unit,1,irow,felem,nelems,nulle,ra(irow),
     &       anynull,status)
            call FTGCVE(unit,2,irow,felem,nelems,nulle,dec(irow),
     &       anynull,status)
            call FTGCVE(unit,3,irow,felem,nelems,nulle,zqso(irow),
     &       anynull,status)
            call FTGCVE(unit,4,irow,felem,nelems,nulle,alpha(irow),
     &       anynull,status)
            call FTGCVE(unit,5,irow,felem,nelems,nulle,vmag(irow),
     &       anynull,status)
            call FTGCVE(unit,6,irow,felem,nelems,nulle,sigblur(irow),
     &       anynull,status)
            call FTGCVE(unit,7,irow,felem,nelems,nulle,s2n(irow),
     &       anynull,status)
            call FTGCVJ(unit,8,irow,felem,nelems,nullj,numlls(irow),
     &       anynull,status)
            write (*,200)irow,ra(irow),dec(irow),zqso(irow),alpha(irow),
     &       numlls(irow)
      end do
 200  format(i2,4x,f10.5,2x,f10.3,2x,f10.5,2x,f6.4,2x,i7)
 300  format(6x,a10,2x,a10,2x,a10,2x,a6,2x,a7)
c Read column names that are relevant for LLS
c      call FTGKNS(unit,'TTYPE',9,20,ttype,nfound,status)
c      write (*,*)'found TTYPE ',nfound
c      call FTGKNE(unit,'blls',1,20,blls,nfound,status)
c      call FTGKNE(unit,'zlls',1,20,zlls,nfound,status)
c------------------------------------------------------------------------------
c Close fits file
      call ftclos(unit,status)
      if (status.eq.0)then 
         print *,status,' File closed'
      else 
         print *,status,' ',errtext
      end if
c------------------------------------------------------------------------------
      RETURN
      END SUBROUTINE readfits
c------------------------------------------------------------------------------
