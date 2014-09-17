c-------------------------------------------------------------------------------
      SUBROUTINE readfits(filename, wstart,wend,dw,nc,nuplim,inoise,
     +                   dvavoid,ra,dec,zqso,alpha,vmag,sigblur,s2n)
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
c------------------------------------------------------------------------------
c GENERAL DECLARATIONS
      IMPLICIT NONE
c Declare variables
      INTEGER :: status,unit,readwrite,blocksize,hdutype,ntable,inoise
      INTEGER,PARAMETER :: nrows=25
      INTEGER :: felem,nelems,nullj,nfound,irow,colnum
      REAL :: nulle
      REAL*8 :: wstart,wend,dw
      REAL*8 :: nc, nuplim,dvavoid
      INTEGER,DIMENSION(nrows):: numlls
      REAL*8,DIMENSION(nrows) :: ra,dec,zqso,alpha,vmag,
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
      call FTGKYD(unit,'wstart',wstart,comment,status)
      call FTGKYD(unit,'wend',wend,comment,status)
      call FTGKYD(unit,'dw',dw,comment,status)
      call FTGKYD(unit,'nc',nc,comment,status)
      call FTGKYD(unit,'nuplim',nuplim,comment,status)
      call FTGKYJ(unit,'inoise',inoise,comment,status)
      call FTGKYD(unit,'dvavoid',dvavoid,comment,status)
 100  format(2x,a10,d9.3,a4)
 125  format(2x,a10,f9.3,a4)
 150  format(2x,a10,i3)
      write (*,125)'wstart =',wstart
      write (*,125)'wend =',wend
      write (*,125)'dw =',dw
      write (*,100)'nc =',nc
      write (*,100)'nuplim =',nuplim
      write (*,150)'inoise =',inoise
      write (*,125)'dvavoid =',dvavoid
c------------------------------------------------------------------------------
c Read contents of 'QSO' (ntable=3) binary table
      ntable=3
      call ftmahd(unit,ntable,hdutype,status)
c      ALLOCATE(ra(nrows),dec(nrows),zqso(nrows),alpha(nrows))
c      ALLOCATE(vmag(nrows),sigblur(nrows))
c      ALLOCATE(numlls(nrows),s2n(nrows))
      if (status.eq.0)then 
         call FTGERR(status, errtext)
         print *,status, errtext
      else 
         call FTGERR(status, errtext)
         print *,status,' ',errtext
      end if
c Read column data, one row at a time, and print them out
      felem=1
      nelems=1
      nulle=0.      
      nullj=0
      write (*,300)'RA','DEC','Z','alpha'
      do irow=1,nrows
            call FTGCVD(unit,1,irow,felem,nelems,nulle,ra(irow),
     &       anynull,status)
            call FTGCVD(unit,2,irow,felem,nelems,nulle,dec(irow),
     &       anynull,status)
            call FTGCVD(unit,3,irow,felem,nelems,nulle,zqso(irow),
     &       anynull,status)
            call FTGCVD(unit,4,irow,felem,nelems,nulle,alpha(irow),
     &       anynull,status)
            call FTGCVD(unit,5,irow,felem,nelems,nulle,vmag(irow),
     &       anynull,status)
            call FTGCVD(unit,6,irow,felem,nelems,nulle,sigblur(irow),
     &       anynull,status)
            call FTGCVD(unit,7,irow,felem,nelems,nulle,s2n(irow),
     &       anynull,status)
            write (*,200)irow,ra(irow),dec(irow),zqso(irow),alpha(irow)
      end do
 200  format(i2,4x,f10.5,2x,f10.3,2x,f10.5,2x,f6.4)
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
         call FTGERR(status, errtext)
         print *,status,' ',errtext
      end if
c------------------------------------------------------------------------------
      RETURN
      END SUBROUTINE readfits
c------------------------------------------------------------------------------
