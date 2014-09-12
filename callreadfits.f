      PROGRAM testmemory
c  PURPOSE: check whether 'call readfits()' returns correct values to the main
c  program, in regards to problems with 'ra' array

      CHARACTER :: filename*20
      INTEGER,PARAMETER :: nrows=25
      INTEGER :: inoise
      INTEGER,DIMENSION(nrows) :: numlls
      REAL :: wstart,wend,dw,nc,nuplim,dvavoid
      REAL,DIMENSION(nrows) :: ra,dec,zqso,alpha,vmag,s2n,sigblur
      REAL*4 :: lambda(262144),flux(262144)
      REAL*4 :: sigma(262144), nnflux(262144)

      filename='sin.fits'
      call readfits(filename, wstart,wend,dw,nc,nuplim,inoise,dvavoid,
     &              ra,dec,zqso,alpha,vmag,sigblur,s2n,numlls)
 100  format(i2,4x,f10.5,2x,f10.5,2x,f8.5,2x,f8.5)
      do i=1,2
         call qsosim9(zqso,alpha,vmag,wstart,wend,dw,nc,nuplim,
     :          sigblur,s2n,inoise,numlls,dvavoid,lambda,flux,sigma,
     :          nnflux)
         write(*,*)zqso
      enddo
      END PROGRAM testmemory
