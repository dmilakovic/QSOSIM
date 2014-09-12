!======================================================================
      PROGRAM main
c  PURPOSE: generate artificial spectra
!======================================================================
      CHARACTER :: infile*20, outfile*20, descriptor*6
      INTEGER,PARAMETER :: nrows=25
      INTEGER,PARAMETER :: npts=1000
      INTEGER :: inoise,npts
      INTEGER,DIMENSION(nrows) :: numlls
      REAL :: wstart,wend,dw,nc,nce,nuplim,nuplime,dvavoid
      REAL,DIMENSION(nrows) :: ra,dec,zqso,alpha,vmag,s2n,sigblur
      REAL*4 :: lambda(262144),flux(262144)
      REAL*4 :: flerr(262144), nnflux(262144)
      real*4, dimension(npts) :: x,y,F,H
      real :: numlin
      EXTERNAL :: qsosim9, spline, readfits, writefits

      infile='sin.fits'
      call readfits(infile, wstart,wend,dw,nc,nuplim,inoise,dvavoid,
     &              ra,dec,zqso,alpha,vmag,sigblur,s2n,numlls)
 100  format(f10.5,2x,f10.5,2x,f8.5,2x,f8.5)
      nc=nc*1E12
      nuplim=nuplim*1E16
      
      call spline(x,y,F,H,numlin)
      do i=1,1
         write (*,*) 'nc: ',nc
         call qsosim9(zqso(i),alpha(i),vmag(i),wstart,wend,dw,nc,
     +          nuplim,sigblur(i),s2n(i),inoise,numlls(i),dvavoid,
     +          npts,lambda,flux,flerr,nnflux)
         write (*,*) 'nc :',nc
c         if (i.eq.5.or.i.eq.7) then
c            write (*,100) lambda,flux,flerr,nnflux
c         end if
         write (descriptor,'(I6.6)') i
         outfile='spec-'//descriptor//'.fits'
         call writefits(outfile,ra(i),dec(i),zqso(i),alpha(i),npts,
     &                     lambda,flux,flerr,nnflux)
         write (*,*)'--------------------------------------------------'
      enddo
      END PROGRAM main
