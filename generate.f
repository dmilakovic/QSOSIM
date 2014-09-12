!======================================================================
      PROGRAM generate
c  PURPOSE: call other subroutines to generate artificial spectra!
c  OUTPUT:  artificial SDSS catalogue
!======================================================================
      IMPLICIT NONE
      CHARACTER :: infile*20, outfile*20, descriptor*6
      INTEGER,PARAMETER :: nrows=25
      INTEGER,PARAMETER :: npoints=2000
      INTEGER :: i,j,inoise, numlin, npts
      INTEGER,DIMENSION(nrows) :: numlls
      REAL :: wstart,wend,dw,nc,nce,nuplim,nuplime,dvavoid
      REAL,DIMENSION(nrows) :: ra,dec,zqso,alpha,vmag,s2n,sigblur
      REAL*8 :: lambda(262144),flux(262144), da4(262144)
      REAL*8 :: flerr(262144), nnflux(262144)
      real*4, dimension(npoints) :: xs,ys,CDDF,H
      EXTERNAL :: qsosim9, spline, readfits, writefits

      infile='sin.fits'
      call readfits(infile, wstart,wend,dw,nc,nuplim,inoise,dvavoid,
     &             ra,dec,zqso,alpha,vmag,sigblur,s2n,numlls)
 100  format(f10.5,2x,f10.5,2x,f8.5,2x,f8.5)
      nc=nc*1E12
      nuplim=nuplim*1E21 
      do i=1,1
         call spline(npoints,3.3,xs,ys,CDDF,numlin)
         call qsosim9(3.3,alpha(i),vmag(i),wstart,wend,dw,
     +          nc,nuplim,sigblur(i),s2n(i),inoise,numlls(i),dvavoid,
     +          npts,lambda,flux,flerr,nnflux,npoints,xs,ys,CDDF,numlin,
     +          da4)
c         do j=1,npoints,10
c            write (*,*) j, xs(j), CDDF(j)
c         end do
         write (descriptor,'(I6.6)') i
         outfile='spec-'//descriptor//'.fits'
         call writefits(outfile,ra(i),dec(i),zqso(i),alpha(i),npts,
     &                     lambda,flux,flerr,nnflux)
         write (*,*)'--------------------------------------------------'
      end do
      
c      call PGBEGIN (0,'/vcps',1,1)
c      call PGSLW(1)
c      call PGENV (wstart-100.,wend+100.,minval(flux),maxval(flux),0,1)
c      call PGLABEL ('lambda','flux','QSO spectrum')
c      call pgline(npts,lambda,flux)
c      call pgsci(2)
c      call pgline(npts,lambda,nnflux)
      END PROGRAM generate
