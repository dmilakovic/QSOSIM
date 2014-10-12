!======================================================================
      PROGRAM generate
c  PURPOSE: call other subroutines to generate artificial spectra!
c  OUTPUT:  artificial SDSS catalogue
!======================================================================
      IMPLICIT NONE
      CHARACTER :: infile*20, outfile*20, descriptor*6
      INTEGER,PARAMETER :: nrows=20
      INTEGER,PARAMETER :: npoints=1000
      INTEGER :: i,j,inoise, numlin, npts, nl
      INTEGER,DIMENSION(nrows) :: numlls
      REAL*8 :: wstart,wend,dw
      REAL*8 :: nc,nuplim,dvavoid
      REAL*8 :: zstart,corr
      REAL*8,DIMENSION(nrows) :: ra,dec,zqso,alpha,vmag,s2n,sigblur
      REAL*8 :: lambda(262144),flux(262144), da4(262144)
      REAL*8 :: flerr(262144), nnflux(262144)
      real*8, dimension(npoints) :: xs,ys,CDDF,H
      REAL*8, dimension(3) :: bigA,gamma,n !change conflict a-bigA in qsosim9
      EXTERNAL :: qsosim9, spline, readfits, writefits, absdist

!====================================================================== 
! ---------------------------------------------------------------------
! DATA USED
! ---------------------------------------------------------------------
      ! from Kim et al 2014
      !     log NHI        log A        log gamma
      !   [13.1,14.0]      1.52+-0.05   1.51+-0.09
      !   [14.0,17.0]      0.72+-0.08   2.16+-0.14
      ! from Unknown
      !   [17.0,22.0]      log(0.175)   1.33

      data bigA/33.11311,5.24807,0.1750/
      data gamma/1.51,2.16,1.33/

! ---------------------------------------------------------------------
! READ INPUT PARAMETERS
! ---------------------------------------------------------------------
      infile='sin.fits'
      write (6,*)'=================================================='
      write (6,*)'                    QSOSIM 10                     '
      write (6,*)'=================================================='
      call readfits(infile,wstart,wend,dw,nc,nuplim,inoise,dvavoid,
     &             ra,dec,zqso,alpha,vmag,sigblur,s2n)
 100  format(f10.5,2x,f10.5,2x,f8.5,2x,f8.5)
! ---------------------------------------------------------------------
! GET THE SPLINE
! ---------------------------------------------------------------------
      call spline(npoints,xs,ys,CDDF)
! ---------------------------------------------------------------------
! GENERATE ARTIFICIAL SPECTRA
! ---------------------------------------------------------------------
      do i=8,8
         zstart=(3650.0/1215.67)-1.
         write (descriptor,'(I6.6)') i
         write (6,*)'=================================================='
         write (6,*)'               Spectrum no. ',descriptor
         write (6,*)'=================================================='
         call power_laws(npoints,zstart,zqso(i),xs,ys,
     +                   bigA,gamma,nl,corr)
c         call absdist(wstart,zqso(i),X)
c         numlin=nint(nl*X)
c         write (6,*) 'Total number of lines = F(NHI) * X(z)'
c         write (6,*) 'Total number of lines = ', numlin
         call qsosim9(zqso(i),alpha(i),vmag(i),wstart,wend,dw,
     +          nc,nuplim,sigblur(i),s2n(i),inoise,dvavoid,
     +          npts,lambda,flux,flerr,nnflux,npoints,xs,ys,CDDF,nl,
     +          bigA,gamma)
c         do j=1,npoints,20
c            write (*,*) j, xs(j), CDDF(j)
c         end do
         outfile='spec-'//descriptor//'.fits'
         call writefits(outfile,ra(i),dec(i),zqso(i),alpha(i),npts,
     &                     lambda,flux,flerr,nnflux)
         write (*,*)'--------------------------------------------------'
      end do
 200  format (4d10.3)
c      do i=1,npts
c      write (6,200) lambda(i), flux(i), flerr(i), nnflux(i)
c      end do
      call PGBEGIN (0,'/xserve',1,1)
      call PGSLW(1)
      call PGENV (3550.,10500.,0.0,real(maxval(flux)),0,1)
      call PGLABEL ('lambda','flux','QSO spectrum')
      call pgline(npts,real(lambda),real(flux))
      call pgsci(2)
      call pgline(npts,real(lambda),real(nnflux))
      call PGEND
! ---------------------------------------------------------------------      
      END PROGRAM generate
