      program plot
! ======================================================================
! plot data from as many sources of log dn/dz and reproduce results for
! quasars from Kim et al. 2014 -- check for validity of reported results
	integer,parameter :: npts=1000
	integer,parameter :: n=18
        character(3) :: label
	real*8, dimension(npts) :: xmesh, ymesh1,ymesh2,ymesh3
	real*8, dimension(npts) :: xs,ys,CDDF
	real*8, dimension(3) :: logA, gamma, naux
	real*8, dimension(n) :: z,dz,logzp1,zstart,zend
	real*8, dimension(n) :: dndz1,dndz2,dndz1_err,dndz2_err
	real*8, dimension(n) :: logdndz1,logdndz2,logdndz3
	real*8, dimension(n,3) :: logdndz_aux
	real*8 :: xstart,xend,dx, ystart,yend,dy
	real*8 :: beta, nl, corrA
        real*4 :: x,y

	data logA/1.52,0.72,-0.7569/
	data gamma/1.51,2.16,1.33/

! ======================================================================
	open (unit=1, file='kim13res.dat')
	do i=1,n
         ! format : <z>, log dn/dz [12.75-14], log dn/dz [14-17], Delta z
         read (1,*,end=100) z(i),zstart(i),zend(i),dz(i),
     +                      dndz1(i),dndz1_err(i),dndz2(i),dndz2_err(i)
         logzp1(i)=dble(alog10(real(z(i)+1.)))
	enddo
 100	close(1)
! ====================================================================== 
	call spline(npts,xs,ys,CDDF,nl)
	do i=1,n
	   call power_laws(npts,zstart(i),zend(i),xs,ys,logA,gamma,corrA,
     +                     naux)
	   logdndz_aux(i,1)=alog10(real(naux(1)/(zend(i)-zstart(i))))
	   logdndz_aux(i,2)=alog10(real(naux(2)/(zend(i)-zstart(i))))
	   logdndz_aux(i,3)=alog10(real(naux(3)/(zend(i)-zstart(i)))) 
	end do
	do i=1,n
	   logdndz1(i)=dble(logdndz_aux(i,1))
	   logdndz2(i)=dble(logdndz_aux(i,2))
	   logdndz3(i)=dble(logdndz_aux(i,3))
        end do
 200    format (f5.3,5x,f5.3)
        write (6,*) 'data     model'
        do i=1,n
           write (6,200) dz(i), zend(i)-zstart(i)
        end do
! ======================================================================
	xstart=0.0
	xend=1.0
	dx=(xend-xstart)/npts
	do i=1,npts
         xmesh(i)  = xstart+dx*dble(i-1)
         ymesh1(i) = alog10(real(corrA)) + gamma(1)*(xmesh(i))
         ymesh2(i) = logA(2) + gamma(2)*(xmesh(i))
	 ymesh3(i) = logA(3) + gamma(3)*(xmesh(i))
	end do
	ymin=min(minval(ymesh1),minval(ymesh2))
	ymax=max(maxval(ymesh1),maxval(ymesh2))
! ======================================================================
!   PLOT DATA + MODELS
! ======================================================================
	call pgbegin (0,'/xserve',1,1)
	call pgslw(2)
	call pgenv(0.35,0.75,0,ymax,0,1)
	call pglabel('log(1+z)','log dn/dz',
     +'Kim et al. 2013 linear regression lines + data')
! ----------------------------------------------------------------------
!       plot models + data
! ----------------------------------------------------------------------
        call pgsci(1)
	call pgline(npts,real(xmesh),real(ymesh1))
	call pgline(npts,real(xmesh),real(ymesh2))
	call pgline(npts,real(xmesh),real(ymesh3))
	call pgpt(n,real(logzp1),real(dndz1),7)
	call pgerry(n,real(logzp1),
     +       real(dndz1+dndz1_err),
     +       real(dndz1-dndz1_err),2.0)
	call pgpt(n,real(logzp1),real(dndz2),6)
	call pgerry(n,real(logzp1),
     +       real(dndz2+dndz2_err),real(dndz2-dndz2_err),2.0)
	call pgsci(2)
        do i=1,n
           write (label,'(i2)') i
           call pgsci(i)
           call pgtext(real(logzp1(i)),0.1,label)
	call pgpt(1,real(logzp1(i)),real(logdndz1(i)),7)
        end do
c	call pgpt(n,real(logzp1),real(logdndz2),6)
	call pgtext(0.38,0.3,'Numerical results (spline+3 power laws)')
	call pgsci(1)
	call pgtext(0.58,0.6,'log N = [12.75,14.00]  ' //char(7))
	call pgtext(0.58,0.5,'log N = [14.00,17.00]  ' //char(6))
	call pgtext(0.38,0.4,'Kim et al 2013 data')
	call pgend

	end program plot
