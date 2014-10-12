      program plot_pl
!=================================================
!plot theoretical power laws against data provided by Kim et al 2013.

      integer,parameter :: npts=1000
      integer,parameter :: n=18
      real*8, dimension(npts) :: xmesh, ymesh1, ymesh2
      real*8, dimension(2) :: A, gamma
      real*8, dimension(n) :: z,dz,dndz1,dndz2, logzp1
      real*8 :: xstart,xend,dx, ystart,yend,dy
      real*8 :: beta
! data set no 1.      
c      data A/1.52,0.72/
c      data gamma/1.51,2.16/
! data set no 2.
c      data A/1.88,1.03/
c      data gamma/0.89,1.61/
! data set no 3.
      data A/1.46,0.03/
      data gamma/1.67,3.40/ 
      
      open (unit=1, file='kim13.dat')
      do i=1,n
         ! format : <z>, log dn/dz [12.75-14], log dn/dz [14-17], Delta z
         read (1,*,end=100) z(i),dndz1(i),dndz2(i),dz(i)
         logzp1(i)=dble(alog10(real(z(i))))
      enddo
 100  close(1)
      
      xstart=0.0
      xend=1.0
      dx=(xend-xstart)/npts
      do i=1,npts
         xmesh(i)  = xstart+dx*dble(i-1)
         ymesh1(i) = A(1) + gamma(1)*(xmesh(i))
         ymesh2(i) = A(2) + gamma(2)*(xmesh(i))
      end do
      ymin=min(minval(ymesh1),minval(ymesh2))
      ymax=max(maxval(ymesh1),maxval(ymesh2))
      call pgbegin (0,'/xserve',1,1)
      call pgenv(0.0,1.0,0,ymax,0,1)
      call pglabel('log(1+z)','log dn/dz',
     + 'Kim et al. 2013 linear regression lines + data')
      call pgline(npts,real(xmesh),real(ymesh1))
      call pgline(npts,real(xmesh),real(ymesh2))
      call pgpt(n,real(logzp1),real(dndz1),2)
      call pgpt(n,real(logzp1),real(dndz2),5)
      call pgtext(0.6,0.6,'log N = [12.75,14.00]  ' //char(2))
      call pgtext(0.6,0.4,'log N = [14.00,17.00]  ' //char(5))
      call pgend

      end program plot_pl
