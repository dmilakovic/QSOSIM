!====================================================================
      subroutine spline(npts,xs,ys,CDDF,nl)
!====================================================================
! Spline interpolation
! Comments: values of function f(x) are calculated in n base points
! then: spline coefficients are computed
!       spline interpolation is computed in 2n-1 points, 
!       a difference sum|f(u)-ispline(u)| 
!====================================================================
      implicit none
      real :: xmin, xmax, ymin, ymax      ! given interval of x()
      integer, parameter :: n=8 ! number of points of spline
      real*8, dimension(n) :: xi,yi,b,c,d
      real*8, dimension(npts) :: xs,xe,ys,ye,z,CDDF,H
      real*8, dimension(npts) :: sumN,sumX,sumW, weight, lines
      integer i,j,numlin,npts,nl
      real*8 ispline
      real*8 dx, dxe, total
      real*8  zstart,zend, dz, bigX
      external func, gamma, gauss16
      real*8 func, gamma, gauss16

c Data points
          data xi/12.0,15.0,17.0,18.0,20.0,21.0,21.5,22.0/
          data yi/-9.71,-14.41,-17.94,-19.39,-21.28,-22.82,-23.95,-25.5/

c Spline coefficients: 
c s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
          data b/-1.09781456,-1.81851768,-1.63473213,-1.18404174,
     +           -1.16128039,-1.96163607,-2.62445641,-3.62053967/
          data c/-0.228617609,-1.16167441E-02,0.103509508,0.347180933,
     +           -0.335800260,-0.464555502,-0.861085415,-1.13108110/
          data d/2.41112094E-02,1.91877093E-02,8.12238082E-02,
     +          -0.113830194,-4.29184213E-02,-0.264353245,-0.179997146,
     +          -0.179997146/
      xmin = 12.0
      xmax = 22.0
!====================================================================
      write (6,*) '--------------------------------------------------'
      write (6,*) 'Started SPLINE!'
      write (6,*) 'Calculating dn/dNHI...'
!====================================================================
!  step 3: interpolation at npts points
c      errav = 0.0
      dx = (xmax-xmin)
      total=0.
!---- calculate normalization constant ------------------------------
      do i=1, npts
         total=total+gauss16(gamma,dble(0.0),dble(i))
      end do
!---- calculate weight of steps -------------------------------------
      do i=1, npts
         weight(i) = 0
         sumW(i) = 0
         do j=1,i
            sumW(i) = sumW(i) + gauss16(gamma,dble(0.0),dble(j))
         end do
         weight(i) =float(npts+1-i+1)/(npts-1)! 1 - sumW(i)/total !
      end do
!---- interpolate spline ys on x points -----------------------------
      do i=1,npts
         xs(i)= xmin + dx*dble(i-1)/(npts-1)
         xe(i) = 10**xs(i)
         ys(i) = ispline(xs(i), xi, yi, b, c, d, n)
         ye(i) = 10**ys(i)             
      end do
 200  format (4d12.4)
 201  format ('           i           sum              F')    
 202  format ('           Average error',f12.5)
!====================================================================
! step 5: calculate CDDF
      total=0
      do i=1,npts
         if (i.eq.npts) then 
            dx=0
            dxe=0
         else 
            dx=xs(i+1)-xs(i)
            dxe=xe(i+1)-xe(i)
         end if
         total=total+ye(i)*dxe  !ys*dx
      end do
      do i=1,npts
         sumN(i)=0
         CDDF(i)=0
         lines(i)=0
         do j=1,i
            if (j.eq.i) then 
               dx=0
               dxe=0
            else 
               dx=xs(j+1)-xs(j)
               dxe=xe(j+1)-xe(j)
            endif
            sumN(i)=sumN(i)+ye(j)*dxe
            lines(i)=lines(i)+ye(j)*dxe
         end do
         CDDF(i)=sumN(i)/total  !sumN(i)
      end do
      
!====================================================================
! step 6: calculate number of lines, n
! n = integrate [f(NHI,z) dX/dz] dNHI dz
! F = integrate [f(NHI,z) dNHI]
! dX/dz = [H0*(1+z)^2]/H(z)
! H(z) = H0 * Sqrt{[OmegaM(1+z)^3 + OmegaRAD(1+z)^4 + OmegaDE]}
!      OmegaM  = 0.27
!      OmegaDE = 0.73
!      OmegaRAD approx 0
! H = integrate [dX/dz] dz
! n(i) = H * sumN(i)
!====================================================================
! calculate 
c      zstart=1.96132994   !zstart=(wstart/1215.67)-1, wstart=3500
c      zend=2.83167458
c      dz=(zend-zstart)/npts
c      total=gauss16(func,zstart,zend)
c      do i=1,npts
c         z(i)=zstart+dz*float(i-1)
c      end do
c      do i=1,npts
c         sumX(i) = gauss16(func,z(1),z(i))
c         H(i) = sumX(i)/total
!         write (*,*) i,sum(i),H(i)
c      end do
!      do i=1,npts
!!         write (*,*) i, x(i), sumN(i), sumX(i)
!      end do
c      bigX=gauss16(func,zstart,zend)
c      write (6,*) bigX
      nl=nint(sumN(npts))!*gauss16(func,zstart,zend))
      write (6,*) 'No. of lines in log NHI range, F(NHI) = ', nl
!====================================================================      
      
      end subroutine spline


!======================================================================
!  Function X(z)
!======================================================================
      function func(x)
      real*8 f, x
      f = ((1+x)**2)*(0.3*((1+x)**3) + 0.7)**(-0.5)
      return
      end function func
!======================================================================
!  Function Gamma(k,x)
!======================================================================
      function gamma(x)
      real*8 f, x, k, th
      k=1.0
      th=0.5
      f = x**(k-1)*exp(-x/th)
      return
      end function gamma

!======================================================================
      function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
      implicit none
      real*8 ispline
      integer n
      real*8  u, x(n), y(n), b(n), c(n), d(n)
      integer i, j, k
      real*8 dx

! if u is ouside the x() interval take a boundary value (left or right)
      if(u <= x(1)) then
         ispline = y(1)
         return
      end if
      if(u >= x(n)) then
         ispline = y(n)
         return
      end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
      i = 1
      j = n+1
      do while (j > i+1)
         k = (i+j)/2
         if(u < x(k)) then
            j=k
         else
            i=k
         end if
      end do
!*
!  evaluate spline interpolation
!*
      dx = u - x(i)
      ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end function ispline

!==========================================================
      Function gauss16(f,a,b)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Gauss 16 points  
! written by: Alex Godunov (October 2009)
! http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch03/gauss.f90
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! OUT:
! gauss8 - Result of integration
!==========================================================
      implicit none
      integer, parameter :: n=8
      real*8 gauss16, f
      real*8 a, b
      real*8 ti(n), ci(n)
      data ti/0.0950125098, 0.2816035507, 0.4580167776, 0.6178762444,   
     &   0.7554044083, 0.8656312023, 0.9445750230, 0.9894009349/ 
      data ci/0.1894506104, 0.1826034150, 0.1691565193, 0.1495959888,
     &   0.1246289712, 0.0951585116, 0.0622535239, 0.0271524594/ 
      real*8 r, m, c
      integer i

      r = 0.0;
      m = (b-a)/2.0;
      c = (b+a)/2.0;

      do i = 1,n 
         r = r + ci(i)*(f(m*(-1.0)*ti(i) + c) + f(m*ti(i) + c))
      end do
      gauss16 = r*m
      return
      end function gauss16
