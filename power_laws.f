! -------------------------------------------------------------------
      subroutine power_laws(npoints,zstart,zqso,xs,ys,
     +                      bigA,gamma,nl,corr)
! -------------------------------------------------------------------
C PURPOSE: prepare constants to be used in calculating the number of 
C          lines and power laws from Kim et al 2013
      implicit none
      integer :: i,index12,index13p1,index14,index17,index22,index12p75
      integer :: nl,npoints
      real*8,dimension(npoints) :: xs,ys
      real*8,dimension(3) :: bigA, gamma, n
      real*8 :: dmin,d12,d12p75,d13p1,d14,d17,d22
      real*8 :: int12to14,int13p1to14,int12p75to14
      real*8 :: n1,n2,n3,total
      real*8 :: z1,z1p1,z2,z2p1,zqso, gp1
      real*8 :: zstart,corr,beta, dx
      real*8 :: s1,s2
! -------------------------------------------------------------------
! CALCULATE INDICES for data manipulation
! -------------------------------------------------------------------     
      dmin=0.1
      do i=1,npoints
         d12=abs(xs(i)-12.0)
         if(d12.lt.dmin)then
            dmin=d12
            index12=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npoints
         d12p75=abs(xs(i)-12.75)
         if(d12p75.lt.dmin)then
            dmin=d12p75
            index12p75=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npoints
         d13p1=abs(xs(i)-13.1)
         if(d13p1.lt.dmin)then
            dmin=d13p1
            index13p1=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npoints
         d14=abs(xs(i)-14.0)
         if(d14.lt.dmin)then
            dmin=d14
            index14=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npoints
         d17=abs(xs(i)-17.0)
         if(d17.lt.dmin)then
            dmin=d17
            index17=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npoints
         d22=abs(xs(i)-22.0)
         if(d22.lt.dmin)then
            dmin=d22
            index22=i
         end if
      end do
! -------------------------------------------------------------------
! INTEGRATION to account for a gap in data from 12 to 13.1
! -------------------------------------------------------------------
      !numerically integrate spline from 12 to 14 and from 13.1 to 14
      !in steps of dx to get the correction factor 'corr'
      int12to14=0
      int12p75to14=0
      do i=index12,index14
         dx=10**xs(i+1)-10**xs(i)             
         int12to14=int12to14+10**ys(i)*dx
      end do
      do i=index12p75,index14
         dx=10**xs(i+1)-10**xs(i)             
         int12p75to14=int12p75to14+10**ys(i)*dx
      end do
      int13p1to14=0
      do i=index13p1,index14
         dx=10**xs(i+1)-10**xs(i)             
         int13p1to14=int13p1to14+10**ys(i)*dx
      end do
      corr=int12p75to14/int13p1to14
c      write (6,*) 'int[12.75->14] = ',int12p75to14
c      write (6,*) 'int[13.10->14] = ',int13p1to14
C     corr - correction factor

! -------------------------------------------------------------------
! DETERMINE THE NUMBER OF LINES IN EACH COLUMN DENSITY BIN         
! -------------------------------------------------------------------
      ! dn/dz = A(1+z)^gamma
      ! n = A/(1+gamma)*[(1+z2)^(1+gamma) - (1+z1)^(1+gamma)]
      ! n1(12->14) = corr*A1/(1+gamma1)*[(1+z2)^(1+gamma1) - (1+z1)^(1+gamma1)]
      ! n2(14->17) = A2/(1+gamma2)*[(1+z2)^(1+gamma2) - (1+z1)^(1+gamma2)]
      ! n3(17->22) = A3/(1+gamma3)*[(1+z2)^(1+gamma3) - (1+z1)^(1+gamma3)]
! -------------------------------------------------------------------
      
      z1=zstart
      z1p1=z1+1.
      z2=zqso
      z2p1=zqso+1.
      beta=abs((ys(index14)-ys(index12))/(xs(index14)-xs(index12)))
c      write (6,*) ((10**12.00)/(10**13.1))**(1.-1.46)
      corr=(int12p75to14/int13p1to14)!*((10**12.00)/(10**13.1))**(1.-1.45)
      write (6,*) 'Numerical correction factor=',corr
! -------------------------------------------------------------------
      total=0.0
      do i=1,3
         gp1=gamma(i)+1.
         if(i.eq.1)then
            n(i) = corr*bigA(i)/(gp1)*((z2p1)**(gp1)-(z1p1)**(gp1))
         else
         n(i) = bigA(i)/(gp1)*((z2p1)**(gp1)-(z1p1)**(gp1))
         end if
         total=total+n(i)
      end do
      nl=nint(total)
      n1=nint(n(1))
      n2=nint(n(2))
      n3=nint(n(3))
c      write (6,*) 'A = corr*bigA/(1+gamma)'
c      write (6,*) 'A (numerical) =',corr*bigA(1)/(gamma(1)+1.)
c      write (6,*) 'gamma = ',gamma(1)
c      write (6,*) 'z1 =',z1
c      write (6,*) 'z2 =',z2
      write (6,*) 'n [12.75->14.00] (numerical) =',nint(n(1))
      write (6,*) 'n [14.00->17.00] (numerical) =',nint(n(2))
      write (6,*) 'n [17.00->22.00] (numerical) =',nint(n(3))
      write (6,*) 'Number of lines =', nl
      end subroutine power_laws
