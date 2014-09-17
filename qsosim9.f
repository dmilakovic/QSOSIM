c ----------------------------------------------------------------------------
c QSOSIM9. Prog. to make synthetic QSO spectrum - John Webb, UNSW, Dec 2013
c Uses VPFIT Voigt profile programs
c Compile using Makefile
c ----------------------------------------------------------------------------
	SUBROUTINE QSOSIM9(zqso,alpha,vmag,wstart,wend,dw,nc,nuplim,
     :          sigblur,s2n,inoise,dvavoid,npts,lambda,flux,
     :          flerr,nnflux,npoints,xs,ys,CDDF,numlin)
	  real*4 wda(262144),danoabs(262144),tau(262144)
	  real*4 da4(262144),da4conv(262144)
	  real*4 da_err4(262144),da4smno(262144)
	  real*8 da(262144),da_err(262144),da_err4mod(262144)
	  real*4 wems(30),relstr(30)
	  real*4 sum,nhi,b,z,g,z1p1,z2p1,beta,x,gp1,w,ff
	  real*4 a12p5,rn,a13p75,a13p1,mbp1,zqsop1,pi
	  real*4 c,d,p,q,r,s,vlight,zleft,zright
	  real*4 epsilon, lognhi, delta,mdelta
	  real*4 nhills(20),blls(20),zlls(20)
	  real*8 zqso,alpha,vmag,wstart,wend,dw,nc,nuplim,sigblur
	  real*8 s2n,dvavoid
	  real*8 lambda(262144),flux(262144)
	  real*8 flerr(262144), nnflux(262144)
          real*8 alm, fik,asm
	  real*4 n1,n2,n3
	  real*4 factor
	  real*8, dimension(npoints) :: xs,ys,CDDF,dummy,nhi4,z4, dummyy
	  integer npts,npoints,idum,numlls,iflag,inoise,index,numlin

	  data pi/3.14159265/
	  data vlight/299792.458/
      external ran3, gasdev3
      common/vpc_ewllns/lbz(2500),lzz(2500),alm(2500),
     :                fik(2500),asm(2500),numel
c ----------------------------------------------------------------------------
c Emission line wavelengths
	  data wems/1026.0,1216.0,1302.0,1335.0,1400.0,1549.0,1640.0,
     +            1858.0,2000.0,2080.0,2140.0,2175.0,2225.0,2326.0,
     +            2423.0,2798.0,2970.0,3130.0,3200.0,3346.0,3426.0,
     +            3727.0,3869.0,3968.0,4068.0,4340.0,4861.0,4959.0,
     +            5007.0,1240.0/
c Approx. relative emission line strengths - from eyeballing a few spectra
	  data relstr/0.093,1.0,0.035,0.025,0.19,0.63,0.18,0.29,0.0049,
     +              0.041,0.0034,0.0076,0.0047,0.06,0.022,0.34,
     +              0.063,0.0073,0.0095,0.0052,0.01,0.0078,0.036,
     +              0.013,0.028,0.13,0.22,0.0093,0.034,0.4/

c Max no. of pixels, hard-coded in array declarations
	  write (6,*) '--------------------------------------------------'
	  write (6,*) 'Started QSOSIM9!'
	  nptsmax=262144
	  nemlines=30

c alpha is QSO spectral index, vmag is the V magnitude
c wstart, wend, dw are the start and end wavelengths and pixel size
c       nc is the column density cut-off (REAL, NOT log10)
	  zqsop1=zqso+1
	  f5550=10**(-(vmag+21.17)/2.5)
	  const=f5550/(1.0/5550.0**(2+alpha))
      sum=0.0
      npts=(wend-wstart)/dw
      if(npts.gt.nptsmax)then
      	write(6,*)' Max. no. of pixels = ',nptsmax
      	stop
      end if
	  do l=1,npts
c Make wavelength scale and underlying power-law QSO continuum
	     wda(l)=wstart+((l-1)*dw)
	     da(l)=const*(1.0/wda(l)**(2+alpha))
	  end do

c Put emission lines in. Guess/approximation for emission line width:
	  fwhm=60.0
	  sigma=fwhm/2.35

	  do m=1,nemlines
        wems(m)=wems(m)*zqsop1
		 do i=1,npts
           x=0.75*relstr(m)*const*(1.0/wems(m)**(2+alpha))
           g=x*exp(-.5*(((wda(i)-wems(m))/sigma)**2))
           da(i)=da(i)+g
         end do
      end do
	  write (6,*) 'Emission lines inputted!'
c Keep unabsorbed QSO spectrum
	  do i=1,npts
	    danoabs(i)=da(i)
	  end do

c Generate optical depth array in preparation for absorption input
	  do i=1,npts
	    tau(i)=0.0
	  end do
         write (*,*) 'Inputting absorption lines ...'
c Next section makes absorption lines
c g is from dn=A(1+z)^g dz, beta is from dn propto N^{-beta} dN.
c a13p75 is the value of A for lines above logN=13.75.
c gp1= gamma+1, mbp1=1-beta, nc is N cutoff, n is total no. of lines
      a13p75=10.0
      beta=1.7
      z1p1=(wstart/1215.67)
      z2p1=zqsop1
      mbp1=-beta+1.0

c Calculate the total no. of lines
      n=numlin
	
c Initialise random numbers
      idum=time()
      do i=1,npoints
	 dummy(i)=dble(i)
	 dummyy(i)=dble((nc**mbp1-(10**xs(i))**mbp1)/(nc)**mbp1)
      end do

c Call Voigt profile generator n times, once for each abs system
	  do i=1,numlin
	     
c Random selection of NHI and redshifts (MODIFIED)
 1	c=ran3(idum)
	index=0
	nhi4(i)=0
	z4(i)=0
	mdelta=1.0
        if(c.lt.1.0) then
! -------------------------------------------------------------------
	  do j=1,npoints
	     delta=abs(c-CDDF(j))
	     if (delta.lt.mdelta) then
		mdelta=delta
		index=j
	     end if
	  end do
        else
          goto 1
        end if
	lognhi=xs(index)
	nhi=10**xs(index)
	nhi4(i)=xs(index)
! --------------------------------------------------------------------
        if ((lognhi.ge.12.0).and.(lognhi.lt.14.0)) then
           g=1.51
           a13p1=10**1.52
	   beta=1.43
	   a=a13p1*((10**12.0)/(10**13.1))**(1.-beta)
	   n1=a/(g+1.)*((z2p1)**(g+1.)-(z1p1)**(g+1.))
        else if ((lognhi.ge.14.0).and.(lognhi.lt.17.0)) then
           g=2.16
	   a=10**0.72
	   n2=a/(g+1.)*((z2p1)**(g+1.)-(z1p1)**(g+1.))
        else if ((lognhi.ge.17.0).and.(lognhi.lt.22.0)) then
	   g=1.33
	   a=0.175
	   n3=a/(g+1.)*((z2p1)**(g+1.)-(z1p1)**(g+1.))
        end if
	d=ran3(idum)	
	gp1=g+1.
        p=z2p1**gp1
        q=z1p1**gp1
        x=alog10(d*(p-q)+q)
        x=x/gp1
        z=10**x -1.0
c	write (*,*) index, c, CDDF(index), mdelta, xs(index), z
	z4(i)=z
! --------------------------------------------------------------------

c b-params.  Guess at sigma and mean of b distribution of 3 and 23 km/s.
        b = 3*gasdev3(idum)+23

c Put the forest lines in, but only if:
c 1. it is not if within the specified avoidance zone of each LLS, and
c 2. if its N(HI) is less than the specified upper limit for forest lines
      iflag=0
c      do j=1,numlls
c      zright=zlls(j) + dvavoid*(1.0+zlls(j))/vlight
c      zleft=zlls(j) - dvavoid*(1.0+zlls(j))/vlight
c      if(z.ge.zleft.and.z.le.zright)iflag=1
c      end do
      if(nhi.ge.nuplim)iflag=1
      
      if(iflag.eq.0)call spvoigt(da,wda,npts,dble(nhi),
     :                           dble(z),dble(b),'H ','I   ')

      end do
	write (6,*) 'Calculating number of lines using three power laws
     & (Kim et al. 2013)'
	  n=nint(n1+n2+n3)
	write (6,*) 'Total no. of lines from power laws = ',n
	write (6,*) 'Calculating number of lines in original QSOSIM9'
	  g=2.0
	  a13p75=10.0
	  beta=1.7
	  gp1=g+1.
	  a=a13p75*((nc)/(10**13.75))**(1.-beta)
	  rn=(a/gp1)*( z2p1**gp1 - z1p1**gp1 )
	  n=nint(rn)
	  write(6,*)' Total no. of lines from QSOSIM9 = ',n
   
c End of loop for forest.  Now put the LLS's in
c      do j=1,numlls
c        call spvoigt (da,wda,npts,dble(nhills(j)),dble(zlls(j)),
c     :                dble(blls(j)),'H ','I   ')
c      end do

c Keep unconvolved real*4 spectrum
      do i=1,npts
        da4(i)=real(da(i))
      end do
Convolution with assumed Gaussian instrumental profile
	   write(6,*) "Convolving..."
c	   call blur(da, npts, sigblur)
	   write (6,*) "Bluring complete..."
c Make real*4 array for pgplot
      do i=1,npts
        da4conv(i)=real(da(i))
      end do

c Make error array
      do i=1,npts
        da_err(i) = da(i)/dble(s2n) + 0.2*dble(danoabs(i))/dble(s2n)
        da_err4(i) = real(da_err(i))
      end do

c Add noise to real*4 convolved spectrum. 2 noise models.
c inoise=0 is constant.  inoise=1 gets worse towards blue. See qsosim9.pdf.
      if(inoise.eq.0)then
      do i=1,npts
        da_err4mod(i) = gasdev3(idum)*da_err4(i)
        da4smno(i) = da4conv(i) + da_err4mod(i)
      end do
      end if
      if(inoise.eq.1)then
      do i=1,npts
        w=(wda(i)-3532)/117
        ff=1.0+exp(-w)
        da_err4mod(i) = ff*gasdev3(idum)*da_err4(i)
        da4smno(i) = da4conv(i) + da_err4mod(i)
      end do
      end if
	
c Plot spectrum
      call PGBEGIN (0,'/vcps',1,3)
      xmin=wstart
      xmax=wend
      ymin=0.0
      ymax=0.0
      do i=1,npts
        if(da(i).gt.ymax)ymax=da(i)
      end do
      ymin=ymax
      do i=1,npts
        if(da(i).lt.ymin)ymin=da(i)
      end do
      ymax=ymax+0.5*ymax
      ymin=ymin-0.5*ymin

c Data to be returned into main program
	
 100	format(f10.5,2x,f10.5,2x,f10.5,2x,f10.5)
      do i=1,npts
	 lambda(i)=wda(i)
	 flux(i)=da4smno(i)!/danoabs(i)
	 flerr(i)=da_err4(i)!/danoabs(i)
	 nnflux(i)=da4conv(i)!/danoabs(i)
c         write(*,100)lambda(i),flux(i),flerr(i),nnflux(i)
      end do


      call PGENV (xmin,xmax,ymin,ymax,0,1)
      call PGLABEL ('Wavelength','f(lambda)','Linear wavelengths')
      call pgline(npts,wda,da4smno)
      call pgsci(5)
      call pgline(npts,wda,da4conv)
      call pgsls(2)
      call pgsci(2)
      call pgline(npts,wda,danoabs)
      call pgsls(1)
      call pgsci(3)
      call pgline(npts,wda,da_err4)

	call pgsci(1)
	!call PGENV(1.0,real(numlin),0.0,1.0,0,1)
	call PGENV(11.5,22.5,0.0,1.0,0,1)
      call PGLABEL('log NHI','CDDF(log NHI)','Column Density Dist')
	call PGLINE(npoints,xs,CDDF)
	call PGTEXT(12.6,0.8,'2')
	call PGTEXT(13.6,0.8,'1')
      call PGPTXT(16.0,0.5,0.0,1.0,'CDDF from Prochaska spline - 1')
	call PGSLS(2)
	call pgline(npoints,xs,dummyy)
      call PGPTXT(16.0,0.4,0.0,1.0,'CDDF from QSOSIM9 (original) - 2')

	call pgsls(1)
	call PGENV(1.0,real(numlin),11.5,22.5,0,1)
	call PGLABEL('Line number','log NHI','Random choice of NHI')
	call PGLINE(numlin,dummy,nhi4)
      call pgend

      return
      end

c Numercial Recipes subroutine
!======================================================================
      FUNCTION ran3(idum)
!======================================================================
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END

c Numercial Recipes subroutine
!======================================================================
      FUNCTION gasdev3(idum)
!======================================================================
      INTEGER idum
      REAL gasdev3
c Uses ran3
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran3
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran3(idum)-1.
        v2=2.*ran3(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev3=v2*fac
        iset=1
      else
        gasdev3=gset
        iset=0
      endif
      return
      END
!======================================================================
c BLUR does Gaussian blurring on array xbuf
	subroutine blur(xbuf,npt,sigma)
!======================================================================
	implicit none
	integer nfilt, npt, i, il, ilow, k
	real*8 xbuf(npt),ybuf(npt),work(512), sum
	real*8 xmns, xmnf, sigma, const, norm

	nfilt=int(10.0*sigma)+1
	if(nfilt.gt.511)stop ' too large a filter'
	if(npt.gt.262144)stop ' too many points in data array'
	if((nfilt/2)*2.eq.nfilt)nfilt=nfilt+1
c *** fill up blurring array
c	const=1.D0/(sqrt(8.D0*atan(1.D0))*sigma)
	const=1.0
	do i=1,nfilt
	   work(i)=const*exp(-0.5D0*(dble(i)-
	1	(dble(nfilt)+1.D0)/2.D0)**2.D0/sigma**2.D0)
	enddo
c *** set first and last edges equal
	il=nfilt/2
	ilow=max0(3,nfilt/4)
	ilow=(ilow/2)*2+1
	sum=0.D0
	do i=1,ilow
 	   sum=sum+xbuf(i)
	enddo
	xmns=sum/dble(ilow)
        sum=0.D0
	do i=1,ilow
 	   sum=sum+xbuf(npt+1-i)
	enddo
	xmnf=sum/dble(ilow)
c *** reflect edges before filtering
	do i=1,il
	   ybuf(i)=2.D0*xmns-xbuf(il+ilow+1-i)
 	   ybuf(npt+i+il)=2.D0*xmnf-xbuf(npt-i-ilow+1)
	enddo
	do i=1,npt
	   ybuf(i+il)=xbuf(i)
	enddo
c *** do blurring
	do k=1,npt
	   sum=0.D0
	   norm=0.D0
	   do i=1,nfilt
	      sum=sum+work(i)*ybuf(i+k-1)
	      norm=norm+work(i)
	   enddo
	   xbuf(k)=sum/norm
	enddo
	return
	end
