c ----------------------------------------------------------------------------
c QSOSIM9. Prog. to make synthetic QSO spectrum - John Webb, UNSW, Dec 2013
c Uses VPFIT Voigt profile programs
c Compile using makeqsosim
c ----------------------------------------------------------------------------
	  real*4 wda(262144),danoabs(262144),tau(262144)
	  real*4 da4(262144),da4conv(262144)
	  real*4 da_err4(262144),da4smno(262144)
	  real*8 da(262144),da_err(262144),da_err4mod(262144)
	  real*4 wems(30),relstr(30)
	  real*4 sum,nhi,b,z,g,z1p1,z2p1,beta,x,gp1,w,ff
	  real*4 a12p5,nc,nuplim,rn,a13p75,mbp1,zqsop1,pi
	  real*4 c,d,p,q,r,s,s2n,dvavoid,vlight,zleft,zright
	  real*4 nhills(20),blls(20),zlls(20)
          real*8 alm, fik,asm,sigblur
	  integer idum,iflag,inoise,i
	  integer,parameter :: nrows=2
	  data pi/3.14159265/
	  data vlight/299792.458/
	  CHARACTER*20 :: filename
	  INTEGER,DIMENSION(nrows) :: numlls
          REAL,DIMENSION(nrows) :: ra,dec,zqso,alpha,vmag,s2n
     &                                     sigblur
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
	  nptsmax=262144
	  nemlines=30

c Read in simulation parameters
	  filename='sin.fits'
	  write(*,400)'reading file:',filename
	  write(*,*)'---------------------------------'
	  write(*,*)ra
	  call readfits(filename, wstart,wend,dw,nc,nuplim,inoise,dvavoid,
     +                 ra,dec,zqso,alpha,vmag,sigblur,numlls)
	  write(*,*)'---------------------------------'
	  do i=1,nrows
	     write(*,500)i,ra(i)
	  enddo
 400	  format(a20,2x,a20)!,2x,f10.3,2x,f10.3,f4.2)
 500	  format(i2,2x,f10.3)
c$$$c Read in any specified LLS
c$$$      do i=1,numlls
c$$$        read(17,*)nhills(i),blls(i),zlls(i)
c$$$c      write(6,*)nhills(i),blls(i),zlls(i)
c$$$      end do
c$$$c alpha is QSO spectral index, vmag is the V magnitude
c$$$c wstart, wend, dw are the start and end wavelengths and pixel size
c$$$c nc is the column density cut-off (REAL, NOT log10)
c$$$
c$$$	  zqsop1=zqso+1
c$$$	  f5550=10**(-(vmag+21.17)/2.5)
c$$$	  const=f5550/(1.0/5550.0**(2+alpha))
c$$$      sum=0.0
c$$$      npts=(wend-wstart)/dw
c$$$      if(npts.gt.nptsmax)then
c$$$      	write(6,*)' Max. no. of pixels = ',nptsmax
c$$$      	stop
c$$$      end if
c$$$
c$$$	  do l=1,npts
c$$$c Make wavelength scale and underlying power-law QSO continuum
c$$$	     wda(l)=wstart+((l-1)*dw)
c$$$	     da(l)=const*(1.0/wda(l)**(2+alpha))
c$$$	  end do
c$$$
c$$$c Put emission lines in. Guess/approximation for emission line width:
c$$$	  fwhm=60.0
c$$$	  sigma=fwhm/2.35
c$$$
c$$$	  do m=1,nemlines
c$$$        wems(m)=wems(m)*zqsop1
c$$$		 do i=1,npts
c$$$           x=0.75*relstr(m)*const*(1.0/wems(m)**(2+alpha))
c$$$           g=x*exp(-.5*(((wda(i)-wems(m))/sigma)**2))
c$$$           da(i)=da(i)+g
c$$$         end do
c$$$      end do
c$$$
c$$$c Keep unabsorbed QSO spectrum
c$$$	  do i=1,npts
c$$$	    danoabs(i)=da(i)
c$$$	  end do
c$$$
c$$$c Generate optical depth array in preparation for absorption input
c$$$	  do i=1,npts
c$$$	    tau(i)=0.0
c$$$	  end do
c$$$
c$$$c Next section makes absorption lines
c$$$c g is from dn=A(1+z)^g dz, beta is from dn propto N^{-beta} dN.
c$$$c a13p75 is the value of A for lines above logN=13.75.
c$$$c gp1= gamma+1, mbp1=1-beta, nc is N cutoff, n is total no. of lines
c$$$      g=2.0
c$$$      a13p75=10.0
c$$$      beta=1.7
c$$$      gp1=g+1.
c$$$      z1p1=(wstart/1215.67)
c$$$      z2p1=zqsop1
c$$$      mbp1=-beta+1.0
c$$$
c$$$c Calculate the total no. of lines
c$$$      a=a13p75*((nc)/(10**13.75))**(1.-beta)
c$$$      rn=(a/gp1)*( z2p1**gp1 - z1p1**gp1 )
c$$$      n=nint(rn)
c$$$      write(6,*)' Total no. of lines = ',n
c$$$
c$$$c Initialise random numbers
c$$$      idum=time()
c$$$
c$$$c Call Voigt profile generator n times, once for each abs system
c$$$	  do i=1,n
c$$$
c$$$c Random selection of redshifts
c$$$        c=ran3(idum)
c$$$        p=z2p1**gp1
c$$$        q=z1p1**gp1
c$$$        x=alog10(c*(p-q)+q)
c$$$        x=x/gp1
c$$$        z=10**x -1.0
c$$$        
c$$$c Random selection of column densities
c$$$ 1      d=ran3(idum)
c$$$        if(d.lt.1.0) then
c$$$          y1=((alog10(1.0-d)) / mbp1) + alog10(nc)
c$$$        else
c$$$          goto 1
c$$$        end if
c$$$        nhi=10**y1
c$$$
c$$$c b-params.  Guess at sigma and mean of b distribution of 3 and 23 km/s.
c$$$        b = 3*gasdev3(idum)+23
c$$$
c$$$c Put the forest lines in, but only if:
c$$$c 1. it is not if within the specified avoidance zone of each LLS, and
c$$$c 2. if its N(HI) is less than the specified upper limit for forest lines
c$$$      iflag=0
c$$$      do j=1,numlls
c$$$      zright=zlls(j) + dvavoid*(1.0+zlls(j))/vlight
c$$$      zleft=zlls(j) - dvavoid*(1.0+zlls(j))/vlight
c$$$      if(z.ge.zleft.and.z.le.zright)iflag=1
c$$$      end do
c$$$      if(nhi.ge.nuplim)iflag=1
c$$$      
c$$$      if(iflag.eq.0)call spvoigt(da,wda,npts,dble(nhi),
c$$$     :                           dble(z),dble(b),'H ','I   ')
c$$$
c$$$      end do
c$$$c End of loop for forest.  Now put the LLS's in
c$$$      do j=1,numlls
c$$$        call spvoigt (da,wda,npts,dble(nhills(j)),dble(zlls(j)),
c$$$     :                dble(blls(j)),'H ','I   ')
c$$$      end do
c$$$
c$$$c Keep unconvolved real*4 spectrum
c$$$      do i=1,npts
c$$$        da4(i)=real(da(i))
c$$$      end do
c$$$Convolution with assumed Gaussian instrumental profile
c$$$	   write(6,*) "Convolving..."
c$$$	   call blur(da, npts, sigblur)
c$$$
c$$$c Make real*4 array for pgplot
c$$$      do i=1,npts
c$$$        da4conv(i)=real(da(i))
c$$$      end do
c$$$
c$$$c Make error array
c$$$      do i=1,npts
c$$$        da_err(i) = da(i)/dble(s2n) + 0.2*dble(danoabs(i))/dble(s2n)
c$$$        da_err4(i) = real(da_err(i))
c$$$      end do
c$$$
c$$$c Add noise to real*4 convolved spectrum. 2 noise models.
c$$$c inoise=0 is constant.  inoise=1 gets worse towards blue. See qsosim9.pdf.
c$$$      if(inoise.eq.0)then
c$$$      do i=1,npts
c$$$        da_err4mod(i) = gasdev3(idum)*da_err4(i)
c$$$        da4smno(i) = da4conv(i) + da_err4mod(i)
c$$$      end do
c$$$      end if
c$$$      if(inoise.eq.1)then
c$$$      do i=1,npts
c$$$        w=(wda(i)-3532)/117
c$$$        ff=1.0+exp(-w)
c$$$        da_err4mod(i) = ff*gasdev3(idum)*da_err4(i)
c$$$        da4smno(i) = da4conv(i) + da_err4mod(i)
c$$$      end do
c$$$      end if
c$$$
c$$$c Plot spectrum
c$$$      call PGBEGIN (0,'?',1,1)
c$$$      xmin=wstart
c$$$      xmax=wend
c$$$      ymin=0.0
c$$$      ymax=0.0
c$$$      do i=1,npts
c$$$        if(da(i).gt.ymax)ymax=da(i)
c$$$      end do
c$$$      ymin=ymax
c$$$      do i=1,npts
c$$$        if(da(i).lt.ymin)ymin=da(i)
c$$$      end do
c$$$      ymax=ymax+0.5*ymax
c$$$      ymin=ymin-0.5*ymin
c$$$
c$$$c Write out normalised spectrum to ascii file
c$$$      open(unit=18,file='spec.dat',status='new')
c$$$      do i=1,npts
c$$$      write(18,*)wda(i),da4smno(i)/danoabs(i),da_err4(i)/danoabs(i),
c$$$     :           da4conv(i)/danoabs(i)
c$$$      end do
c$$$
c$$$      call PGENV (xmin,xmax,ymin,ymax,0,0)
c$$$      call PGLABEL ('Wavelength','f(lambda)','Linear wavelengths')
c$$$      call pgline(npts,wda,da4smno)
c$$$      call pgsci(5)
c$$$      call pgline(npts,wda,da4conv)
c$$$      call pgsls(2)
c$$$      call pgsci(2)
c$$$      call pgline(npts,wda,danoabs)
c$$$      call pgsls(1)
c$$$      call pgsci(3)
c$$$      call pgline(npts,wda,da_err4)
c$$$      call pgend
c$$$
c$$$      stop
c$$$      end
c$$$
c$$$c Numercial Recipes subroutine
c$$$      FUNCTION ran3(idum)
c$$$      INTEGER idum
c$$$      INTEGER MBIG,MSEED,MZ
c$$$C     REAL MBIG,MSEED,MZ
c$$$      REAL ran3,FAC
c$$$      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
c$$$C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
c$$$      INTEGER i,iff,ii,inext,inextp,k
c$$$      INTEGER mj,mk,ma(55)
c$$$C     REAL mj,mk,ma(55)
c$$$      SAVE iff,inext,inextp,ma
c$$$      DATA iff /0/
c$$$      if(idum.lt.0.or.iff.eq.0)then
c$$$        iff=1
c$$$        mj=MSEED-iabs(idum)
c$$$        mj=mod(mj,MBIG)
c$$$        ma(55)=mj
c$$$        mk=1
c$$$        do 11 i=1,54
c$$$          ii=mod(21*i,55)
c$$$          ma(ii)=mk
c$$$          mk=mj-mk
c$$$          if(mk.lt.MZ)mk=mk+MBIG
c$$$          mj=ma(ii)
c$$$11      continue
c$$$        do 13 k=1,4
c$$$          do 12 i=1,55
c$$$            ma(i)=ma(i)-ma(1+mod(i+30,55))
c$$$            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
c$$$12        continue
c$$$13      continue
c$$$        inext=0
c$$$        inextp=31
c$$$        idum=1
c$$$      endif
c$$$      inext=inext+1
c$$$      if(inext.eq.56)inext=1
c$$$      inextp=inextp+1
c$$$      if(inextp.eq.56)inextp=1
c$$$      mj=ma(inext)-ma(inextp)
c$$$      if(mj.lt.MZ)mj=mj+MBIG
c$$$      ma(inext)=mj
c$$$      ran3=mj*FAC
c$$$      return
c$$$      END
c$$$
c$$$c Numercial Recipes subroutine
c$$$      FUNCTION gasdev3(idum)
c$$$      INTEGER idum
c$$$      REAL gasdev3
c$$$c Uses ran3
c$$$      INTEGER iset
c$$$      REAL fac,gset,rsq,v1,v2,ran3
c$$$      SAVE iset,gset
c$$$      DATA iset/0/
c$$$      if (iset.eq.0) then
c$$$1       v1=2.*ran3(idum)-1.
c$$$        v2=2.*ran3(idum)-1.
c$$$        rsq=v1**2+v2**2
c$$$        if(rsq.ge.1..or.rsq.eq.0.)goto 1
c$$$        fac=sqrt(-2.*log(rsq)/rsq)
c$$$        gset=v1*fac
c$$$        gasdev3=v2*fac
c$$$        iset=1
c$$$      else
c$$$        gasdev3=gset
c$$$        iset=0
c$$$      endif
c$$$      return
c$$$      END
c$$$
c$$$c BLUR does Gaussian blurring on array xbuf
c$$$	subroutine blur(xbuf,npt,sigma)
c$$$	implicit none
c$$$	integer nfilt, npt, i, il, ilow, k
c$$$	real*8 xbuf(npt),ybuf(npt),work(512), sum
c$$$	real*8 xmns, xmnf, sigma, const, norm
c$$$
c$$$	nfilt=int(10.0*sigma)+1
c$$$	if(nfilt.gt.511)stop ' too large a filter'
c$$$	if(npt.gt.262144)stop ' too many points in data array'
c$$$	if((nfilt/2)*2.eq.nfilt)nfilt=nfilt+1
c$$$c *** fill up blurring array
c$$$c	const=1.D0/(sqrt(8.D0*atan(1.D0))*sigma)
c$$$	const=1.0
c$$$	do i=1,nfilt
c$$$	   work(i)=const*exp(-0.5D0*(dble(i)-
c$$$	1	(dble(nfilt)+1.D0)/2.D0)**2.D0/sigma**2.D0)
c$$$	enddo
c$$$c *** set first and last edges equal
c$$$	il=nfilt/2
c$$$	ilow=max0(3,nfilt/4)
c$$$	ilow=(ilow/2)*2+1
c$$$	sum=0.D0
c$$$	do i=1,ilow
c$$$ 	   sum=sum+xbuf(i)
c$$$	enddo
c$$$	xmns=sum/dble(ilow)
c$$$        sum=0.D0
c$$$	do i=1,ilow
c$$$ 	   sum=sum+xbuf(npt+1-i)
c$$$	enddo
c$$$	xmnf=sum/dble(ilow)
c$$$c *** reflect edges before filtering
c$$$	do i=1,il
c$$$	   ybuf(i)=2.D0*xmns-xbuf(il+ilow+1-i)
c$$$ 	   ybuf(npt+i+il)=2.D0*xmnf-xbuf(npt-i-ilow+1)
c$$$	enddo
c$$$	do i=1,npt
c$$$	   ybuf(i+il)=xbuf(i)
c$$$	enddo
c$$$c *** do blurring
c$$$	do k=1,npt
c$$$	   sum=0.D0
c$$$	   norm=0.D0
c$$$	   do i=1,nfilt
c$$$	      sum=sum+work(i)*ybuf(i+k-1)
c$$$	      norm=norm+work(i)
c$$$	   enddo
c$$$	   xbuf(k)=sum/norm
c$$$	enddo
	return
	end
