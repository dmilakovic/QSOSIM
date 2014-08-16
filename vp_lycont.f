      subroutine vp_lycont(cold,flx,wav,nw,zdp1)
*
	implicit none
*	include 'vp_sizes.f'
	double precision cold,zdp1
        integer i,na,nb,nlast,nlp,nup,nw
	real flx(nw), wav(nw)
*
	
	double precision fxx,fxy,temp1,temp2,tzero,x
	double precision wvd,dxd,dfd,dem,temp
*
*	Minimum wavelength Lyman line (for interpolation to Lyman cont.)
	double precision wlsmin,vblstar,collsmin
	common/vpc_lycont/wlsmin,vblstar,collsmin
*	control for verbose output:
	
c
c	subroutine to take care of Lyman limit absorption
c
	write(6,*) ' Lyman limit abs included', cold, zdp1-1D0

*	copy coeffts to single array
*	call vp_cfcopy(ichunk,1)
*	wavelengths in observers frame
	wvd=911.7536d0*zdp1
        temp = 1.0d20
        do i=1,nw
           if ( dabs(dble(wav(i)) - wvd) .lt. temp) then
              temp = dabs(dble(wav(i)) - wvd)
              nup = i
           endif
        enddo

*	Lyman limit channel (which may be negative!)
*	You may wish to include:	if(nup.le.0) return

	if(nup.gt.0) then
*	  set continuum absorption if range includes it:
	  tzero=cold*6.30d-18
	  do i=1,nup
	    x=dble(wav(i))/(zdp1*911.7536d0)
	    x=x*x*x
	    flx(i)=flx(i)*real(exp(-tzero*x))
	  end do
	end if
	if(nup.lt.nw) then
c	  set up interpolation limits
c	  find first channel at local minimum

*	  step 1: find channel corresponding to lowest wavelength hydrogen line
*
*	  rfc 27.6.95: new wavelength routines used
	  wvd=wlsmin*zdp1
          temp = 1.0d20
          do i=1,nw
             if ( dabs(dble(wav(i)) - wvd) .lt. temp) then
                temp = dabs(dble(wav(i)) - wvd)
                nlast = i
             endif
          enddo

*	  step 2: search up for local flux minimum in Lyman series spectrum
c	  start at int(ced)
	  if(nlast.gt.nw) then
*	    lowest line outside region, so just extrapolate the cubic
*	    since there is not much else you can do..
	    nlast=nw
	    x=dble(wav(nlast))/(zdp1*911.7536d0)
	    x=x*x*x
	    flx(nlast)=real(exp(-tzero*x))
	  end if
1	  nlp=nlast+1
	  if(nlp.gt.nw) then
*	    write(6,*) 'WARNING: High order Lyman series minimum flux'
*	    write(6,*) '         is outside wavelength range for spectral'
*	    write(6,*) '         region -- try expanding the region'
	    nlp=nw
	    goto 2
	  end if
	  if(flx(nlp).gt.flx(nlast)) goto 2
	  nlast=nlp
	  goto 1
2	  continue
c	  linearly interpolate in wavelength from nup to nlast
	  na=nup+1
	  nb=nlast-1
	  if(nb.lt.na) return
*	  interpolate linearly:
	  if(nup.gt.0) then
	    fxy=dble(flx(nup))
	   else
	    x=dble(wav(nup))/(zdp1*911.7536d0)
	    x=x*x*x
	    fxy=exp(-tzero*x)
	  end if
	  fxx=dble(flx(nlast))
*	  interpolation uses wavelength ratios, so no need for zdp1 correction
	  dfd=dble(wav(nup))
	  dxd=dble(wav(nlast))
	  dem=dxd-dfd
	  if(na.lt.1) then
	    na=1
	  end if    
	  do i=na,nb
	    wvd=dble(wav(i))
	    temp1=(wvd-dfd)/dem
	    temp2=(dxd-wvd)/dem
	    flx(i)=real(fxx*temp1+fxy*temp2)
	  end do
	end if
	return
	end
