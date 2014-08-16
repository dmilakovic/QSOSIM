      subroutine spvoigt(contin,wav,npts,col,z,b,atom,ion)
*     PVOIGT  computes Voigt profile
*
*     Input:
*     flx(npts)         continuum(emission) or anything(absn)
*     contin(npts)	continuum level
*     npts	integer	array size for data
*     col	dble	column density (may be log) (and scaled!!)
*     z	dble	redshift
*     b	dble	Doppler parameter (may be log) (in km/s?)
*     atom  ch*2	element
*     ion  ch*4	ion
*     ichunk		chunk number
*     Output:
*     flx(npts)		absorption on continuum
*
*      include 'vp_sizes.f'
*     Wavelength scale need not be linear
      real flx(npts),wav(npts), expf
      double precision col,z,b, contin(npts)
*     
*     LOCAL variables
      double precision v,a,cne
      double precision zdp1,vely,veld,binvel
      double precision bl,cns,dxd,wlod,whid
      double precision wvd
      double precision wavlow,wavhigh,vxbl
      double precision dtemp,wtemp,temp
      double precision xtxx,cold
      double precision ww, wv
*     real variables for working with flux values:
      real binlength,e,s,subtau,sumtot
*     
*     FUNCTIONS:
      double precision voigt
*

*     ion is the 'level' i.e. ionization stage


      integer nexpd, chand, i, jmid, j, jstep
   
      integer numel 
      character*2 atom,lbz	
      character*4 ion,lzz
*     nhmc: These are atom.dat things, keep for now...
      double precision alm,fik,asm 
      common/vpc_ewllns/lbz(2500),lzz(2500),alm(2500),
     :     fik(2500),asm(2500),numel
*	parameters for smoothing
*	indvar=1 for log N, 0 linear
      integer indvar

*	Lyman series minimum wavelength in table
      double precision wlsmin,vblstar,collsmin
      	common/vpc_lycont/wlsmin,vblstar,collsmin
*     Column density for inclusion of ion in ALL regions
      double precision fcollallzn
      double precision clvaldrop

c      write(*,*)"   col, z, b", col, z, b

*     set column densities used to be linear
      indvar=0
*     Smallest column density used (log)      
      clvaldrop=8.0
*     No idea
      vblstar=3.0
*     minimum col den for LL inclusion
      collsmin=1.0d15   
*     minimum line strength for inclusion
      fcollallzn=4.162d18
*
*
*	set velocity parameter for internal use
*	  straight wavelengths used in program
      vely=b
*	antilog col if necessary, but not if added component
      if(indvar.eq.1.) then
         xtxx=col
*	  mim@phys.unsw.edu.au: 17.10.01 Old default of 10.0
*	  replaced by clvaldrop
         if(xtxx.gt.34.0.or.xtxx.lt.clvaldrop) then
            xtxx=clvaldrop
*           col=clvaldrop*scalelog ! not used subsequently
         end if
         cold=10.0d0**xtxx
      else
         cold=col
      end if
      zdp1=1.0d0+z
*	set up parameters
      nexpd=5
*

      do i=1,npts
*	    rfc 24.2.95: replace by unit continuum and multiply later,
*	    so that can interpolate from continuous absorption to
*	    highest series line:
         flx(i)=1.0
      end do
 
*
*     nhmc: ewred reads in atom.dat.  m is the number of transitions listed
*     in atom.dat
      if(numel.le.0) call ewred(0)
*	smallest rest wavelength for later use
      wavlow=dble(wav(1))/zdp1
*	upper rest wavelength:
      wavhigh=dble(wav(npts))/zdp1
      if(vely.eq.0.0d0.or.vely.eq.veld) goto 302
      veld=vely
*     
 302  vxbl=dble(wav(npts/2))/zdp1
*	want at least nexpd (default 5) b-parameters per sub-bin
      binvel=2.997924562d5*(wavhigh-wavlow)/(vxbl*dble(npts-1))
*	velocity units cm/s internal, while km/s used!
*      vxbl=veld*vxbl/1.0d5  N: this doesn't seem to be used again, so
*                               commenting out.
*
*	now find the lines and slot them in
*
      do k = 1, numel
*         ajc 8-apr-92  include high column density lines always
*	  [rfc 18-Nov-97: Was as well as continuum altering thing - for 
*         that order (otherwise, this routine not called) - but now uses
*	  redshift as well]
         if ( ( atom .eq. lbz( k ) .and. ion .eq. lzz( k ) ) .and.
     :        ( (alm(k).ge.wavlow.and.alm(k).le.wavhigh) .or.
     :        cold*fik(k).gt. fcollallzn )) then 
c            WRITE(*,*) "Match!", lbz(k), lzz(k)
*         for each bin, starting from the line center
	    wvd=alm(k)*zdp1
*     Need to find array number with wavelength closest to wavelength wvd...
            temp = 1.0d20
            do i=1,npts
               if ( dabs(dble(wav(i)) - wvd) .lt. temp) then
                  temp = dabs(dble(wav(i)) - wvd)
                  chand = i
               endif
            enddo
c            write(*,*) "Line in channel no.",chand

	    jmid=chand
	    tauchan=1.0e25
	    jstep=-1
	    j=jmid
*
            do while (tauchan.gt.1.0e-6.and.j.le.npts)
               if(jstep.eq.-1) then
                  if(j.gt.1) then
                     j=j-1
                  else
                     tauchan=1.0e-7
                     if(jmid.le.0) jmid=1
                     goto 901
                  end if
	       else
                  if(j.lt.npts) then
                     j=j+1
                  else
                     tauchan=1.0e-7
                     goto 901
                  end if
               end if
               stau=10.0
*

*               only want to be here if normal ion 
*               line parameters
               wv=alm(k)*1.0d-8
               bl=veld*wv/2.997924562d5
               a=asm(k)*wv*wv/(3.7699d11*bl)
               cns=wv*wv*fik(k)/(bl*2.00213d12)
               cne=cold*cns
*                set bin length so that nexpd Doppler parameters per bin
               dtemp=veld/binvel
               if(dtemp.gt.dble(nexpd)) then
                  binlength=1.0
                  nexbin=1
               else
                  nexbin=int(dble(nexpd)*binvel/veld+0.5d0)
                  if(nexbin/2*2.eq.nexbin) nexbin=nexbin+1
                  binlength=1.0/real(nexbin)
               end if
           
*               high and low wavelength points (interpolate linearly)
               wlod=(dble(wav(j))+dble(wav(j-1)))/2.0 * 1.0d-8/zdp1
               whid=(dble(wav(j))+dble(wav(j+1)))/2.0 * 1.0d-8/zdp1
*               sub-divide
               dxd=(whid-wlod)/dble(nexbin)
               sumtot=0.0
               stau=0.0
               do klv=1,nexbin
                  fluxcur=flx(j)
                  ww=wlod+(dble(klv)-0.5d0)*dxd
                  v=wv*wv*(1.0d0/ww-1.0d0/wv)/bl
*                   absorb...
*                   this seems to be voigt profile...
                  subtau=real(cne*voigt( v, a ))
                  e = expf(-subtau)
                  if(e.gt.1e-30.and.fluxcur.gt.1e-30)then
                     if ( log10( e ) + log10( fluxcur ) 
     :                    .lt. -30.0 ) then
                        s = 0.0
                     else
                        s=fluxcur*e
                     end if 
                  else
                     s = 0.0
                  end if
                     
*                 end of v.p.
*	          add as simple sub-chunks
                  stau=subtau+stau
                  sumtot=sumtot+s
               end do
               flx(j)=sumtot*binlength
               stau=stau*binlength
               
               tauchan=stau
 901           continue
               if(jstep.lt.0.and.(tauchan.lt.1.0e-6.or.j.le.1)) then
                  jstep=1
                  j=jmid-1
                  tauchan=1.0e25
               end if
	    end do
         end if
      end do
*
*	rfc 24.2.95:
*
      
*	  add Lyman continuum absorption if hydrogen, and if in the
*	  right wavelength region:
      if(atom.eq.'H '.and.ion.eq.'I   ') then
*	    hydrogen, so check if need Lyman limit:
*	    wavelength limit is set by Doppler parameter (vblstar=3 default)
         wtemp=wlsmin*(1.0d0+vblstar/2.997924562d5)
*	    column density limit is logN(HI)=15.0 by default 
*	    (i.e. collsmin=1.0d15)
         if(wavlow.le.wtemp.and.dble(cold).gt.collsmin) then
            call vp_lycont(cold,flx,wav,npts,zdp1)
*              if(verbose) write(6,*) 'vp_lycont exit OK'
         end if
      end if
*
*     multiply by the continuum which came in, to get result:
      do i=1,npts
c         write(*,*) "flx, contin", flx(i), contin(i)
         contin(i)=contin(i)*dble(flx(i))
      end do
     
*     
      return
      end


      function expf(x)
      if(x.lt.-80.0) goto 1
      expf=exp(x)
 2    return
 1    expf=0.0
      goto 2
      end
