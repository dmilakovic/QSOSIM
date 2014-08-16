      subroutine ewred(ind)
*	reads the atomic data from a file into common arrays
*	IN:
*	ind	integer	If zero, use environment variable to 
*			determine data source, if non-zero, prompt for
*			filename.
*
*	include 'vp_sizes.f'
*
	integer ind
	logical loldstyle
	character*64 filnm
        character*200 atomdr
*
        integer i, lastchpos
	double precision tmass
*       character input and separation variables
        character*60 cvstr(15)
        real rvstr(15)
        integer ivstr(15)
        character*132 inchstr
        common/sepspace/inchstr,rvstr,ivstr,cvstr,nvstr
*
        integer numel
	character*2 lbz
	character*4 lzz
	double precision alm,fik,asm 
	common/vpc_ewllns/lbz(2500),lzz(2500),alm(2500),
     :                fik(2500),asm(2500),numel
*	atomic mass table
	character*2 lbm
	double precision amass
	common/vpc_atmass/lbm(2500),amass(2500),mmass
*	Minimum wavelength Lyman line (for interpolation to Lyman cont.)
	double precision wlsmin,vblstar,collsmin
	common/vpc_lycont/wlsmin,vblstar,collsmin
* 	get directory for atomic data from logical ATOMDIR or use bob's
*	
	if(ind.gt.0) goto 990
        atomdr=' '
        call getenv ( 'ATOMDIR', atomdr )
        if ( atomdr .eq. ' ' ) then
           atomdr = './atom.dat'
        end if
	open(unit=18,file=atomdr,status='old',err=99)
	goto 11
99	i=lastchpos(atomdr)
        if(atomdr(i:i).eq.'/') then
          atomdr=atomdr(1:i)//'atom.dat'
         else
          atomdr=atomdr(1:i)//'/atom.dat'
        end if
        open(unit=18,file=atomdr,status='old',err=990)
        goto 11
990	write(6,*) ' Atomic data filename'
        read(5,'(a64)',end=98) filnm
        if(filnm(1:4).eq.'    ') then
          write(6,*) ' Uses Ly-a only'
          numel=1
          lbz(1)='H '
          lzz(1)='I   '
          alm(1)=1215.67d0
          fik(1)=0.4162d0
          asm(1)=6.265d8
          goto 98
        end if
        open(unit=18,file=filnm,status='old',err=99)
c	print name of file
11      continue
        inquire( 18, name=filnm )
        write ( *, * ) ' Using data from : '//filnm
  	numel=0
	mmass=0
	loldstyle=.true.
1	numel=numel+1
	tmass=0.0d0
	read(18,'(a)',end=97) inchstr
	if(loldstyle) then
           read(inchstr,100,iostat=ierr) lbz(numel),lzz(numel),
     1          alm(numel),fik(numel),asm(numel)
 100       format(a2,a4,f8.3,f8.6,e7.3)
*     if first pass, check if it is old style input
           if(numel.eq.1) then
*     fixed format error => not fixed format!
              if(ierr.ne.0) then
                 loldstyle=.false.
              else
                 if(lzz(1).ne.'I'.and.lzz(1).ne.'V'.and.
     1                lzz(1).ne.'X'.and.lbz(1)(1:1).ne.'>') then
                    loldstyle=.false.
                 else
                    if(alm(1).le.0.0d0.or.fik(1).le.0.0d0.or.
     1                   asm(1).le.0.0d0) then
                       loldstyle=.false.
                    end if
                 end if
              end if
              if(loldstyle) then
                 write(6,*) ' Using old style atomic data file'
                 write(6,*) ' '
              else
*     separate inchstr string (in common) to get variables
               call vp_atomsep(lbz(numel),lzz(numel),alm(numel),
     :                fik(numel),asm(numel),tmass)
                end if
           end if
        else
*	  free format atomic data table, with atomic masses
*	  separate inchstr string (in common) to get variables
           call vp_atomsep(lbz(numel),lzz(numel),alm(numel),
     :          fik(numel),asm(numel),tmass)
  	end if
*
*	set up mass table if a mass is there
	if(tmass.gt.0.0d0) then
*	  check have not already got a mass
           if(mmass.gt.0) then
              i=1
              do while (lbz(numel).ne.lbm(i).and.i.le.mmass)
                 i=i+1
              end do
*	    i>mmass => not already done
              if(i.gt.mmass) then
                 mmass=mmass+1
                 lbm(mmass)=lbz(numel)
                 amass(mmass)=tmass
              end if
	   else
              mmass=mmass+1
              lbm(mmass)=lbz(numel)
              amass(mmass)=tmass
           end if
	end if
*
       
	if(alm(numel).gt.0.0d0.and.numel.lt.2500) goto 1
	if(alm(numel).gt.0.0d0) then
*	  must have run out of array space
           write(6,'('' WARNING: ATOMIC DATA TABLE FULL'')')
           write(6,'(''  last record: '',a2,a4,f9.3,f10.6,1pe10.2)')
     1          lbz(numel),lzz(numel),alm(numel),fik(numel),asm(numel)
*	  add 1 to m so it can be subtracted later.
           numel=numel+1
	end if
 97     numel=numel-1
	rewind 18
	close(unit=18)
c	find shortest wavelength in Lyman series used
c	for interpolation to Lyman limit if needed
 98     wlsmin=1d20
	if(numel.ge.1) then
           do j=1,numel
              if(lbz(j).eq.'H ') then
                 if(alm(j).lt.wlsmin.and.asm(j).gt.2.0D0) wlsmin=alm(j)
              end if
           end do
	end if
	return
	end

      subroutine vp_atomsep(lbz,lzz,alm,fik,asm,amass)
	character*2 lbz
	character*4 lzz
	double precision alm,fik,asm,amass
	integer nvstr
        character*60 cvstr(15)
        real rvstr(15)
        integer ivstr(15)
        character*132 inchstr
        common/sepspace/inchstr,rvstr,ivstr,cvstr,nvstr
	double precision dvstr(15)
	common/vpc_dsepspace/dvstr


	call dsepvar(inchstr,6,dvstr,ivstr,cvstr,nvstr)
*	WAS: if(cvstr(2)(1:1).eq.'I'.or.cvstr(2)(1:1).eq.'V'.or.
*     :     cvstr(2)(1:1).eq.'X'.or.cvstr(2)(1:1).eq.'J') then
*	Replaced by any capital letter:
	  if(ichar(cvstr(2)(1:1)).ge.ichar('A').and.
     :       ichar(cvstr(2)(1:1)).le.ichar('Z')) then
	  lbz=cvstr(1)(1:2)
	  lzz=cvstr(2)(1:4)
	  alm=dvstr(3)
	  fik=dvstr(4)
	  asm=dvstr(5)
	  amass=dvstr(6)
	 else
*	  element/ion is in first character string, so separate it
*	  WAS: if(cvstr(1)(2:2).eq.'I'.or.cvstr(1)(2:2).eq.'V'.or.
*     :       cvstr(1)(2:2).eq.'X'.or.cvstr(1)(1:1).eq.'J') then
*	  Replaced by any CAPITAL letter
	  if(ichar(cvstr(1)(2:2)).ge.ichar('A').and.
     :       ichar(cvstr(1)(2:2)).le.ichar('Z')) then
	    lbz=cvstr(1)(1:1)//' '
	    lzz=cvstr(1)(2:5)
	   else
	    lbz=cvstr(1)(1:2)
	    lzz=cvstr(1)(3:6)
	  end if
	  alm=dvstr(2)
	  fik=dvstr(3)
	  asm=dvstr(4)
	  amass=dvstr(5)
	end if
	return
	end

      integer function lastchpos(chstr)
*	finds the position of the last non-blank character in
*	the character string chstr
      character*(*) chstr
      lastchpos=len(chstr)
      do while(lastchpos.gt.1.and.
     :     chstr(lastchpos:lastchpos).eq.' ')
         lastchpos=lastchpos-1
      end do
      if(lastchpos.eq.1.and.chstr(1:1).eq.' ') then
         lastchpos=0
      end if
      return
      end
