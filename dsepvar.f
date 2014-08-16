      SUBROUTINE DSEPVAR(CHST,NVAR,RV,IV,CV,NV)
c
c	Separate a string of variables, separated by blanks, commas or &,
c	into variable arrays. ' or " are used to delineate character
c	variables which contain any special characters (space, comma,
c	' or "), so ' and " cannot appear in the same string. Successive
c	commas are treated as zero or blank, as appropriate. Variables
c	which do not translate to real or integer are set to zero.
c	Note - spaces after a character before commas are treated
c	as another separator, so
c		1 ,,2     becomes 1  0  0  2
c	and
c	        1,, 2        "    1  0  2
c
c	input:	chst	character variable containing data
c		nvar	number of variables expected
c
c	output:	cv	character variables
c		rv	double precision variables
c		iv	integer variables
c		nv	actual number of variables found
c
c
	character*132 cx
	character*(*) cv(*),chst
	character*1 chdg,ctab
	character*3 cdum
	dimension iv(*)
	double precision rv(*)
	double precision xxx
*	tab character
	ctab='	'
c	single quote fix
	cdum='x''x'
	chdg=cdum(2:2)
c
c
c	string lengths
	lstr=len(chst)
	if(lstr.le.0) then
c	  no data at all
	  cv(1)=' '
	  rv(1)=0.0d0
	  iv(1)=0
	  nv=0
	  return
	end if
c	lv=len(cv(1))
c
	kk=1
c
c	find first non-blank character
	k1=1
1	continue
	do while (chst(k1:k1).eq.' ')
	  if(k1.ge.lstr) goto 3
	  k1=k1+1
	end do
c
c	check that it is not a special character
	if(chst(k1:k1).eq.chdg.or.chst(k1:k1).eq.'"') then
c
c	  character string in quotes - search for matching quote
	  k2=k1+1
	  do while (chst(k2:k2).ne.chst(k1:k1).and.
     1       k2.le.lstr)
	    k2=k2+1
	  end do
	  if(k2.le.k1+1) then
	    cv(kk)=' '
	   else
	    cv(kk)=chst(k1+1:k2-1)
	  end if
	  rv(kk)=0.0d0
	  iv(kk)=0
	  k1=k2+1
	  kk=kk+1
c	  if separator added, go one more step
	  if(chst(k1:k1).eq.','.or.chst(k1:k1).eq.'&'.or.
     :       chst(k1:k1).eq.ctab) k1=k1+1
	 else
	  if(chst(k1:k1).eq.','.or.chst(k1:k1).eq.'&'.or.
     :       chst(k1:k1).eq.ctab) then
c	    comma separator, straight after previous separator
	    cv(kk)=' '
	    rv(kk)=0.0d0
	    iv(kk)=0
	    kk=kk+1
	    k1=k1+1
	   else
c
c	    some other character - search for the next 
c	    space or comma or '&'
	    k2=k1+1
	    do while (chst(k2:k2).ne.' '.and.chst(k2:k2).ne.','
     1       .and.chst(k2:k2).ne.'&'.and.chst(k2:k2).ne.ctab
     :       .and.k2.le.lstr)
	      k2=k2+1
	    end do
	    k2m=k2-1
	    cv(kk)=chst(k1:k2m)
	    kvar=k2m-k1
c	    check if a reasonable integer
	    if(kvar.gt.9) goto 91
	    cx=' '
	    if(kvar.lt.20) then
	      cx(20-kvar:20)=chst(k1:k2m)
	    end if
c	    read as an integer
	    read(cx,'(i20)',err=91)  iv(kk)
	    rv(kk)=iv(kk)
	    goto 92
c
c	    integer failed - try as a real
91	    cx=' '
	    if(kvar.lt.30) then
	      cx(30-kvar:30)=chst(k1:k2m)
	    end if
	    read(cx,'(f30.0)',err=93) rv(kk)
	    xxx=rv(kk)
	    if(xxx.lt.0.0) xxx=-xxx
	    if(xxx.le.2.0d9) then
	      iv(kk)=int(rv(kk))
	     else
	      iv(kk)=0
	    end if
	    goto 92
c
c	    character string only thing that works
93	    rv(kk)=0.0d0
	    iv(kk)=0
92	    k1=k2+1
	    kk=kk+1
	  end if
	end if
	if(k1.le.lstr.and.kk.le.nvar) goto 1
c
c	end of data
3	nv=kk-1
c
c	zero or blank off nv+1 to nvar
	if(nv.lt.nvar) then
	  do i=kk,nvar
	    cv(i)=' '
	    rv(i)=0.0d0
	    iv(i)=0
	  end do
	end if
	return
	end
