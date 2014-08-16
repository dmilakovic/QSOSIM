      subroutine f13_read(inchstr,n,b,z,atom,ion)
      implicit none
      character*132 inchstr
      real*4 n,b,z
      character*2 atom
      character*4 ion
c      integer i
      character*60 cvstr(15)
      real*8 dvstr(15)
      integer ivstr(15), nvstr
      
c      write(*,*) "inchstr inside f13_read |", inchstr,"|"
      
      call dsepvar(inchstr,5,dvstr,ivstr,cvstr,nvstr)
      
c      do i=1,nvstr
c         write(*,*)"cvstr",i," is |", cvstr(i),"|"
c         write(*,*)"dvstr",i," is ",dvstr(i)
c      enddo

      if(ichar(cvstr(1)(2:2)).ge.ichar('a').and.
     :     ichar(cvstr(1)(2:2)).le.ichar('z')) then
         atom = cvstr(1)(1:2)
         ion = cvstr(1)(3:6)
      
         n = real(10.0D0**dvstr(2))
         z = real(dvstr(3))
         b = real(dvstr(4))
      else
         atom = cvstr(1)(1:2)
         ion = cvstr(2)(1:4)
      
         n = real(10.0D0**dvstr(3))
         z = real(dvstr(4))
         b = real(dvstr(5))
      endif

      return
      end
