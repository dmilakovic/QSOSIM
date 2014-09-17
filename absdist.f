!======================================================================
	subroutine absdist(wstart,zqso,bigX)
!======================================================================	
	IMPLICIT NONE
	real*8 :: wstart, bigX
	real*8 :: zstart, zend, zqso
	real*8 f, gauss16
	external f, gauss16
	
	write (6,*) '--------------------------------------------------'
        write (6,*) 'Calculating absorption distance!'
	zstart=(wstart/1215.67) - 1.   !zstart=1.96, wstart=3650
	zend=zqso
	bigX = gauss16(f,zstart,zend)
	write (6,*) 'Quasar redshift = ', zqso
	write (6,*) 'Total absorption distance, X(z) = ', bigX
	end subroutine absdist
!======================================================================
!  Function X(z)
!======================================================================
      function f(x)
      real*8 f, x
      f = ((1+x)**2)*(0.3*((1+x)**3) + 0.7)**(-0.5)
      return
      end function f
