      PROGRAM readfile 
      IMPLICIT NONE 
      REAL, DIMENSION(:), ALLOCATABLE :: mydata 
      INTEGER, PARAMETER :: maxrecs = 133334
      INTEGER :: J, NR, ios 
      CHARACTER(LEN=100) :: inputfile 
      CHARACTER(LEN=1) :: junk 
      write(*,*) 'Enter name of file to read in...'
      inputfile='spec.dat' 
!     Determine total number of lines in file 
      NR = 0
      OPEN(UNIT=1,FILE=inputfile) 
      DO J=1,maxrecs 
         READ(1,*,IOSTAT=ios) junk 
         IF (ios /= 0) EXIT
         IF (J == maxrecs) THEN
            write(*,*) 'Error: Maximum number of records exceeded...'
            write(*,*) 'Exiting program now...'
            STOP 
         ENDIF 
         NR = NR + 1
      ENDDO 
      REWIND(1) 
!Now we can allocate data variables 
      ALLOCATE(mydata(NR)) 
!Now read data into mydata 
      DO J=1,NR 
         READ(1,*) mydata(J) 
         print *,mydata(J)
      ENDDO 
      CLOSE(1) 
      END PROGRAM readfile
