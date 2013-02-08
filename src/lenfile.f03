MODULE lenfile

CONTAINS


INTEGER FUNCTION len_file(filename)
! Purpose: determine the number of non commented lines in a given file
 IMPLICIT NONE
! Data dictionary
 CHARACTER(len=100), INTENT(IN) :: filename
 CHARACTER(len=100) :: buffer
 CHARACTER(len=1) :: dummy
 INTEGER :: pos
 INTEGER, PARAMETER :: fh = 16
 INTEGER :: ierror = 0
 INTEGER :: line = 0

 OPEN(fh, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
 WRITE(*,*) 'ierror: ', ierror

 DO WHILE (ierror == 0)
   READ(fh, '(A)', IOSTAT=ierror) buffer
   IF (ierror == 0) THEN
! Find first white space
! If the first character is #, then don't consider this line
     pos = scan(buffer, '    ')
     dummy = TRIM(buffer(1:pos))
     
     IF (LGE(dummy,'#').AND.LLE(dummy,'#')) THEN
     ELSE
       line = line + 1
     END IF
   END IF
 END DO

 CLOSE(fh)

 len_file = line

END FUNCTION len_file

END MODULE lenfile
