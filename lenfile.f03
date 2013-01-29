MODULE lenfile

CONTAINS

INTEGER FUNCTION len_file(filename)
 IMPLICIT NONE
! Data dictionary
 CHARACTER(len=100), INTENT(IN) :: filename
 CHARACTER(len=100) :: buffer
 INTEGER, PARAMETER :: fh = 16
 INTEGER :: ierror = 0
 INTEGER :: line = 0

 OPEN(fh, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

 DO WHILE (ierror == 0)
   READ(fh, '(A)', IOSTAT=ierror) buffer
   IF (ierror == 0) THEN
     line = line + 1
   END IF
 END DO

 len_file = line

END FUNCTION len_file

END MODULE lenfile
