MODULE factorial
CONTAINS

RECURSIVE FUNCTION fact(n) RESULT(answer)
! Purpose: to calculate the factorial function
!          | n(n-1)!      n >= 1
!     n! = |
!          | 1            n = 0

 IMPLICIT NONE
! Data dictionary
 REAL, INTENT(IN) :: n !factorial of this number
 REAL :: answer        !result

 IF (n >= 1 ) THEN
   answer = n * fact(n-1)
 ELSE
   answer = 1
 END IF 

END FUNCTION fact

END MODULE factorial
