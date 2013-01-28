MODULE wigner3j
! Purpose: to calculate the Wigner 3j symbol for any given combination
! of j1, m1, j2, m2, j3, m3
! ( j1  j2  j3)
! (-m1  m2  m3)
! find Edmonds74 for formula
! be aware of the MINUS sign
use factorial

CONTAINS
! Test set of j and m values
! REAL :: j1 =  1.0
! REAL :: j2 =  1.0
! REAL :: j3 =  0.0
! REAL :: m1 =  1.0
! REAL :: m2 = -1.0
! REAL :: m3 =  0.0

REAL FUNCTION eval_wigner3j(j1,j2,j3,m1,m2,m3)
 IMPLICIT NONE
 
! Data dictionary
 REAL, INTENT(IN) :: j1,j2,j3,m1,m2,m3
 REAL :: fac

 fac = fact(j1+j2)
 

END FUNCTION eval_wigner3j

END MODULE wigner3j
