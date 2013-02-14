MODULE wigner3j
! Purpose: to calculate the Wigner 3j symbol for any given combination
! of j1, m1, j2, m2, j3, m3
! (j1  j2  j3)
! (m1  m2  m3)
! using one of the formulas from Edmonds74
! be aware of the MINUS sign, which is introduced for m1 to m1m !!!
! 
use factorial

CONTAINS

REAL FUNCTION eval_wigner3j(j1orig,j2orig,j3orig,m1orig,m2orig,m3orig)
 IMPLICIT NONE
 
! Data dictionary
 REAL,INTENT(IN) :: j1orig,j2orig,j3orig,m1orig,m2orig,m3orig
 REAL :: j1,j2,j3,m1,m2,m3
 REAL :: m1m
 REAL :: fac
 REAL :: zaehler1,zaehler2,zaehler3
 REAL :: j1add,j1diff,j2add,j2diff,j3add,j3diff
 REAL :: nenner1,nenner2,nenner3,nenner4,nenner5,nenner6
 REAL :: prefactor
 REAL :: adddiff, root
 REAL :: nenner,sum_nenner=0
 REAL :: vc_coef,expon
 INTEGER :: k,kmax !counter
 REAL :: l !real of k
 REAL :: temp

 sum_nenner= 0 

! assign the original values to dummy values, which can be changed
 j1 = j1orig
 j2 = j2orig
 j3 = j3orig
 m1 = m1orig
 m2 = m2orig
 m3 = m3orig

! Change input of 3j symbol to m1m, which will be used
! for the evaluation of the 3j symbol
 m1m = -m1
! m1m = m1


 shuffle:IF (j1-j2-m3 < 0) THEN
   temp = j1
   j1   = j2
   j2   = j3
   j3   = temp
   temp = -m1m
   m1m  = -m2
   m2   = m3
   m3   = temp

 END IF shuffle

! WRITE(*,*) 'm1m= ', m1m
! WRITE(*,*) 'm2+m3= ', m2+m3

! define variables
 zaehler1 = j2+j3-j1
 zaehler2 = j1+j2-j3
 zaehler3 = j1+j3-j2
 j1add    = j1+m1m
 j1diff   = j1-m1m
 j2add    = j2+m2
 j2diff   = j2-m2
 j3add    = j3+m3
 j3diff   = j3-m3
 nenner1  = j1+j2+j3+1.0
 nenner2  = j2+j3-j1
 nenner3  = j2-m2
 nenner4  = j3+m3
 nenner5  = j1-j3+m2
 nenner6  = j1-j2-m3
 expon    = INT(j2-j3+m1m)

!   WRITE(*,*) nenner2
!   WRITE(*,*) nenner3
!   WRITE(*,*) nenner4
!   WRITE(*,*) nenner5
!   WRITE(*,*) nenner6

! Evaluate delta function
! WRITE(*,*) 'm1m = ', m1m
! WRITE(*,*) 'm2+m3 = ', m2+m3
 delta:IF ( m1m == m2 + m3) THEN

   IF (nenner2 > nenner3) THEN
     IF (nenner3 > nenner4) THEN
       kmax = INT(nenner4)
     ELSE
       kmax = INT(nenner3)
     END IF
   ELSE
     IF (nenner2 > nenner4) THEN
       kmax = INT(nenner4)
     ELSE
       kmax = INT(nenner2)
     END IF
   END IF

!   WRITE(*,*) 'kmax= ', kmax

   prefactor = SQRT(fact(zaehler1)*fact(zaehler2)*fact(zaehler3)*(2*j1+1)/fact(nenner1))
   adddiff   = fact(j1add)*fact(j1diff)*fact(j2add)*fact(j2diff)*fact(j3add)*fact(j3diff)
   root      = SQRT(adddiff)
   
   DO k=0,kmax
     l = REAL(k)
     nenner = fact(l)*fact(nenner2-l)*fact(nenner3-l)*fact(nenner4-l)*fact(nenner5+l)&
             &*fact(nenner6+l)
     sum_nenner = sum_nenner + (-1.0)**k / nenner
!     WRITE(*,*) 'Nenner in loop', nenner
!     WRITE(*,*) 'sum in loop = ', sum_nenner
   END DO


!   WRITE(*,*) 'prefactor= ', prefactor
!   WRITE(*,*) 'adddiff= ', adddiff
!   WRITE(*,*) 'root= ', root
!   WRITE(*,*) 'nenner= ', nenner
!   WRITE(*,*) 'sum_nenner= ', sum_nenner
  
   vc_coef = prefactor*root*sum_nenner

 ELSE delta
   vc_coef = 0.0 
 END IF delta

! WRITE(*,*) 'vc_coef= ', vc_coef
 eval_wigner3j = (-1)**expon/SQRT(2*j1+1) * vc_coef


END FUNCTION eval_wigner3j

END MODULE wigner3j
