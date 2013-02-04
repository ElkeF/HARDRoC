MODULE geometry_fun

CONTAINS

REAL FUNCTION dist (xyz_in, xyz_fin)
! Purpose: Evaluate the distance between two atoms defined
! by the 3D coordinates in incoord(i,1:3) and fin2coord(j,1:3)  

  IMPLICIT NONE

! Data dictionary
  REAL, DIMENSION(1,3),INTENT(IN) :: xyz_in, xyz_fin
  REAL :: diff, square, sum_squares
  INTEGER :: i,j

  sum_squares = 0.0

  DO i=1,3
    diff = xyz_in(1,i) - xyz_fin(1,i)
    square = diff**2
    sum_squares = sum_squares + square
!    WRITE(of,*) 'Diff: ', diff
  ENDDO

  dist = SQRT(sum_squares)

END FUNCTION dist



INTEGER FUNCTION no_ind_entries(array,length)
! Purpose: To determine the number of different entries in an array.

  IMPLICIT NONE

! Data dictionary
  INTEGER, INTENT(IN) :: length
  INTEGER :: i,j !counter
  REAL :: thresh = 0.001
  REAL :: temp
  REAL, DIMENSION(length), INTENT(IN) :: array
  REAL, DIMENSION(length) :: temparray

  no_ind_entries = 0

  temparray = array

  DO i=1,length
    IF (ABS(temparray(i))>thresh) THEN
      temp = temparray(i)
      temparray(i) = 0.0
      no_ind_entries = no_ind_entries + 1

      DO j=i+1, length
        IF (ABS(temparray(j)-temp) < thresh) THEN
          temparray(j) = 0.0
        END IF
      END DO
    END IF
  END DO

END FUNCTION no_ind_entries

END MODULE geometry_fun
