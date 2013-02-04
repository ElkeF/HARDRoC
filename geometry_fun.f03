MODULE geometry_fun

CONTAINS

FUNCTION eval_com2(xyz1,xyz2)
! Purpose: to calculate the center of mass of two atoms
! the coordinates are red in as row vectors and the atomic masses
! are taken from the module atomic_masses. The atom types are stored
! in the module control.

  use control
  use atomic_masses

  IMPLICIT NONE

! Data dictionary: input parameters
  REAL, DIMENSION(1,3), INTENT(IN) :: xyz1,xyz2

! Data dictionary: Output array
  REAL, DIMENSION(1,3) :: eval_com2

! Temporary variables
  CHARACTER(len=2), DIMENSION(2) :: atom_types
  REAL, DIMENSION(2) :: masses
  REAL :: sum_m=0.0
  REAL,DIMENSION(1,3) :: sum_m_r=0.0

! Counters
  INTEGER :: i

  atom_types(1) = in_atom_type
  atom_types(2) = fin_atom_type1


  DO i=1,2
    SELECT CASE(atom_types(i))
      CASE ('Ar')
        masses(i) = m_Ar
      CASE ('Xe')
        masses(i) = m_Xe
      CASE default
        WRITE(of,*) 'Invalid atom type.'
    END SELECT
    
    sum_m   = sum_m + masses(i)

  END DO

  sum_m_r   = masses(1)*xyz1 + xyz2*masses(2)

  eval_com2 = sum_m_r / sum_m

END FUNCTION eval_com2






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
