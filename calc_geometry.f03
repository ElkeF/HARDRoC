MODULE calc_geometry


CONTAINS

SUBROUTINE calc_distances(incoord,fin2coord,distances,number_of_in,&
                         &number_of_fin2,no_pairs,no_dist)
use control
use geometry_fun

IMPLICIT NONE

! Data dictionary
  INTEGER :: ierror = 0, row=1

  ! Arrays of coordinates and distances
  INTEGER :: i, j !counters
  INTEGER :: number_of_in, number_of_fin2, no_pairs
  INTEGER, INTENT(OUT) :: no_dist
  REAL, DIMENSION(number_of_in,3), INTENT(IN) :: incoord
  REAL, DIMENSION(number_of_fin2,3), INTENT(IN) :: fin2coord
  REAL, DIMENSION(no_pairs), INTENT(OUT) :: distances 
  REAL, DIMENSION(1,3) :: xyz_in, xyz_fin ! help arrays to evaluate function dist

  DO i=1,number_of_in
    DO j=1,number_of_fin2
    xyz_in = incoord(i:i,1:3)
    xyz_fin = fin2coord(j:j,1:3)
    distances(row) = dist(xyz_in,xyz_fin)
!    WRITE (of,*) "Distance ", row, "is ", distances(row)
    row = row+1
    END DO
  END DO

  no_dist = no_ind_entries(distances,no_pairs)
  WRITE(of,*) 'Number of independent distances: ', no_dist

END SUBROUTINE calc_distances




SUBROUTINE create_dist_stat(distances,dist_stat,no_pairs,no_dist)
! Purpose: To take the array distances and to write an new array
!          dist_stat containing how often a distance is invoked
!          and the distance itself

  use control
  use geometry_fun

  IMPLICIT NONE
  
! Data dictionary: exchanged with the main program
  INTEGER, INTENT(IN) :: no_pairs, no_dist
  REAL, DIMENSION(no_pairs) :: distances
  REAL, DIMENSION(no_dist,2), INTENT(OUT) :: dist_stat

! Data dictionary: temporary variables
  INTEGER :: i,j,row=0 !counters
  REAL :: temp, ne_dist
  REAL :: thresh = 0.001

  
  DO i=1,no_pairs
    IF (ABS(distances(i)) > thresh) THEN
      row = row + 1
      temp = distances(i)
      distances(i) = 0.0
      ne_dist = 1
      
      DO j=i+1, no_pairs
        IF (ABS(distances(j)-temp) < thresh) THEN
          distances(j) = 0.0
          ne_dist = ne_dist + 1.0
        END IF
      END DO
      
      dist_stat(row,1) = ne_dist
      dist_stat(row,2) = temp
    END IF
  END DO

END SUBROUTINE create_dist_stat



END MODULE calc_geometry
