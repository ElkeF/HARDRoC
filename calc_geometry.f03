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




SUBROUTINE calc_triples(incoord,fin1coord,fin2coord,jacobi3,number_of_in&
                      &,number_of_fin1,number_of_fin2,no_triples)

  use control
  use geometry_fun

  IMPLICIT NONE

!Data dictionary: Input coordinates and output
  INTEGER, INTENT(IN) :: number_of_in,number_of_fin1,number_of_fin2,no_triples
  REAL, INTENT(IN), DIMENSION(number_of_in,3) :: incoord
  REAL, INTENT(IN), DIMENSION(number_of_fin1,3) :: fin1coord
  REAL, INTENT(IN), DIMENSION(number_of_fin2,3) :: fin2coord
  REAL, DIMENSION(no_triples,4) :: jacobi3

! Data dictionary: counters
  INTEGER :: i,j,k,l,row=1

! Temporary arrays and output variables
  REAL, DIMENSION(1,3) :: xyz_in,xyz_fin1,xyz_fin2
  REAL, DIMENSION(1,3) :: COM
  REAL, DIMENSION(1,3) :: a_vec,b_vec
  REAL :: distance
  REAL :: thresh = 1E-3
  REAL :: R,Q,theta,Coulomb_dist

  WRITE(of,*) 'Parameters of triples'

! Assign the temp arrays and calculate the geometry parameters
! for all possible combinations
  DO i=1,number_of_in
    xyz_in = incoord(i:i,1:3)
    DO j=1,number_of_fin1
      xyz_fin1 = fin1coord(j:j,1:3)
      DO k=1,number_of_fin2
        xyz_fin2 = fin2coord(k:k,1:3)

! Only evaluate parameters of fin1coord and fin2coord are not the same
        distance = dist(xyz_fin1,xyz_fin2)
!        WRITE(of,*) 'Distance between final ', distance
        not_same:IF (distance > thresh) THEN

          COM = eval_com2(xyz_in,xyz_fin1)
!          WRITE(of,*) 'COM= ', COM
!          WRITE(of,*) 'xyz_in= ', xyz_in
!          WRITE(of,*) 'xyz_fin1= ', xyz_fin1
!          WRITE(of,*) 'xyz_fin2= ', xyz_fin2

          a_vec = xyz_in - COM
          b_vec = xyz_fin2 - COM

          Q = dist(xyz_in,xyz_fin1)
          R = dist(COM,xyz_fin2)
          theta = scalar_prod_row(a_vec,b_vec,3)/(norm_row(a_vec,3)*norm_row(b_vec,3))
          Coulomb_dist = dist(xyz_fin1,xyz_fin2)


          jacobi3(row,1) = Q
          jacobi3(row,2) = R
          jacobi3(row,3) = theta
          jacobi3(row,4) = Coulomb_dist


!          WRITE(of,130) (jacobi3(row,l),l=1,4)
          130 FORMAT(' ',4F8.3)

          row = row + 1

        END IF not_same
      END DO
    END DO
  END DO  



END SUBROUTINE calc_triples


END MODULE calc_geometry
