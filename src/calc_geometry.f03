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



SUBROUTINE calc_dist_eff(incoord,fin2coord,distances,number_of_in,&
                         &number_of_fin2,no_pairs,no_dist,smallest_R_open)
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
  REAL, INTENT(IN) :: smallest_R_open
  REAL :: eff !ICD efficiency

  INTEGER :: no_no_decay, partners ! counters

  no_no_decay = 0

  DO i=1,number_of_in
    partners = 0
    DO j=1,number_of_fin2
      xyz_in = incoord(i:i,1:3)
      xyz_fin = fin2coord(j:j,1:3)
      distances(row) = dist(xyz_in,xyz_fin)
!      WRITE (of,*) "Distance ", row, "is ", distances(row)
      IF (distances(row) > smallest_R_open) THEN
        partners = partners + 1
      END IF
      row = row+1
    END DO
    IF (partners == 0) THEN
      no_no_decay = no_no_decay + 1
    END IF
  END DO

  eff = 1 - (no_no_decay / number_of_in)

  WRITE(of,*)
  WRITE(of,*) '--------------------------------'
  WRITE(of,*) 'ICD efficiency = ', eff
  WRITE(of,*) '--------------------------------'

  no_dist = no_ind_entries(distances,no_pairs)
  WRITE(of,*) ''
  WRITE(of,*) 'Number of independent distances: ', no_dist
  WRITE(of,*) ''

END SUBROUTINE calc_dist_eff




SUBROUTINE create_dist_temp(distances,dist_temp,no_pairs,no_dist,no_dist_thresh)
! Purpose: To take the array distances and to write an new array
!          dist_stat containing how often a distance is invoked
!          and the distance itself

  use control
  use geometry_fun

  IMPLICIT NONE
  
! Data dictionary: exchanged with the main program
  INTEGER, INTENT(IN) :: no_pairs, no_dist
  REAL, DIMENSION(no_pairs) :: distances
  REAL, DIMENSION(no_dist,2), INTENT(OUT) :: dist_temp
  INTEGER             :: no_dist_thresh

! Data dictionary: temporary variables
  INTEGER :: i,j,row=0 !counters
  REAL :: temp, ne_dist
  REAL :: thresh = 0.001
  INTEGER :: ierror = 0

  
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
      
      dist_temp(row,1) = ne_dist
      dist_temp(row,2) = temp
    END IF
  END DO

  no_dist_thresh = 0

  DO i=1, no_dist
    IF (dist_temp(i,2) < Rmax) THEN
      no_dist_thresh = no_dist_thresh + 1
    END IF
  END DO

  

  OPEN(dist_outf, FILE='dist_stat', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
  WRITE(dist_outf,*) '#  R        number'
  WRITE(dist_outf,'(1X,F7.3,3X,F5.1)') (dist_temp(i,2),dist_temp(i,1), i=1,no_dist)

END SUBROUTINE create_dist_temp


SUBROUTINE create_dist_stat(dist_temp,dist_stat,no_dist,no_dist_thresh)
! Purpose: Reduce the number of calculated distances according to the threshold
!          given in the input file

  use control

!Data dictionary
  INTEGER :: no_dist, no_dist_thresh
  REAL,INTENT(IN),DIMENSION(no_dist,2)         :: dist_temp
  REAL,INTENT(OUT),DIMENSION(no_dist_thresh,2) :: dist_stat

! Data dictionary: Local Variables
  INTEGER :: i
  INTEGER :: row = 0

  DO i = 1, no_dist
    IF (dist_temp(i,2) < Rmax) THEN
      row = row + 1
      dist_stat(row,1) = dist_temp(i,1)
      dist_stat(row,2) = dist_temp(i,2)
    END IF
  END DO

END SUBROUTINE create_dist_stat



SUBROUTINE calc_triples(incoord,fin1coord,fin2coord,jacobi3,number_of_in&
                      &,number_of_fin1,number_of_fin2,no_triples)

  use control
  use geometry_fun

  IMPLICIT NONE

!Data dictionary: Input coordinates and output
  INTEGER, INTENT(IN) :: number_of_in,number_of_fin1,number_of_fin2
  INTEGER, INTENT(INOUT) :: no_triples
  REAL, INTENT(IN), DIMENSION(number_of_in,3) :: incoord
  REAL, INTENT(IN), DIMENSION(number_of_fin1,3) :: fin1coord
  REAL, INTENT(IN), DIMENSION(number_of_fin2,3) :: fin2coord

! Data dictionary: counters
  INTEGER :: i,j,k,l,row=0

! Temporary arrays and output variables
  REAL, DIMENSION(1,3) :: xyz_in,xyz_fin1,xyz_fin2
  REAL, DIMENSION(1,3) :: COM
  REAL, DIMENSION(1,3) :: a_vec,b_vec
  REAL :: distance
  REAL :: thresh = 1E-3
  REAL :: R,Q,theta,Coulomb_dist
  REAL :: argument
  REAL, ALLOCATABLE, DIMENSION(:,:) :: temparray

! Data dictionary: Output
  REAL, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: jacobi3


  ALLOCATE(temparray(no_triples,4))

!  WRITE(of,*) 'Parameters of triples'

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

          Q = dist(xyz_in,xyz_fin1)
          Qselect:IF (Q <= Qmax) THEN
  
  ! COM as reference
  !------------------------------------------------------
  !          COM = eval_com2(xyz_in,xyz_fin1) !checked
  
  !          WRITE(of,*) 'COM= ', COM
  !          WRITE(of,*) 'xyz_in= ', xyz_in
  !          WRITE(of,*) 'xyz_fin1= ', xyz_fin1
  !          WRITE(of,*) 'xyz_fin2= ', xyz_fin2
  
  !          a_vec = xyz_in - COM ! checked
  !          b_vec = xyz_fin2 - COM ! checked
  !------------------------------------------------------
  
  ! in as reference
  !------------------------------------------------------
            !a_vec = xyz_in ! checked
            a_vec = xyz_fin1 - xyz_in
            b_vec = xyz_fin2 - xyz_in ! checked
  !------------------------------------------------------
  
  
            argument = scalar_prod_row(a_vec,b_vec,3)/(norm_row(a_vec,3)*norm_row(b_vec,3))
            argument = ANINT(argument*10000)/10000
  
  !          WRITE(*,*) 'a_vec: ', a_vec
  !          WRITE(*,*) 'b_vec: ', b_vec
  
  !          WRITE(*,*) 'norm a_vec: ', norm_row(a_vec,3) !checked
  !          WRITE(*,*) 'norm b_vec: ', norm_row(b_vec,3) !checked
  !          WRITE(*,*) 'Scalar product a_vec,b_vec: ', scalar_prod_row(a_vec,b_vec,3)! checked
            WRITE(*,*) 'argument', argument
  
  
  !          R = dist(COM,xyz_fin2)
            R = dist(xyz_in,xyz_fin2)
            theta = ACOS(argument)
            Coulomb_dist = dist(xyz_fin1,xyz_fin2)
  
  

            row = row + 1

            temparray(row,1) = Q
            temparray(row,2) = R
            temparray(row,3) = theta
            temparray(row,4) = Coulomb_dist


!          WRITE(*,130) (jacobi3(row,l),l=1,4)
            130 FORMAT(' ',4F8.3)


          END IF Qselect

        END IF not_same
      END DO
    END DO
  END DO  

  ALLOCATE(jacobi3(row,4))

  jacobi3(1:row,1:4) = temparray(1:row,1:4)
  no_triples = row
  DEALLOCATE(temparray)

END SUBROUTINE calc_triples



END MODULE calc_geometry
