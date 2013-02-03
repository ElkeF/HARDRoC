PROGRAM hardroc
! This is the program HARDRoC -- Hunting Asymptotic
! Relativistic Decay Rates of Clusters. The purpose of this
! program is, to read in xyz files of cluster structures,
! to calculate statistics of pair and triple parameters,
! read experimental parameters for the chosen systems and
! to calculate the partial and total decay rates of the cluster.

use control
use input_routines
use icd_calc
use lenfile
use wigner3j
use physical_constants
IMPLICIT NONE

!Data dictionary: Parameters for Inputfile reading
CHARACTER(len=100) :: ctrl_file !Filename of the controlfile
CHARACTER(len=100) :: xyz_file  !Filename of the xyz-file
CHARACTER(len=100) :: icd_channel_file  !Filename of the xyz-file

!Data dictionary: wirting output
CHARACTER(len=106)  :: nmout != 'output'


INTEGER :: number_of_in, number_of_fin1, number_of_fin2
INTEGER :: allocstatin=0, allocstatfin1=0, allocstatfin2=0
INTEGER :: ierror = 0 ! used for outputfile
INTEGER :: no_pairs, no_dist
INTEGER :: no_channels ! How many lines the file channels has


!Test
REAL :: factest

! Allocatable array which will be allocated during the
! calculation
REAL, ALLOCATABLE, DIMENSION(:,:) :: incoord !coordinates of initially ionized
REAL, ALLOCATABLE, DIMENSION(:,:) :: fin1coord ! coordinates of fin1
REAL, ALLOCATABLE, DIMENSION(:,:) :: fin2coord !coords of fin2
REAL, ALLOCATABLE, DIMENSION(:)   :: distances !array of distances for ICD
REAL, ALLOCATABLE, DIMENSION(:,:) :: dist_stat !number distance | distance
REAL, ALLOCATABLE, DIMENSION(:,:) :: channels  ! array specifying the parameters of the channels


! The control file name is the first command-line argument
CALL get_command_argument(1, ctrl_file)
CALL get_command_argument(2, xyz_file)
CALL get_command_argument(3, icd_channel_file)

!Set name of the outputfile
WRITE(nmout,*) TRIM(xyz_file)//'.out'

OPEN(of, FILE=TRIM(ADJUSTL(nmout)), STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)

! Call a subroutine to parse the control file
CALL read_control_file(ctrl_file)

! Call a subroutine to determine the length of coord arrays
CALL len_atarray(xyz_file,number_of_in,number_of_fin1,number_of_fin2)

! Allocate arrays for coordinates
ALLOCATE(incoord(number_of_in,3), STAT=allocstatin)
ALLOCATE(fin1coord(number_of_fin1,3), STAT=allocstatfin1)
ALLOCATE(fin2coord(number_of_fin2,3), STAT=allocstatfin2)

WRITE(of,*) 'Incoord allocated: ', ALLOCATED(incoord)
WRITE(of,*) 'Fin1coord allocated: ', ALLOCATED(fin1coord)
WRITE(of,*) 'Fin2coord allocated: ', ALLOCATED(fin2coord)

! Call a subroutine to read in the xyz_file
CALL read_xyz_file(xyz_file,incoord,fin1coord,fin2coord,&
     &             number_of_in,number_of_fin1,number_of_fin2)

! Write this to ouput file in the end, skipping now
WRITE(of,120) incoord
120 FORMAT(' ',3F8.3)

!Display the physical constants used
!CALL write_phys_const()

! At this point the coordinates and the input parameters have been
! read into the programme and we can now start to calculate the distances
! and internal coordinates of the pairs and triples
pairs:IF (do_pairs) THEN

  no_channels = len_file(icd_channel_file)
  WRITE(of,*) 'No channels= ', no_channels

  ! Allocate the memory for the channels array
  ALLOCATE(channels(no_channels,15))

  CALL read_icd_channels(icd_channel_file,channels,no_channels)

  no_pairs = number_of_in*number_of_fin2
  ALLOCATE(distances(no_pairs))

  CALL calc_distances(incoord,fin2coord,distances,number_of_in,&
                     &number_of_fin2,no_pairs,no_dist)


! Create array dist_stat with the format with the number of equal
! distances ne_dist and the corresponding distance.
! The rank 2 array is totally real, because arrays are not allowed
! to contain different types
  ALLOCATE(dist_stat(no_dist,2))
  CALL create_dist_stat(distances,dist_stat,no_pairs,no_dist)

  WRITE(of,121) dist_stat
  121 FORMAT (' ',F5.1,F8.3)

  DEALLOCATE(distances)

  CALL calc_icd_gamma(channels,dist_stat,no_channels,no_dist,number_of_in)
  DEALLOCATE(channels)

END IF pairs




triples:IF (do_triples) THEN
END IF triples

END PROGRAM hardroc











SUBROUTINE calc_distances(incoord,fin2coord,distances,number_of_in,&
                         &number_of_fin2,no_pairs,no_dist)
use control

IMPLICIT NONE

! Data dictionary
  INTEGER :: ierror = 0, row=1

  ! Arrays of coordinates and distances
  INTEGER :: i, j !counters
  INTEGER :: number_of_in, number_of_fin2, no_pairs
  INTEGER :: no_ind_entries
  INTEGER, INTENT(OUT) :: no_dist
  REAL, DIMENSION(number_of_in,3), INTENT(IN) :: incoord
  REAL, DIMENSION(number_of_fin2,3), INTENT(IN) :: fin2coord
  REAL, DIMENSION(no_pairs), INTENT(OUT) :: distances 
  REAL, DIMENSION(1,3) :: xyz_in, xyz_fin ! help arrays to evaluate function dist
  REAL :: dist

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

