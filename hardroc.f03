PROGRAM hardroc
! This is the program HARDRoC -- Hunting Asymptotic
! Relativistic Decay Rates of Clusters. The purpose of this
! program is, to read in xyz files of cluster structures,
! to calculate statistics of pair and triple parameters,
! read experimental parameters for the chosen systems and
! to calculate the partial and total decay rates of the cluster.

use input_parameters
IMPLICIT NONE

!Data dictionary: Parameters for Inputfile reading
CHARACTER(len=100) :: ctrl_file !Filename of the controlfile
CHARACTER(len=100) :: xyz_file  !Filename of the xyz-file

INTEGER :: number_of_in, number_of_fin1, number_of_fin2
INTEGER :: allocstatin=0, allocstatfin1=0, allocstatfin2=0
INTEGER :: no_pairs, no_dist

REAL, ALLOCATABLE, DIMENSION(:,:) :: incoord !coordinates of initially ionized
REAL, ALLOCATABLE, DIMENSION(:,:) :: fin1coord ! coordinates of fin1
REAL, ALLOCATABLE, DIMENSION(:,:) :: fin2coord !coords of fin2
REAL, ALLOCATABLE, DIMENSION(:)   :: distances !array of distances for ICD
REAL, ALLOCATABLE, DIMENSION(:,:) :: dist_stat !number distance | distance

! The control file name is the first command-line argument

CALL get_command_argument(1, ctrl_file)
CALL get_command_argument(2, xyz_file)

! Call a subroutine to parse the control file
CALL read_control_file(ctrl_file)

! Call a subroutine to determine the length of coord arrays
CALL len_atarray(xyz_file,number_of_in,number_of_fin1,number_of_fin2)

! Allocate arrays for coordinates
ALLOCATE(incoord(number_of_in,3), STAT=allocstatin)
ALLOCATE(fin1coord(number_of_fin1,3), STAT=allocstatfin1)
ALLOCATE(fin2coord(number_of_fin2,3), STAT=allocstatfin2)

WRITE(*,*) 'Incoord allocated: ', ALLOCATED(incoord)
WRITE(*,*) 'Fin1coord allocated: ', ALLOCATED(fin1coord)
WRITE(*,*) 'Fin2coord allocated: ', ALLOCATED(fin2coord)

! Call a subroutine to read in the xyz_file
CALL read_xyz_file(xyz_file,incoord,fin1coord,fin2coord,&
     &             number_of_in,number_of_fin1,number_of_fin2)

! Write this to ouput file in the end, skipping now
!WRITE(*,120) incoord
!120 FORMAT(' ',3F8.3)

! At this point the coordinates and the input parameters have been
! read into the programme and we can now start to calculate the distances
! and internal coordinates of the pairs and triples
pairs:IF (do_pairs) THEN
  no_pairs = number_of_in*number_of_fin1
  ALLOCATE(distances(no_pairs))

  CALL calc_distances(incoord,fin1coord,distances,number_of_in,&
                     &number_of_fin1,no_pairs,no_dist)
  

! Create array dist_stat with the format with the number of equal
! distances ne_dist and the corresponding distance.
! The rank 2 array is totally real, because arrays are not allowed
! to contain different types
  ALLOCATE(dist_stat(no_dist,2))
  CALL create_dist_stat(distances,dist_stat,no_pairs,no_dist)

  WRITE(*,121) dist_stat
  121 FORMAT (' ',F5.1,F8.3)

END IF pairs

triples:IF (do_triples) THEN
END IF triples

END PROGRAM hardroc




SUBROUTINE read_control_file(ctrl_file)
use input_parameters
IMPLICIT NONE
! Input related variables
  CHARACTER(len=100) :: ctrl_file
  CHARACTER(len=100) :: buffer, label
  INTEGER :: pos
  INTEGER, PARAMETER :: fh = 15
  INTEGER :: ierror = 0
  INTEGER :: line = 0

  ! Control file variables
  integer, dimension(5) :: vector

  OPEN(fh, FILE=ctrl_file, STATUS='OLD', ACTION='READ',IOSTAT=ierror)

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

  DO WHILE (ierror == 0)
     READ(fh, '(A)', IOSTAT=ierror) buffer
     IF (ierror == 0) THEN
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, '    ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        SELECT CASE (label)
        CASE ('pairs')
           READ(buffer, *, IOSTAT=ierror) do_pairs
           PRINT *, 'Read do_pairs: ', do_pairs
        CASE ('triples')
           READ(buffer, *, IOSTAT=ierror) do_triples
           PRINT *, 'Read do_triples: ', do_triples
        CASE ('atin')
           READ(buffer, *, IOSTAT=ierror) in_atom_type
           PRINT *, 'Read in_atom_type: ', in_atom_type
        CASE ('atfin1')
           READ(buffer, *, IOSTAT=ierror) fin_atom_type1
           PRINT *, 'Read fin_atom_type1: ', fin_atom_type1
        CASE ('atfin2')
           READ(buffer, *, IOSTAT=ierror) fin_atom_type2
           PRINT *, 'Read fin_atom_type2: ', fin_atom_type2
        CASE ('vector')
           READ(buffer, *, IOSTAT=ierror) vector
           PRINT *, 'Read vector: ', vector
        CASE default
           PRINT *, 'Skipping invalid label at line', line
        END SELECT
     END IF
  END DO
END SUBROUTINE read_control_file


SUBROUTINE len_atarray(xyz_file,number_of_in,number_of_fin1,number_of_fin2)
use input_parameters
  IMPLICIT NONE

! Data dictionary
  CHARACTER(len=100), INTENT(IN) :: xyz_file
  CHARACTER(len=100) :: buffer
  CHARACTER(len=2)   :: at !atom type, requieres no space in front of atomtype in input file
  INTEGER :: pos
  INTEGER, PARAMETER :: fh = 16
  INTEGER :: ierror = 0
  INTEGER :: line = 0

  ! Control file variables
  INTEGER,INTENT(OUT) :: number_of_in, number_of_fin1, number_of_fin2

! Find out how many initial and final state atoms are given in order
! to allocate arrays of appropriate length
  number_of_in = 0
  number_of_fin1 = 0
  number_of_fin2 = 0

  OPEN(fh, FILE=xyz_file, STATUS='OLD', ACTION='READ',IOSTAT=ierror)

  DO WHILE (ierror == 0)
    READ(fh, '(A)', IOSTAT=ierror) buffer
     IF (ierror == 0) THEN
        line = line + 1

        ! Find the first instance of whitespace.  Split atomtype in at from coords
        pos = scan(buffer, '    ')
        at = TRIM(buffer(1:pos))

! Find out, how many atoms we have of each type
! These will give us the length of the coordinate arrays
! in the following.

        IF (LGE(at,in_atom_type).AND.LLE(at,in_atom_type)) THEN
          number_of_in = number_of_in +1
!          WRITE(*,*) 'Number of initially ionized atoms: ', number_of_in
        END IF
        IF (LGE(at,fin_atom_type1).AND.LLE(at,fin_atom_type1)) THEN
          number_of_fin1 = number_of_fin1 +1
!          WRITE(*,*) 'Number of atoms ionized in the final state1', number_of_fin1
        END IF
        IF (LGE(at,fin_atom_type2).AND.LLE(at,fin_atom_type2)) THEN
          number_of_fin2 = number_of_fin2 +1
!          WRITE(*,*) 'Number of atoms ionized in the final state2', number_of_fin2
        END IF
     END IF
  END DO

  CLOSE(UNIT=fh)

  WRITE(*,*) '#initially  #final1     #final2'
  WRITE(*,11) number_of_in, number_of_fin1, number_of_fin2
  11 FORMAT (' ', 3I12)

END SUBROUTINE len_atarray



SUBROUTINE read_xyz_file(xyz_file,incoord,fin1coord,fin2coord,&
           &              number_of_in,number_of_fin1,number_of_fin2)
use input_parameters
  IMPLICIT NONE

! Data dictionary
  CHARACTER(len=100), INTENT(IN) :: xyz_file
  CHARACTER(len=100) :: buffer
  CHARACTER(len=2)   :: at !atom type, requieres no space in front of atomtype in input file
  INTEGER :: pos
  INTEGER, PARAMETER :: fh = 16
  INTEGER :: ierror = 0
  INTEGER :: line = 0
  INTEGER :: allocstatin=0, allocstatfin1=0, allocstatfin2=0
  INTEGER :: xyz !Count coordinates x=1, y=2, z=3

  ! Control file variables
  INTEGER :: number_of_in, number_of_fin1, number_of_fin2
  REAL, DIMENSION(number_of_in,3) :: incoord
  REAL, DIMENSION(number_of_fin1,3) :: fin1coord
  REAL, DIMENSION(number_of_fin2,3) :: fin2coord

  OPEN(fh, FILE=xyz_file, STATUS='OLD', ACTION='READ',IOSTAT=ierror)

! Reactivate number_of_in as counter
  number_of_in = 0
  number_of_fin1 = 0
  number_of_fin2 = 0

  DO WHILE (ierror == 0)

     READ(fh, '(A)', IOSTAT=ierror) buffer
     IF (ierror == 0) THEN
        line = line + 1

        ! Find the first instance of whitespace.  Split atom type and coords.
        pos = scan(buffer, '    ')
        at = TRIM(buffer(1:pos))
        buffer = buffer(pos+1:)

!Check the correct atomtype, skip for now
!        WRITE(*,*) 'at: ', at

! Write the coordinates into the arrays incoord, fin1coord, fin2coord
        IF (LGE(at,in_atom_type).AND.LLE(at,in_atom_type)) THEN
          number_of_in = number_of_in +1
          READ(buffer, *, IOSTAT=ierror) (incoord(number_of_in,xyz), xyz=1,3)
!          WRITE(*,*) 'Coordinates in incoord: ', (incoord(number_of_in,xyz), xyz=1,3)
        END IF
        IF (LGE(at,fin_atom_type1).AND.LLE(at,fin_atom_type1)) THEN
          number_of_fin1 = number_of_fin1 +1
          READ(buffer, *, IOSTAT=ierror) (fin1coord(number_of_fin1,xyz), xyz=1,3)
!          WRITE(*,*) 'Coordinates in fin1coord: ', (fin1coord(number_of_fin1,xyz), xyz=1,3)
        END IF
        IF (LGE(at,fin_atom_type2).AND.LLE(at,fin_atom_type2)) THEN
          number_of_fin2 = number_of_fin2 +1
          READ(buffer, *, IOSTAT=ierror) (fin2coord(number_of_fin2,xyz), xyz=1,3)
!          WRITE(*,*) 'Coordinates in fin2coord: ', (fin2coord(number_of_fin2,xyz), xyz=1,3)
        END IF

     END IF
  END DO
END SUBROUTINE read_xyz_file




SUBROUTINE calc_distances(incoord,fin1coord,distances,number_of_in,&
                         &number_of_fin1,no_pairs,no_dist)
IMPLICIT NONE

! Data dictionary
  INTEGER :: ierror = 0, row=1

  ! Arrays of coordinates and distances
  INTEGER :: i, j !counters
  INTEGER :: number_of_in, number_of_fin1, no_pairs
  INTEGER :: no_ind_entries
  INTEGER, INTENT(OUT) :: no_dist
  REAL, DIMENSION(number_of_in,3), INTENT(IN) :: incoord
  REAL, DIMENSION(number_of_fin1,3), INTENT(IN) :: fin1coord
  REAL, DIMENSION(no_pairs), INTENT(OUT) :: distances 
  REAL, DIMENSION(1,3) :: xyz_in, xyz_fin ! help arrays to evaluate function dist
  REAL :: dist

  DO i=1,number_of_in
    DO j=1,number_of_fin1
    xyz_in = incoord(i:i,1:3)
    xyz_fin = fin1coord(j:j,1:3)
    distances(row) = dist(xyz_in,xyz_fin)
!    WRITE (*,*) "Distance ", row, "is ", distances(row)
    row = row+1
    END DO
  END DO

  no_dist = no_ind_entries(distances,no_pairs)
  WRITE(*,*) 'Number of independent distances: ', no_dist

END SUBROUTINE calc_distances




SUBROUTINE create_dist_stat(distances,dist_stat,no_pairs,no_dist)
! Purpose: To take the array distances and to write an new array
!          dist_stat containing how often a distance is invoked
!          and the distance itself

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
! by the 3D coordinates in incoord(i,1:3) and fin1coord(j,1:3)  

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
!    WRITE(*,*) 'Diff: ', diff
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
