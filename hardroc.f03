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
REAL, ALLOCATABLE, DIMENSION(:,:) :: incoord
REAL, ALLOCATABLE, DIMENSION(:,:) :: fin1coord
REAL, ALLOCATABLE, DIMENSION(:,:) :: fin2coord

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

        WRITE(*,*) 'at: ', at

! Write the coordinates into the arrays incoord, fin1coord, fin2coord
        IF (LGE(at,in_atom_type).AND.LLE(at,in_atom_type)) THEN
          number_of_in = number_of_in +1
          READ(buffer, *, IOSTAT=ierror) (incoord(number_of_in,xyz), xyz=1,3)
          WRITE(*,*) 'Coordinates in incoord: ', (incoord(number_of_in,xyz), xyz=1,3)
        END IF
        IF (LGE(at,fin_atom_type1).AND.LLE(at,fin_atom_type1)) THEN
          number_of_fin1 = number_of_fin1 +1
          READ(buffer, *, IOSTAT=ierror) (fin1coord(number_of_fin1,xyz), xyz=1,3)
          WRITE(*,*) 'Coordinates in fin1coord: ', (fin1coord(number_of_fin1,xyz), xyz=1,3)
        END IF
        IF (LGE(at,fin_atom_type2).AND.LLE(at,fin_atom_type2)) THEN
          number_of_fin2 = number_of_fin2 +1
          READ(buffer, *, IOSTAT=ierror) (fin2coord(number_of_fin2,xyz), xyz=1,3)
          WRITE(*,*) 'Coordinates in fin2coord: ', (fin2coord(number_of_fin2,xyz), xyz=1,3)
        END IF

     END IF
  END DO
END SUBROUTINE read_xyz_file
