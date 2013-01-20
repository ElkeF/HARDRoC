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
  character(len=100) :: ctrl_file !Filename of the controlfile
  character(len=100) :: xyz_file  !Filename of the xyz-file

! The control file name is the first command-line argument

call get_command_argument(1, ctrl_file)
call get_command_argument(2, xyz_file)

! Call a subroutine to parse the control file
call read_control_file(ctrl_file)

! Call a subroutine to read in the xyz_file
call read_xyz_file(xyz_file)

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
        CASE ('vector')
           READ(buffer, *, IOSTAT=ierror) vector
           PRINT *, 'Read vector: ', vector
        CASE default
           PRINT *, 'Skipping invalid label at line', line
        END SELECT
     END IF
  END DO
END SUBROUTINE read_control_file



SUBROUTINE read_xyz_file(xyz_file)
  IMPLICIT NONE

! Data dictionary
  CHARACTER(len=100), INTENT(IN) :: xyz_file
  CHARACTER(len=100) :: buffer, label
  INTEGER :: pos
  INTEGER, PARAMETER :: fh = 16
  INTEGER :: ierror = 0
  INTEGER :: line = 0

  ! Control file variables
  INTEGER                        :: number_of_atoms
  REAL, dimension(5) :: vector

  OPEN(fh, FILE=xyz_file, STATUS='OLD', ACTION='READ',IOSTAT=ierror)

  ! ierror is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

  DO WHILE (ierror == 0)
     READ(fh,*,IOSTAT=ierror) number_of_atoms
     IF (ierror == 0) THEN
     WRITE(*,*) number_of_atoms
     END IF

     READ(fh, '(A)', IOSTAT=ierror) buffer
     IF (ierror == 0) THEN
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, '    ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        SELECT CASE (label)
!        CASE ('pairs')
!           READ(buffer, *, IOSTAT=ierror) do_pairs
!           PRINT *, 'Read do_pairs: ', do_pairs
        CASE ('vector')
           READ(buffer, *, IOSTAT=ierror) vector
           PRINT *, 'Read vector: ', vector
        CASE default
           PRINT *, 'Skipping invalid label at line', line
        END SELECT
     END IF
  END DO
END SUBROUTINE read_xyz_file
