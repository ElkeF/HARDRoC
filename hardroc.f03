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

! Exit the program if neither pairs not triples are to be calculated
!Funktioniert nicht wie gedacht
!IF (do_pairs.EQV..FALSE. .AND. do_triples.EQV..FALSE.) THEN
!  WRITE(*,*) 'Please choose either pairs or triples to be calculated.'
!  WRITE(*,*) 'Program terminates now'
!END IF

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



SUBROUTINE read_xyz_file(xyz_file)
use input_parameters
  IMPLICIT NONE

! Data dictionary
  CHARACTER(len=100), INTENT(IN) :: xyz_file
  CHARACTER(len=100) :: buffer, label
  CHARACTER(len=2)   :: at !atom type, requieres no space in front of atomtype in input file
  INTEGER :: pos
  INTEGER, PARAMETER :: fh = 16
  INTEGER :: ierror = 0
  INTEGER :: line = 0

  ! Control file variables
  INTEGER :: number_of_atoms=0, number_of_in=0, number_of_fin1=0, number_of_fin2=0
  REAL, ALLOCATABLE, DIMENSION(:,:) :: inarray
  REAL, ALLOCATABLE, DIMENSION(:,:) :: fin1array
  REAL, ALLOCATABLE, DIMENSION(:,:) :: fin2array
  INTEGER, DIMENSION(5) :: vector

! Find out how many initial and final state atoms are given in order
! to allocate arrays of appropriate length

  OPEN(fh, FILE=xyz_file, STATUS='OLD', ACTION='READ',IOSTAT=ierror)

  DO WHILE (ierror == 0)
    READ(fh, '(A)', IOSTAT=ierror) buffer
     IF (ierror == 0) THEN
        line = line + 1

        ! Find the first instance of whitespace.  Split atomtype in label from coords
        pos = scan(buffer, '    ')
        at = TRIM(buffer(1:pos))

! Check Character variables
!        WRITE(*,*) 'at: ', at
!        WRITE(*,*) 'Length of at: ', LEN(at)
!        WRITE(*,*) 'Atom type: ', in_atom_type
!        WRITE(*,*) 'Length atom type: ', LEN(in_atom_type)

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

  WRITE(*,*) '#initially  #final1     #final2'
  WRITE(*,10) number_of_in, number_of_fin1, number_of_fin2
  10 FORMAT (' ', 3I12)

  REWIND(UNIT=fh)

  ! ierror is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

  DO WHILE (ierror == 0)
!     READ(fh,*,IOSTAT=ierror) number_of_atoms
!     IF (ierror == 0) THEN
!     WRITE(*,*) number_of_atoms
!     END IF

     READ(fh, '(A)', IOSTAT=ierror) buffer
     IF (ierror == 0) THEN
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, '    ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        SELECT CASE (label)
        CASE ('vector')
           READ(buffer, *, IOSTAT=ierror) vector
           PRINT *, 'Read vector: ', vector
        CASE default
           PRINT *, 'Skipping invalid label at line', line
        END SELECT
     END IF
  END DO
END SUBROUTINE read_xyz_file
