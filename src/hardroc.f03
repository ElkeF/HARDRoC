PROGRAM hardroc
! This is the program HARDRoC -- Hunting Asymptotic
! Relativistic Decay Rates of Clusters. The purpose of this
! program is, to read in xyz files of cluster structures,
! to calculate statistics of pair and triple parameters,
! read experimental parameters for the chosen systems and
! to calculate the partial and total decay rates of the cluster.

use control
use input_routines
use calc_geometry
use icd_calc
use etmd_calc
use lenfile
use array_operations
!use wigner3j
!use physical_constants
IMPLICIT NONE

!Data dictionary: Parameters for Inputfile reading
CHARACTER(len=100) :: ctrl_file !Filename of the controlfile
CHARACTER(len=100) :: xyz_file  !Filename of the xyz-file
CHARACTER(len=100) :: icd_channel_file  !Filename
CHARACTER(len=100) :: etmd_channel_file  !Filename


!Data dictionary: wirting output
CHARACTER(len=106)  :: nmout != 'output'


INTEGER :: number_of_in, number_of_fin1, number_of_fin2
INTEGER :: allocstatin=0, allocstatfin1=0, allocstatfin2=0
INTEGER :: ierror = 0 ! used for outputfile
INTEGER :: no_pairs, no_triples, no_dist, no_ind_triples
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
REAL, ALLOCATABLE, DIMENSION(:,:) :: channels  ! array specifying the parameters of the icd channels
REAL, ALLOCATABLE, DIMENSION(:,:) :: jacobi3 ! array to hold jacobi coords of all triples
REAL, ALLOCATABLE, DIMENSION(:,:) :: triple_parameters !no Q R theta Coulomb_dist


! The control file name is the first command-line argument
CALL get_command_argument(1, ctrl_file)
CALL get_command_argument(2, xyz_file)
CALL get_command_argument(3, icd_channel_file)
CALL get_command_argument(4, etmd_channel_file)


!Set name of the outputfile
WRITE(nmout,*) TRIM(xyz_file)//'.out'

OPEN(of, FILE=TRIM(ADJUSTL(nmout)), STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)

! Test the correct filename
WRITE(of,*) 'icd_file: ', icd_channel_file
WRITE(of,*) 'etmd_file: ', etmd_channel_file


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
!WRITE(of,120) incoord ! Take care wrong order in print out! todo
120 FORMAT(' ',3F8.3)

!Display the physical constants used
!CALL write_phys_const()

! At this point the coordinates and the input parameters have been
! read into the programme and we can now start to calculate the distances
! and internal coordinates of the pairs and triples
pairs:IF (do_pairs) THEN

  WRITE(of,*) ''
  WRITE(of,*) ''
  WRITE(of,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  WRITE(of,*) '%%%%%%%%%%%%%%%%% ICD %%%%%%%%%%%%%%%%%%%'
  WRITE(of,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  WRITE(of,*) ''

  no_channels = len_file(icd_channel_file)
  WRITE(of,*) 'No ICD channels= ', no_channels

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

!  WRITE(of,121) dist_stat
!  121 FORMAT (' ',F5.1,F8.3)

  DEALLOCATE(distances)

  CALL calc_icd_gamma(channels,dist_stat,no_channels,no_dist,number_of_in)
  DEALLOCATE(channels)
  DEALLOCATE(dist_stat)

END IF pairs


triples:IF (do_triples) THEN
  
  WRITE(of,*) ''
  WRITE(of,*) ''
  WRITE(of,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  WRITE(of,*) '%%%%%%%%%%%%%%%%% ETMD %%%%%%%%%%%%%%%%%%%'
  WRITE(of,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  WRITE(of,*) ''

  no_channels = len_file(etmd_channel_file)
  WRITE(of,*) 'Number of ETMD channels', no_channels

! Allocate the array for the channels array
  ALLOCATE(channels(no_channels,12))

  CALL read_etmd_channels(etmd_channel_file,channels,no_channels)

  no_triples = number_of_in*number_of_fin1*number_of_fin2

  ALLOCATE(jacobi3(no_triples,4))

  CALL calc_triples(incoord,fin1coord,fin2coord,jacobi3,number_of_in&
                  &,number_of_fin2,number_of_fin2,no_triples)


! Determine the number of independent entries in the array jacobi3
  no_ind_triples = diff_rows(jacobi3,no_triples,4)

  WRITE(of,*) 'Number of triples: ', no_triples
  WRITE(of,*) 'Number of different triples: ', no_ind_triples

! Allocate the array of the statistical triple parameters
  ALLOCATE(triple_parameters(no_ind_triples,5))

  CALL shrink_rows(jacobi3,triple_parameters,no_triples,no_ind_triples,4)

! Now we have the much smaller array containing the parameters of the triples
! and we deallocate the huge array
  DEALLOCATE(jacobi3)

! We have read in everything, now we need to get the corresponding values
! and calculate Gamma
  CALL calc_etmd_gamma(channels,triple_parameters,no_channels,&
                       no_ind_triples,number_of_in)

! Clean up rest of triple arrays
  DEALLOCATE(triple_parameters)
  DEALLOCATE(channels)

END IF triples

END PROGRAM hardroc











