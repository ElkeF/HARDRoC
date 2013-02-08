MODULE etmd_calc

CONTAINS

SUBROUTINE calc_etmd_gamma(channels,triple_parameters,no_channels,&
                           no_ind_triples,number_of_in)
! Purpose: to assign the parameters to the corresponding variables
! and calculate the ETMD decay rates

  use physical_constants
  use control
  use channel_char

  IMPLICIT NONE

! Data dictionary: Input parameters
  INTEGER, INTENT(IN) :: no_channels, no_ind_triples, number_of_in
  REAL, DIMENSION(no_channels,11), INTENT(IN) :: channels
  REAL, DIMENSION(no_ind_triples,5), INTENT(IN) :: triple_parameters

! Data dictionary: Parameters of the possible triples
  INTEGER :: neq_pairs !number of equivalent pairs, to be read in from triple_parameters(:,1)
  REAL :: Q_angstrom !read from triple_paramters(:,2)
  REAL :: Q_bohr
  REAL :: R_angstrom !to be read in from triple_parameters(:,3)
  REAL :: R_bohr !R_angstrom converted to bohr for calculation
  REAL :: theta !read from triple_parameters(:,4)
  REAL :: R_Coulomb_angstrom !read from triple_parameters(:,5)
  REAL :: R_Coulomb_bohr

! Data dictionary: experimental input parameters
  REAL :: trdm !transition dipole moment, to be calculated from fit in trdm_library
  REAL :: omega_vp !virtual photon energy
  REAL :: omega_vp_ev !in ev
  REAL :: sigma !ionization cross section, interpolated from lib
  REAL :: sigma_au ! in atomic units

! Data dictionary: Energies
  REAL :: E_in, E_fin1, E_fin2 !shifted energies
  REAL :: E_Coulomb !Coulomb energy
  REAL :: E_sec !Secondary electron produced in the ICD/ETMD process

! Data dictionary: Counters
  INTEGER :: ichannel
  INTEGER :: itriple

! Data dictionary: outputfile variables
  CHARACTER(len=15) :: specfile
  INTEGER :: ierror=0


  each_channel:DO ichannel=1,no_channels

  ! Assign channel parameters to values in module channel_char
    J_A         = channels(ichannel,1) / 2.0
    M_A         = channels(ichannel,2) / 2.0
    SIP_in      = channels(ichannel,3)
    shift_in    = channels(ichannel,4)

    J_D         = channels(ichannel,5) / 2.0
    M_D         = channels(ichannel,6) / 2.0
    SIP_fin1    = channels(ichannel,7)
    shift_fin1  = channels(ichannel,8)

    j_Bp        = channels(ichannel,9) / 2.0
    SIP_fin2    = channels(ichannel,10)
    shift_fin2  = channels(ichannel,11)


! Energies used for comparison and calculation
    E_in        = SIP_in + shift_in
    E_fin1      = SIP_fin1 + shift_fin1
    E_fin2      = SIP_fin2 + shift_fin2
    omega_vp_ev = E_in - E_fin1
    omega_vp    = (E_in - E_fin1) * ev_to_hartree


! Determine the ionization cross section



  END DO each_channel

END SUBROUTINE calc_etmd_gamma

END MODULE etmd_calc
