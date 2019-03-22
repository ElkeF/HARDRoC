MODULE efficiency

CONTAINS

  SUBROUTINE channel_opening_distances(no_channels,channels,&
                                       smallest_R_open)
! Purpose: Calculate the channel opening distances for all specified
!          ICD channels, print them and determine the smallest channel
!          opening distance.

    use physical_constants
    use control

    IMPLICIT NONE

! Data dictionary: Communicated data
    INTEGER, INTENT(IN)                         :: no_channels
    REAL, DIMENSION(no_channels,14), INTENT(IN) :: channels
    REAL, INTENT(OUT)                           :: smallest_R_open

! Data dictionary: Local data
    REAL, DIMENSION(no_channels)      :: SIP_in, SIP_fin1, SIP_fin2
    REAL, DIMENSION(no_channels)      :: shift_in, shift_fin1, shift_fin2
    INTEGER                           :: ichannel,i
    REAL, DIMENSION(no_channels)      :: R_angstrom
    REAL                              :: R_bohr
    REAL                              :: E_au, E_eV
    REAL                              :: E_in, E_fin1, E_fin2

! Get the ionization energies
    IF (do_fit) THEN
      DO ichannel = 1, no_channels
        SIP_in(ichannel)     = channels(ichannel,3)
        SIP_fin1(ichannel)   = channels(ichannel,6)
        SIP_fin2(ichannel)   = channels(ichannel,9)
        shift_in(ichannel)   = channels(ichannel,4)
        shift_fin1(ichannel) = channels(ichannel,7)
        shift_fin2(ichannel) = channels(ichannel,10)
      END DO
    ELSE
      DO ichannel = 1, no_channels
        SIP_in(ichannel)     = channels(ichannel,3)
        SIP_fin1(ichannel)   = channels(ichannel,7)
        SIP_fin2(ichannel)   = channels(ichannel,12)
        shift_in(ichannel)   = channels(ichannel,4)
        shift_fin1(ichannel) = channels(ichannel,8)
        shift_fin2(ichannel) = channels(ichannel,13)
      END DO
    END IF

    DO ichannel=1, no_channels
      E_in    = SIP_in(ichannel) + shift_in(ichannel)
      E_fin1  = SIP_fin1(ichannel) + shift_fin1(ichannel)
      E_fin2  = SIP_fin2(ichannel) + shift_fin2(ichannel)
      E_eV    = E_in - E_fin1 - E_fin2
      E_au    = E_eV * ev_to_hartree
      R_bohr  = 1/E_au
      R_angstrom(ichannel) = R_bohr * bohr_to_angstrom
    END DO

    IF (do_RICD) THEN
      smallest_R_open = 0.0
    ELSE
      smallest_R_open = MINVAL(R_angstrom)
    END IF

    WRITE(of,*) '' 
    WRITE(of,*) 'Channel   channel opening distance'
    WRITE(of,*) '----------------------------------'
    132 FORMAT(' ',I3,10X,F8.3)
    WRITE(of,132) (ichannel,R_angstrom(ichannel), ichannel=1,no_channels)
    WRITE(of,*) '' 
    WRITE(of,*) 'Smallest channel opening distance: ', smallest_R_open
    WRITE(of,*) '' 


  END SUBROUTINE channel_opening_distances

END MODULE efficiency
