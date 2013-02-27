MODULE etmd_calc

CONTAINS

SUBROUTINE calc_etmd_gamma(channels,triple_parameters,no_channels,&
                           no_ind_triples,number_of_in)
! Purpose: to assign the parameters to the corresponding variables
! and calculate the ETMD decay rates

  use physical_constants
  use control
  use channel_char
  use select_sigma_para
  use select_trdm_para

  IMPLICIT NONE

! Data dictionary: Input parameters
  INTEGER, INTENT(IN) :: no_channels, no_ind_triples, number_of_in
  REAL, DIMENSION(no_channels,12), INTENT(IN) :: channels
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
  REAL :: sigma_abs
  REAL :: sigma_rel
  REAL :: sigma !ionization cross section, interpolated from lib
  REAL :: sigma_au ! in atomic units
  REAL :: trdm_x !transition dipole moment
  REAL :: trdm_z !transition dipole moment
  REAL :: trdm_x_au, trdm_z_au

! Data dictionary: Energies
  REAL :: E_in, E_fin1, E_fin2 !shifted energies
  REAL :: E_Coulomb !Coulomb energy
  REAL :: E_sec !Secondary electron produced in the ICD/ETMD process

! Data dictionary: parameters used for calculation of the trdm
! at a specific distance Q and ionization cross section for a given
! energy of the virtual photon
  REAL :: factor,alpha !trdm
  CHARACTER(len=1) :: dir
  REAL :: quad,lin,oneover !sigma
  REAL :: const !both

! Data dictionary: Counters
  INTEGER :: ichannel
  INTEGER :: itriple
  INTEGER :: iM_D !counter for loop over all M_D

! Data dictionary: Results
  REAL :: gamma_all_triples
  REAL :: gamma_b_all_M
  REAL :: gamma_b
  REAL :: gamma_b_triples
  REAL :: gamma_b_all_triples
  

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
    sigmarel    = channels(ichannel,12)


! Write the characteristics of each channel
    WRITE(of,*) ''
    WRITE(of,*) ''
    WRITE(of,210) 'Processing ETMD channel ', ichannel
    210 FORMAT (' ',A23,I3)
    WRITE(of,*) ''
    WRITE(of,*) 'Characterized by'
    WRITE(of,220) 'J_A = ',INT(2*J_A),'M_A = ',INT(2*M_A),"J_D = ",INT(2*J_D),&
                  "M_D = ",INT(2*M_D),"j_B' = ",INT(2*j_Bp)
    220 FORMAT (' ',5(A7,I2,4X))
    WRITE(of,*) ''
    


! Energies used for comparison and calculation
    E_in        = SIP_in + shift_in
    E_fin1      = SIP_fin1 + shift_fin1
    E_fin2      = SIP_fin2 + shift_fin2
    omega_vp_ev = E_in - E_fin1
    omega_vp    = (E_in - E_fin1) * ev_to_hartree


! Determine the ionization cross section
    CALL select_sigma_fit_para(quad,lin,const,oneover)
    sigma_abs = quad*omega_vp_ev**2 + lin*omega_vp_ev + const + oneover/omega_vp_ev
    sigma     = sigma_abs / (1 + 1/sigmarel)
    sigma_au  = sigma * megabarn_to_sqmeter * meter_to_bohr**2


! Special case if M_Ap = 88 calculate gamma for all possible values of M_A'
    calc_them_all:IF (INT(2*M_D) == 88) THEN

      gamma_all_triples = 0.0


      WRITE(specfile, '(A5,2(I1,A1),A4,I1)') 'ETMD_', INT(2*J_A), '_', INT(2*J_D)&
                                          & ,'_', 'all_', INT(2*j_Bp)

! Open the specfile for output
      OPEN(ETMD_outf,FILE=TRIM(ADJUSTL(specfile)),STATUS='UNKNOWN',ACTION='WRITE'&
          &,IOSTAT=ierror)


! Test whether this channel makes sense at all
      channel_sense_all:IF (E_in - E_fin1 - E_fin2 > 0) THEN

        WRITE(of,*) 'Processing all values of M_Donor'
        WRITE(of,*) 'All decay rates are given in eV'
        WRITE(of,*) ''
        WRITE(of,230) 'J_A','M_A',"J_A'","j_B'",'no triples','Q [AA]','R [AA]',&
                      'theta','E_ETMD','Gamma one',&
                      'Gamma all'
        230 FORMAT (' ',4(1X,A4),5(2X,A9),2(4X,A9))
        WRITE(of,*) '-----------------------------------------------------------------------------------------------------'

        gamma_b_triples = 0

! Now start for each and every triple
        DO itriple=1,no_ind_triples


! Find values for given triple
          neq_pairs  = INT(triple_parameters(itriple,1))
          Q_angstrom = triple_parameters(itriple,2)
          R_angstrom = triple_parameters(itriple,3)
          R_bohr     = R_angstrom * angstrom_to_bohr
          theta      = triple_parameters(itriple,4)
          R_Coulomb_angstrom = triple_parameters(itriple,5)
          R_Coulomb_bohr  = R_Coulomb_angstrom * angstrom_to_bohr

          E_Coulomb  = 1/R_Coulomb_bohr * hartree_to_ev
          E_sec      = E_in - E_fin1 - E_fin2 - E_Coulomb

!          WRITE(of,*) 'Q_angstrom = ', Q_angstrom
            !WRITE(of,*) 'R_angstrom = ', R_angstrom
!            WRITE(of,*) 'R_bohr = ', R_bohr

          gamma_b_all_M = 0
          gamma_b_all_triples = 0

          all_M_Ap:DO iM_D=(INT(2*(M_A-1))),(INT(2*(-M_A+1)))

            M_D = M_A + iM_D
! Set the M_D dependent variables
! Transition dipole moments :-)
            CALL select_trdm_fit_para(factor,alpha,const,dir)
          
            IF (dir == 'x') THEN
              trdm_x = factor * EXP(-alpha*Q_angstrom) + const
              trdm_z = 0.0
            ELSE IF (dir == 'z') THEN
              trdm_x = 0.0
              trdm_z = factor * EXP(-alpha*Q_angstrom) + const
            END IF 

!            WRITE(of,*) 'trdm_x = ', trdm_x
!            WRITE(of,*) 'trdm_z = ', trdm_z
!            WRITE(of,*) 'theta = ', theta
!            WRITE(of,*) 'sin   = ', SIN(theta)
!            WRITE(of,*) 'sin2  = ', SIN(theta)**2

!            WRITE(of,*) 'omega_vp = ', omega_vp
!            WRITE(of,*) 'sigma_au = ', sigma_au

! Verfahre nur weiter, wenn die Sekundaerenergie >=0 ist
            IF (E_sec >= 0.0) THEN
!              WRITE(of,*) 'E_sec= ', E_sec
              gamma_b = 2*pi/R_bohr**6  & !check
                        *(2*(trdm_x**2 * (1+COS(theta)**2) + trdm_z**2 * SIN(theta)**2) &
                        +4*(trdm_x**2 * SIN(theta)**2 + trdm_z**2 * COS(theta)**2)) &
                        * c_au *sigma_au/(4*pi**2 * omega_vp) * hartree_to_ev&
                        / number_of_in !Normalize to one ionization
!              WRITE(of,*) 'Gamma beta = ', gamma_b
!              gamma_b_triples = neq_pairs * gamma_b
              gamma_b_all_M = gamma_b_all_M + gamma_b

            END IF
          END DO all_M_Ap
          
          gamma_b_triples = neq_pairs * gamma_b_all_M

!Write result to output file
          WRITE(of,240) INT(2*J_A), INT(2*M_A), INT(2*J_D),&
                        INT(2*j_Bp), INT(neq_pairs), Q_angstrom, R_angstrom,&
                        theta, E_sec, gamma_b_all_M, gamma_b_triples
          240 FORMAT (' ',4(1X,I4),4X,I5,6X,4(F7.3,4X),2(ES9.3,4X))

!Write result to specfile
          WRITE(ETMD_outf,141) E_sec, gamma_b_triples
          141 FORMAT (' ',F12.4,ES15.5)


          gamma_all_triples = gamma_all_triples + gamma_b_triples

        END DO

        WRITE(of,*) '-----------------------------------------------------------------------------------------------------'
        WRITE(of,250) gamma_all_triples
        250 FORMAT (' ',92X,ES9.3)

      END IF channel_sense_all

      CLOSE(ETMD_outf)



! If we only want to calculate a certain value of M_A', then proceed
! with the following
    ELSE calc_them_all


! Open the specfile for output
      IF (INT(M_D) >= 0) THEN
        WRITE(specfile, '(A5,3(I1,A1),I1)') 'ETMD_', INT(2*J_A), '_', INT(2*J_D)&
                                          & ,'_', INT(2*M_D),'_', INT(2*j_Bp)
      ELSE
        WRITE(specfile, '(A5,2(I1,A1),I2,A1,I1)') 'ETMD_', INT(2*J_A), '_', INT(2*J_D)&
                                                  & ,'_', INT(2*M_Ap),'_', INT(2*j_D)
      END IF
      OPEN(ETMD_outf,FILE=TRIM(ADJUSTL(specfile)),STATUS='UNKNOWN',ACTION='WRITE'&
          &,IOSTAT=ierror)


! Test whether this channel makes sense at all
      channel_sense:IF (E_in - E_fin1 - E_fin2 > 0) THEN

      WRITE(of,430) 'J_A','M_A',"J_A'","M_A'","j_B'",'no triples','Q [AA]','R [AA]',&
                    'theta','E_ETMD','Gamma one',&
                    'Gamma all'
      430 FORMAT (' ',5(1X,A4),5(2X,A9),2(4X,A9))
      WRITE(of,*) '----------------------------------------------------------------------------------------------------------'


        DO itriple=1,no_ind_triples

! Find values for given triple
          neq_pairs  = INT(triple_parameters(itriple,1))
          Q_angstrom = triple_parameters(itriple,2)
          R_angstrom = triple_parameters(itriple,3)
          R_bohr     = R_angstrom * angstrom_to_bohr
          theta      = triple_parameters(itriple,4)
          R_Coulomb_angstrom = triple_parameters(itriple,2)
          R_Coulomb_bohr  = R_Coulomb_angstrom * angstrom_to_bohr

          E_Coulomb  = 1/R_Coulomb_bohr * hartree_to_ev
          E_sec      = E_in - E_fin1 - E_fin2 - E_Coulomb

!Set the M_Ap dependent variables
! Transition dipole moments :-)
          CALL select_trdm_fit_para(factor,alpha,const,dir)
          
          IF (dir == 'x') THEN
            trdm_x = factor * EXP(-alpha*Q_angstrom) + const
            trdm_z = 0.0
          ELSE IF (dir == 'z') THEN
            trdm_x = 0
            trdm_z = factor * EXP(-alpha*Q_angstrom) + const
          END IF 

!          WRITE(of,*) 'trdm_x', trdm_x
!          WRITE(of,*) 'trdm_z', trdm_z
!          WRITE(of,*) 'sigma_au', sigma_au
!          WRITE(of,*) 'omega_vp', omega_vp
!          WRITE(of,*) 'c_au', c_au
!          WRITE(of,*) 'cos theta', COS(theta)
!          WRITE(of,*) 'sin theta', SIN(theta)
  

! Verfahre nur weiter, wenn die Sekundaerenergie >=0 ist
          IF (E_sec >= 0.0) THEN
!            WRITE(of,*) 'E_sec= ', E_sec
              gamma_b = 2*pi/R_bohr**6  & !check
                        *(2*(trdm_x**2 * (1+COS(theta)**2) + trdm_z**2 * SIN(theta)**2) &
                        +4*(trdm_x**2 * SIN(theta)**2 + trdm_z**2 * COS(theta)**2)) &
                        * c_au *sigma_au/(4*pi**2 * omega_vp) * hartree_to_ev&
                        / number_of_in !Normalize to one ionization
!              WRITE(of,*) 'Gamma beta = ', gamma_b
              gamma_b_triples = neq_pairs * gamma_b
              gamma_b_all_triples = gamma_b_all_triples + gamma_b_triples


!Write result to output file
          WRITE(of,440) INT(2*J_A), INT(2*M_A), INT(2*J_D), INT(2*M_D),&
                        INT(2*j_Bp), INT(neq_pairs), Q_angstrom, R_angstrom,&
                        theta, E_sec, gamma_b, gamma_b_triples
          440 FORMAT (' ',5(1X,I4),4X,I5,6X,4(F7.3,4X),2(ES9.3,4X))

!Write result to specfile
          WRITE(ETMD_outf,141) E_sec, gamma_b_triples
!            141 FORMAT (' ',F12.4,ES15.5) is already defined in all
          END IF
        END DO

      WRITE(of,*) '----------------------------------------------------------------------------------------------------------'
      WRITE(of,450) gamma_b_all_triples
      450 FORMAT (' ',97X,ES9.3)     

      END IF channel_sense

      CLOSE(ETMD_outf)

   END IF calc_them_all

  END DO each_channel

END SUBROUTINE calc_etmd_gamma

END MODULE etmd_calc
