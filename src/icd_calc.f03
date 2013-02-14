MODULE icd_calc

CONTAINS

SUBROUTINE calc_icd_gamma(channels,dist_stat,no_channels,no_dist,number_of_in)
!Purpose: to calculate the partial and total ICD Gammas

!Modules
  use control
  use channel_char
  use select_sigma_para
  use physical_constants
  use wigner3j

  IMPLICIT NONE

!Data dictionary: arrays
  INTEGER, INTENT(IN) :: no_channels,no_dist,number_of_in
  REAL, DIMENSION(no_channels,14), INTENT(IN) :: channels
  REAL, DIMENSION(no_dist,2) :: dist_stat

! Data dictionary: distances and number of equivalent pairs
  INTEGER :: neq_pairs !number of equivalent pairs, to be read in from dist_stat(:,1)
  REAL :: R_angstrom !to be read in from dist_stat(:,2)
  REAL :: R_bohr !R_angstrom converted to bohr for calculation
  
! Data dictionary: actual variables in calculation in correct units
  REAL :: wigner !will hold the alue of the wigner3j symbol
  REAL :: one = 1.0! for wigner
  INTEGER :: B_MAMAp !value of B depending on M_Ap-M_A
  REAL :: tau !Radiative lifetime of the channel
  REAL :: tau_au !in atomic units
  REAL :: omega_vp !energy of the virtual photon transferred
  REAL :: omega_vp_ev
  REAL :: sigma ! ionization cross section
  REAL :: sigma_au !in atomic units

!Data dictionary
  INTEGER :: ichannel !counter of the channels
  INTEGER :: idist ! counter of the collected pairs
  REAL :: E_in, E_fin1, E_fin2 !shifted energies
  REAL :: E_Coulomb !Coulomb energy
  REAL :: E_sec !Secondary electron produced in the ICD/ETMD process
  REAL :: gamma_b !Decay rate of channel beta
  REAL :: gamma_b_pairs !Decay rate considering all equivalent pairs
  REAL :: gamma_b_all_pairs !Sum of all gamma_b_pairs
  REAL :: gamma_b_all_M

! Data dictionary: outputfile variables
  CHARACTER(len=15) :: specfile
  INTEGER :: ierror=0

! Data dictionary: fit parameters for ionization cross section
  REAL :: quad,lin,const,oneover

! Special for all variables
  INTEGER :: iM_Ap !counter for loop over all M_A'
  REAL :: gamma_all !will hold the result for all M_A' values of one distance


  each_channel:DO ichannel=1,no_channels

    gamma_all = 0.0


! Assign channel parameters to values in module channel_char
    J_A         = channels(ichannel,1) / 2.0
    M_A         = channels(ichannel,2) / 2.0
    SIP_in      = channels(ichannel,3)
    shift_in    = channels(ichannel,4)

    J_Ap        = channels(ichannel,5) / 2.0
    M_Ap        = channels(ichannel,6) / 2.0
    SIP_fin1    = channels(ichannel,7)
    shift_fin1  = channels(ichannel,8)
    tottau      = channels(ichannel,9)
    taurel      = channels(ichannel,10)

    j_Bp        = channels(ichannel,11) / 2.0
    SIP_fin2    = channels(ichannel,12)
    shift_fin2  = channels(ichannel,13)
!    sigmaabs    = channels(ichannel,14)
    sigmarel    = channels(ichannel,14)



! Write the characteristics of each channel
    WRITE(of,*) ''
    WRITE(of,*) ''
    WRITE(of,210) 'Processing ICD channel ', ichannel
    210 FORMAT (' ',A23,I3)
    WRITE(of,*) ''
    WRITE(of,*) 'Characterized by'
    WRITE(of,220) 'J_A = ',INT(2*J_A),'M_A = ',INT(2*M_A),"J_A' = ",INT(2*J_Ap),&
                  "M_A' = ",INT(2*M_Ap),"j_B' = ",INT(2*j_Bp)
    220 FORMAT (' ',5(A7,I2,4X))
    WRITE(of,*) ''
    
  

! Energies used for comparison and calculation
    E_in        = SIP_in + shift_in
    E_fin1      = SIP_fin1 + shift_fin1
    E_fin2      = SIP_fin2 + shift_fin2
    omega_vp_ev = E_in - E_fin1
    omega_vp    = (E_in - E_fin1) * ev_to_hartree

!Determine the ionization cross section
    CALL select_sigma_fit_para(quad,lin,const,oneover)
    sigmaabs = quad*omega_vp_ev**2 + lin*omega_vp_ev + const + oneover/omega_vp_ev
    sigma     = sigmaabs / (1 + 1/sigmarel)
    sigma_au  = sigma * megabarn_to_sqmeter * meter_to_bohr**2

    tau      = tottau/(1 + 1/taurel)
    tau_au   = tau * second_to_atu
!    WRITE(of,*) 'tau_au= ', tau_au



! Special case if M_Ap = 88 calculate gamma for all possible values of M_A'
    calc_them_all:IF (INT(2*M_Ap) == 88) THEN


      WRITE(of,*) 'Processing all values of M_A prime'
!      WRITE(of,*) ''
      WRITE(of,*) 'All decay rates are given in eV'
      WRITE(of,*) ''
      WRITE(of,230) 'J_A','M_A',"J_A'","j_B'",'no pairs','R [AA]','Gamma one',&
                    'Gamma all'
      230 FORMAT (' ',4(1X,A4),2(2X,A9),2(4X,A9))
      WRITE(of,*) '---------------------------------------------------------------------'


      WRITE(specfile, '(A4,2(I1,A1),A4,I1)') 'ICD_', INT(2*J_A), '_', INT(2*J_Ap)&
                                          & ,'_', 'all_', INT(2*j_Bp)

! Open the specfile for output
      OPEN(ICD_outf,FILE=TRIM(ADJUSTL(specfile)),STATUS='UNKNOWN',ACTION='WRITE'&
          &,IOSTAT=ierror)


! Test whether this channel makes sense at all
      channel_sense_all:IF (E_in - E_fin1 - E_fin2 > 0) THEN
  

! Each and every pair.
         DO idist=1,no_dist

           gamma_b_all_M     = 0.0
           gamma_b_all_pairs = 0.0

! Find values for given pair
           neq_pairs  = INT(dist_stat(idist,1))
           R_angstrom = dist_stat(idist,2)
           R_bohr     = R_angstrom * angstrom_to_bohr
     
           E_Coulomb  = 1/R_bohr * hartree_to_ev
           E_sec      = E_in - E_fin1 - E_fin2 - E_Coulomb


           all_M_Ap:DO iM_Ap=(INT(2*(M_A-1))),(INT(2*(-M_A+1)))

             M_Ap = M_A + iM_Ap

! Set the M_A# dependent variables
             SELECT CASE (INT(M_Ap - M_A))
               CASE (-1,1)
                 B_MAMAp = 1
               CASE (0)
                 B_MAMAp = -2
               CASE DEFAULT
                 B_MAMAp = 0
             END SELECT
!             WRITE(of,*) 'B_MAMAp= ', B_MAMAp

             wigner   = eval_wigner3j(J_Ap,one,J_A,-M_Ap,M_Ap-M_A,M_A)
!             WRITE(of,*) 'Wigner3j symbol: ', wigner


! Verfahre nur weiter, wenn die Sekundaerenergie >=0 ist
            IF (E_sec >= 0.0) THEN
!              WRITE(of,*) 'E_sec= ', E_sec
              gamma_b = 2*pi/R_bohr**6 * B_MAMAp**2 * wigner**2 *(2*J_A+1)&
                      & *3*c_au**4 *sigma_au/(16*pi**2 * omega_vp**4 * tau_au) * hartree_to_ev&
                      & / number_of_in !Normalize to one ionization
!              WRITE(of,*) 'Gamma beta = ', gamma_b
              gamma_b_all_M = gamma_b_all_M + gamma_b
    
            END IF
          END DO all_M_Ap

          gamma_b_pairs = neq_pairs * gamma_b_all_M
          
!Write result to output file
          WRITE(of,240) INT(2*J_A), INT(2*M_A), INT(2*J_Ap),&
                        INT(2*j_Bp), INT(neq_pairs), R_angstrom,&
                        gamma_b_all_M, gamma_b_pairs
          240 FORMAT (' ',4(1X,I4),4X,I5,6X,F7.3,4X,2(ES9.3,4X))
! Write results to specfile
          WRITE(ICD_outf,141) E_sec, gamma_b_all_M
          141 FORMAT (' ',F12.4,ES15.5)

          gamma_all = gamma_all + gamma_b_pairs

        END DO 

        WRITE(of,*) '---------------------------------------------------------------------'
        WRITE(of,250) gamma_all
        250 FORMAT (' ',53X,ES15.3)

      END IF channel_sense_all

      CLOSE(ICD_outf)





    ELSE calc_them_all

      WRITE(of,*) 'Processing values of M_A prime separately'
      WRITE(of,*) 'All decay rates are given in eV'
      WRITE(of,*) ''
      WRITE(of,410) 'J_A','M_A',"J_A'","M_A'","j_B'",'no pairs','R [AA]','Gamma one',&
                    'Gamma all'
      410 FORMAT (' ',5(1X,A4),2(2X,A9),2(4X,A9))
      WRITE(of,*) '-------------------------------------------------------------------------'

!Set the M_Ap dependent variables
      SELECT CASE (INT(M_Ap - M_A))
        CASE (-1,1)
          B_MAMAp = 1
        CASE (0)
          B_MAMAp = -2
        CASE DEFAULT
          B_MAMAp = 0
      END SELECT
!      WRITE(of,*) 'B_MAMAp= ', B_MAMAp

      wigner   = eval_wigner3j(J_Ap,one,J_A,-M_Ap,M_Ap-M_A,M_A)
      WRITE(of,*) 'Wigner3j symbol: ', wigner

! Open the specfile for output
      IF (INT(M_Ap) >= 0) THEN
        WRITE(specfile, '(A4,3(I1,A1),I1)') 'ICD_', INT(2*J_A), '_', INT(2*J_Ap)&
                                          & ,'_', INT(2*M_Ap),'_', INT(2*j_Bp)
      ELSE
        WRITE(specfile, '(A4,2(I1,A1),I2,A1,I1)') 'ICD_', INT(2*J_A), '_', INT(2*J_Ap)&
                                                  & ,'_', INT(2*M_Ap),'_', INT(2*j_Bp)
      END IF
      OPEN(ICD_outf,FILE=TRIM(ADJUSTL(specfile)),STATUS='UNKNOWN',ACTION='WRITE'&
          &,IOSTAT=ierror)


! Test whether this channel makes sense at all
      channel_sense:IF (E_in - E_fin1 - E_fin2 > 0) THEN
  
!        WRITE(of,*) 'Processing ICD channel ', ichannel

        DO idist=1,no_dist

! Find values for given pair
          neq_pairs  = INT(dist_stat(idist,1))
          R_angstrom = dist_stat(idist,2)
!          WRITE(of,*) 'R_ang = ', R_angstrom
          R_bohr     = R_angstrom * angstrom_to_bohr
    
          E_Coulomb  = 1/R_bohr * hartree_to_ev
          E_sec      = E_in - E_fin1 - E_fin2 - E_Coulomb

! Verfahre nur weiter, wenn die Sekundaerenergie >=0 ist
          IF (E_sec >= 0.0) THEN
!            WRITE(of,*) 'E_sec= ', E_sec
            gamma_b = 2*pi/R_bohr**6 * B_MAMAp**2 * wigner**2 *(2*J_A+1)&
                    & *3*c_au**4 *sigma_au/(16*pi**2 * omega_vp**4 * tau_au) * hartree_to_ev&
                    & / number_of_in !Normalize to one ionization
!            WRITE(of,*) 'Gamma beta = ', gamma_b
            gamma_b_pairs = neq_pairs * gamma_b
            gamma_b_all_pairs = gamma_b_all_pairs + gamma_b_pairs

! Write summary to output file
            WRITE(of,420) INT(2*J_A), INT(2*M_A), INT(2*J_Ap),INT(2*M_Ap),&
                          INT(2*j_Bp), INT(neq_pairs), R_angstrom,&
                          gamma_b, gamma_b_pairs
            420 FORMAT (' ',5(1X,I4),4X,I5,6X,F7.3,4X,2(ES9.3,4X))

!Write result to specfile
            WRITE(ICD_outf,141) E_sec, gamma_b_pairs
!            141 FORMAT (' ',F12.4,ES15.5) is already defined in all
          END IF
        END DO

      WRITE(of,*) '-------------------------------------------------------------------------'
      WRITE(of,430) gamma_b_all_pairs
      430 FORMAT (' ',64X,ES9.3)

      END IF channel_sense

      CLOSE(ICD_outf)

   END IF calc_them_all

  END DO each_channel

!    WRITE(of,*) 'Sum of all Gammas for this channel and geometry ', gamma_b_all_pairs


END SUBROUTINE calc_icd_gamma

END MODULE icd_calc
