MODULE channel_char

IMPLICIT NONE
! all variables characterising the different ICD, ETMD
! channels
! The default parameters are set to 8.8, which will make it easy to detect
! errors due to wrong input reading

! Initial state properties
REAL :: J_A = 8.8 !Total angular momentum of the initial state, alternatively L_A
REAL :: M_A = 8.8 !Projection of J_A
REAL :: SIP_in = 0.0 !SIP of initial state
REAL :: shift_in = 0.0 !shift due to the surrounding of SIP(in)

! Final state properties
REAL :: J_Ap = 8.8 !Total angular momentum of fin1
REAL :: M_Ap = 8.8 !Projection of J_Ap
REAL :: SIP_fin1 = 0.0 !SIP of final state 1
REAL :: shift_fin1 = 0.0 !shift of fin1 due to surroundings
REAL :: tottau = 8.8 !total tau
REAL :: taurel = 8.8 !relation of taus, tau_lowerSIP/tau_higherSIP

REAL :: j_Bp = 8.8 !Total angular momentum of atom fin2
REAL :: SIP_fin2 = 0.0 !SIP of final state 2
REAL :: shift_fin2 = 0.0 !Shift due to surroundings
REAL :: sigmaabs = 0.0 !ionization cross section of atom B
REAL :: sigmarel = 0.0 !relation of sigma, sigma_lowerSIP/sigma_higherSIP

! ETMD properties
REAL :: J_D = 8.8 ! Total angular momentum of the donor atom in the e transition
REAL :: M_D = 8.8 ! Projection of J_D

! Coefficients for fitting procedure
REAL(8)  :: c_pre = 0.0
REAL(8)  :: c_exp = 0.0
REAL(8)  :: c_6   = 0.0

END MODULE
