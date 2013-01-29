MODULE control
!Default values of the variables
LOGICAL      :: do_pairs = .FALSE.
LOGICAL      :: do_triples = .FALSE.
CHARACTER(2) :: in_atom_type = '', fin_atom_type1 = '', fin_atom_type2 = ''
!REAL         :: C = 300.0      !Speed of light
!REAL         :: SIPIN = 0.0    !Single ionization potential of the initial state
!REAL         :: SIP1 = 0.0     !First single ionization potential of the final state
!REAL         :: SIP2 = 0.0     !Second single ionization potential of the final state
!REAL         :: tottau = 0.0   !Total photonemission lifetime
!REAL         :: relatau = 0.0  !Relation between two taus
!REAL         :: trdm = 0.0     !Transition dipole moment for the electron transfer
END MODULE control
