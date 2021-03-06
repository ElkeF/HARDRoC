MODULE control
!Default values of the variables
LOGICAL      :: do_pairs = .FALSE.
LOGICAL      :: do_fit   = .FALSE.
LOGICAL      :: do_eff   = .FALSE.
LOGICAL      :: do_triples = .FALSE.
LOGICAL      :: do_RICD = .FALSE.
CHARACTER(2) :: in_atom_type = '', fin_atom_type1 = '', fin_atom_type2 = ''

! Parameters for reducing numerical effort
REAL :: Qmax = 1000.0 !Maximum distance of Q which will be taken into account
REAL :: Rmax = 1000.0! Maximum ditance for R in ICD

!Output file variables
! Data dictionary: io number start of output files
!CHARACTER(len=30) :: nmout !Name of the output file
INTEGER,PARAMETER :: of = 30
INTEGER :: ICD_outf = of+1
INTEGER :: ETMD_outf = of+20
INTEGER :: dist_outf = of+30
END MODULE control
