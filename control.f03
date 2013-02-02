MODULE control
!Default values of the variables
LOGICAL      :: do_pairs = .FALSE.
LOGICAL      :: do_triples = .FALSE.
CHARACTER(2) :: in_atom_type = '', fin_atom_type1 = '', fin_atom_type2 = ''

!Output file variables
! Data dictionary: io number start of output files
!CHARACTER(len=30) :: nmout !Name of the output file
INTEGER :: of = 30
END MODULE control
