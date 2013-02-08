MODULE select_trdm_para

CONTAINS

SUBROUTINE select_trdm_fit_para(factor,alpha,const,dir)
! Purpose: select the correct parameters of an exponential fit of
! the transition dipole moments of molecules.

  use control
  use channel_char
  use trdm_library

  IMPLICIT NONE

! Data dictionary: output parameters
  REAL, INTENT(OUT) :: factor,alpha,const
  CHARACTER(len=1), INTENT(OUT) :: dir !in which direction

! Data dictionary: 

  inatom: SELECT CASE (in_atom_type)
  
    CASE('Ar')
   
      IF ((INT(2*J_A) == 1).AND.(INT(2*M_A) == 1)) THEN

        donor:SELECT CASE (fin_atom_type1)

          CASE('Xe')
            SELECT CASE (INT(2*J_D))
              CASE(3)
                SELECT CASE(INT(2*M_D))
                  CASE(3)
                    factor = ArXe_1_1_3_3_factor
                    alpha  = ArXe_1_1_3_3_alpha
                    const  = ArXe_1_1_3_3_const
                    dir    = 'x'
                  CASE(1)
                    factor = ArXe_1_1_3_1_factor
                    alpha  = ArXe_1_1_3_1_alpha
                    const  = ArXe_1_1_3_1_const
                    dir    = 'z'
                  CASE(-1)
                    factor = ArXe_1_1_3_m1_factor
                    alpha  = ArXe_1_1_3_m1_alpha
                    const  = ArXe_1_1_3_m1_const
                    dir    = 'x'
                  CASE DEFAULT
                    WRITE(of,*) 'No TRDM provided for this M_D.'
                    factor = 0
                    alpha  = 0
                    const  = 0
                END SELECT
              CASE(1)
                SELECT CASE(INT(2*M_D))
                  CASE(1)
                    factor = ArXe_1_1_1_1_factor
                    alpha  = ArXe_1_1_1_1_alpha
                    const  = ArXe_1_1_1_1_const
                    dir    = 'z'
                  CASE(-1)
                    factor = ArXe_1_1_1_m1_factor
                    alpha  = ArXe_1_1_1_m1_alpha
                    const  = ArXe_1_1_1_m1_const
                    dir    = 'x'
                  CASE DEFAULT
                    WRITE(of,*) 'No TRDM provided for this M_D.'
                    factor = 0
                    alpha  = 0
                    const  = 0
                END SELECT
              CASE DEFAULT
                WRITE(of,*) 'No TRDM provided for this value of J_D.'
            END SELECT
            
          CASE DEFAULT
            WRITE(of,*) 'No TRDM provided for this donor atom type'

        END SELECT donor

      END IF

    CASE DEFAULT
      WRITE(of,*) 'No TRDM provided for this atom type'

  END SELECT inatom

END SUBROUTINE select_trdm_fit_para

END MODULE select_trdm_para
