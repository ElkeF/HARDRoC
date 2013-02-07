MODULE array_operations

CONTAINS

SUBROUTINE shrink_rows(inarray,outarray,rows_in,no_out,columns)
! Purpose: to get some statistics on the differnt type of entries
! and to reduce the number of lines without loosing information

  use control
  use geometry_fun

  IMPLICIT NONE

! Data dictionary: Input parameters
  INTEGER, INTENT(IN) :: rows_in, columns, no_out
  REAL, DIMENSION(rows_in,columns), INTENT(IN) :: inarray

! Data dictionary: Output parameters
  REAL, DIMENSION(no_out,columns+1), INTENT(OUT) :: outarray

! Data dictionary: Temporary variables
  INTEGER :: row, ne_entries
  REAL, DIMENSION(rows_in,columns) :: temparray
  REAL, DIMENSION(1,columns) :: line1, line2, diff
  REAL :: thresh = 1E-4
  REAL :: norm

! Data dictionary: Counters
  INTEGER :: i,j

  temparray = inarray
  row = 0

  dorows:DO i=1,rows_in

    line1 = temparray(i:i,1:columns)

    IF (norm_row(line1,columns) > thresh) THEN
      row = row + 1 
      ne_entries = 1
!    WRITE(of,*) 'Line1: ', line1
!    WRITE(of,*) 'row = ', row
!    WRITE(of,*) 'ne_entries= ', ne_entries

! Zero out all equivalent entries
      DO j=i+1,rows_in
      
        line2 = temparray(j:j,1:columns)
!        WRITE(of,*) 'Line 2: ', line2

        IF (norm_row(line2,columns) > thresh) THEN

          diff  = line1-line2
          norm  = norm_row(diff,columns)
!          WRITE(of,*) 'diff ', diff
!          WRITE(of,*) 'line2 ', line2
!          WRITE(of,*) 'Difference in "position": ', norm

          IF (norm < thresh) THEN
    
            temparray(j:j,1:columns) = 0.0
            ne_entries = ne_entries + 1
!            WRITE(of,*) 'norm ', norm
!            WRITE(of,*) 'ne_entries ', ne_entries

          END IF
        END IF
      END DO
! Store result in outarray
      outarray(row,1) = ne_entries
      outarray(row,2:columns+1) = line1(1,1:columns)

      WRITE(of,*) 'outarray: ', outarray(row,1:columns+1)
    END IF
  END DO dorows

 

END SUBROUTINE shrink_rows





INTEGER FUNCTION diff_rows(array,rows,columns)
! Purpose: to determine the number of different rows in an array

  use geometry_fun

  IMPLICIT NONE

! Data dictionary: Input data
  INTEGER, INTENT(IN) :: rows,columns
  REAL, DIMENSION(rows,columns), INTENT(IN) :: array

! Data dictionary: Temporary variables
  REAL, DIMENSION(rows,columns) :: temparray
  REAL, DIMENSION(1,columns) :: line1, line2, diff
  REAL :: thresh = 1E-4
  REAL :: norm

! Data dictionary: Counters
  INTEGER :: i,j

  temparray = array
  diff_rows = 0

  dorows:DO i=1,rows

    line1 = temparray(i:i,1:columns)

    IF (norm_row(line1,columns) > thresh) THEN
      diff_rows = diff_rows + 1 ! we have one more possibility

! Zero out all equivalent entries
      DO j=i+1,rows
      
        line2 = temparray(j:j,1:columns)

        IF (norm_row(line2,columns) > thresh) THEN

          diff  = line1-line2
          norm  = norm_row(diff,columns)

          IF (norm < thresh) THEN
     
            temparray(j:j,1:columns) = 0.0

          END IF
        END IF
      END DO
    END IF
  END DO dorows

END FUNCTION



END MODULE array_operations