      !Author : Louisda16th a.k.a Ashwith J. Rego
      !These set of subroutines find the largest eigenvalue and eigenmatrix of the matrix.
      !The algorithm is based on Rayleigh's power method
      !Please note that the subroutine used to multiply the two matrices is not general.
      !Also note that the number of iterations must be specified
      SUBROUTINE FINDEIGEN (Matrix, n, x, EigenVal, steps)
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: n, steps  !n = order of matrix, steps = number of iterations
          REAL, INTENT(IN), DIMENSION(n,n) :: Matrix(n,n)  !Input Matrix
          REAL, INTENT(INOUT), DIMENSION(n) :: x !Eigenvector
          REAL, INTENT(INOUT) :: EigenVal !Eigenvalue
          INTEGER :: i, j
           
          x  = 1 !Initialize eigen vector to any value.
           
          DO i = 1, steps
              CALL MULMATRIX (Matrix, x, n)       !Multiply input matrix by eigenvector
              CALL FINDLARGEST (x, n, EigenVal)   !Find eigenvalue
              IF(EigenVal == 0) EXIT         
              DO j = 1, n                        !Find eigenvector
                  x(j) = x(j)/EigenVal
              END DO 
          END DO
       
      END SUBROUTINE FINDEIGEN
       
      SUBROUTINE MULMATRIX (a, b, n)
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: n !matrix size
          REAL, INTENT(IN), DIMENSION(n,n) :: a  !Matrix of order > 1
          REAL, INTENT(INOUT), DIMENSION(n) :: b !1x1 matrix
           
          INTEGER i, j
          REAL, DIMENSION(n) :: temp !temporary matrix
       
          temp = 0
           
          !These two loops to the multiplication
          DO i = 1, n
              DO j = 1, n
                  temp(i) = temp(i) + a(i,j)*b(j)
              END DO
          END DO
          b = temp
       
      END SUBROUTINE MULMATRIX
       
      SUBROUTINE FINDLARGEST (x, n, l)
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: n
          REAL, INTENT(IN), DIMENSION(n) :: x
          REAL, INTENT(INOUT) :: l !Largest value
           
          INTEGER :: i
          !Algorithm is easy
          !Let the largest number be the first one.
          !If you find a number larger than it, store this number and then continue
          l = ABS(x(1))
          DO i = 2, n
              IF (ABS(x(i)) > l) l = ABS(x(i))
          END DO
               
      END SUBROUTINE FINDLARGEST
