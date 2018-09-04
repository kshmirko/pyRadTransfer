      FUNCTION TRAPZ(X, Y, N)
C Интегрирование функции, заданной векторами X, Y по методу трапеций
CF2PY INTENT(IN)  X, Y, N
      INTEGER N
      REAL*8  X(N), Y(N), TRAPZ
      
      TRAPZ=SUM((Y(1+1:N-0) + Y(1+0:N-1))*(X(1+1:N-0) - X(1+0:N-1)))/2.0
      
      END 