      SUBROUTINE MAKE_LAYER_FILE(LAYF, SCATF)
C Создает файл с оаписанием атмосферы - входной файл для программы
C PolRadTran.
      IMPLICIT NONE
      CHARACTER*64  LAYF, SCATF
Cf2PY INTETN(IN)  LAYF, SCATF

      OPEN(UNIT=200, FILE=LAYF, STATUS='UNKNOWN')
      WRITE(200, "(2F7.2,F7.3,A,A)") 1.0, 0.0, 0.0, '   ', 
     .                                "'"//TRIM(SCATF)//"'"
      WRITE(200, "(2F7.2,F7.3,A,A)") 0.0, 0.0, 0.0, '   ', "'         '"  
      CLOSE(200)

      END