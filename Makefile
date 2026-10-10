# Makefile de la factorización LU paralela  (make / make clean)
# Se compila sin optimización, igual que en el TP1, para medir el algoritmo y no el compilador.

CC     = gcc
MPICC  = mpicc
CFLAGS =
LIBS   = -lm        # necesario por sqrt()/fabs() de math.h

all: lu_serial lu_mpi

lu_serial: lu_serial.c
	$(CC) $(CFLAGS) lu_serial.c -o lu_serial $(LIBS)

lu_mpi: lu_mpi.c
	$(MPICC) $(CFLAGS) lu_mpi.c -o lu_mpi $(LIBS)

clean:
	rm -f lu_serial lu_mpi

.PHONY: all clean
