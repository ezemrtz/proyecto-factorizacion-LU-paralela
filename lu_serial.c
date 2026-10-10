// FACTORIZACIÓN LU POR FILAS CON PIVOTEO PARCIAL - VERSIÓN SERIAL

#include <stdio.h>
#include <stdlib.h>
#include <math.h>   
#include <time.h>   

#define SEMILLA 42  

// ============================== 1. FUNCIONES AUXILIARES ==============================

// Carga la matriz inicial con valores aleatorios en [-1, 1]
void inicializarMatriz(int n, double A[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
        }
    }
}

// Carga el vector de términos independientes 'b' con valores aleatorios en [-1, 1]
void inicializarVector(int n, double b[n]) {
    for (int i = 0; i < n; i++) {
        b[i] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
    }
}

// Muestra la matriz compacta formateada en pantalla (solo para N chico)
void imprimirMatriz(int n, double A[n][n]) {
    printf("Matriz compacta A (L debajo de la diagonal, U en la diagonal y arriba):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%9.4f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Busca en la columna 'k' el elemento con mayor valor absoluto (pivoteo parcial)
int buscarPivote(int n, double A[n][n], int k) {
    int fila_max = k;
    double max_val = fabs(A[k][k]);
    for (int i = k + 1; i < n; i++) {
        if (fabs(A[i][k]) > max_val) {
            max_val = fabs(A[i][k]);
            fila_max = i;
        }
    }
    return fila_max;
}

// Intercambia dos filas completas de la matriz
void intercambiarFilas(int n, double A[n][n], int f1, int f2) {
    for (int j = 0; j < n; j++) {
        double temp = A[f1][j];
        A[f1][j] = A[f2][j];
        A[f2][j] = temp;
    }
}

// Realiza el cálculo matemático de una fila:
// Guarda el multiplicador in-place y actualiza los elementos restantes
void actualizarFila(int n, double fila_actual[n], double fila_pivote[n], int k) {
    fila_actual[k] /= fila_pivote[k];
    for (int j = k + 1; j < n; j++) {
        fila_actual[j] -= fila_actual[k] * fila_pivote[j];
    }
}

// Resuelve A x = b a partir de la LU compacta y el vector de pivotes.
void resolverSistema(int n, double A[n][n], int pivotes[n], double b[n], double x[n]) {
    double *y = malloc(n * sizeof(double));

    // Aplica a 'b' los mismos intercambios de filas que tuvo A: P b
    for (int k = 0; k < n - 1; k++) {
        if (pivotes[k] != k) {
            double temp = b[k];
            b[k] = b[pivotes[k]];
            b[pivotes[k]] = temp;
        }
    }

    // Sustitución hacia adelante: L y = P b  
    for (int i = 0; i < n; i++) {
        double suma = b[i];
        for (int j = 0; j < i; j++) {
            suma -= A[i][j] * y[j];
        }
        y[i] = suma;
    }

    // Sustitución hacia atrás: U x = y
    for (int i = n - 1; i >= 0; i--) {
        double suma = y[i];
        for (int j = i + 1; j < n; j++) {
            suma -= A[i][j] * x[j];
        }
        x[i] = suma / A[i][i];
    }

    free(y);
}

// Calcula el error residual ||A x - b|| sobre la matriz y el vector ORIGINALES. Devuelve el residual absoluto y, por puntero, el relativo.
double calcularResidual(int n, double A[n][n], double b[n], double x[n], double *residual_rel) {
    double norma_r = 0.0, norma_b = 0.0;
    for (int i = 0; i < n; i++) {
        double fila = 0.0;
        for (int j = 0; j < n; j++) {
            fila += A[i][j] * x[j];
        }
        double r_i = fila - b[i];
        norma_r += r_i * r_i;
        norma_b += b[i] * b[i];
    }
    norma_r = sqrt(norma_r);
    *residual_rel = norma_r / sqrt(norma_b);
    return norma_r;
}

// ============================== 2. FUNCIÓN PRINCIPAL ==============================

int main(int argc, char *argv[]) {
    struct timespec t_inicio, t_fin;

    // Dimensión de la matriz: argv[1] (por defecto 4)
    int n = (argc > 1) ? atoi(argv[1]) : 4;

    // Reserva dinámica
    double (*A)[n]      = malloc(sizeof(double[n][n])); // se sobrescribe con la LU
    double (*A_orig)[n] = malloc(sizeof(double[n][n])); // copia para el residual
    double *b      = malloc(n * sizeof(double));        // término independiente original
    double *pb     = malloc(n * sizeof(double));        // copia de b para resolver la cual se permuta
    double *x      = malloc(n * sizeof(double));        // solución
    int    *pivotes = malloc(n * sizeof(int));          // registro de intercambios de filas

    // Datos reproducibles y copias de respaldo
    srand(SEMILLA);
    inicializarMatriz(n, A);
    inicializarVector(n, b);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) A_orig[i][j] = A[i][j];
        pb[i] = b[i];
    }

    // ------------- FACTORIZACIÓN LU CON PIVOTEO PARCIAL -------------
    clock_gettime(CLOCK_MONOTONIC, &t_inicio);

    for (int k = 0; k < n - 1; k++) {
        // Pivoteo parcial
        int fila_max = buscarPivote(n, A, k);
        pivotes[k] = fila_max;
        if (fila_max != k) {
            intercambiarFilas(n, A, k, fila_max);
        }
        // Actualiza todas las filas por debajo de la diagonal
        for (int i = k + 1; i < n; i++) {
            actualizarFila(n, A[i], A[k], k);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t_fin);
    double tiempo = (t_fin.tv_sec - t_inicio.tv_sec)
                  + (t_fin.tv_nsec - t_inicio.tv_nsec) / 1e9;

    // ------------- RESOLUCIÓN DEL SISTEMA Y ERROR RESIDUAL -------------
    resolverSistema(n, A, pivotes, pb, x);
    double residual_rel;
    double residual_abs = calcularResidual(n, A_orig, b, x, &residual_rel);

    // Salida legible (N chico) o compacta (N grande)
    if (n <= 8) {
        imprimirMatriz(n, A);
        printf("Solución x:\n");
        for (int i = 0; i < n; i++) printf("  x[%d] = %12.6f\n", i, x[i]);
        printf("\n");
    }
    printf("N = %d | tiempo = %.9f s | ||Ax-b|| = %.6e | relativo = %.6e\n",
           n, tiempo, residual_abs, residual_rel);

    free(A); free(A_orig); free(b); free(pb); free(x); free(pivotes);
    return 0;
}
