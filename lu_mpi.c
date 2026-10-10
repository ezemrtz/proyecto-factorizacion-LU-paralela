// FACTORIZACIÓN LU POR FILAS CON PIVOTEO PARCIAL - VERSIÓN MPI

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define SEMILLA 42

// ============================== 1. DESCRIPTOR DE DISTRIBUCIÓN ==============================
// modo = 0 -> cíclico;  modo = 1 -> bloques
typedef struct {
    int n, p, rank, modo;
    int filas_por_bloque;   
} dist_t;

// Proceso dueño de la fila global g
static inline int dueno(const dist_t *d, int g) {
    return d->modo ? g / d->filas_por_bloque : g % d->p;
}
// Índice local de la fila global g
static inline int indice_local(const dist_t *d, int g) {
    return d->modo ? g - (g / d->filas_por_bloque) * d->filas_por_bloque : g / d->p;
}
// Índice global de la fila local il de este proceso
static inline int indice_global(const dist_t *d, int il) {
    return d->modo ? d->rank * d->filas_por_bloque + il : d->rank + il * d->p;
}

// ============================== 2. FUNCIONES AUXILIARES ==============================

// Carga la matriz inicial con valores aleatorios en [-1, 1]
void inicializarMatriz(int n, double A[n][n]) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
}

// Carga el vector de términos independientes 'b' con valores aleatorios en [-1, 1]
void inicializarVector(int n, double b[n]) {
    for (int i = 0; i < n; i++)
        b[i] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
}

// Muestra la matriz compacta formateada en pantalla (solo para N chico)
void imprimirMatriz(int n, double A[n][n]) {
    printf("Matriz compacta A (L debajo de la diagonal, U en la diagonal y arriba):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) printf("%9.4f ", A[i][j]);
        printf("\n");
    }
    printf("\n");
}

// Eliminación de una fila: guarda el multiplicador in-place y actualiza el resto
void actualizarFila(int n, double fila_actual[n], double fila_pivote[n], int k) {
    fila_actual[k] /= fila_pivote[k];
    for (int j = k + 1; j < n; j++)
        fila_actual[j] -= fila_actual[k] * fila_pivote[j];
}

// Resuelve A x = b a partir de la LU compacta y el vector de pivotes.
void resolverSistema(int n, double A[n][n], int pivotes[n], double b[n], double x[n]) {
    double *y = malloc(n * sizeof(double));
    for (int k = 0; k < n - 1; k++)                 // aplica P a b
        if (pivotes[k] != k) {
            double t = b[k]; b[k] = b[pivotes[k]]; b[pivotes[k]] = t;
        }
    for (int i = 0; i < n; i++) {                   // L y = P b 
        double suma = b[i];
        for (int j = 0; j < i; j++) suma -= A[i][j] * y[j];
        y[i] = suma;
    }
    for (int i = n - 1; i >= 0; i--) {              // U x = y 
        double suma = y[i];
        for (int j = i + 1; j < n; j++) suma -= A[i][j] * x[j];
        x[i] = suma / A[i][i];
    }
    free(y);
}

// Error residual ||A x - b|| sobre la matriz y el vector ORIGINALES
double calcularResidual(int n, double A[n][n], double b[n], double x[n], double *residual_rel) {
    double norma_r = 0.0, norma_b = 0.0;
    for (int i = 0; i < n; i++) {
        double fila = 0.0;
        for (int j = 0; j < n; j++) fila += A[i][j] * x[j];
        double r_i = fila - b[i];
        norma_r += r_i * r_i;
        norma_b += b[i] * b[i];
    }
    norma_r = sqrt(norma_r);
    *residual_rel = norma_r / sqrt(norma_b);
    return norma_r;
}

// ============================== 3. FUNCIÓN PRINCIPAL ==============================

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Status estado;
    double t_inicio, t_fin, ts;
    struct { double val; int idx; } pivote_local, pivote_global;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Argumentos: N, modo (0 cíclico / 1 bloques), desglose (0 / 1)
    int n        = (argc > 1) ? atoi(argv[1]) : 4;
    int modo     = (argc > 2) ? atoi(argv[2]) : 0;
    int desglose = (argc > 3) ? atoi(argv[3]) : 0;

    dist_t d = { n, size, rank, modo, (n + size - 1) / size };

    // Cantidad de filas locales de este proceso
    int filas_locales = 0;
    for (int g = 0; g < n; g++) if (dueno(&d, g) == rank) filas_locales++;

    // Cada proceso guarda SOLO sus filas
    int maxfilas = d.filas_por_bloque;
    double (*local)[n]  = malloc(sizeof(double[maxfilas > 0 ? maxfilas : 1][n]));
    double *fila_pivote = malloc(n * sizeof(double));

    // Solo el proceso 0 reserva la matriz completa y los vectores
    double (*A)[n] = NULL, (*A_orig)[n] = NULL;
    double *b = NULL, *pb = NULL, *x = NULL;
    int *pivotes = NULL;
    if (rank == 0) {
        A = malloc(sizeof(double[n][n]));  A_orig = malloc(sizeof(double[n][n]));
        b = malloc(n*sizeof(double)); pb = malloc(n*sizeof(double)); x = malloc(n*sizeof(double));
        pivotes = malloc(n * sizeof(int));
        srand(SEMILLA);
        inicializarMatriz(n, A);
        inicializarVector(n, b);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) A_orig[i][j] = A[i][j];
            pb[i] = b[i];
        }
    }

    // ---- Reparto de filas (fuera de la medición) ----
    for (int g = 0; g < n; g++) {
        int due = dueno(&d, g);
        if (rank == 0) {
            if (due == 0) for (int j = 0; j < n; j++) local[indice_local(&d, g)][j] = A[g][j];
            else          MPI_Send(A[g], n, MPI_DOUBLE, due, 40, MPI_COMM_WORLD);
        } else if (rank == due) {
            MPI_Recv(local[indice_local(&d, g)], n, MPI_DOUBLE, 0, 40, MPI_COMM_WORLD, &estado);
        }
    }

    // Contadores del desglose de tiempos
    double c_cmp = 0.0, c_com = 0.0, c_idl = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    t_inicio = MPI_Wtime();

    // ------------- BUCLE PRINCIPAL DE FACTORIZACIÓN LU -------------
    for (int k = 0; k < n - 1; k++) {

        // --- 1. Pivoteo parcial: máximo local (CÓMPUTO) ---
        ts = MPI_Wtime();
        pivote_local.val = -1.0; pivote_local.idx = k;
        for (int r = 0; r < filas_locales; r++) {
            int g = indice_global(&d, r);
            if (g < k) continue;
            double v = fabs(local[r][k]);
            if (v > pivote_local.val) { pivote_local.val = v; pivote_local.idx = g; }
        }
        c_cmp += MPI_Wtime() - ts;

        // --- 2. Reducción global del pivote (COMUNICACIÓN) ---
        if (desglose) { ts = MPI_Wtime(); MPI_Barrier(MPI_COMM_WORLD); c_idl += MPI_Wtime() - ts; }
        ts = MPI_Wtime();
        MPI_Allreduce(&pivote_local, &pivote_global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        c_com += MPI_Wtime() - ts;
        int fila_max = pivote_global.idx;
        if (rank == 0) pivotes[k] = fila_max;

        // --- 3. Intercambio de las filas k y fila_max (COMUNICACIÓN) ---
        ts = MPI_Wtime();
        if (fila_max != k) {
            int due_k = dueno(&d, k), due_max = dueno(&d, fila_max);
            if (due_k == due_max) {
                if (rank == due_k) {
                    int lk = indice_local(&d, k), lm = indice_local(&d, fila_max);
                    for (int j = 0; j < n; j++) { double t = local[lk][j]; local[lk][j] = local[lm][j]; local[lm][j] = t; }
                }
            } else if (rank == due_k) {
                MPI_Sendrecv_replace(local[indice_local(&d, k)], n, MPI_DOUBLE, due_max, 50, due_max, 50, MPI_COMM_WORLD, &estado);
            } else if (rank == due_max) {
                MPI_Sendrecv_replace(local[indice_local(&d, fila_max)], n, MPI_DOUBLE, due_k, 50, due_k, 50, MPI_COMM_WORLD, &estado);
            }
        }
        c_com += MPI_Wtime() - ts;

        // --- 4. Difusión de la fila pivote (COMUNICACIÓN) ---
        int due_k = dueno(&d, k);
        if (rank == due_k) for (int j = 0; j < n; j++) fila_pivote[j] = local[indice_local(&d, k)][j];
        if (desglose) { ts = MPI_Wtime(); MPI_Barrier(MPI_COMM_WORLD); c_idl += MPI_Wtime() - ts; }
        ts = MPI_Wtime();
        MPI_Bcast(fila_pivote, n, MPI_DOUBLE, due_k, MPI_COMM_WORLD);
        c_com += MPI_Wtime() - ts;

        // --- 5. Actualización de las filas locales g > k (CÓMPUTO) ---
        ts = MPI_Wtime();
        for (int r = 0; r < filas_locales; r++) {
            int g = indice_global(&d, r);
            if (g <= k) continue;
            actualizarFila(n, local[r], fila_pivote, k);
        }
        c_cmp += MPI_Wtime() - ts;
    }

    t_fin = MPI_Wtime();

    // ---- Reunión de filas en el proceso 0 (fuera de la medición) ----
    for (int g = 0; g < n; g++) {
        int due = dueno(&d, g);
        if (rank == 0) {
            if (due == 0) for (int j = 0; j < n; j++) A[g][j] = local[indice_local(&d, g)][j];
            else          MPI_Recv(A[g], n, MPI_DOUBLE, due, 60, MPI_COMM_WORLD, &estado);
        } else if (rank == due) {
            MPI_Send(local[indice_local(&d, g)], n, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD);
        }
    }

    // ---- Promedio de los tiempos por componente entre procesos ----
    double s_cmp = 0, s_com = 0, s_idl = 0;
    MPI_Reduce(&c_cmp, &s_cmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&c_com, &s_com, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&c_idl, &s_idl, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // ------------- RESOLUCIÓN, RESIDUAL Y SALIDA (solo el proceso 0) -------------
    if (rank == 0) {
        resolverSistema(n, A, pivotes, pb, x);
        double residual_rel;
        double residual_abs = calcularResidual(n, A_orig, b, x, &residual_rel);

        if (n <= 8) imprimirMatriz(n, A);

        // tiempo  residual_abs  residual_rel  [t_comp t_comm t_idle]  (promedios)
        printf("%.9f %.6e %.6e", t_fin - t_inicio, residual_abs, residual_rel);
        if (desglose)
            printf(" %.9f %.9f %.9f", s_cmp / size, s_com / size, s_idl / size);
        printf("\n");

        free(A); free(A_orig); free(b); free(pb); free(x); free(pivotes);
    }

    free(local);
    free(fila_pivote);
    MPI_Finalize();
    return 0;
}
