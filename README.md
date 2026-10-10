# Proyecto 6 — Factorización LU Paralela (A = LU)

Factorización LU con pivoteo parcial (eliminación gaussiana), en versión **serial**
y **paralela con MPI**, con medición de **speedup** y del **error residual ‖Ax − b‖**.

- **Integrantes:** Naessens, Maia Guadalupe · Martínez, Ezequiel — Ingeniería en Informática
- **Materia:** Métodos Numéricos II — Computación Científica · FACET, UNT · 2026

---

## Consigna

> **Proyecto 6: Factorización LU Paralela (A = LU).** La factorización LU descompone
> una matriz cuadrada invertible *A* en el producto de una matriz triangular inferior
> *L* y una matriz triangular superior *U*, constituyendo la base del cálculo de
> sistemas de ecuaciones lineales *Ax = b*. Investigue la implementación paralela de
> la eliminación gaussiana. Analice cómo escala el tiempo de cómputo y la comunicación
> al incrementar N, evaluando el *speedup* y la precisión numérica del error residual
> ‖Ax − b‖.

---

## Estructura del proyecto

| Archivo / carpeta | Para qué |
|---|---|
| `lu_serial.c` | Versión **serial** de referencia (1 proceso, sin MPI): factoriza, resuelve *Ax = b* y calcula el residual. |
| `lu_mpi.c` | Versión **paralela con MPI** (distribución cíclica por filas). Soporta además reparto por bloques y la medición separada de cómputo / comunicación / idle. |
| `medir.sh` | Script de **mediciones**. Subcomandos: `pc` y `cluster` (speedup y eficiencia) y `desglose` (cómputo/comunicación/idle). Compila, corre varias veces, toma la mediana y genera un CSV. |
| `graficar.py` | Genera las **figuras** a partir de un CSV (detecta el tipo por sus columnas). |
| `Makefile` | Compila ambos programas (opcional; `medir.sh` ya compila). |
| `resultados/Cluster/`, `resultados/pc/` | CSV y figuras de cada equipo. |
| `informe/` | Informe final en Word (`.docx`). |

---

## Requisitos

- Linux (o WSL en Windows), con **OpenMPI** (`mpicc`, `mpirun`) y `gcc`.
- **Python 3** con `matplotlib` (`pip install matplotlib`) para las figuras.

---

## Cómo correr

**Compilar** (o dejar que lo haga `medir.sh`):

```bash
make
```

**Una corrida suelta** (la dimensión N va como argumento):

```bash
./lu_serial 1000
mpirun -n 4 ./lu_mpi 1000
```

**Mediciones** (generan un CSV) y **figuras**:

```bash
./medir.sh pc                                 # notebook: N=100..800, p=1..6
./medir.sh cluster                            # cluster:  N=4000,8000,12000, p hasta 48
./medir.sh desglose 3000 "2 4 8 16 32 48"     # cómputo/comunicación/idle (cíclico vs bloques)
python3 graficar.py resultados/Cluster/resultados.csv
```

Opciones por variable de entorno: `REP=20 ./medir.sh pc` (repeticiones por punto) y
`MAXP=48 ./medir.sh cluster` (límite de procesos).

El CSV de speedup tiene columnas `n, p, mediana_s, speedup, eficiencia, residual_abs,
residual_rel`. El tiempo de referencia T₁ para el speedup es el de `lu_serial` (serial
puro). Las figuras se guardan junto al CSV, con el nombre del CSV como prefijo (así
distintas corridas no se pisan).
