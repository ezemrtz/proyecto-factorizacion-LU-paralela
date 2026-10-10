#!/bin/bash
# Mediciones de la factorización LU paralela.
#   ./medir.sh pc  [salida.csv]            Speedup/eficiencia, preset notebook (p 1..6)
#   ./medir.sh cluster [salida.csv]        Speedup/eficiencia, preset cluster (p hasta 48)
#   ./medir.sh desglose <N> "<lista_p>"    Cómputo/comunicación/idle, cíclico vs bloques
# Variables de entorno: REP=<n> repeticiones (por defecto 5), MAXP=<n> límite de procesos.
# No se corren más procesos que núcleos físicos: sin oversubscribe el speedup es real.

SUB=${1:-pc}
REP=${REP:-5}

# --- Compilación ---
echo "Compilando lu_mpi.c y lu_serial.c ..."
mpicc lu_mpi.c    -o lu_mpi    -lm || { echo "ERROR al compilar lu_mpi.c";    exit 1; }
gcc   lu_serial.c -o lu_serial -lm || { echo "ERROR al compilar lu_serial.c"; exit 1; }

# --- Límite de procesos = núcleos físicos ---
CORES=$(lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l)
if ! [ "$CORES" -ge 1 ] 2>/dev/null; then CORES=$(nproc); fi
LIMITE=${MAXP:-$CORES}
echo "Núcleos físicos ~ $CORES   (límite de procesos = $LIMITE)"

# =========================================================================
# SUBCOMANDO: desglose (cómputo / comunicación / idle, cíclico vs bloques)
# =========================================================================
if [ "$SUB" = "desglose" ]; then
    N=${2:-1500}
    PS=${3:-"1 2 4 8 16 24 32 48"}
    OUT="desglose_N${N}.csv"
    echo "Desglose  N=$N  REP=$REP"
    echo "n,modo,p,tiempo,comp,comm,idle" > "$OUT"
    for modo in 0 1; do
        nombre=$([ "$modo" -eq 0 ] && echo ciclico || echo bloques)
        for p in $PS; do
            if [ "$p" -gt "$LIMITE" ]; then echo "  (salteo p=$p > $LIMITE)"; continue; fi
            # desglose=1 imprime:  tiempo res rel comp comm idle  -> promedio de REP corridas
            L=$(for r in $(seq 1 $REP); do mpirun -n $p ./lu_mpi $N $modo 1; done \
                | awk '{t+=$1;c+=$4;m+=$5;i+=$6;k++} END{printf "%.6f,%.6f,%.6f,%.6f",t/k,c/k,m/k,i/k}')
            echo "$N,$nombre,$p,$L" >> "$OUT"
            echo "  N=$N $nombre p=$p -> $L"
        done
    done
    echo "Listo. Resultados en $OUT   (figuras: python3 graficar.py $OUT)"
    exit 0
fi

# =========================================================================
# SUBCOMANDO: pc | cluster  (speedup y eficiencia)
# =========================================================================
PERFIL=$SUB
OUT=${2:-resultados.csv}
if [ "$PERFIL" = "cluster" ]; then
    # N usados en el informe.
    NS="4000 8000 12000"; PS="1 2 4 8 16 24 32 48"
else
    NS="100 200 400 600 800"; PS="1 2 3 4 5 6"
fi
echo "Perfil=$PERFIL  REP=$REP"
echo "N = $NS"
echo "p = $PS"
echo "n,p,mediana_s,speedup,eficiencia,residual_abs,residual_rel" > "$OUT"
for N in $NS; do
    # T1 = mediana del tiempo de la versión SERIAL PURA (lu_serial), no del MPI con
    # p=1: usar el MPI con un proceso inflaría el speedup. 
    T1=$(for r in $(seq 1 $REP); do ./lu_serial $N | sed -n 's/.*tiempo = \([0-9.]*\).*/\1/p'; done \
         | sort -g | awk '{t[NR]=$1} END{print (NR%2)?t[(NR+1)/2]:(t[NR/2]+t[NR/2+1])/2}')
    echo "  T1 (serial puro) n=$N -> $T1 s"
    for p in $PS; do
        if [ "$p" -gt "$LIMITE" ]; then echo "  (salteo p=$p: supera el límite de $LIMITE)"; continue; fi
        echo "  [corriendo] n=$N p=$p  ($REP repeticiones)..."
        L=$(for r in $(seq 1 $REP); do mpirun -n $p ./lu_mpi $N; done | sort -g \
            | awk -v p=$p -v t1="$T1" '{t[NR]=$1; ra=$2; rr=$3} END{
                  m=(NR%2)?t[(NR+1)/2]:(t[NR/2]+t[NR/2+1])/2;
                  s=t1/m;
                  printf "%.9f,%.3f,%.1f,%s,%s", m, s, 100*s/p, ra, rr }')
        echo "$N,$p,$L" >> "$OUT"
        echo "  n=$N p=$p -> $L"
    done
done
echo "Listo. Resultados en $OUT   (figuras: python3 graficar.py $OUT)"
