#!/usr/bin/env python3
# Genera las figuras de un CSV de medir.sh (detecta el tipo por sus columnas).
# Uso:  python3 graficar.py <archivo.csv>

import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

archivo = sys.argv[1] if len(sys.argv) > 1 else "resultados.csv"
SALIDA = os.path.dirname(os.path.abspath(archivo))
STEM = os.path.splitext(os.path.basename(archivo))[0]   # nombre del CSV sin extensión
def ruta(nombre):
    # Prefijo = nombre del CSV, así distintas corridas no se pisan entre sí
    return os.path.join(SALIDA, f"{STEM}_{nombre}")

with open(archivo) as f:
    filas = list(csv.DictReader(f))
cols = list(filas[0].keys()) if filas else []

# =========================================================================
# CASO 1: desglose (barras apiladas cómputo / comunicación / idle)
# =========================================================================
if "idle" in cols:
    datos = {"ciclico": {}, "bloques": {}}
    N = None
    for r in filas:
        N = int(float(r["n"]))
        datos[r["modo"]][int(float(r["p"]))] = (float(r["comp"]), float(r["comm"]), float(r["idle"]))
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    for col, modo in enumerate(["ciclico", "bloques"]):
        d = datos[modo]; ps = sorted(d)
        comp = [d[p][0] for p in ps]; comm = [d[p][1] for p in ps]; idle = [d[p][2] for p in ps]
        xs = [str(p) for p in ps]; base = [c + m for c, m in zip(comp, comm)]
        ax[col].bar(xs, comp, label="cómputo", color="#4C78A8")
        ax[col].bar(xs, comm, bottom=comp, label="comunicación", color="#F58518")
        ax[col].bar(xs, idle, bottom=base, label="idle (desbalance)", color="#E45756")
        ax[col].set_title(f"Reparto {modo}"); ax[col].set_xlabel("procesos (p)"); ax[col].grid(axis="y", alpha=0.3)
    ax[0].set_ylabel("tiempo (s)"); ax[0].legend()
    fig.suptitle(f"Composición del tiempo de factorización (N={N})", y=1.02)
    fig.savefig(ruta("apilado.png"), dpi=150, bbox_inches="tight")
    print("Figura generada:", ruta("apilado.png"))
    sys.exit(0)

# =========================================================================
# CASO 2: speedup / eficiencia / tiempo / residual
# =========================================================================
datos = [{
    "n": int(float(r["n"])), "p": int(float(r["p"])), "t": float(r["mediana_s"]),
    "sp": float(r["speedup"]), "ef": float(r["eficiencia"]), "res": float(r["residual_abs"]),
} for r in filas]
Ns = sorted({d["n"] for d in datos}); Ps = sorted({d["p"] for d in datos})

def por_n(n, campo):
    fs = sorted((d for d in datos if d["n"] == n), key=lambda d: d["p"]); return [d["p"] for d in fs], [d[campo] for d in fs]
def por_p(p, campo):
    fs = sorted((d for d in datos if d["p"] == p), key=lambda d: d["n"]); return [d["n"] for d in fs], [d[campo] for d in fs]

# Speedup
plt.figure()
for n in Ns:
    xs, ys = por_n(n, "sp"); plt.plot(xs, ys, marker="o", label=f"N = {n}")
plt.plot(Ps, Ps, "k--", linewidth=1, label="ideal (S = p)")
plt.xlabel("Cantidad de procesos (p)"); plt.ylabel("Speedup  $S_p = T_1 / T_p$")
plt.title("Speedup de la factorización LU"); plt.legend(); plt.grid(True, alpha=0.3)
plt.savefig(ruta("speedup_vs_p.png"), dpi=150, bbox_inches="tight")

# Eficiencia
plt.figure()
for n in Ns:
    xs, ys = por_n(n, "ef"); plt.plot(xs, ys, marker="s", label=f"N = {n}")
plt.axhline(100, color="k", linestyle="--", linewidth=1, label="ideal (100 %)")
plt.xlabel("Cantidad de procesos (p)"); plt.ylabel("Eficiencia  $E_p$ (%)")
plt.title("Eficiencia de la factorización LU"); plt.legend(); plt.grid(True, alpha=0.3)
plt.savefig(ruta("eficiencia_vs_p.png"), dpi=150, bbox_inches="tight")

# Tiempo vs N (log-log)
plt.figure()
for p in Ps:
    xs, ys = por_p(p, "t"); plt.plot(xs, ys, marker="o", label=f"p = {p}")
plt.xscale("log"); plt.yscale("log")
plt.xlabel("Dimensión N"); plt.ylabel("Tiempo de factorización (s)")
plt.title("Escalabilidad del tiempo de cómputo"); plt.legend(); plt.grid(True, which="both", alpha=0.3)
plt.savefig(ruta("tiempo_vs_N.png"), dpi=150, bbox_inches="tight")

# Residual vs N (semilog)
p0 = Ps[0]; xs, ys = por_p(p0, "res")
plt.figure()
plt.plot(xs, ys, marker="d", color="C3"); plt.yscale("log")
plt.xlabel("Dimensión N"); plt.ylabel("Error residual  $\\|Ax - b\\|_2$")
plt.title("Precisión numérica del residual"); plt.grid(True, which="both", alpha=0.3)
plt.savefig(ruta("residual_vs_N.png"), dpi=150, bbox_inches="tight")

print("Figuras generadas en:", SALIDA, "con prefijo", STEM + "_")
print(f"  - {STEM}_speedup_vs_p.png   - {STEM}_eficiencia_vs_p.png")
print(f"  - {STEM}_tiempo_vs_N.png    - {STEM}_residual_vs_N.png")
