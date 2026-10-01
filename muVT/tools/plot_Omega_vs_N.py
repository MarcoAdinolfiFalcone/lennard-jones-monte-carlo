#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_two_col(path):
    x, y = [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            x.append(int(float(parts[0])))
            y.append(float(parts[1]))
    return np.array(x, dtype=int), np.array(y, dtype=float)

def read_bad_windows(path):
    bad = set()
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            try:
                bad.add(int(float(parts[0])))
            except:
                pass
    return bad

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_Omega_vs_N.py <run_dir> [Tstar]")
        raise SystemExit(1)

    run_dir = sys.argv[1]
    Tstar = float(sys.argv[2]) if len(sys.argv) >= 3 else 1.4

    # usa il file con gauge esplicito se presente
    cand = ["logW_muRef.dat", "logW_from_windows.dat"]
    for name in cand:
        p = os.path.join(run_dir, name)
        if os.path.isfile(p):
            in_path = p
            src_name = name
            break
    else:
        raise FileNotFoundError(f"Missing {cand} in {run_dir}")

    N, logW = read_two_col(in_path)

    # ordina
    idx = np.argsort(N)
    N = N[idx]
    logW = logW[idx]

    # DEFINIZIONE: Omega*(N) = -T* ln W(N)
    Omega_star = -Tstar * logW

    # bad windows (se presenti)
    bad = set()
    for name in ["bad_windows.dat", "bad_windows_from_windows.dat"]:
        p = os.path.join(run_dir, name)
        if os.path.isfile(p):
            bad = read_bad_windows(p)
            break

    plt.figure()
    plt.plot(N, Omega_star, marker="o", linestyle="None", markersize=3)

    if bad:
        mask = np.array([n in bad for n in N], dtype=bool)
        if np.any(mask):
            plt.plot(
                N[mask], Omega_star[mask],
                marker="x", linestyle="None",
                markersize=6, mew=1.5,
                label="poorly sampled window"
            )
            plt.legend(loc="best")

    plt.xlabel(r"$N$")
    plt.ylabel(r"$\Omega^*(N) = -T^* \ln W(N;\mu_{\mathrm{ref}}^*)$")
    plt.title(rf"SUS reconstruction at $T^*={Tstar}$, gauge fixed by $\ln W(0)=0$")
    plt.grid(True, alpha=0.3)

    out_png = os.path.join(run_dir, "OmegaStar_vs_N_muRef.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print("Saved:", out_png)

if __name__ == "__main__":
    main()
