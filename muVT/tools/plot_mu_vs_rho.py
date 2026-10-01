#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_eos(path):
    mu, rho = [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            # columns: mu  rho  P  <N>  lnXi
            mu.append(float(parts[0]))
            rho.append(float(parts[1]))
    return np.array(mu, dtype=float), np.array(rho, dtype=float)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_mu_vs_rho.py <run_dir>")
        raise SystemExit(1)

    run_dir = sys.argv[1]
    eos_path = os.path.join(run_dir, "EOS_mu_rho_P.dat")
    if not os.path.isfile(eos_path):
        raise FileNotFoundError(f"Missing: {eos_path}")

    mu, rho = read_eos(eos_path)

    # ordina per densità
    idx = np.argsort(rho)
    rho = rho[idx]
    mu  = mu[idx]

    plt.figure()
    plt.scatter(rho, mu, s=14)   # <<< SOLO PUNTI

    plt.xlabel(r"$\rho^*=\langle N\rangle/V^*$")
    plt.ylabel(r"$\mu^*$")
    plt.title(r"$\mu^*(\rho^*)$ from $\mu VT$+SUS reweighting ($T^*=1.4$, $V^*=343$)")
    plt.grid(True, alpha=0.3)

    out_png = os.path.join(run_dir, "mu_vs_rho.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print("Saved:", out_png)

if __name__ == "__main__":
    main()
