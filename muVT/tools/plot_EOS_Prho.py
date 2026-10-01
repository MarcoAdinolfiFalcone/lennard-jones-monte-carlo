#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_eos(path):
    """
    Read EOS_mu_rho_P.dat
    Expected columns (whitespace):
      mu   rho(<N>/V)   P   <N>   lnXi
    Returns: rho, P
    """
    rho, P = [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            rho.append(float(parts[1]))
            P.append(float(parts[2]))
    return np.array(rho, dtype=float), np.array(P, dtype=float)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_EOS_Prho.py <run_dir> [Tstar]")
        raise SystemExit(1)

    run_dir = sys.argv[1]
    Tstar = float(sys.argv[2]) if len(sys.argv) >= 3 else 1.4

    eos_path = os.path.join(run_dir, "EOS_mu_rho_P.dat")
    if not os.path.isfile(eos_path):
        raise FileNotFoundError(f"Missing: {eos_path}")

    rho, P = read_eos(eos_path)

    # sort by rho for a clean curve
    idx = np.argsort(rho)
    rho = rho[idx]
    P   = P[idx]

    # Ideal gas line: P* = rho* T*
    P_id = Tstar * rho

    plt.figure()

    # Data: ONLY POINTS (no connecting line)
    plt.scatter(rho, P, s=14, label=r"$\mu VT$ + SUS")

    # Ideal gas reference
    plt.plot(rho, P_id, linestyle="--", linewidth=1.5,
             label=rf"Ideal gas: $P^*=T^*\rho^*$ (here $T^*={Tstar}$)")

    plt.xlabel(r"$\rho^*=\langle N\rangle/V^*$")
    plt.ylabel(r"$P^*$")
    plt.title(rf"Equation of state from $\mu VT$ + SUS ($T^*={Tstar}$, $V^*=343$)")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best")

    out_png = os.path.join(run_dir, "EOS_SUS_Prho_with_ideal.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print("Saved:", out_png)

if __name__ == "__main__":
    main()
