#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt


def load_sus_eos(path):
    """
    EOS from SUS (muVT) via reweighting
    expected columns (whitespace separated):
      mu   rho   P   <N>   lnXi
    with possible comment lines starting with '#'
    """
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    mu = data[:, 0]
    rho = data[:, 1]
    P = data[:, 2]
    return mu, rho, P


def load_npt_eos(path):
    """
    EOS from NPT (point 1)
    expected columns (whitespace separated):
      P   rho   U_over_N   V   L
    with possible comment lines starting with '#'
    """
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    P = data[:, 0]
    rho = data[:, 1]
    return rho, P


def main():
    # ====== EDIT THESE PATHS ======
    sus_file = "nmax=130/eos_sus_muVT_nmax130.dat"   # muVT+SUS output (mu rho P <N> lnXi)
    npt_file = "nmax=130/eos_NPT.dat"   # NPT output (P rho U/N V L)

    # Reduced temperature (for ideal gas reference line)
    Tstar = 1.4

    # ====== LOAD DATA ======
    mu_sus, rho_sus, P_sus = load_sus_eos(sus_file)
    rho_npt, P_npt = load_npt_eos(npt_file)

    # ====== PLOT ======
    plt.figure()
    plt.scatter(rho_sus, P_sus, marker="o", s=18, label="µVT + SUS (reweighting)")
    plt.scatter(rho_npt, P_npt, marker="s", s=35, label="NPT (discrete runs)")

    # Ideal gas reference
    rho_min = min(np.min(rho_sus), np.min(rho_npt))
    rho_max = max(np.max(rho_sus), np.max(rho_npt))
    rho_line = np.linspace(max(0.0, rho_min), rho_max, 300)
    P_ideal = Tstar * rho_line
    plt.plot(rho_line, P_ideal, label=f"Ideal gas: P*=ρ*T* (T*={Tstar})")

    plt.xlabel(r"$\rho^*$")
    plt.ylabel(r"$P^*$")
    plt.title("Equation of State comparison (LJ reduced units)")
    plt.grid(True)
    plt.legend()

    # Save + show
    out_png = "EOS_compare.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print(f"Saved: {out_png}")
    plt.show()


if __name__ == "__main__":
    main()
