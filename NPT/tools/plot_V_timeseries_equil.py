import os
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# PATHS
# ==========================================
BASE_DIR = "../data"
EOS_FILE = os.path.join(BASE_DIR, "EOS_NPT_T1p4.dat")

OUT_DIR  = os.path.join(BASE_DIR, "figures")
OUT_PNG  = os.path.join(OUT_DIR, "EOS_NPT_T1p4_P_vs_rho.png")

# ==========================================
# PARAMETERS
# ==========================================
TSTAR = 1.4   # reduced temperature

# ==========================================
# MAIN
# ==========================================
def main():

    if not os.path.exists(EOS_FILE):
        raise FileNotFoundError(f"Missing file: {EOS_FILE}")

    os.makedirs(OUT_DIR, exist_ok=True)

    # Columns in EOS file (from your C++ header):
    # 0: P
    # 1: <rho>
    data = np.loadtxt(EOS_FILE, comments="#")

    P   = data[:, 0]
    rho = data[:, 1]

    # Sort by density (ONLY for nicer visualization)
    idx = np.argsort(rho)
    rho = rho[idx]
    P   = P[idx]

    # Ideal gas reference
    rho_line = np.linspace(rho.min(), rho.max(), 300)
    P_ideal  = TSTAR * rho_line

    # ======================================
    # PLOT
    # ======================================
    plt.figure(figsize=(6,4))

    # NPT EOS: discrete points only
    plt.scatter(rho, P, s=60, zorder=3, label="NPT MC")

    # Ideal gas
    plt.plot(rho_line, P_ideal, linestyle="--", zorder=2, label="Ideal gas")

    plt.xlabel(r"$\rho^*$")
    plt.ylabel(r"$P^*$")
    plt.title(rf"NPT equation of state at $T^*={TSTAR}$")

    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=250)
    plt.close()

    print("EOS NPT plot saved in:")
    print(" ", OUT_PNG)

# ==========================================
if __name__ == "__main__":
    main()
