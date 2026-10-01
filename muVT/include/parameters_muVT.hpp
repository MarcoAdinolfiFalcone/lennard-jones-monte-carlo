#ifndef PARAMETERS_MUVT_HPP
#define PARAMETERS_MUVT_HPP

#include <cmath>
#include <algorithm>
#include "pvector.hpp"

// ======================================================
//  MuVTParams — parametri per LJ in ensemble µV
// ======================================================
struct MuVTParams
{
    // ---------- Thermo ----------
    double T  = 1.4;
    double mu = -3.5;

    // ---------- LJ ----------
    double epsilon = 1.0;
    double sigma   = 1.0;
    double rc      = 2.5;

    // ---------- Box / init ----------
    // Se V_box > 0, impone il volume. Altrimenti usa rho_init.
    double V_box   = 0.0;
    double rho_init = 0.3;
    int    N_init  = 50;

    // L (cubo)
    pvec3d L = {0.0, 0.0, 0.0};

    // N corrente (gestito dal system)
    int Np = 0;

    // ---------- MC ----------
    long totsteps    = 400000;
    long eqstps      = 100000;
    long savemeasure = 50;
    long outstps     = 10000;

    double deltra = 0.2;

    // probabilità assolute (somma=1)
    double p_tra = 0.50;
    double p_ins = 0.25;
    double p_del = 0.25;

    // ---------- RNG ----------
    int seed = 12345; // se <0 -> random seed (gestito nel system)

    // ---------- activity ----------
    double z = 0.0;   // z = exp(beta*mu) (assumo lambda^3=1)

    // ---------------------------
    // helper
    // ---------------------------
    double beta() const { return 1.0 / T; }

    double volume() const { return L(0) * L(1) * L(2); }

    double density() const
    {
        const double V = volume();
        return (V > 0.0) ? (double(Np) / V) : 0.0;
    }

    void normalize_move_probs()
    {
        const double s = p_tra + p_ins + p_del;
        if (s <= 0.0) {
            p_tra = 1.0; p_ins = 0.0; p_del = 0.0;
            return;
        }
        p_tra /= s;
        p_ins /= s;
        p_del /= s;
    }

    void set_box_from_init()
    {
        double V = 0.0;

        if (V_box > 0.0) {
            V = V_box;
        } else {
            // V = N/rho
            V = (rho_init > 0.0) ? (double(N_init) / rho_init) : 1.0;
        }

        const double Lc = std::cbrt(V);
        L = {Lc, Lc, Lc};
    }

    void update_activity()
    {
        z = std::exp(beta() * mu);
    }
};

#endif // PARAMETERS_MUVT_HPP
