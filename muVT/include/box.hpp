// box.hpp
#ifndef BOX_HPP
#define BOX_HPP

#include "pvector.hpp"

struct Box {
    double L;  // lato della scatola cubica

    Box() : L(0.0) {}
    explicit Box(double L_) : L(L_) {}

    double volume() const {
        return L * L * L;
    }

    double density(int N) const {
        return static_cast<double>(N) / volume();
    }

    // Applica PBC riportando r in [0, L) in ogni direzione
    void apply_pbc(pvector<double,3>& r) const {
        for (int k = 0; k < 3; ++k) {
            while (r[k] >= L) r[k] -= L;
            while (r[k] <  0.0) r[k] += L;
        }
    }

    // Vettore distanza a immagine minima: ri - rj, corretto con PBC
    pvector<double,3> delta_pbc(const pvector<double,3>& ri,
                                const pvector<double,3>& rj) const
    {
        pvector<double,3> dr = ri - rj;
        for (int k = 0; k < 3; ++k) {
            if (dr[k] >  0.5 * L) dr[k] -= L;
            if (dr[k] < -0.5 * L) dr[k] += L;
        }
        return dr;
    }

    // Scala la scatola (serve per mosse di volume in NPT)
    void scale(double factor) {
        L *= factor;
    }
};

#endif // BOX_HPP
