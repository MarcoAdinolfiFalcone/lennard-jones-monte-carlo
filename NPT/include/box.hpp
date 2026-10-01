// ======================================================
// box.hpp — Box cubico con PBC + Minimum Image Convention
// (compatibile col tuo codice + stile prof)
// ======================================================

#ifndef BOX_HPP
#define BOX_HPP

#include "pvector.hpp"
#include <cmath>   // floor

struct Box
{
    double L;  // lato della scatola cubica

    Box() : L(0.0) {}
    explicit Box(double L_) : L(L_) {}

    double volume() const { return L * L * L; }

    double density(int N) const
    {
        return static_cast<double>(N) / volume();
    }

    // ---- stile prof: box dimensions come vettore (Lx,Ly,Lz) ----
    pvec3d Lvec() const
    {
        return pvec3d{L, L, L};
    }

    // --------------------------------------------------
    // Applica PBC: riporta r in [0, L)
    // Versione robusta: usa floor invece dei while
    // --------------------------------------------------
    void apply_pbc(pvec3d& r) const
    {
        for (int k = 0; k < 3; ++k) {
            // r = r - L * floor(r/L)
            // garantisce r in [0,L) anche se r è fuori di molte box-length
            r[k] -= L * std::floor(r[k] / L);

            // per sicurezza numerica (caso limite r==L)
            if (r[k] >= L) r[k] -= L;
            if (r[k] < 0.0) r[k] += L;
        }
    }

    // --------------------------------------------------
    // Minimum Image Convention: dr = ri - rj corretto
    //   dr = dr - L * rint(dr/L)
    // dove rint è arrotondamento al più vicino intero.
    // --------------------------------------------------
    pvec3d delta_pbc(const pvec3d& ri, const pvec3d& rj) const
    {
        pvec3d dr = ri - rj;
        pvec3d Lv = Lvec();

        // richiede che pvector supporti divcw/mulcw e rint(pvec)
        // (nel tuo progetto li hai già in pvector.hpp)
        dr = dr - Lv.mulcw(rint(dr.divcw(Lv)));

        return dr;
    }

    // Scala la scatola (utile se vuoi, anche se nel tuo NPT setti L direttamente)
    void scale(double factor)
    {
        L *= factor;
    }
};

#endif // BOX_HPP
