// ======================================================
// particle_npt.hpp  — particella minimale per MC NPT
// ======================================================

#ifndef PARTICLE_NPT_HPP
#define PARTICLE_NPT_HPP

#include "pvector.hpp"

struct Particle
{
    pvec3d r;      // posizione corrente
    pvec3d rold;   // posizione salvata (per reject)

    Particle() : r{0.0, 0.0, 0.0}, rold{0.0, 0.0, 0.0} {}

    // salva lo stato (prima di una mossa proposta)
    void store()
    {
        rold = r;
    }

    // ripristina lo stato (se rifiuti la mossa)
    void restore()
    {
        r = rold;
    }
    void tra_move(const pvec3d& dr)
    {
        r += dr;
    }
};

#endif // PARTICLE_NPT_HPP
