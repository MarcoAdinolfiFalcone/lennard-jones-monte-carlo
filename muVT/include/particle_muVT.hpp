#ifndef PARTICLE_MUVT_HPP
#define PARTICLE_MUVT_HPP

#include "pvector.hpp"

// Particella usata nell'ensemble µVT
// Per ora contiene solo la posizione.
// In futuro si può estendere (es: velocità, tipo, ID, ecc.)

struct Particle {
    pvector<double,3> r;   // posizione nello spazio 3D
};

#endif // PARTICLE_MUVT_HPP
