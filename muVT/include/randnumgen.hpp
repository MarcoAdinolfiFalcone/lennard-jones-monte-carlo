// ======================================================
//  randnumgen.hpp  — RNG globale (C++17)
// ======================================================
#ifndef RANDNUMGEN_HPP
#define RANDNUMGEN_HPP

#include <random>

class Random {
private:
    std::mt19937 engine;
    std::uniform_real_distribution<double> dist;

public:
    Random() : engine(std::random_device{}()), dist(0.0, 1.0) {}

    void seed(unsigned int s) { engine.seed(s); }
    void rseed() { engine.seed(std::random_device{}()); }

    double ranf() { return dist(engine); } // in [0,1)
};

// C++17: variabile inline header-only
inline Random rng;

#endif // RANDNUMGEN_HPP
