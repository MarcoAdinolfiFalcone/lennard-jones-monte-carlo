#ifndef _RANDNUMGEN_HPP_
#define _RANDNUMGEN_HPP_

#include <random>

// RNG semplice con interfaccia ranf()
class Random {
private:
    std::mt19937 engine;
    std::uniform_real_distribution<double> dist;

public:
    Random() : engine(std::random_device{}()), dist(0.0, 1.0) {}

    void seed(unsigned int s) { engine.seed(s); }

    double ranf() { return dist(engine); }
};

// definizione globale (C++17)
inline Random rng;

#endif