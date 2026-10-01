#ifndef _PVECTOR_
#define _PVECTOR_

#include <initializer_list>
#include <iostream>
#include <cmath>
#include <string>
#include <type_traits>     // std::conditional, std::is_integral
#include <cstdlib>         // std::exit
#include "randnumgen.hpp"

// =====================
//  WRAPPER wrapt
// =====================
template <typename dstntype>
class wrapt
{
    dstntype wt;

public:
    template <typename srcntype>
    wrapt(const srcntype t)
    {
        wt = dstntype(t);
    }

    operator dstntype() const
    {
        return wt;
    }
};

// =====================
//  CLASS pvector
// =====================
template <typename ntype, int NT>
class pvector
{
    ntype v[NT];   // membro privato

public:
    // ----------------------
    //  COSTRUTTORI
    // ----------------------
    pvector()
    {
        for (int i = 0; i < NT; i++) v[i] = 0;
    }

    pvector(std::initializer_list<wrapt<ntype>> list)
    {
        int c = 0;
        for (ntype el : list)
        {
            if (c < NT) v[c] = el;
            c++;
        }
        for (; c < NT; c++) v[c] = ntype(0);
    }

    ~pvector() {}

    // ----------------------
    //  METODO show()
    // ----------------------
    void show(std::string s = "") const
    {
        std::cout << s << "(";
        for (int i = 0; i < NT; i++)
        {
            std::cout << v[i];
            if (i < NT - 1) std::cout << ",";
        }
        std::cout << ")\n";
    }

    // ----------------------
    //  ASSEGNAMENTO =
    // ----------------------
    pvector& operator=(const pvector& v2)
    {
        for (int i = 0; i < NT; i++) v[i] = v2.v[i];
        return *this;
    }

    // ----------------------
    //  OPERATORI +  -
    // ----------------------
    pvector operator+(const pvector& v2) const
    {
        pvector vs;
        for (int i = 0; i < NT; i++) vs.v[i] = v[i] + v2.v[i];
        return vs;
    }

    pvector operator-(const pvector& v2) const
    {
        pvector vs;
        for (int i = 0; i < NT; i++) vs.v[i] = v[i] - v2.v[i];
        return vs;
    }

    pvector sum(const pvector& v2) const
    {
        return (*this) + v2;
    }

    // ----------------------
    //  GET / SET
    // ----------------------
    ntype get(int i) const { return v[i]; }
    ntype set(int i, ntype val) { return v[i] = val; }

    // ----------------------
    //  += e -=
    // ----------------------
    pvector& operator+=(const pvector& v2)
    {
        for (int i = 0; i < NT; i++) v[i] += v2.v[i];
        return *this;
    }

    pvector& operator-=(const pvector& v2)
    {
        for (int i = 0; i < NT; i++) v[i] -= v2.v[i];
        return *this;
    }

    // ----------------------
    //  ACCESSO CON ()
    // ----------------------
    ntype operator()(int idx) const
    {
        return v[idx];
    }

    ntype& operator()(int idx)
    {
        return v[idx];
    }

    // ----------------------
    //  ACCESSO CON []
    // ----------------------
    ntype& operator[](int idx)
    {
        return v[idx];
    }

    const ntype& operator[](int idx) const
    {
        return v[idx];
    }

    // ----------------------
    //  DOT PRODUCT
    // ----------------------
    ntype operator*(const pvector& vec) const
    {
        ntype sp = 0;
        for (int i = 0; i < NT; i++) sp += v[i] * vec.v[i];
        return sp;
    }

    // ----------------------
    //  NORM
    // ----------------------
    typename std::conditional<std::is_integral<ntype>::value, double, ntype>::type
    norm(void) const
    {
        return std::sqrt((*this) * (*this));
    }

    // ----------------------
    //  VETTORE * SCALARE
    // ----------------------
    pvector operator*(ntype s) const
    {
        pvector vt;
        for (int i = 0; i < NT; i++) vt.v[i] = v[i] * s;
        return vt;
    }

    pvector operator/(ntype s) const
    {
        pvector vt;
        for (int i = 0; i < NT; i++) vt.v[i] = v[i] / s;
        return vt;
    }

    pvector& operator*=(ntype s)
    {
        for (int i = 0; i < NT; i++) v[i] *= s;
        return *this;
    }

    pvector& operator/=(ntype s)
    {
        for (int i = 0; i < NT; i++) v[i] /= s;
        return *this;
    }

    // ----------------------
    //  == e !=
    // ----------------------
    bool operator==(const pvector& vec) const
    {
        for (int i = 0; i < NT; i++)
            if (v[i] != vec.v[i]) return false;
        return true;
    }

    bool operator!=(const pvector& vec) const
    {
        return !((*this) == vec);
    }

    // ----------------------
    //  CROSS PRODUCT (solo NT=3)
    // ----------------------
    pvector operator^(const pvector& vec) const
    {
        if (NT != 3)
        {
            std::cout << "Cross product not defined\n";
            std::exit(1);
        }

        pvector vt;
        vt.v[0] = v[1] * vec.v[2] - v[2] * vec.v[1];
        vt.v[1] = v[2] * vec.v[0] - v[0] * vec.v[2];
        vt.v[2] = v[0] * vec.v[1] - v[1] * vec.v[0];
        return vt;
    }

    // ----------------------
    //  FRIEND scalar * vector
    // ----------------------
    friend pvector operator*(ntype s, const pvector& vec)
    {
        pvector vt;
        for (int i = 0; i < NT; i++) vt.v[i] = s * vec.v[i];
        return vt;
    }

    // ----------------------
    //  RINT() (metodo membro)  <-- reso const
    // ----------------------
    pvector rint() const
    {
        pvector vt;
        for (int i = 0; i < NT; i++) vt.v[i] = std::rint(v[i]);
        return vt;
    }

    // ----------------------
    //  <<
    // ----------------------
    friend std::ostream& operator<<(std::ostream& os, const pvector& vec)
    {
        os << "(";
        for (int i = 0; i < NT; i++)
        {
            os << vec.v[i];
            if (i < NT - 1) os << ",";
        }
        os << ")";
        return os;
    }

    // ----------------------
    //  COMPONENT-WISE OPERATIONS  <-- rese const
    // ----------------------
    pvector mulcw(const pvector& vec) const
    {
        pvector vt;
        for (int i = 0; i < NT; i++) vt(i) = v[i] * vec(i);
        return vt;
    }

    pvector divcw(const pvector& vec) const
    {
        pvector vt;
        for (int i = 0; i < NT; i++) vt(i) = v[i] / vec(i);
        return vt;
    }

    // ----------------------
    //  RANDOM
    // ----------------------
    pvector& random(const ntype& L)
    {
        for (int i = 0; i < NT; i++)
            v[i] = (rng.ranf() - 0.5) * L;
        return *this;
    }

    void random_orient(void)
    {
        ntype rS, S, V1, V2;

        if (NT != 3)
        {
            std::cout << "[random_orient] Only 3D vectors are supported\n";
            return;
        }

        do {
            V1 = 2.0 * rng.ranf() - 1.0;
            V2 = 2.0 * rng.ranf() - 1.0;
            S = V1 * V1 + V2 * V2;
        } while (S >= 1.0);

        rS = std::sqrt(1.0 - S);
        (*this) = { 2.0 * rS * V1, 2.0 * rS * V2, 1.0 - 2.0 * S };
    }
};

// ===============================================
//  TEMPLATE LIBERO: rint(vec)
// ===============================================
template<typename ntype, int NT>
pvector<ntype, NT> rint(const pvector<ntype, NT>& vec)
{
    pvector<ntype, NT> vt;
    for (int i = 0; i < NT; i++) vt(i) = std::rint(vec(i));
    return vt;
}

// alias utile
using pvec3d = pvector<double, 3>;

#endif
