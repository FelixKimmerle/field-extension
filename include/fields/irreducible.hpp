#pragma once

#include "polynomial.hpp"
#include "rational.hpp"
#include "real.hpp"

namespace fields
{

enum class Irreducibility
{
    Irreducible,
    Reducible,
    // no algorithm available for this field (e.g. Q with degree >= 4)
    Unknown,
};

namespace detail
{

// Ben-Or: f is irreducible over F_q iff gcd(f, X^(q^i) - X) = 1 for 1 <= i <= deg(f) / 2,
// since X^(q^i) - X is the product of all monic irreducibles with degree dividing i
template <FiniteField K> constexpr Irreducibility ben_or(const Polynomial<K>& f)
{
    constexpr uint64_t q = FieldTraits<K>::cardinality;
    const Polynomial<K> x = Polynomial<K>::x();
    Polynomial<K> h = x % f;
    for (int64_t i = 1; i <= f.degree() / 2; i++)
    {
        // h = X^(q^i) mod f
        h = Polynomial<K>::pow_mod(h, q, f);
        if (Polynomial<K>::gcd(f, h - x).degree() != 0)
        {
            return Irreducibility::Reducible;
        }
    }
    return Irreducibility::Irreducible;
}

constexpr std::vector<int64_t> divisors(int64_t n)
{
    n = n < 0 ? -n : n;
    std::vector<int64_t> result;
    for (int64_t d = 1; d <= n / d; d++)
    {
        if (n % d == 0)
        {
            result.push_back(d);
            if (d != n / d)
            {
                result.push_back(n / d);
            }
        }
    }
    return result;
}

constexpr int64_t gcd(int64_t a, int64_t b)
{
    while (b != 0)
    {
        int64_t t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

// rational root test (isk 2.12): after clearing denominators every root u/v
// in lowest terms has u | a_0 and v | a_n
constexpr bool has_rational_root(const Polynomial<Rational>& f)
{
    if (f.coefficient(0) == Rational(0))
    {
        return true;
    }
    int64_t lcm = 1;
    for (const Rational& c : f.coefficients())
    {
        lcm = lcm / gcd(lcm, c.den) * c.den;
    }
    const Rational a0 = f.coefficient(0) * Rational(lcm);
    const Rational an = f.leading_coefficient() * Rational(lcm);
    for (int64_t u : divisors(a0.num))
    {
        for (int64_t v : divisors(an.num))
        {
            if (f(Rational(u, v)) == Rational(0) || f(Rational(-u, v)) == Rational(0))
            {
                return true;
            }
        }
    }
    return false;
}

} // namespace detail

// Decides whether f is irreducible over K (isk 2.5), if there is a method for K.
// Irreducible means deg(f) >= 1 and f has no factorization into factors of degree >= 1.
template <Field K> constexpr Irreducibility check_irreducible(const Polynomial<K>& polynomial)
{
    const Polynomial<K> f = polynomial.monic();
    const int64_t n = f.degree();
    if (n < 1)
    {
        return Irreducibility::Reducible;
    }
    if (n == 1)
    {
        return Irreducibility::Irreducible;
    }

    if constexpr (FieldTraits<K>::algebraically_closed)
    {
        // every polynomial splits into linear factors (isk 2.12), e.g. over C
        return Irreducibility::Reducible;
    }
    else if constexpr (is_finite_field_v<K>)
    {
        return detail::ben_or(f);
    }
    else if constexpr (std::same_as<K, double>)
    {
        // over R only linear factors and quadratics without real roots are irreducible,
        // since C = R[X]/(X^2 + 1) is algebraically closed
        if (n > 2)
        {
            return Irreducibility::Reducible;
        }
        const double p = f.coefficient(1);
        const double q = f.coefficient(0);
        return p * p - 4 * q < 0 ? Irreducibility::Irreducible : Irreducibility::Reducible;
    }
    else if constexpr (std::same_as<K, Rational>)
    {
        // degree 2 and 3: irreducible iff there is no root (isk 2.5)
        if (n <= 3)
        {
            return detail::has_rational_root(f) ? Irreducibility::Reducible : Irreducibility::Irreducible;
        }
        return Irreducibility::Unknown;
    }
    else
    {
        return Irreducibility::Unknown;
    }
}

template <Field K> constexpr bool is_irreducible(const Polynomial<K>& f)
{
    return check_irreducible(f) == Irreducibility::Irreducible;
}

} // namespace fields
