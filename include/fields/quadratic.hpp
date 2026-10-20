#pragma once

#include "extension.hpp"

#include <cmath>
#include <numbers>

namespace fields
{

template <typename E>
concept QuadraticExtension = Field<E> && requires {
    typename E::Base;
    { E::degree } -> std::convertible_to<std::size_t>;
} && (E::degree == 2);

// R[X]/(X^2 + 1) (isk 2.7, 2.12)
template <typename E>
concept ComplexField =
    QuadraticExtension<E> && std::same_as<typename E::Base, double> && (E::modulus == std::array<double, 3>{1, 0, 1});

// Conjugation of a quadratic extension: the other root of f = X^2 + a_1 X + a_0 is
// x' = -a_1 - x, so a + bx -> (a - b a_1) - bx.
// For C this is a + bi -> a - bi, for F9 a + bι -> a - bι, for F4 it is the frobenius.
template <QuadraticExtension E> constexpr E conjugate(const E& z)
{
    using K = typename E::Base;
    const K a = z.coeffs[0];
    const K b = z.coeffs[1];
    return E::make(a - b * E::modulus[1], -b);
}

// z * conj(z), lies in K; 1/z = conj(z) / norm(z) (isk 2.7, 2.10)
template <QuadraticExtension E> constexpr typename E::Base norm(const E& z)
{
    return (z * conjugate(z)).coeffs[0];
}

// z + conj(z), lies in K
template <QuadraticExtension E> constexpr typename E::Base trace(const E& z)
{
    return (z + conjugate(z)).coeffs[0];
}

template <ComplexField E> constexpr double re(const E& z)
{
    return z.coeffs[0];
}

template <ComplexField E> constexpr double im(const E& z)
{
    return z.coeffs[1];
}

// |a + bi| = sqrt(a^2 + b^2)
template <ComplexField E> double abs(const E& z)
{
    return std::hypot(z.coeffs[0], z.coeffs[1]);
}

// the angle phi in [0, 2 pi) with z = |z| (cos(phi) + i sin(phi)) (isk 2.12)
template <ComplexField E> double arg(const E& z)
{
    if (z.is_zero())
    {
        throw std::domain_error("0 has no argument");
    }
    double phi = std::atan2(z.coeffs[1], z.coeffs[0]);
    return phi < 0 ? phi + 2 * std::numbers::pi : phi;
}

// r (cos(phi) + i sin(phi))
template <ComplexField E> E from_polar(double r, double phi)
{
    return E::make(r * std::cos(phi), r * std::sin(phi));
}

} // namespace fields
