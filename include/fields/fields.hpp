#pragma once

// Fields after "Mathematik für Informatiker" (M. Künzer), chapter 2

#include "extension.hpp"
#include "field.hpp"
#include "galois.hpp"
#include "irreducible.hpp"
#include "polynomial.hpp"
#include "prime_field.hpp"
#include "quadratic.hpp"
#include "rational.hpp"
#include "real.hpp"

namespace fields
{

using Q = Rational;
using R = Real;

// isk 2.7: C = R[X]/(X^2 + 1)R[X], i^2 = -1
using Complex = Extension<Real, poly<Real>(1, 0, 1), "i">;

using F2 = Fp<2>;
using F3 = Fp<3>;
using F5 = Fp<5>;
using F7 = Fp<7>;

// isk 2.8: F4 = F2[X]/(X^2 + X + 1), α^2 = α + 1
using F4 = Extension<F2, poly<F2>(1, 1, 1), "α">;

// isk 2.9: F8 = F2[X]/(X^3 + X + 1), β^3 = β + 1
using F8 = Extension<F2, poly<F2>(1, 1, 0, 1), "β">;

// isk 2.10: F9 = F3[X]/(X^2 + 1), ι^2 = -1
using F9 = Extension<F3, poly<F3>(1, 0, 1), "ι">;

} // namespace fields
