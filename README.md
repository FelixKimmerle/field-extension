# fields

Header-only C++20 library for fields, based on chapter 2 of *Mathematik für Informatiker*
(M. Künzer). Its core is the construction from §2.6: for a field K and an irreducible
f ∈ K[X], the quotient K[X]/f(X)K[X] is a field. ℂ, F4, F8, F9, and every GF(pⁿ) come
out of this one construction.

The lecture page with the script is
[Mathematik für Informatiker (Künzer, Uni Stuttgart)](https://pnp.mathematik.uni-stuttgart.de/lexmath/kuenzer/info1920/).
That's where the mathematical foundations for everything here come from. Section numbers
like §2.6 in the code (`isk 2.6`) refer to this script.

```cpp
#include <fields/fields.hpp>
using namespace fields;

// isk 2.7 - 2.10 (predefined in fields.hpp)
using Complex = Extension<Real, poly<Real>(1, 0, 1), "i">;   // R[X]/(X^2 + 1)
using F8      = Extension<F2,   poly<F2>(1, 1, 0, 1), "β">;  // F2[X]/(X^3 + X + 1)

GF<9>   x;   // like FiniteField.create(9): F3[X]/(X^2 + 1)
GF<256> y;   // F2[X]/(X^8 + X^4 + X^3 + X + 1)

// towers: the base field can be an extension itself
using F16 = Extension<F4, poly<F4>(F4::generator(), 1, 1), "γ">;  // F4[X]/(X^2 + X + α)

const F8 b = F8::generator();
pow(b, 7) == F8(1);         // also works in static_assert
1 / (b + 1);                // β + β^2
print_multiplication_table<F4>(std::cout);
```

## Overview

| header | contents |
|---|---|
| `field.hpp` | `Field` / `FiniteField` concepts, `FieldTraits`, `pow`, `inverse` |
| `prime_field.hpp` | `Fp<p>` = ℤ/pℤ (§2.1, §2.3) |
| `real.hpp`, `rational.hpp` | ℝ (`double`), exact ℚ (`Rational`) |
| `polynomial.hpp` | `Polynomial<K>`: arithmetic, division with remainder, gcd, extended gcd, `pow_mod`, evaluation |
| `irreducible.hpp` | `check_irreducible`: Ben-Or for finite fields, degree/discriminant for ℝ, rational root test for ℚ up to degree 3 |
| `extension.hpp` | `Extension<K, f, name>` = K[X]/f(X) (§2.6) |
| `galois.hpp` | `GF<q>`, `all_elements`, `multiplicative_order`, `primitive_element`, `frobenius`, `irreducible_polynomials`, addition and multiplication tables |
| `quadratic.hpp` | `conjugate`, `norm`, `trace` for degree-2 extensions; `re`, `im`, `abs`, `arg`, `from_polar` for ℂ (§2.12) |

Notes:
- Coefficients are listed in ascending order (`poly<F2>(1, 1, 0, 1)` = 1 + X + X³), the same as in `Polynomial`.
- The modulus is a template argument. Where possible it is checked for irreducibility at
  compile time (all finite fields, ℝ, ℂ, ℚ up to degree 3), so a modulus like
  `Extension<F2, poly<F2>(1, 0, 1)>` fails to compile.
- `GF<q>` picks the first monic irreducible polynomial in lexicographic order. That gives
  the script's moduli for F4, F8, and F9, so `GF<8, "β">` is the same type as `F8`.
- Elements of F_p print as 0…p−1 (so −1 in F9 prints as `2`).

## Build

```sh
cmake -S . -B build
cmake --build build
./build/tests
./build/demo
```
