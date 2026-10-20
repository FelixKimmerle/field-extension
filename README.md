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
| `eigen.hpp` | optional Eigen adapter: matrices over finite fields and ℚ, see below |

Notes:
- Coefficients are listed in ascending order (`poly<F2>(1, 1, 0, 1)` = 1 + X + X³), the same as in `Polynomial`.
- The modulus is a template argument. Where possible it is checked for irreducibility at
  compile time (all finite fields, ℝ, ℂ, ℚ up to degree 3), so a modulus like
  `Extension<F2, poly<F2>(1, 0, 1)>` fails to compile.
- `GF<q>` picks the first monic irreducible polynomial in lexicographic order. That gives
  the script's moduli for F4, F8, and F9, so `GF<8, "β">` is the same type as `F8`.
- Elements of F_p print as 0…p−1 (so −1 in F9 prints as `2`).

## Linear algebra with Eigen

`eigen.hpp` is an optional adapter that turns the exact fields into
[Eigen](https://eigen.tuxfamily.org) scalars. With it you can solve linear systems and
invert matrices over F_p, F9, GF(q) or ℚ, like the exercises in the lecture, without doing
Gauss elimination by hand. It needs Eigen 3.4 or newer (tested with 3.4 and 5.0.1); the
CMake build downloads Eigen for you (see [Build](#build)).

```cpp
#include <fields/eigen.hpp>
using namespace fields;

const F9 i = F9::generator();
Matrix<F9> a(3, 3);
a << 1,     i, 0,
     1 + i, 2, i,
     0,     1, 1 - i;
Vector<F9> b(3);
b << 1, i, 2;

a.determinant();                 // 2 + ι
auto inv = invert(a);            // std::optional<Matrix<F9>>
auto x = solve(a, b);            // std::optional<Vector<F9>>: (2, ι, 2ι)
print_matrix(std::cout, *inv);
```

`examples/linear_algebra.cpp` has more worked examples (F9, a singular matrix over F5,
the Hilbert matrix over ℚ, a 200×200 matrix over GF(256)).

### Supported scalar types

| works | does not work |
|---|---|
| `Fp<p>` (F2, F3, F5, …) | `Real` (`double`): use Eigen's normal `double` matrices |
| finite extensions: F4, F8, F9, `GF<q>`, towers like F16 over F4 | `Complex`: use `std::complex<double>` (gets SIMD and BLAS) |
| `Rational` (ℚ) and extensions of ℚ, e.g. ℚ(√2) | |

The concept `fields::ExactField<K>` says whether a type is supported. Using an
unsupported type, or mixing two fields in one expression (e.g. an F5 matrix times an F9
vector), is a compile error.

### Functions

| function | result |
|---|---|
| `Matrix<K>`, `Vector<K>` | aliases for dynamically sized `Eigen::Matrix<K, Dynamic, Dynamic>` / `Eigen::Matrix<K, Dynamic, 1>` |
| `invert(A)` | `std::optional<Matrix<K>>`: A⁻¹, or `nullopt` if A is singular or not square |
| `solve(A, b)` | `std::optional`: one x with A·x = b, or `nullopt` if there is none. `b` may have several columns. If the solution is not unique, all solutions are x + `kernel(A)` |
| `rank(A)` | rank as `Eigen::Index` |
| `kernel(A)` | basis of {x : A·x = 0} as columns; a matrix with 0 columns if the kernel is {0} |
| `print_matrix(os, A)` | prints A with aligned columns |

The rest of Eigen works directly on these matrices: `+ - *`, `transpose()`, blocks,
`determinant()`, `inverse()`, `partialPivLu()`, `fullPivLu()` (`image()`, `isInvertible()`,
…), for both fixed-size (`Eigen::Matrix<F7, 4, 4>`) and dynamic matrices. The helpers
above are thin wrappers around `fullPivLu()` that report "no solution" and "singular" as
`nullopt` instead of returning garbage.

Don't use anything that needs square roots, norms or a real order: QR (Householder),
SVD, LLT/LDLT, eigenvalue solvers, `norm()`, `normalized()`. Over a finite field these
make no sense, and the adapter's `sqrt` throws `std::domain_error`.

### How it works

- **No BLAS:** BLAS/MKL only exist for `float`, `double` and `std::complex`. For any other
  scalar type Eigen uses its own generic C++ kernels (plain loops, no SIMD), and those
  only need `+ - * /`.
- **What Eigen is told about the types:** a `Eigen::NumTraits<K>` specialization with
  `epsilon() = 0` (exact arithmetic, so pivots are compared against 0 instead of a
  tolerance), plus `abs`, `abs2`, `conj`, `real`, `imag` and `sqrt`, found via ADL.
- **Pivoting:** Eigen picks the pivot with the largest `abs()`. That only matters for
  floating point stability; over an exact field any nonzero pivot is fine. So `abs(x)`
  is 0 or 1 (for ℚ it's the real |x|), and `eigen.hpp` adds an arbitrary total order with 0
  as the smallest element, because Eigen compares pivots with `<`. This order exists only
  for Eigen. It is not a field ordering (finite fields have none), and it does not change
  how elements print.
- **Exact checks:** since there is no rounding, `solve` checks `A·x == b` exactly and
  `rank`/`kernel` are exact. With `double`, these would need a tolerance.

### Limitations

- **ℚ:** Gauss elimination makes fractions grow, and `Rational` stores 64-bit numerators
  and denominators. Large or badly conditioned systems can throw `std::overflow_error`.
- **Speed:** elements are computed coefficient by coefficient, without lookup tables. In
  a Release build, inverting a random 200×200 matrix over GF(256) takes about 1 s; a Debug
  build is about 30× slower.
- **Printing:** Eigen's `operator<<` works, but it counts bytes, so columns with names
  like `ι` are misaligned. Use `print_matrix`.

### Using it in your own CMake project

```cmake
add_subdirectory(field-extension)   # or FetchContent this repository
target_link_libraries(my_app PRIVATE fields Eigen3::Eigen)
```

`Eigen3::Eigen` exists when `FIELDS_WITH_EIGEN` is `ON` (the default). If your project
already provides Eigen, set `-DFIELDS_WITH_EIGEN=OFF` and link your own Eigen target.

## Build

```sh
cmake -S . -B build
cmake --build build
./build/tests
./build/demo
```

Eigen 5.0.1 is downloaded with CMake's FetchContent when you configure. It is controlled
by the option `FIELDS_WITH_EIGEN` (default `ON`) and builds two more targets:

```sh
./build/eigen_tests
./build/linear_algebra   # worked examples over F9, F5, ℚ, GF(256)
```

To build without Eigen (nothing is downloaded, and `eigen.hpp` can't be used):

```sh
cmake -S . -B build -DFIELDS_WITH_EIGEN=OFF
```

The rest of the library has no dependencies. It needs a C++20 compiler and is tested with
clang 23 and gcc 16.
