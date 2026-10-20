// Walks through the examples of chapter 2 of the script (isk)

#include <fields/fields.hpp>

#include <iostream>
#include <numbers>

using namespace fields;

void complex_numbers()
{
    std::cout << "=== " << field_name<Complex>() << " (isk 2.7) ===\n";
    const Complex i = Complex::generator();
    const Complex z = 2 + 5 * i;
    const Complex w = 4 - 3 * i;
    std::cout << "(" << z << ") + (" << w << ") = " << z + w << '\n';
    std::cout << "(" << z << ") * (" << w << ") = " << z * w << '\n';
    std::cout << "i^2 = " << i * i << '\n';
    std::cout << "1 / (" << z << ") = " << 1 / z << "  (= (2 - 5i) / 29)\n";
    std::cout << "conj(" << z << ") = " << conjugate(z) << ", norm = " << norm(z) << '\n';

    const Complex v = -1 + i;
    std::cout << "|" << v << "| = " << abs(v) << ", arg = " << arg(v) / std::numbers::pi
              << " pi\n";
    std::cout << "(" << v << ")^2 = " << pow(v, 2) << '\n';
    std::cout << "from_polar(1, pi/2) = " << from_polar<Complex>(1, std::numbers::pi / 2) << "\n\n";
}

void field_f4()
{
    std::cout << "=== F4 = " << field_name<F4>() << " (isk 2.8) ===\n";
    print_addition_table<F4>(std::cout);
    std::cout << '\n';
    print_multiplication_table<F4>(std::cout);
    const F4 a = F4::generator();
    std::cout << "1/(1 + α) = " << 1 / (1 + a) << ", 1/α = " << 1 / a << "\n\n";
}

void field_f8()
{
    std::cout << "=== F8 = " << field_name<F8>() << " (isk 2.9) ===\n";
    const F8 b = F8::generator();
    for (int k = 0; k <= 7; k++)
    {
        std::cout << "β^" << k << " = " << pow(b, k) << '\n';
    }
    std::cout << "β^-1 = " << pow(b, -1) << ", (β^2)^-1 = " << 1 / (b * b)
              << ", (β + 1)^-1 = " << 1 / (b + 1) << "\n\n";
}

void field_f9()
{
    std::cout << "=== F9 = " << field_name<F9>() << " (isk 2.10) ===\n";
    const F9 i = F9::generator();
    std::cout << "ι^2 = " << i * i << '\n';
    std::cout << "(1 + ι)^3 = " << pow(1 + i, 3) << '\n';
    std::cout << "1 / (1 + ι) = " << 1 / (1 + i) << '\n';
    std::cout << "order of ι = " << multiplicative_order(i)
              << ", primitive element: " << primitive_element<F9>() << "\n\n";
}

void galois_fields()
{
    std::cout << "=== GF(q), like FiniteField.create(q) ===\n";
    std::cout << "GF(7)   = " << field_name<GF<7>>() << '\n';
    std::cout << "GF(16)  = " << field_name<GF<16>>() << '\n';
    std::cout << "GF(25)  = " << field_name<GF<25>>() << '\n';
    std::cout << "GF(27)  = " << field_name<GF<27>>() << '\n';
    std::cout << "GF(256) = " << field_name<GF<256>>() << '\n';

    // F16 as a tower over F4, the modulus has a coefficient α in F4
    using F16 = Extension<F4, poly<F4>(F4::generator(), 1, 1), "γ">;
    const F16 g = F16::generator();
    std::cout << "F16 = " << field_name<F16>() << ": γ^2 = " << g * g
              << ", 1/γ = " << 1 / g << ", |F16| = " << F16::cardinality << '\n';

    std::cout << "irreducible polynomials of degree 4 over F2:";
    for (const auto& f : irreducible_polynomials<F2>(4))
    {
        std::cout << "  " << f;
    }
    std::cout << "\n\n";
}

void rationals_and_polynomials()
{
    std::cout << "=== Polynomials (isk 2.4, 2.5) ===\n";
    const Polynomial<Real> f{1, 1, 0, 0, 0, 0, 1};
    const Polynomial<Real> g{2, -1, 1};
    const auto [h, r] = f.divmod(g);
    std::cout << f << " = (" << g << ") * (" << h << ") + (" << r << ")\n";

    // Q(√2) = Q[X]/(X^2 - 2)
    using Q2 = Extension<Q, poly<Q>(-2, 0, 1), "√2">;
    const Q2 s = Q2::generator();
    std::cout << "in " << field_name<Q2>() << ": 1 / (1 + √2) = " << 1 / (1 + s)
              << ", 1 / (3 + 2/7 √2) = " << 1 / (3 + Q(2, 7) * s) << '\n';
}

int main()
{
    complex_numbers();
    field_f4();
    field_f8();
    field_f9();
    galois_fields();
    rationals_and_polynomials();
}
