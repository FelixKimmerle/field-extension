#include <fields/fields.hpp>

#include <cmath>
#include <iostream>
#include <numbers>

using namespace fields;

static int failures = 0;

#define CHECK(cond)                                                                      \
    do                                                                                   \
    {                                                                                    \
        if (!(cond))                                                                     \
        {                                                                                \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n";   \
            failures++;                                                                  \
        }                                                                                \
    } while (0)

// ----------------------------------------------------------------------------------
// compile time checks
// ----------------------------------------------------------------------------------

static_assert(Field<Real> && Field<Q> && Field<F2> && Field<F4> && Field<Complex>);
static_assert(FiniteField<F9> && !FiniteField<Complex>);
static_assert(F4::cardinality == 4 && F8::cardinality == 8 && F9::cardinality == 9);
static_assert(F9::characteristic == 3 && Complex::characteristic == 0);
static_assert(Complex::algebraically_closed && !FieldTraits<Real>::algebraically_closed);

// GF finds the moduli of the script
static_assert(std::same_as<GF<2>, F2>);
static_assert(std::same_as<GF<4, "α">, F4>);
static_assert(std::same_as<GF<8, "β">, F8>);
static_assert(std::same_as<GF<9, "ι">, F9>);

// arithmetic works in constant expressions
static_assert(pow(F8::generator(), 7) == F8(1));
static_assert(Complex::generator() * Complex::generator() == Complex(-1));
static_assert(F3(2).inverse() == F3(2));

// irreducibility (isk 2.5)
static_assert(is_irreducible(Polynomial<Real>{1, 0, 1}));
static_assert(!is_irreducible(Polynomial<Real>{-2, 0, 1}));
static_assert(!is_irreducible(Polynomial<Real>{1, 0, 2, 0, 1})); // (X^2 + 1)^2
static_assert(is_irreducible(Polynomial<F2>{1, 1, 1}));
static_assert(!is_irreducible(Polynomial<F2>{1, 0, 1})); // (X + 1)^2
static_assert(is_irreducible(Polynomial<F3>{1, 2, 0, 1})); // X^3 - X + 1 (isk 2.12)
static_assert(is_irreducible(Polynomial<Q>{-2, 0, 1}));
static_assert(!is_irreducible(Polynomial<Q>{-4, 11, 2, 3})); // root 1/3 (isk 2.12)
static_assert(check_irreducible(Polynomial<Q>{1, 0, 0, 0, 1}) == Irreducibility::Unknown);
static_assert(!is_irreducible(Polynomial<Complex>{1, 0, 1}));
static_assert(prime_power(9) == std::pair<uint64_t, uint64_t>{3, 2});
static_assert(prime_power(12).first == 0);

// ----------------------------------------------------------------------------------
// helpers
// ----------------------------------------------------------------------------------

// all field axioms of isk 2.3, checked exhaustively
template <FiniteField K> void check_axioms()
{
    const auto elements = all_elements<K>();
    CHECK(elements.size() == FieldTraits<K>::cardinality);
    CHECK(K(0) != K(1));
    for (const K& x : elements)
    {
        CHECK(K::from_index(x.index()) == x);
        CHECK(x + K(0) == x);
        CHECK(x * K(1) == x);
        CHECK(x + (-x) == K(0));
        CHECK(K(static_cast<int64_t>(FieldTraits<K>::characteristic)) * x == K(0));
        if (x != K(0))
        {
            CHECK(x * inverse(x) == K(1));
            // Lagrange: x^(q-1) = 1
            CHECK(pow(x, static_cast<int64_t>(FieldTraits<K>::cardinality - 1)) == K(1));
        }
        for (const K& y : elements)
        {
            CHECK(x + y == y + x);
            CHECK(x * y == y * x);
            if (FieldTraits<K>::cardinality <= 32)
            {
                for (const K& z : elements)
                {
                    CHECK((x + y) + z == x + (y + z));
                    CHECK((x * y) * z == x * (y * z));
                    CHECK(x * (y + z) == x * y + x * z);
                }
            }
        }
    }
}

bool close(const Complex& a, const Complex& b)
{
    return abs(a - b) < 1e-12;
}

template <typename T> std::string str(const T& x)
{
    return detail::to_string(x);
}

// ----------------------------------------------------------------------------------
// tests
// ----------------------------------------------------------------------------------

void test_polynomials()
{
    using P = Polynomial<Real>;

    // isk 2.4
    CHECK((P{7, 2, 1} * P{3, 5}) == (P{21, 41, 13, 5}));

    // isk 2.5: X^6 + X + 1 = (X^2 - X + 2)(X^4 + X^3 - X^2 - 3X - 1) + (6X + 3)
    auto [q, r] = P{1, 1, 0, 0, 0, 0, 1}.divmod(P{2, -1, 1});
    CHECK(q == (P{-1, -3, -1, 1, 1}));
    CHECK(r == (P{3, 6}));

    // the cases of the old main.cpp
    CHECK((P{2, 3, 1} / P{1, 1}) == (P{2, 1}));
    CHECK((P{2, 3, 1} % P{1, 1}).is_zero());
    CHECK((P{5, 7, 7, 2} / P{2, 1}) == (P{1, 3, 2}));
    CHECK((P{5, 7, 7, 2} % P{2, 1}) == P{3});
    CHECK((P{1, 1} / P{1, 0, 1}).is_zero());
    CHECK((P{1, 1} % P{1, 0, 1}) == (P{1, 1}));
    CHECK((P{6, 2, 4} / P{2}) == (P{3, 1, 2}));
    CHECK((P{0, 1, 1, 1} / P{0, 1, 1}) == (P{0, 1}));
    CHECK((P{0, 1, 1, 1} % P{0, 1, 1}) == (P{0, 1}));

    CHECK(P().degree() == -1);
    CHECK((P{1, 2, 3})(2) == 17);
    CHECK((P{1, 2, 3}).derivative() == (P{2, 6}));
    CHECK(P::gcd(P{-1, 0, 1}, P{1, 2, 1}) == (P{1, 1}));

    auto [g, s, t] = Polynomial<Q>::extended_gcd(Polynomial<Q>{-1, 0, 1}, Polynomial<Q>{1, 2, 1});
    CHECK(g == (Polynomial<Q>{1, 1}));
    CHECK((s * Polynomial<Q>{-1, 0, 1} + t * Polynomial<Q>{1, 2, 1}) == g);

    CHECK(str(P{1, 1, 0, 1}) == "X^3 + X + 1");
    CHECK(str(P{3, -2, 0, -1}) == "-X^3 - 2X + 3");
    CHECK(str(Polynomial<F2>{1, 1, 1}) == "X^2 + X + 1");
    CHECK(str(P()) == "0");
}

void test_prime_fields()
{
    check_axioms<F2>();
    check_axioms<F3>();
    check_axioms<F7>();
    check_axioms<Fp<31>>();
    CHECK(F3(7) == F3(1));
    CHECK(F3(-1) == F3(2));
    CHECK(F7(3) / F7(5) * F7(5) == F7(3));
}

void test_rationals()
{
    CHECK(Q(2, 4) == Q(1, 2));
    CHECK(Q(1, -2) == Q(-1, 2));
    CHECK(Q(1, 2) + Q(1, 3) == Q(5, 6));
    CHECK(Q(1, 2) < Q(2, 3));
    CHECK(Q(3, 4) / Q(3, 4) == Q(1));
    CHECK(str(Q(-3, 6)) == "-1/2");

    using Q2 = Extension<Q, poly<Q>(-2, 0, 1), "√2">;
    const Q2 s = Q2::generator();
    CHECK(s * s == Q2(2));
    CHECK((1 + s) * (1 - s) == Q2(-1));
    CHECK(1 / (1 + s) == -1 + s);
    const Q2 z = Q(3, 5) + Q(-7, 2) * s;
    CHECK(z * inverse(z) == Q2(1));
    CHECK(norm(z) == Q(9, 25) - 2 * Q(49, 4));
}

void test_complex()
{
    const Complex i = Complex::generator();
    // isk 2.7
    CHECK((2 + 5 * i) + (4 - 3 * i) == 6 + 2 * i);
    CHECK((2 + 5 * i) * (4 - 3 * i) == 23 + 14 * i);
    CHECK(close(1 / (2 + 5 * i), Complex::make(2.0 / 29, -5.0 / 29)));
    // isk 2.12
    CHECK(conjugate(2 + 5 * i) == 2 - 5 * i);
    CHECK(norm(2 + 5 * i) == 29);
    CHECK(abs(3 + 4 * i) == 5);
    CHECK(std::abs(arg(-1 + i) - 3 * std::numbers::pi / 4) < 1e-12);
    CHECK(pow(-1 + i, 2) == -2 * i);
    CHECK(close(from_polar<Complex>(std::sqrt(2), 3 * std::numbers::pi / 4), -1 + i));
    CHECK(re(2 + 5 * i) == 2 && im(2 + 5 * i) == 5);
    CHECK(Complex(2.5).in_base_field());

    CHECK(str(23 + 14 * i) == "23 + 14i");
    CHECK(str(2 - 5 * i) == "2 - 5i");
    CHECK(str(-i) == "-i");
    CHECK(str(Complex()) == "0");

    // division by zero
    bool thrown = false;
    try
    {
        (void)(1 / Complex());
    }
    catch (const std::domain_error&)
    {
        thrown = true;
    }
    CHECK(thrown);
}

void test_f4()
{
    check_axioms<F4>();
    const F4 a = F4::generator();
    // isk 2.8
    CHECK(a * a == a + 1);
    CHECK((1 + a) + a == F4(1));
    CHECK((1 + a) * a == F4(1));
    CHECK(1 / (1 + a) == a);
    CHECK(1 / a == a + 1);
    CHECK(conjugate(a) == frobenius(a));
    CHECK(str(1 + a) == "1 + α");
}

void test_f8()
{
    check_axioms<F8>();
    const F8 b = F8::generator();
    // isk 2.9
    CHECK(pow(b, 3) == b + 1);
    CHECK(pow(b, 4) == b * b + b);
    CHECK(pow(b, 5) == b * b + b + 1);
    CHECK(pow(b, 6) == b * b + 1);
    CHECK(pow(b, 7) == F8(1));
    CHECK(inverse(b) == b * b + 1);
    CHECK(inverse(b * b) == b * b + b + 1);
    CHECK(inverse(b + 1) == b * b + b);
    CHECK(is_primitive(b));
    CHECK(str(1 + b + b * b) == "1 + β + β^2");
}

void test_f9()
{
    check_axioms<F9>();
    const F9 i = F9::generator();
    // isk 2.10
    CHECK(i * i == F9(-1));
    CHECK(pow(1 + i, 3) == 1 - i);
    CHECK(1 / (1 + i) == -1 + i);
    CHECK(multiplicative_order(i) == 4);
    CHECK(is_primitive(1 + i));
    CHECK(norm(1 + i) == F3(2));
    CHECK(str(2 + 2 * i) == "2 + 2ι");
}

void test_galois_fields()
{
    check_axioms<GF<16>>();
    check_axioms<GF<25>>();
    check_axioms<GF<27>>();
    check_axioms<GF<32>>();
    check_axioms<GF<49>>();
    check_axioms<GF<64>>();
    check_axioms<GF<125>>();
    check_axioms<GF<256>>();
    check_axioms<GF<3 * 3 * 3 * 3 * 3>>();

    // number of monic irreducible polynomials over F2 of degree n (necklace polynomial)
    const std::size_t counts[] = {0, 2, 1, 2, 3, 6, 9, 18, 30};
    for (std::size_t n = 1; n <= 8; n++)
    {
        CHECK(irreducible_polynomials<F2>(n).size() == counts[n]);
    }
    CHECK(irreducible_polynomials<F3>(2).size() == 3);
    CHECK(irreducible_polynomials<F4>(2).size() == 6);

    // every element of F_q is a root of X^q - X
    for (const auto& x : all_elements<GF<27>>())
    {
        CHECK(pow(x, 27) == x);
    }
}

void test_towers()
{
    // F16 = F4[X]/(X^2 + X + α)
    using F16 = Extension<F4, poly<F4>(F4::generator(), 1, 1), "γ">;
    static_assert(F16::cardinality == 16 && F16::characteristic == 2);
    check_axioms<F16>();
    const F16 g = F16::generator();
    const F4 a = F4::generator();
    CHECK(g * g == g + a);
    CHECK(str(g * g * g) == "α + (1 + α)γ");

    // F81 = F9[X]/(X^2 - (1 + ι)), 1 + ι is no square in F9 since it is primitive
    using F81 = Extension<F9, poly<F9>(-(1 + F9::generator()), 0, 1), "κ">;
    check_axioms<F81>();
    CHECK(multiplicative_order(F81::generator()) == 16);
}

int main()
{
    test_polynomials();
    test_prime_fields();
    test_rationals();
    test_complex();
    test_f4();
    test_f8();
    test_f9();
    test_galois_fields();
    test_towers();

    if (failures == 0)
    {
        std::cout << "all tests passed\n";
        return 0;
    }
    std::cout << failures << " checks failed\n";
    return 1;
}
