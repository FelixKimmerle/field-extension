#pragma once

#include "irreducible.hpp"

#include <array>
#include <tuple>
#include <type_traits>

namespace fields
{

// coefficient array for the modulus of an extension, ascending like Polynomial:
// poly<F2>(1, 1, 0, 1) is 1 + X + X^3
template <Field K, typename... Cs> constexpr std::array<K, sizeof...(Cs)> poly(const Cs&... cs)
{
    return {K(cs)...};
}

namespace detail
{

template <typename A, typename K> struct is_coefficient_array : std::false_type
{
};

template <typename K, std::size_t N> struct is_coefficient_array<std::array<K, N>, K> : std::bool_constant<(N >= 2)>
{
};

template <Field K, std::size_t N> constexpr std::array<K, N> make_monic(std::array<K, N> f)
{
    const K inv_lead = K(1) / f[N - 1];
    for (K& c : f)
    {
        c = c * inv_lead;
    }
    return f;
}

} // namespace detail

// The field L = K[X]/f(X)K[X] for an irreducible polynomial f of degree n (isk 2.6).
//
// Modulus are the coefficients a_0, ..., a_n of f (see poly), Name is the name of the
// adjoined element x := [X]_f. Every element is uniquely b_0 + b_1 x + ... + b_(n-1) x^(n-1)
// with b_j in K, and f(x) = 0, i.e. x^n = -a_(n-1) x^(n-1) - ... - a_0.
//
// Examples (isk 2.7 - 2.10):
//   C  = Extension<Real, poly<Real>(1, 0, 1), "i">     i^2 = -1
//   F4 = Extension<F2, poly<F2>(1, 1, 1), "α">         α^2 = α + 1
//   F8 = Extension<F2, poly<F2>(1, 1, 0, 1), "β">      β^3 = β + 1
//   F9 = Extension<F3, poly<F3>(1, 0, 1), "ι">         ι^2 = -1
//
// Irreducibility of f is checked at compile time where possible (finite fields, R, C,
// Q up to degree 3); for other base fields it is the responsibility of the user.
// Since K may itself be an extension, towers like F16 = F4[X]/(X^2 + X + α) work too.
template <Field K, auto Modulus, FixedString Name = "x">
    requires detail::is_coefficient_array<std::remove_cvref_t<decltype(Modulus)>, K>::value
struct Extension
{
    using Base = K;

    static constexpr std::size_t degree = std::tuple_size_v<std::remove_cvref_t<decltype(Modulus)>> - 1;

    static_assert(Modulus[degree] != K(0), "the leading coefficient of the modulus must not be 0");

    // f normalized to be monic
    static constexpr std::array<K, degree + 1> modulus = detail::make_monic(Modulus);

    static_assert(check_irreducible(Polynomial<K>(modulus)) != Irreducibility::Reducible,
                  "the modulus has to be irreducible, otherwise K[X]/f(X)K[X] has zero divisors "
                  "and is no field (isk 2.6)");

    static constexpr std::string_view generator_name = Name.view();

    static constexpr uint64_t characteristic = FieldTraits<K>::characteristic;
    // |L| = |K|^n (isk 2.6)
    static constexpr uint64_t cardinality = detail::checked_pow(FieldTraits<K>::cardinality, degree);
    // every quadratic extension of R is C
    static constexpr bool algebraically_closed = std::same_as<K, double> && degree == 2;

    // b_0, ..., b_(n-1)
    // public only so Extension is a structural type and can be the coefficient of another
    // modulus; prefer the accessors
    std::array<K, degree> coeffs{};

    constexpr Extension() = default;

    // n * 1
    template <std::integral I> constexpr Extension(I n)
    {
        coeffs[0] = K(n);
    }

    // K is a subfield of L (isk 2.11): a = a + 0x + ... + 0x^(n-1)
    constexpr Extension(const K& a)
    {
        coeffs[0] = a;
    }

    constexpr explicit Extension(const std::array<K, degree>& c) : coeffs(c)
    {
    }

    // b_0 + b_1 x + ..., missing coefficients are 0
    template <typename... Cs>
        requires(sizeof...(Cs) <= degree)
    static constexpr Extension make(const Cs&... cs)
    {
        return Extension(std::array<K, degree>{K(cs)...});
    }

    // the adjoined root x = [X]_f
    static constexpr Extension generator()
    {
        return from_polynomial(Polynomial<K>::x());
    }

    static constexpr Polynomial<K> modulus_polynomial()
    {
        return Polynomial<K>(modulus);
    }

    // [p(X)]_f, reduced by polynomial division by f
    static constexpr Extension from_polynomial(const Polynomial<K>& p)
    {
        const Polynomial<K> r = p % modulus_polynomial();
        Extension e;
        for (std::size_t j = 0; j < degree; j++)
        {
            e.coeffs[j] = r.coefficient(j);
        }
        return e;
    }

    // the representative of degree < n
    constexpr Polynomial<K> to_polynomial() const
    {
        return Polynomial<K>(coeffs);
    }

    constexpr const K& operator[](std::size_t j) const
    {
        return coeffs[j];
    }

    constexpr bool is_zero() const
    {
        return *this == Extension();
    }

    // lies in the subfield K
    constexpr bool in_base_field() const
    {
        for (std::size_t j = 1; j < degree; j++)
        {
            if (coeffs[j] != K(0))
            {
                return false;
            }
        }
        return true;
    }

    static std::string field_name()
    {
        // towers like (F_2[X]/(X^2 + X + 1))[X]/(X^2 + X + α)
        std::string base = FieldTraits<K>::name();
        if (base.find('[') != std::string::npos)
        {
            base = "(" + base + ")";
        }
        return base + "[X]/(" + detail::to_string(modulus_polynomial()) + ")";
    }

    constexpr Extension inverse() const
    {
        if (is_zero())
        {
            throw std::domain_error("division by zero");
        }
        // s * f + t * u = gcd(f, u) = 1 since f is irreducible, so t * u = 1 in K[X]/f
        const auto [g, s, t] = Polynomial<K>::extended_gcd(modulus_polynomial(), to_polynomial());
        if (g.degree() != 0)
        {
            throw std::domain_error("the modulus is reducible, this element is a zero divisor");
        }
        return from_polynomial(t);
    }

    // numbering of the elements 0, ..., |L| - 1 by the digits b_j (base |K|),
    // like numberToElement in the java version
    constexpr uint64_t index() const
        requires is_finite_field_v<K>
    {
        uint64_t result = 0;
        for (std::size_t j = degree; j > 0; j--)
        {
            result = result * FieldTraits<K>::cardinality + coeffs[j - 1].index();
        }
        return result;
    }

    static constexpr Extension from_index(uint64_t i)
        requires is_finite_field_v<K>
    {
        Extension e;
        for (std::size_t j = 0; j < degree; j++)
        {
            e.coeffs[j] = K::from_index(i % FieldTraits<K>::cardinality);
            i /= FieldTraits<K>::cardinality;
        }
        return e;
    }

    friend constexpr Extension operator+(Extension a, const Extension& b)
    {
        for (std::size_t j = 0; j < degree; j++)
        {
            a.coeffs[j] = a.coeffs[j] + b.coeffs[j];
        }
        return a;
    }

    friend constexpr Extension operator-(Extension a)
    {
        for (K& c : a.coeffs)
        {
            c = -c;
        }
        return a;
    }

    friend constexpr Extension operator-(Extension a, const Extension& b)
    {
        for (std::size_t j = 0; j < degree; j++)
        {
            a.coeffs[j] = a.coeffs[j] - b.coeffs[j];
        }
        return a;
    }

    friend constexpr Extension operator*(const Extension& a, const Extension& b)
    {
        // product in K[X] ...
        std::array<K, 2 * degree - 1> t{};
        for (std::size_t i = 0; i < degree; i++)
        {
            for (std::size_t j = 0; j < degree; j++)
            {
                t[i + j] = t[i + j] + a.coeffs[i] * b.coeffs[j];
            }
        }
        // ... reduced from the top with x^n = -a_(n-1) x^(n-1) - ... - a_0
        for (std::size_t k = 2 * degree - 1; k-- > degree;)
        {
            const K c = t[k];
            for (std::size_t j = 0; j < degree; j++)
            {
                t[k - degree + j] = t[k - degree + j] - c * modulus[j];
            }
        }
        Extension r;
        std::copy_n(t.begin(), degree, r.coeffs.begin());
        return r;
    }

    friend constexpr Extension operator/(const Extension& a, const Extension& b)
    {
        return a * b.inverse();
    }

    constexpr Extension& operator+=(const Extension& b)
    {
        return *this = *this + b;
    }
    constexpr Extension& operator-=(const Extension& b)
    {
        return *this = *this - b;
    }
    constexpr Extension& operator*=(const Extension& b)
    {
        return *this = *this * b;
    }
    constexpr Extension& operator/=(const Extension& b)
    {
        return *this = *this / b;
    }

    friend constexpr bool operator==(const Extension&, const Extension&) = default;

    // prints it like the script: 1 + β + β^2
    friend std::ostream& operator<<(std::ostream& os, const Extension& e)
    {
        std::vector<std::pair<K, std::size_t>> terms;
        for (std::size_t j = 0; j < degree; j++)
        {
            terms.emplace_back(e.coeffs[j], j);
        }
        detail::print_terms(os, terms, generator_name);
        return os;
    }
};

} // namespace fields
