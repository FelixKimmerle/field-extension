#pragma once

#include "extension.hpp"
#include "prime_field.hpp"

#include <functional>
#include <iomanip>

namespace fields
{

// q = p^n -> {p, n}, {0, 0} if q is no prime power
constexpr std::pair<uint64_t, uint64_t> prime_power(uint64_t q)
{
    if (q < 2)
    {
        return {0, 0};
    }
    uint64_t p = q;
    for (uint64_t d = 2; d <= q / d; d++)
    {
        if (q % d == 0)
        {
            p = d;
            break;
        }
    }
    uint64_t n = 0;
    while (q % p == 0)
    {
        q /= p;
        n++;
    }
    if (q != 1)
    {
        return {0, 0};
    }
    return {p, n};
}

constexpr std::vector<uint64_t> prime_factors(uint64_t n)
{
    std::vector<uint64_t> result;
    for (uint64_t d = 2; d <= n / d; d++)
    {
        if (n % d == 0)
        {
            result.push_back(d);
            while (n % d == 0)
            {
                n /= d;
            }
        }
    }
    if (n > 1)
    {
        result.push_back(n);
    }
    return result;
}

// The first monic irreducible polynomial of degree N over F_P, where X^N + a_(N-1) X^(N-1) +
// ... + a_0 is ordered by the number a_0 + a_1 P + ... + a_(N-1) P^(N-1).
// This yields exactly the moduli of the script: X^2 + X + 1 (F4), X^3 + X + 1 (F8), X^2 + 1 (F9).
template <uint64_t P, std::size_t N> constexpr std::array<Fp<P>, N + 1> first_irreducible()
{
    for (uint64_t i = 0;; i++)
    {
        std::array<Fp<P>, N + 1> f{};
        uint64_t m = i;
        for (std::size_t k = 0; k < N; k++)
        {
            f[k] = Fp<P>(m % P);
            m /= P;
        }
        f[N] = Fp<P>(1);
        // divisible by X
        if (N > 1 && f[0] == Fp<P>(0))
        {
            continue;
        }
        if (check_irreducible(Polynomial<Fp<P>>(f)) == Irreducibility::Irreducible)
        {
            return f;
        }
    }
}

namespace detail
{

template <uint64_t Q, FixedString Name, uint64_t P = prime_power(Q).first, uint64_t N = prime_power(Q).second>
struct GaloisField
{
    static_assert(P != 0, "a finite field with q elements exists iff q is a prime power");
    using type = Extension<Fp<P>, first_irreducible<P, N>(), Name>;
};

template <uint64_t Q, FixedString Name, uint64_t P> struct GaloisField<Q, Name, P, 1>
{
    using type = Fp<P>;
};

} // namespace detail

// The finite field with q = p^n elements, like FiniteField.create(q) in the java version:
// F_p for n = 1, otherwise F_p[X]/f(X) with f = first_irreducible<p, n>().
// Name is the name of the adjoined element, GF<8, "β"> is the same type as F8.
template <uint64_t Q, FixedString Name = "a"> using GF = typename detail::GaloisField<Q, Name>::type;

// all elements, ordered by index()
template <FiniteField K> std::vector<K> all_elements()
{
    std::vector<K> result;
    result.reserve(FieldTraits<K>::cardinality);
    for (uint64_t i = 0; i < FieldTraits<K>::cardinality; i++)
    {
        result.push_back(K::from_index(i));
    }
    return result;
}

// smallest k >= 1 with x^k = 1, a divisor of q - 1
template <FiniteField K> constexpr uint64_t multiplicative_order(const K& x)
{
    if (x == K(0))
    {
        throw std::domain_error("0 has no multiplicative order");
    }
    uint64_t order = FieldTraits<K>::cardinality - 1;
    for (uint64_t r : prime_factors(order))
    {
        while (order % r == 0 && pow(x, static_cast<int64_t>(order / r)) == K(1))
        {
            order /= r;
        }
    }
    return order;
}

// x generates the cyclic group K \ {0}
template <FiniteField K> constexpr bool is_primitive(const K& x)
{
    return x != K(0) && multiplicative_order(x) == FieldTraits<K>::cardinality - 1;
}

// the primitive element with the smallest index
template <FiniteField K> constexpr K primitive_element()
{
    for (uint64_t i = 1;; i++)
    {
        if (is_primitive(K::from_index(i)))
        {
            return K::from_index(i);
        }
    }
}

// x -> x^p, an automorphism of K fixing F_p
template <FiniteField K> constexpr K frobenius(const K& x)
{
    return pow(x, static_cast<int64_t>(FieldTraits<K>::characteristic));
}

// all monic irreducible polynomials of degree n over K
template <FiniteField K> std::vector<Polynomial<K>> irreducible_polynomials(std::size_t n)
{
    constexpr uint64_t q = FieldTraits<K>::cardinality;
    const uint64_t count = detail::checked_pow(q, n);
    std::vector<Polynomial<K>> result;
    for (uint64_t i = 0; i < count; i++)
    {
        std::vector<K> f(n + 1);
        uint64_t m = i;
        for (std::size_t k = 0; k < n; k++)
        {
            f[k] = K::from_index(m % q);
            m /= q;
        }
        f[n] = K(1);
        Polynomial<K> p(std::move(f));
        if (is_irreducible(p))
        {
            result.push_back(std::move(p));
        }
    }
    return result;
}

namespace detail
{

template <FiniteField K>
void print_table(std::ostream& os, std::string_view symbol, const std::function<K(const K&, const K&)>& op)
{
    const std::vector<K> elements = all_elements<K>();
    std::vector<std::string> names;
    std::size_t width = display_width(symbol);
    for (const K& e : elements)
    {
        names.push_back(to_string(e));
        width = std::max(width, display_width(names.back()));
    }
    auto cell = [&](std::string_view s) { os << std::string(width - display_width(s) + 1, ' ') << s << ' '; };

    cell(symbol);
    os << '|';
    for (const std::string& name : names)
    {
        cell(name);
    }
    os << '\n' << std::string(width + 2, '-') << '+' << std::string((width + 2) * names.size(), '-') << '\n';
    for (const K& a : elements)
    {
        cell(to_string(a));
        os << '|';
        for (const K& b : elements)
        {
            cell(to_string(op(a, b)));
        }
        os << '\n';
    }
}

} // namespace detail

// Additionstafel like in isk 2.1 and 2.8
template <FiniteField K> void print_addition_table(std::ostream& os)
{
    detail::print_table<K>(os, "(+)", [](const K& a, const K& b) { return a + b; });
}

// Multiplikationstafel like in isk 2.1 and 2.8
template <FiniteField K> void print_multiplication_table(std::ostream& os)
{
    detail::print_table<K>(os, "(·)", [](const K& a, const K& b) { return a * b; });
}

} // namespace fields
