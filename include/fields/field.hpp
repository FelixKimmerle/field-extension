#pragma once

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace fields
{

// compile time string, used to give the adjoined element of an extension a name
// (i for C, α for F4, β for F8, ι for F9 like in the script)
template <std::size_t N> struct FixedString
{
    char data[N]{};

    constexpr FixedString(const char (&s)[N])
    {
        std::copy_n(s, N, data);
    }

    constexpr std::string_view view() const
    {
        return {data, N - 1};
    }
};

// Information about a field that is not part of its elements.
// Field types provide it as static members, fundamental types (double) specialize it.
template <typename K> struct FieldTraits
{
};

template <typename K>
    requires requires {
        K::characteristic;
        K::cardinality;
        K::algebraically_closed;
    }
struct FieldTraits<K>
{
    // 0 for fields of characteristic 0 (Q, R, C)
    static constexpr uint64_t characteristic = K::characteristic;
    // number of elements, 0 for infinite fields
    static constexpr uint64_t cardinality = K::cardinality;
    // every polynomial of degree >= 1 has a root (isk 2.12)
    static constexpr bool algebraically_closed = K::algebraically_closed;

    static std::string name()
    {
        return K::field_name();
    }
};

// A field in the sense of isk 2.3: a set with +, -, *, / (not by 0) and 0 != 1.
// K(n) has to be n * 1_K for integers n.
template <typename K>
concept Field = std::regular<K> && std::constructible_from<K, int> && requires(const K a, const K b, std::ostream& os) {
    { a + b } -> std::convertible_to<K>;
    { a - b } -> std::convertible_to<K>;
    { a * b } -> std::convertible_to<K>;
    { a / b } -> std::convertible_to<K>;
    { -a } -> std::convertible_to<K>;
    { os << a } -> std::convertible_to<std::ostream&>;
    { FieldTraits<K>::characteristic } -> std::convertible_to<uint64_t>;
    { FieldTraits<K>::cardinality } -> std::convertible_to<uint64_t>;
    { FieldTraits<K>::algebraically_closed } -> std::convertible_to<bool>;
};

template <typename K> inline constexpr bool is_finite_field_v = FieldTraits<K>::cardinality != 0;

template <typename K>
concept FiniteField = Field<K> && is_finite_field_v<K>;

template <Field K> std::string field_name()
{
    return FieldTraits<K>::name();
}

template <Field K> constexpr K inverse(const K& x)
{
    return K(1) / x;
}

// x^e via square and multiply, negative exponents use the inverse
template <Field K> constexpr K pow(K base, int64_t exponent)
{
    if (exponent < 0)
    {
        base = inverse(base);
        exponent = -exponent;
    }
    K result(1);
    while (exponent > 0)
    {
        if (exponent % 2 == 1)
        {
            result = result * base;
        }
        base = base * base;
        exponent /= 2;
    }
    return result;
}

namespace detail
{

// R and Q have a meaningful sign; finite fields and extensions never do,
// even if they get an ordering (eigen.hpp adds one for pivoting)
template <typename K>
inline constexpr bool has_sign_v = std::totally_ordered<K> && FieldTraits<K>::characteristic == 0 &&
                                   !requires { typename K::Base; };

constexpr uint64_t checked_pow(uint64_t base, uint64_t exponent)
{
    if (base == 0)
    {
        return 0;
    }
    uint64_t result = 1;
    for (uint64_t i = 0; i < exponent; i++)
    {
        if (result > std::numeric_limits<uint64_t>::max() / base)
        {
            throw std::overflow_error("cardinality does not fit into 64 bits");
        }
        result *= base;
    }
    return result;
}

template <typename T> std::string to_string(const T& x)
{
    std::ostringstream os;
    os << x;
    return os.str();
}

// number of displayed characters of an utf-8 string
inline std::size_t display_width(std::string_view s)
{
    return std::count_if(s.begin(), s.end(), [](char c) { return (c & 0xC0) != 0x80; });
}

// prints the sum of c * var^k for the given (c, k) pairs in the given order
template <Field K>
void print_terms(std::ostream& os, const std::vector<std::pair<K, std::size_t>>& terms, std::string_view var)
{
    bool first = true;
    for (const auto& [c, k] : terms)
    {
        if (c == K(0))
        {
            continue;
        }

        K coeff = c;
        bool negative = false;
        if constexpr (has_sign_v<K>)
        {
            if (coeff < K(0))
            {
                negative = true;
                coeff = -coeff;
            }
        }
        if (first)
        {
            os << (negative ? "-" : "");
        }
        else
        {
            os << (negative ? " - " : " + ");
        }
        first = false;

        if (k == 0)
        {
            os << coeff;
            continue;
        }
        if (coeff != K(1))
        {
            // coefficients of nested extensions like "1 + α" need parentheses
            std::string s = to_string(coeff);
            if (s.find_first_of(" /") != std::string::npos)
            {
                os << '(' << s << ')';
            }
            else
            {
                os << s;
            }
        }
        os << var;
        if (k > 1)
        {
            os << '^' << k;
        }
    }
    if (first)
    {
        os << '0';
    }
}

} // namespace detail

} // namespace fields
