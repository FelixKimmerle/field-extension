#pragma once

#include "field.hpp"

#include <compare>

namespace fields
{

namespace detail
{
// __extension__ silences -Wpedantic, 128 bit integers are a gcc/clang extension
__extension__ using int128 = __int128;
} // namespace detail

// Q, exact fractions num / den
// Intermediate results are computed with 128 bits, results that do not fit
// into 64 bits throw std::overflow_error.
struct Rational
{
    static constexpr uint64_t characteristic = 0;
    static constexpr uint64_t cardinality = 0;
    static constexpr bool algebraically_closed = false;

    // invariant: den > 0 and gcd(num, den) = 1
    // public only so Rational is a structural type, never write them directly
    int64_t num = 0;
    int64_t den = 1;

    constexpr Rational() = default;

    template <std::integral I> constexpr Rational(I n) : num(static_cast<int64_t>(n))
    {
    }

    template <std::integral I, std::integral J> constexpr Rational(I n, J d)
    {
        *this = make(n, d);
    }

    static std::string field_name()
    {
        return "Q";
    }

    constexpr Rational inverse() const
    {
        if (num == 0)
        {
            throw std::domain_error("division by zero");
        }
        return make(den, num);
    }

    friend constexpr Rational operator+(const Rational& a, const Rational& b)
    {
        return make(static_cast<detail::int128>(a.num) * b.den + static_cast<detail::int128>(b.num) * a.den,
                    static_cast<detail::int128>(a.den) * b.den);
    }

    friend constexpr Rational operator-(const Rational& a)
    {
        return make(-static_cast<detail::int128>(a.num), a.den);
    }

    friend constexpr Rational operator-(const Rational& a, const Rational& b)
    {
        return a + (-b);
    }

    friend constexpr Rational operator*(const Rational& a, const Rational& b)
    {
        return make(static_cast<detail::int128>(a.num) * b.num, static_cast<detail::int128>(a.den) * b.den);
    }

    friend constexpr Rational operator/(const Rational& a, const Rational& b)
    {
        return a * b.inverse();
    }

    constexpr Rational& operator+=(const Rational& b)
    {
        return *this = *this + b;
    }
    constexpr Rational& operator-=(const Rational& b)
    {
        return *this = *this - b;
    }
    constexpr Rational& operator*=(const Rational& b)
    {
        return *this = *this * b;
    }
    constexpr Rational& operator/=(const Rational& b)
    {
        return *this = *this / b;
    }

    friend constexpr bool operator==(const Rational&, const Rational&) = default;

    friend constexpr std::strong_ordering operator<=>(const Rational& a, const Rational& b)
    {
        return static_cast<detail::int128>(a.num) * b.den <=> static_cast<detail::int128>(b.num) * a.den;
    }

    friend std::ostream& operator<<(std::ostream& os, const Rational& a)
    {
        os << a.num;
        if (a.den != 1)
        {
            os << '/' << a.den;
        }
        return os;
    }

    explicit operator double() const
    {
        return static_cast<double>(num) / static_cast<double>(den);
    }

  private:
    static constexpr detail::int128 gcd(detail::int128 a, detail::int128 b)
    {
        a = a < 0 ? -a : a;
        b = b < 0 ? -b : b;
        while (b != 0)
        {
            detail::int128 t = a % b;
            a = b;
            b = t;
        }
        return a;
    }

    static constexpr Rational make(detail::int128 n, detail::int128 d)
    {
        if (d == 0)
        {
            throw std::domain_error("division by zero");
        }
        if (d < 0)
        {
            n = -n;
            d = -d;
        }
        detail::int128 g = gcd(n, d);
        n /= g;
        d /= g;
        constexpr detail::int128 max = std::numeric_limits<int64_t>::max();
        if (n > max || n < -max || d > max)
        {
            throw std::overflow_error("rational number does not fit into 64 bits");
        }
        Rational r;
        r.num = static_cast<int64_t>(n);
        r.den = static_cast<int64_t>(d);
        return r;
    }
};

} // namespace fields
