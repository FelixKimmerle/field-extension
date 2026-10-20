#pragma once

#include "field.hpp"

namespace fields
{

constexpr bool is_prime(uint64_t n)
{
    if (n < 2)
    {
        return false;
    }
    for (uint64_t d = 2; d <= n / d; d++)
    {
        if (n % d == 0)
        {
            return false;
        }
    }
    return true;
}

// F_p = Z/pZ, the residue classes modulo a prime p (isk 2.1, 2.3)
template <uint64_t P> struct Fp
{
    static_assert(is_prime(P), "Z/nZ is a field iff n is prime (isk 2.3)");
    static_assert(P < (uint64_t(1) << 32), "products of representatives must fit into 64 bits");

    static constexpr uint64_t characteristic = P;
    static constexpr uint64_t cardinality = P;
    static constexpr bool algebraically_closed = false;

    // representative in [0, P)
    // public only so Fp is a structural type and can be used in template arguments,
    // never write it directly
    uint64_t value = 0;

    constexpr Fp() = default;

    // [n]_p
    template <std::integral I> constexpr Fp(I n) : value(reduce(n))
    {
    }

    static std::string field_name()
    {
        return "F_" + std::to_string(P);
    }

    // numbering of the elements 0, ..., P - 1
    constexpr uint64_t index() const
    {
        return value;
    }

    static constexpr Fp from_index(uint64_t i)
    {
        return Fp(i);
    }

    constexpr Fp inverse() const
    {
        if (value == 0)
        {
            throw std::domain_error("division by zero");
        }
        // extended euclidean algorithm (isk 1.9): t * value + s * P = gcd = 1
        int64_t r0 = P, r1 = value;
        int64_t t0 = 0, t1 = 1;
        while (r1 != 0)
        {
            int64_t q = r0 / r1;
            int64_t r2 = r0 - q * r1;
            r0 = r1;
            r1 = r2;
            int64_t t2 = t0 - q * t1;
            t0 = t1;
            t1 = t2;
        }
        return Fp(t0);
    }

    friend constexpr Fp operator+(Fp a, const Fp& b)
    {
        a.value += b.value;
        if (a.value >= P)
        {
            a.value -= P;
        }
        return a;
    }

    friend constexpr Fp operator-(const Fp& a)
    {
        Fp r;
        r.value = a.value == 0 ? 0 : P - a.value;
        return r;
    }

    friend constexpr Fp operator-(const Fp& a, const Fp& b)
    {
        return a + (-b);
    }

    friend constexpr Fp operator*(const Fp& a, const Fp& b)
    {
        return Fp(a.value * b.value);
    }

    friend constexpr Fp operator/(const Fp& a, const Fp& b)
    {
        return a * b.inverse();
    }

    constexpr Fp& operator+=(const Fp& b)
    {
        return *this = *this + b;
    }
    constexpr Fp& operator-=(const Fp& b)
    {
        return *this = *this - b;
    }
    constexpr Fp& operator*=(const Fp& b)
    {
        return *this = *this * b;
    }
    constexpr Fp& operator/=(const Fp& b)
    {
        return *this = *this / b;
    }

    friend constexpr bool operator==(const Fp&, const Fp&) = default;

    friend std::ostream& operator<<(std::ostream& os, const Fp& a)
    {
        return os << a.value;
    }

  private:
    template <std::integral I> static constexpr uint64_t reduce(I n)
    {
        if constexpr (std::is_signed_v<I>)
        {
            int64_t r = static_cast<int64_t>(n) % static_cast<int64_t>(P);
            return static_cast<uint64_t>(r < 0 ? r + static_cast<int64_t>(P) : r);
        }
        else
        {
            return static_cast<uint64_t>(n) % P;
        }
    }
};

} // namespace fields
