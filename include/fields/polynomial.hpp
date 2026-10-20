#pragma once

#include "field.hpp"

#include <array>
#include <initializer_list>
#include <vector>

namespace fields
{

// Polynomial in K[X] over a field K (isk 2.4, 2.5)
// It has the invariant that there are no leading zeros
template <Field K> class Polynomial
{
  private:
    // little endian storage: coeffs[j] is the coefficient of X^j
    std::vector<K> coeffs;

    // remove "leading" zeros
    // since we use little endian storage we remove zeros from the back
    constexpr void leading_zero_suppression()
    {
        while (!coeffs.empty() && coeffs.back() == K(0))
        {
            coeffs.pop_back();
        }
    }

  public:
    struct DivisionResult
    {
        Polynomial quotient;
        Polynomial remainder;
    };

    // s * a + t * b = gcd
    struct ExtendedGcdResult
    {
        Polynomial gcd;
        Polynomial s;
        Polynomial t;
    };

    constexpr Polynomial() = default;

    constexpr Polynomial(std::vector<K> coeffs) : coeffs(std::move(coeffs))
    {
        leading_zero_suppression();
    }

    constexpr Polynomial(std::initializer_list<K> init) : coeffs(init)
    {
        leading_zero_suppression();
    }

    template <std::size_t N>
    constexpr explicit Polynomial(const std::array<K, N>& init) : coeffs(init.begin(), init.end())
    {
        leading_zero_suppression();
    }

    // c * X^k
    static constexpr Polynomial monomial(const K& c, std::size_t k)
    {
        std::vector<K> result(k + 1, K(0));
        result[k] = c;
        return Polynomial(std::move(result));
    }

    static constexpr Polynomial constant(const K& c)
    {
        return monomial(c, 0);
    }

    // the formal variable X
    static constexpr Polynomial x()
    {
        return monomial(K(1), 1);
    }

    // check if it is the zero polynomial
    constexpr bool is_zero() const
    {
        return coeffs.empty();
    }

    // returns the degree (highest exponent) of the polynomial
    // if the polynomial is the 0 polynomial it returns -1 (stands for -infinity, isk 2.5)
    constexpr int64_t degree() const
    {
        return static_cast<int64_t>(coeffs.size()) - 1;
    }

    // coefficient of X^j, also defined beyond the degree
    constexpr K coefficient(std::size_t j) const
    {
        return j < coeffs.size() ? coeffs[j] : K(0);
    }

    constexpr const std::vector<K>& coefficients() const
    {
        return coeffs;
    }

    constexpr K leading_coefficient() const
    {
        return is_zero() ? K(0) : coeffs.back();
    }

    constexpr bool is_monic() const
    {
        return !is_zero() && coeffs.back() == K(1);
    }

    // divide by the leading coefficient
    constexpr Polynomial monic() const
    {
        if (is_zero())
        {
            return *this;
        }
        return *this * (K(1) / coeffs.back());
    }

    // polynomial function x -> f(x) (isk 2.4), evaluated with horner's method
    constexpr K operator()(const K& x) const
    {
        K result(0);
        for (std::size_t i = coeffs.size(); i > 0; i--)
        {
            result = result * x + coeffs[i - 1];
        }
        return result;
    }

    // formal derivative
    constexpr Polynomial derivative() const
    {
        if (coeffs.size() <= 1)
        {
            return Polynomial();
        }
        std::vector<K> result(coeffs.size() - 1);
        for (std::size_t j = 1; j < coeffs.size(); j++)
        {
            result[j - 1] = K(static_cast<int64_t>(j)) * coeffs[j];
        }
        return Polynomial(std::move(result));
    }

    friend constexpr Polynomial operator+(const Polynomial& a, const Polynomial& b)
    {
        std::size_t n = std::max(a.coeffs.size(), b.coeffs.size());
        std::vector<K> result(n);
        for (std::size_t i = 0; i < n; i++)
        {
            result[i] = a.coefficient(i) + b.coefficient(i);
        }
        return Polynomial(std::move(result));
    }

    friend constexpr Polynomial operator-(const Polynomial& a)
    {
        std::vector<K> result(a.coeffs.size());
        for (std::size_t i = 0; i < result.size(); i++)
        {
            result[i] = -a.coeffs[i];
        }
        return Polynomial(std::move(result));
    }

    friend constexpr Polynomial operator-(const Polynomial& a, const Polynomial& b)
    {
        return a + (-b);
    }

    // isk 2.4: the coefficient at X^l is the sum of a_j * b_(l-j)
    friend constexpr Polynomial operator*(const Polynomial& a, const Polynomial& b)
    {
        if (a.is_zero() || b.is_zero())
        {
            return Polynomial();
        }
        std::vector<K> result(a.coeffs.size() + b.coeffs.size() - 1, K(0));
        for (std::size_t i = 0; i < a.coeffs.size(); i++)
        {
            for (std::size_t j = 0; j < b.coeffs.size(); j++)
            {
                result[i + j] = result[i + j] + a.coeffs[i] * b.coeffs[j];
            }
        }
        return Polynomial(std::move(result));
    }

    friend constexpr Polynomial operator*(const Polynomial& a, const K& c)
    {
        std::vector<K> result(a.coeffs.size());
        for (std::size_t i = 0; i < result.size(); i++)
        {
            result[i] = a.coeffs[i] * c;
        }
        return Polynomial(std::move(result));
    }

    friend constexpr Polynomial operator*(const K& c, const Polynomial& a)
    {
        return a * c;
    }

    // polynomial division (isk 2.5): *this = divisor * quotient + remainder
    // with deg(remainder) < deg(divisor)
    constexpr DivisionResult divmod(const Polynomial& divisor) const
    {
        if (divisor.is_zero())
        {
            throw std::domain_error("division by zero");
        }

        // deg(divisor) > deg(*this) -> quotient = 0, remainder = *this
        if (degree() < divisor.degree())
        {
            return {Polynomial(), *this};
        }

        const std::size_t e = static_cast<std::size_t>(divisor.degree());
        const std::size_t n = static_cast<std::size_t>(degree()) - e + 1;
        std::vector<K> quotient(n, K(0));
        std::vector<K> rest = coeffs;

        const K inv_lead = K(1) / divisor.coeffs.back();
        for (std::size_t i = n; i > 0; i--)
        {
            // the leading term of the remainder is at (i - 1) + deg(divisor)
            K q = rest[(i - 1) + e] * inv_lead;
            quotient[i - 1] = q;
            for (std::size_t j = 0; j < e; j++)
            {
                rest[(i - 1) + j] = rest[(i - 1) + j] - q * divisor.coeffs[j];
            }
            // vanishes by construction, set it explicitly so rounding (R) can not leave garbage
            rest[(i - 1) + e] = K(0);
        }
        return {Polynomial(std::move(quotient)), Polynomial(std::move(rest))};
    }

    friend constexpr Polynomial operator/(const Polynomial& a, const Polynomial& b)
    {
        return a.divmod(b).quotient;
    }

    friend constexpr Polynomial operator%(const Polynomial& a, const Polynomial& b)
    {
        return a.divmod(b).remainder;
    }

    constexpr Polynomial& operator+=(const Polynomial& b)
    {
        return *this = *this + b;
    }
    constexpr Polynomial& operator-=(const Polynomial& b)
    {
        return *this = *this - b;
    }
    constexpr Polynomial& operator*=(const Polynomial& b)
    {
        return *this = *this * b;
    }

    friend constexpr bool operator==(const Polynomial&, const Polynomial&) = default;

    // f | g (isk 2.6)
    constexpr bool divides(const Polynomial& g) const
    {
        return (g % *this).is_zero();
    }

    // Returns the monic greatest common divisor of two polynomials (0 if both are 0)
    static constexpr Polynomial gcd(Polynomial a, Polynomial b)
    {
        // Continue while b is not the zero polynomial
        while (!b.is_zero())
        {
            Polynomial temp = b;
            b = a % b;
            a = std::move(temp);
        }
        return a.monic();
    }

    // extended euclidean algorithm in K[X], gcd is monic
    static constexpr ExtendedGcdResult extended_gcd(const Polynomial& a, const Polynomial& b)
    {
        Polynomial r0 = a, r1 = b;
        Polynomial s0 = constant(K(1)), s1;
        Polynomial t0, t1 = constant(K(1));
        while (!r1.is_zero())
        {
            auto [q, r] = r0.divmod(r1);
            r0 = std::move(r1);
            r1 = std::move(r);
            Polynomial s2 = s0 - q * s1;
            s0 = std::move(s1);
            s1 = std::move(s2);
            Polynomial t2 = t0 - q * t1;
            t0 = std::move(t1);
            t1 = std::move(t2);
        }
        if (r0.is_zero())
        {
            return {r0, s0, t0};
        }
        const K inv_lead = K(1) / r0.leading_coefficient();
        return {r0 * inv_lead, s0 * inv_lead, t0 * inv_lead};
    }

    // Computes (base^exp) % mod_poly
    static constexpr Polynomial pow_mod(Polynomial base, uint64_t exp, const Polynomial& mod_poly)
    {
        Polynomial result = constant(K(1)) % mod_poly;
        base = base % mod_poly;

        while (exp > 0)
        {
            if (exp % 2 == 1)
            {
                result = (result * base) % mod_poly;
            }
            base = (base * base) % mod_poly;
            exp /= 2;
        }
        return result;
    }

    // prints it like the script: X^3 + X + 1
    friend std::ostream& operator<<(std::ostream& os, const Polynomial& polynomial)
    {
        std::vector<std::pair<K, std::size_t>> terms;
        for (std::size_t i = polynomial.coeffs.size(); i > 0; i--)
        {
            terms.emplace_back(polynomial.coeffs[i - 1], i - 1);
        }
        detail::print_terms(os, terms, "X");
        return os;
    }
};

} // namespace fields
