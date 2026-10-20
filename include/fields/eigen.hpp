#pragma once

// Optional adapter that makes the exact fields usable as Eigen scalars:
// finite fields (Fp, F4, F9, GF<q>, towers), Q and extensions of Q like Q(√2).
//
// Eigen uses its generic scalar kernels for these types (BLAS/MKL only exist for
// float, double and std::complex). Supported are everything built on +, -, *, / and
// the LU decompositions: products, determinant, inverse, solve, rank, kernel, image.
// Algorithms that need square roots or a real order (QR, SVD, LLT/LDLT, eigenvalues)
// make no sense over these fields.
//
// Pivoting: Eigen picks the pivot with the largest abs(). That only matters for floating
// point stability, over an exact field any nonzero pivot is fine. So abs(x) is 0 or 1
// (|x| for Q), and the ordering defined below is an arbitrary total order with 0 minimal.
// It exists for Eigen only, it is no field ordering (finite fields have none).

#include "fields.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include <optional>

namespace fields
{

namespace detail
{

template <typename K> struct is_exact : std::false_type
{
};

template <uint64_t P> struct is_exact<Fp<P>> : std::true_type
{
};

template <> struct is_exact<Rational> : std::true_type
{
};

template <typename K, auto M, FixedString N>
struct is_exact<Extension<K, M, N>> : is_exact<K>
{
};

// lexicographic, highest coefficient first
template <typename K> constexpr std::strong_ordering eigen_compare(const K& a, const K& b)
{
    if constexpr (std::same_as<K, Rational>)
    {
        return a <=> b;
    }
    else if constexpr (requires { a.value; })
    {
        return a.value <=> b.value;
    }
    else
    {
        for (std::size_t j = K::degree; j > 0; j--)
        {
            const std::strong_ordering c = eigen_compare(a.coeffs[j - 1], b.coeffs[j - 1]);
            if (c != 0)
            {
                return c;
            }
        }
        return std::strong_ordering::equal;
    }
}

} // namespace detail

// fields with exact arithmetic: no rounding, x == y is meaningful
// (the trait is checked first, checking Field on Eigen's own types would recurse into Eigen)
template <typename K>
concept ExactField = detail::is_exact<K>::value && Field<K>;

// the exact fields without a natural order, these get the pivoting order
template <typename K>
concept EigenOrderedField = ExactField<K> && !std::same_as<K, Rational>;

template <EigenOrderedField K> constexpr bool operator<(const K& a, const K& b)
{
    return detail::eigen_compare(a, b) < 0;
}
template <EigenOrderedField K> constexpr bool operator>(const K& a, const K& b)
{
    return b < a;
}
template <EigenOrderedField K> constexpr bool operator<=(const K& a, const K& b)
{
    return !(b < a);
}
template <EigenOrderedField K> constexpr bool operator>=(const K& a, const K& b)
{
    return !(a < b);
}

// math functions Eigen calls (found via ADL)
template <ExactField K> constexpr K abs(const K& x)
{
    if constexpr (std::same_as<K, Rational>)
    {
        return x < K(0) ? -x : x;
    }
    else
    {
        return x == K(0) ? K(0) : K(1);
    }
}
template <ExactField K> constexpr K abs2(const K& x)
{
    return abs(x) * abs(x);
}
template <ExactField K> constexpr K conj(const K& x)
{
    return x;
}
template <ExactField K> constexpr K real(const K& x)
{
    return x;
}
template <ExactField K> constexpr K imag(const K&)
{
    return K(0);
}
template <ExactField K> K sqrt(const K&)
{
    throw std::domain_error("sqrt is not available over " + field_name<K>() +
                            ", use the LU based algorithms");
}

template <ExactField K> using Matrix = Eigen::Matrix<K, Eigen::Dynamic, Eigen::Dynamic>;
template <ExactField K> using Vector = Eigen::Matrix<K, Eigen::Dynamic, 1>;

// A x = b, nullopt if there is no solution; one solution if there are many (see kernel)
template <typename Derived, typename DerivedB>
auto solve(const Eigen::MatrixBase<Derived>& a, const Eigen::MatrixBase<DerivedB>& b)
    -> std::optional<Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, DerivedB::ColsAtCompileTime>>
{
    using K = typename Derived::Scalar;
    const Matrix<K> am = a;
    Eigen::Matrix<K, Eigen::Dynamic, DerivedB::ColsAtCompileTime> x = am.fullPivLu().solve(b);
    // the arithmetic is exact, so the check is too
    if (am * x != b)
    {
        return std::nullopt;
    }
    return x;
}

// A^-1, nullopt if A is singular
template <typename Derived>
auto invert(const Eigen::MatrixBase<Derived>& a) -> std::optional<Matrix<typename Derived::Scalar>>
{
    using K = typename Derived::Scalar;
    if (a.rows() != a.cols())
    {
        return std::nullopt;
    }
    const Eigen::FullPivLU<Matrix<K>> lu(a);
    if (!lu.isInvertible())
    {
        return std::nullopt;
    }
    return Matrix<K>(lu.inverse());
}

template <typename Derived> Eigen::Index rank(const Eigen::MatrixBase<Derived>& a)
{
    return Matrix<typename Derived::Scalar>(a).fullPivLu().rank();
}

// basis of { x : A x = 0 } as columns, no columns if the kernel is 0
template <typename Derived> auto kernel(const Eigen::MatrixBase<Derived>& a)
{
    using K = typename Derived::Scalar;
    const Eigen::FullPivLU<Matrix<K>> lu(a);
    if (lu.rank() == a.cols())
    {
        return Matrix<K>(a.cols(), 0);
    }
    return Matrix<K>(lu.kernel());
}

// prints a matrix with aligned columns (Eigen's operator<< miscounts utf-8 names like ι)
template <typename Derived> void print_matrix(std::ostream& os, const Eigen::MatrixBase<Derived>& m)
{
    std::vector<std::string> cells;
    std::vector<std::size_t> widths(m.cols(), 0);
    for (Eigen::Index r = 0; r < m.rows(); r++)
    {
        for (Eigen::Index c = 0; c < m.cols(); c++)
        {
            cells.push_back(detail::to_string(m(r, c)));
            widths[c] = std::max(widths[c], detail::display_width(cells.back()));
        }
    }
    for (Eigen::Index r = 0; r < m.rows(); r++)
    {
        os << '|';
        for (Eigen::Index c = 0; c < m.cols(); c++)
        {
            const std::string& s = cells[r * m.cols() + c];
            os << std::string(widths[c] - detail::display_width(s) + 1, ' ') << s << ' ';
        }
        os << "|\n";
    }
}

} // namespace fields

namespace Eigen
{

template <fields::ExactField K> struct NumTraits<K> : GenericNumTraits<K>
{
    using Real = K;
    using NonInteger = K;
    using Literal = K;
    using Nested = K;

    enum
    {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = std::same_as<K, fields::Rational>,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 3,
        MulCost = 3,
    };

    // exact arithmetic: pivots are compared against 0, not against a tolerance
    static constexpr K epsilon()
    {
        return K(0);
    }
    static constexpr K dummy_precision()
    {
        return K(0);
    }
    static constexpr int digits10()
    {
        return 0;
    }
    static constexpr int max_digits10()
    {
        return 0;
    }
    static constexpr int digits()
    {
        return 0;
    }
};

} // namespace Eigen
