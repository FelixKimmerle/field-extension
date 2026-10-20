#include <fields/eigen.hpp>

#include <iostream>
#include <random>

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

static_assert(ExactField<F9> && ExactField<Q> && ExactField<GF<256>>);
static_assert(!ExactField<Real> && !ExactField<Complex>);

template <FiniteField K> Matrix<K> random_matrix(Eigen::Index rows, Eigen::Index cols, std::mt19937& rng)
{
    std::uniform_int_distribution<uint64_t> dist(0, FieldTraits<K>::cardinality - 1);
    Matrix<K> m(rows, cols);
    for (Eigen::Index r = 0; r < rows; r++)
    {
        for (Eigen::Index c = 0; c < cols; c++)
        {
            m(r, c) = K::from_index(dist(rng));
        }
    }
    return m;
}

// invert, solve, rank and kernel of random matrices, including singular ones
template <FiniteField K> void check_random(Eigen::Index n, std::mt19937& rng)
{
    const Matrix<K> id = Matrix<K>::Identity(n, n);
    for (int round = 0; round < 5; round++)
    {
        const Matrix<K> a = random_matrix<K>(n, n, rng);
        const auto inv = invert(a);
        CHECK(inv.has_value() == (a.determinant() != K(0)));
        if (inv)
        {
            CHECK(a * *inv == id);
            CHECK(*inv * a == id);
            CHECK(rank(a) == n);
            const Vector<K> b = random_matrix<K>(n, 1, rng);
            const auto x = solve(a, b);
            CHECK(x.has_value() && a * *x == b);
        }

        // rank deficient: product of n x (n - 2) and (n - 2) x n
        const Matrix<K> s = random_matrix<K>(n, n - 2, rng) * random_matrix<K>(n - 2, n, rng);
        CHECK(s.determinant() == K(0));
        CHECK(!invert(s).has_value());
        const Matrix<K> k = kernel(s);
        CHECK(rank(s) + k.cols() == n);
        CHECK(s * k == Matrix<K>::Zero(n, k.cols()));
    }
}

void test_kuenzer_style()
{
    // systems over F5
    Matrix<F5> a(3, 3);
    a << 1, 2, 3, 4, 0, 1, 2, 2, 1;
    Vector<F5> b(3);
    b << 1, 2, 3;
    const auto x = solve(a, b);
    CHECK(x.has_value() && a * *x == b);
    // det = 10 = 0 in F5, and b is not in the image (checked by brute force over all 125 x)
    a(2, 2) = 2;
    CHECK(a.determinant() == F5(0));
    CHECK(!solve(a, b).has_value());
    // but a consistent right hand side has a solution plus the kernel
    const Vector<F5> consistent = a * Vector<F5>::Ones(3);
    const auto y = solve(a, consistent);
    CHECK(y.has_value() && a * *y == consistent);
    CHECK(kernel(a).cols() == 1);

    // over F9 with ι
    const F9 i = F9::generator();
    Matrix<F9> m(3, 3);
    m << 1, i, 0, 1 + i, 2, i, 0, 1, 1 - i;
    CHECK(m.determinant() == 2 + i);
    const auto mi = invert(m);
    CHECK(mi.has_value() && m * *mi == Matrix<F9>::Identity(3, 3));

    // no solution: second row is a multiple of the first, b is not
    Matrix<F9> s(2, 2);
    s << 1, i, i, -1;
    Vector<F9> c(2);
    c << 1, 1;
    CHECK(!solve(s, c).has_value());
    CHECK(rank(s) == 1);

    // printing is unaffected by the pivoting order
    CHECK(detail::to_string(2 - i) == "2 + 2ι");
}

void test_rationals()
{
    // hilbert matrix, its inverse has integer entries
    const int n = 5;
    Matrix<Q> h(n, n);
    for (int r = 0; r < n; r++)
    {
        for (int c = 0; c < n; c++)
        {
            h(r, c) = Q(1, r + c + 1);
        }
    }
    const auto hi = invert(h);
    CHECK(hi.has_value());
    CHECK(*hi * h == Matrix<Q>::Identity(n, n));
    CHECK((*hi)(0, 0) == Q(25));
    CHECK((*hi)(4, 4) == Q(44100));
    CHECK(h.determinant() == Q(1, 266716800000));

    // Q(√2)
    using Q2 = Extension<Q, poly<Q>(-2, 0, 1), "√2">;
    const Q2 w = Q2::generator();
    Matrix<Q2> m(2, 2);
    m << 1, w, w, 3;
    CHECK(m.determinant() == Q2(1));
    CHECK(invert(m).has_value() && m * *invert(m) == Matrix<Q2>::Identity(2, 2));
}

int main()
{
    std::mt19937 rng(42);

    test_kuenzer_style();
    test_rationals();

    check_random<F2>(10, rng);
    check_random<F3>(7, rng);
    check_random<F4>(8, rng);
    check_random<F9>(6, rng);
    check_random<GF<256>>(40, rng); // large enough for Eigen's blocked LU
    check_random<Extension<F4, poly<F4>(F4::generator(), 1, 1), "γ">>(12, rng);

    // fixed size matrices take other code paths (closed formulas for det/inverse)
    Eigen::Matrix<F7, 4, 4> f = Eigen::Matrix<F7, 4, 4>::Identity();
    f(0, 3) = 3;
    f(2, 1) = 5;
    CHECK(f.determinant() == F7(1));
    CHECK(f * f.inverse() == (Eigen::Matrix<F7, 4, 4>::Identity()));

    if (failures == 0)
    {
        std::cout << "all eigen tests passed\n";
        return 0;
    }
    std::cout << failures << " checks failed\n";
    return 1;
}
