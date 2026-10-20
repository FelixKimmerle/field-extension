// Inverting matrices and solving systems over finite fields and Q with Eigen

#include <fields/eigen.hpp>

#include <chrono>
#include <iostream>
#include <random>

using namespace fields;

int main()
{
    // an exercise style system over F9
    const F9 i = F9::generator();
    Matrix<F9> a(3, 3);
    a << 1, i, 0, 1 + i, 2, i, 0, 1, 1 - i;
    Vector<F9> b(3);
    b << 1, i, 2;

    std::cout << "A over " << field_name<F9>() << " =\n";
    print_matrix(std::cout, a);
    std::cout << "det(A) = " << a.determinant() << "\n\nA^-1 =\n";
    print_matrix(std::cout, *invert(a));
    std::cout << "\nsolution of A x = (1, ι, 2):\n";
    print_matrix(std::cout, *solve(a, b));

    // singular matrix over F5: rank and kernel
    Matrix<F5> s(3, 3);
    s << 1, 2, 3, 2, 4, 1, 3, 1, 4;
    std::cout << "\nS over F_5 =\n";
    print_matrix(std::cout, s);
    std::cout << "rank(S) = " << rank(s) << ", kernel basis:\n";
    print_matrix(std::cout, kernel(s));

    // exact over Q
    Matrix<Q> h(3, 3);
    h << 1, Q(1, 2), Q(1, 3), Q(1, 2), Q(1, 3), Q(1, 4), Q(1, 3), Q(1, 4), Q(1, 5);
    std::cout << "\nhilbert matrix over Q, inverse:\n";
    print_matrix(std::cout, *invert(h));

    // something no one wants to do by hand
    using K = GF<256>;
    const int n = 200;
    std::mt19937 rng(1);
    std::uniform_int_distribution<uint64_t> dist(0, 255);
    Matrix<K> big(n, n);
    for (int r = 0; r < n; r++)
    {
        for (int c = 0; c < n; c++)
        {
            big(r, c) = K::from_index(dist(rng));
        }
    }
    const auto start = std::chrono::steady_clock::now();
    const auto inv = invert(big);
    const auto ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start);
    std::cout << "\ninverted a random " << n << "x" << n << " matrix over GF(256) in " << ms.count()
              << " ms, A * A^-1 == I: " << std::boolalpha
              << (inv && big * *inv == Matrix<K>::Identity(n, n)) << '\n';
}
