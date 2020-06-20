#pragma once
#include "base_meta_func.hpp"

namespace my {
/*!
 * integrral on simplex 2D
 */
template <size_t alpha, size_t beta>
struct integral_2d
{
    static constexpr size_t value = factorial<alpha + beta + 2>::value / factorial<alpha>::value / factorial<beta>::value;
};
/*!
 * struct Pow
 */
    template <size_t N, class T>
    constexpr T pow(const T& x)
    {
        return N > 1 ? x*pow<(N-1)*(N > 1)>(x)
                     : N < 0 ? T(1)/pow<(-N)*(N < 0)>(x)
                             : N == 1 ? x
                                      : T(1);
    }

/*!
 * struct meta_loop
 */
    template <size_t N, size_t I, class Closure>
    typename boost::enable_if_t<(I == N)> is_meta_loop(Closure& closure) {}

    template <size_t N, size_t I, class Closure>
    typename boost::enable_if_t<(I < N)> is_meta_loop(Closure& closure) {
        closure.template apply<I>();
        is_meta_loop<N, I + 1>(closure);
    }
    template <size_t N, class Closure>
    void meta_loop(Closure& closure) {
        is_meta_loop<N, 0>(closure);
    }
/*!
 * struct abstract_sum
 */
    template <typename Closure>
    struct abstr_sum_closure
    {
        typedef typename Closure::value_type value_type;
        Closure& closure;
        value_type result;
        abstr_sum_closure(Closure& closure) : closure(closure), result(value_type()) {}

        template<size_t I>
        void apply() {
            result += closure.template value<I>();
        }
    };
    template <size_t N, class Closure>
    typename Closure::value_type abstr_sum(Closure& closure) {
        abstr_sum_closure<Closure> my_closure(closure);
        meta_loop<N>(my_closure);
        return my_closure.result;
    }

    template <size_t N, typename T>
    struct binom_closure
    {
        typedef T value_type;
    private:
        const value_type a, b;
    public:
        binom_closure(value_type a, value_type b) : a(a), b(b) {}
        template<size_t K>
        value_type value() const {
            return binom_C<N, K>::value * pow<K>(a) * pow<N - K>(b) / integral_2d<K, N - K>::value;
        }
    };

    template <size_t N, typename T>
    struct binom_s {
        T value(T a, T b) const {
            binom_closure<N, T> closure(a, b);
            return abstr_sum<N + 1>(closure);
        }
    };
} //close namespace my

/*!
 * Реализация последовательности частичных рядов от интеграла Бинома
 */

template<typename T = double>
using abstr_sums_binom = boost::mpl::vector<my::binom_s<0, T>, my::binom_s<1, T>, my::binom_s<2, T>,
                            my::binom_s<3, T>, my::binom_s<4, T>, my::binom_s<5, T>, my::binom_s<6, T>, my::binom_s<7, T>, my::binom_s<8, T>, my::binom_s<20, T>> ;
template<size_t N>
struct abstr_sums {
    BOOST_STATIC_ASSERT((N < 10));
    using type = typename boost::mpl::at_c<abstr_sums_binom<>, N>::type;
};
/*!
 * Метод Галеркина
 */
namespace met_Galerkin {

#define BASE_FUNCTION_COUNT 10

/*!
 * integrral on simplex 3D
 */
template <size_t alpha, size_t beta, size_t gamma>
struct integral_3d
{
    static constexpr size_t value = factorial<alpha + beta + gamma + 3>::value / factorial<alpha>::value / factorial<beta>::value / factorial<gamma>::value;
};

template <size_t N, size_t K, typename T1, typename Closure, typename T2>
typename boost::enable_if_t<(K > 0), T1> partial_sum(const T2& args) {
    return partial_sum<N, K - 1, T1, Closure, T2>(args) + Closure::template value<N, K>(args);
}

template <size_t N, size_t K, typename T1, typename Closure, typename T2>
typename boost::enable_if_t<(K == 0), T1> partial_sum(const T2& args) {
    return Closure::template value<N, K>(args);
}

template <size_t ALPHA, size_t BETA, size_t GAMMA>
struct base_func
{
    static constexpr size_t alpha = ALPHA;
    static constexpr size_t beta  = BETA;
    static constexpr size_t gamma = GAMMA;
};
struct bf_1 : base_func<0, 0, 0>{};
struct bf_x : base_func<1, 0, 0>{};
struct bf_y : base_func<0, 1, 0>{};
struct bf_z : base_func<0, 0, 1>{};
struct bf_x2 : base_func<2, 0, 0>{};
struct bf_y2 : base_func<0, 2, 0>{};
struct bf_z2 : base_func<0, 0, 2>{};
struct bf_xy : base_func<1, 1, 0>{};
struct bf_xz : base_func<1, 0, 1>{};
struct bf_yz : base_func<0, 1, 1>{};

typedef boost::mpl::vector<bf_1, bf_x, bf_y, bf_z, bf_x2, bf_y2, bf_z2, bf_xy, bf_xz, bf_yz> bf_t;

template <typename T>
struct abc_struct
{
private:
    const T *a, *b, *c;
public:
    abc_struct(const T* a_, const T* b_, const T* c_) : a(a_), b(b_), c(c_) {}

    template<size_t step>
    typename boost::enable_if_t<step == 1, const T*> get() const { return a; }

    template<size_t step>
    typename boost::enable_if_t<step == 2, const T*> get() const { return b; }

    template<size_t step>
    typename boost::enable_if_t<step == 3, const T*> get() const { return c; }
};
///three inserts summ in integral_3D

template <size_t step, typename next_step, size_t PX, size_t PY, size_t PZ>
struct simplex_integral_sum3
{
    template <size_t J, size_t K, typename T>
    static T value(const abc_struct<T>& args) {
        return binom_C<J, K>::value * my::pow<K>(args.template get<step>()[0]) *
                my::pow<J - K>(args.template get<step>()[1]) * next_step::template value<PX + K, PY + J - K, PZ>(args);
    }
};
template <size_t step, typename next_step, size_t PX, size_t PY, size_t PZ>
struct simplex_integral_sum2
{
    template <size_t I, size_t J, typename T>
    static T value(const abc_struct<T>& args) {
        return binom_C<I, J>::value * my::pow<I - J>(args.template get<step>()[2]) *
                partial_sum<J, J, T, simplex_integral_sum3<step, next_step, PX, PY, PZ + I - J>>(args);
    }
};
template <size_t step, typename next_step, size_t PX, size_t PY, size_t PZ>
struct simplex_integral_sum1
{
    template <size_t N, size_t I, typename T>
    static T value(const abc_struct<T>& args) {
        return binom_C<N, I>::value * my::pow<N - I>(args.template get<step>()[3]) *
                partial_sum<I, I, T, simplex_integral_sum2<step, next_step, PX, PY, PZ>>(args);
    }
};
struct simplex_integral
{
    template <size_t PX, size_t PY, size_t PZ, typename T>
    static T value(const abc_struct<T>& ) { return T(1) / integral_3d<PX, PY, PZ>::value; }
};
template <size_t GAMMA>
struct simplex_integral_step_3
{
    template <size_t PX, size_t PY, size_t PZ, typename T>
    static T value(const abc_struct<T>& args) {
        return partial_sum<GAMMA, GAMMA, T, simplex_integral_sum1<3, simplex_integral, PX, PY, PZ>>(args);
    }
};
template <size_t BETA, size_t GAMMA>
struct simplex_integral_step_2
{
    template <size_t PX, size_t PY, size_t PZ, typename T>
    static T value(const abc_struct<T>& args) {
        return partial_sum<BETA, BETA, T, simplex_integral_sum1<3, simplex_integral_step_3<GAMMA>, PX, PY, PZ>>(args);
    }
};
template<size_t ALPHA, size_t BETA, size_t GAMMA, typename T>
T simplex_integral_polinom3(const T a[4], const T b[4], const T c[4]) {
    return partial_sum<ALPHA, ALPHA, T, simplex_integral_sum1<1, simplex_integral_step_2<BETA, GAMMA>, 0, 0, 0>>(abc_struct<T>(a, b, c));
}

template <size_t i, typename Matrix, typename Vector>
struct CALC_A
{
private:
    Matrix& A;
    const Vector& dx;
    const double J, *a, *b, *c;
public:
    typedef typename boost::mpl::at_c<bf_t, i>::type type_i;

    CALC_A(Matrix& A, double J, const double* a, const double* b, const double* c, const Vector& dx):
            A(A), J(J), a(a), b(b), c(c), dx(dx) {}

    template <size_t j>
    void apply() {
        typedef typename boost::mpl::at_c<bf_t, j>::type type_j;
        constexpr size_t alpha = type_i::alpha + type_j::alpha;
        constexpr size_t beta  = type_i::beta + type_j::beta;
        constexpr size_t gamma = type_i::gamma + type_j::gamma;
        A(i, j) = J * simplex_integral_polinom3<alpha, beta, gamma>(a, b, c) / my::pow<alpha>(dx.value_x()) / my::pow<beta>(dx.value_y()) / my::pow<gamma>(dx.value_z());
    }
};

template <typename Matrix, typename Vector>
struct CALC_A_i
{
private:
    Matrix& A;
    const Vector& dx;
    const double J, *a, *b, *c;
public:
    CALC_A_i(Matrix& A, double J, const double* a, const double* b, const double* c, const Vector& dx):
            A(A), J(J), a(a), b(b), c(c), dx(dx) {}

    template <size_t i>
    void apply() {
        CALC_A<i, Matrix, Vector> closure(A, J, a, b, c, dx);
        my::meta_loop<i + 1>(closure);
    }
};
template <typename Matrix, typename Vector>
void calc_matrix(Matrix& A, double J, const double* a, const double* b, const double* c, const Vector& dx) {
    CALC_A_i<Matrix, Vector> closure(A, J, a, b, c, dx);
    my::meta_loop<BASE_FUNCTION_COUNT>(closure);
}
}//namespace Метод Галеркина

