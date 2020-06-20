#pragma once
#include "base_meta_func.hpp"
///Simpson Rule

namespace Simpson {
/*!
 * struct meta_loop
 */
    template <size_t N, size_t I, class Closure>
    typename boost::enable_if_t<(I == N)> is_meta_loop(Closure& closure) {}

    template <size_t N, size_t I, class Closure>
    typename boost::enable_if_t<(I < N)> is_meta_loop(Closure& closure) {
        closure.template apply<I>();
        is_meta_loop<N + 1, I + 2>(closure);
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
/*!
 *
 */
template <size_t N, typename Func>
struct integral_Simpson2
{
public:
    typedef typename Func::value_type value_type;
    value_type h;
private:
    Func& func;
    value_type x0, xn;
public:
    integral_Simpson2(Func& func, const double x0, const double xn) : func(func), x0(x0), xn(xn) {
        h = (xn - x0) / (N);
    }
    template<size_t j>
    value_type value() {
        return 4 * func(x0 + j*h) + 2 * func(x0 + (j - 1)*h);
    }
    template<>
    value_type value<0>() {
        return func(x0) + func(xn) + 4 * func(x0 + h);
    }
    template<>
    value_type value<1>() {
        return value_type(0);
    }
    template<>
    value_type value<2>() {
        return value_type(0);
    }
    template<>
    value_type value<N + 1>() {
        return value_type(0);
    }
};

template <size_t N, typename Func>
struct Simpson2
{
    typedef typename Func::value_type value_type;
    value_type value(Func& func, value_type a, value_type b) const {
        integral_Simpson2<N, Func> closure(func, a, b);
        return abstr_sum<N>(closure) * closure.h / 3;
    }
};
template <typename T>
struct Func
{
    typedef T value_type;
    T operator()(T x) const {
        return 1/(1 + x*x);
    }
};
}
