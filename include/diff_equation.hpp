#pragma once

#include "base_meta_func.hpp"

#include <boost/array.hpp>
#include <boost/mpl/divides.hpp>
#include <boost/mpl/multiplies.hpp>
#include <boost/mpl/minus.hpp>
#include <boost/type_traits/enable_if.hpp>

template<size_t N>
struct int_constant : expression<int_constant<N>> 
{
    static constexpr size_t value = N;
    typedef int_constant<0> diff_type;
    typedef int_constant<0> type;

    diff_type diff() const {
        return diff_type();
    }

    template<typename T>
    size_t operator() (const T& x) const {
        return value;
    }
};

template<typename T> 
struct is_int_constant : boost::mpl::false_ {};

template<size_t N> 
struct is_int_constant<int_constant<N>> : boost::mpl::true_ {};

template<typename T> 
struct is_constant_value : boost::mpl::integral_c<size_t, 0> {};

template<size_t N> 
struct is_constant_value<int_constant<N>> : boost::mpl::integral_c<size_t, N> {};

template<typename T>
struct  scalar : expression<scalar<T>> {
    typedef T value_type;
    typedef int_constant<0> diff_type;
    const value_type value;

    scalar(const value_type& value) : value(value) {}
    diff_type diff() const {
        return diff_type();
    }
    template<typename E>
    value_type operator() (const E& x) const {
        return value;
    }
};

template<typename E>
struct is_scalar : boost::mpl::false_ {};

template<typename T>
struct is_scalar<scalar<T>> : boost::mpl::true_ {};

template<typename T>
scalar<T> _ (const T& val) {
    return scalar<T>(val);
}


template<typename E> struct negate_expression;

template<typename E> 
struct negate_expression_type {
    typedef typename boost::mpl::if_c<is_int_constant<E>::value, int_constant<- is_constant_value<E>::value>, 
                            typename boost::mpl::if_c<is_scalar<E>::value, E, negate_expression<E>>::type>::type type;
};

template<typename T>
struct negate_expression : expression<negate_expression<T>> {
    typedef typename negate_expression_type<typename T::diff_type>::type diff_type;

    const T& e;
    negate_expression(const expression<T>& e) : e(e.self()) {}

    diff_type diff() const {
        return -e.diff();
    }

    template<typename E>
    E operator() (const E& x) const {
        return -e(x);
    }
}; 

template<typename T>
negate_expression<T> operator - (const expression<T>& e ) {
    return negate_expression<T>(e);
}
template<size_t N>
int_constant<- N> operator - (const int_constant<N>& ) {
    return int_constant<- N>();
}
template<typename T>
scalar<T> operator - (const scalar<T>& e ) {
    return scalar<T>(-e.value);
}
template<typename T1, char OP, typename T2> struct operator_expression;
template<typename T1, char OP, typename T2>
using operator_type_plus_ = typename boost::mpl::if_c<is_int_constant<T1>::value || is_int_constant<T2>::value,
                                                int_constant<is_constant_value<T1>::value + is_constant_value<T2>::value>,
                            typename boost::mpl::if_c<is_scalar<T1>::value || is_scalar<T2>::value,
                                                boost::is_same<T1, T2>, operator_expression<T1, OP, T2>>::type>::type;
template<typename T1, char OP, typename T2>
using operator_type_minus_ = typename boost::mpl::if_c<is_int_constant<T1>::value || is_int_constant<T2>::value,
                                                int_constant<is_constant_value<T1>::value - is_constant_value<T2>::value>,
                            typename boost::mpl::if_c<is_scalar<T1>::value || is_scalar<T2>::value,
                                                boost::is_same<T1, T2>, operator_expression<T1, OP, T2>>::type>::type;

template<typename T1, char OP, typename T2>
struct operator_expression_type {
    typedef typename boost::mpl::if_c<OP == '+', operator_type_plus_<T1, OP, T2>,
                     boost::enable_if_<OP == '-', operator_type_minus_<T1, OP, T2>>>::type type;
};

template<char OP, typename T1, typename T2>
typename boost::enable_if_t<OP == '+', typename operator_expression_type<typename T1::diff_type, OP, typename T2::diff_type>::type>
operator_expression_diff(const T1& x1, const T2& x2) {
    return x1.diff() + x2.diff();
}
template<char OP, typename T1, typename T2>
typename boost::enable_if_t<OP == '-', typename operator_expression_type<typename T1::diff_type, OP, typename T2::diff_type>::type>
operator_expression_diff(const T1& x1, const T2& x2) {
    return x1.diff() - x2.diff();
}
template<typename T1, char OP, typename T2>
struct operator_expression : expression<operator_expression<T1, OP, T2>> {
    typedef typename operator_expression_type<typename T1::diff_type, OP, typename T2::diff_type>::type diff_type;

    const T1 x1;
    const T2 x2;

    operator_expression(const expression<T1>& x1, const expression<T2>& x2) : x1(x1.self()), x2(x2.self()) {}

    diff_type diff() const {
        return operator_expression_diff<OP>(x1, x2);
    }
    template<typename T>
    typename boost::enable_if_t<OP == '+', T> operator() (const T&x) const {
        return x1(x) + x2(x);
    }
    template<typename T>
    typename boost::enable_if_t<OP == '-', T> operator() (const T&x) const {
        return x1(x) - x2(x);
    }
};

template<typename T1, typename T2>
operator_expression<T1, '+', T2> operator + (const expression<T1>& x1, const expression<T2>& x2) {
    return operator_expression<T1, '+', T2>(x1, x2);
}
template<typename T>
const T& operator + (const expression<T>& x, const int_constant<0>&) {
    return x.self();
}
template<typename T>
const T& operator + (const int_constant<0>&, const expression<T>& x) {
    return x.self();
}

template<typename T1, typename T2>
operator_expression<T1, '-', T2> operator - (const expression<T1>& x1, const expression<T2>& x2) {
    return operator_expression<T1, '-', T2>(x1, x2);
}
template<typename T>
const T& operator - (const expression<T>& x, const int_constant<0>&) {
    return x.self();
}
template<typename T>
const T& operator - (const int_constant<0>&, const expression<T>& x) {
    return x.self();
}
template<size_t N1, size_t N2>
int_constant<N1 + N2> operator + (const int_constant<N1>&, const int_constant<N2>&) {
    return int_constant<N1 + N2>();
}
template<size_t N1>
int_constant<N1> operator + (const int_constant<N1>&, const int_constant<0>&) {
    return int_constant<N1>();
}
template<size_t N1>
int_constant<N1> operator + (const int_constant<0>&, const int_constant<N1>&) {
    return int_constant<N1>();
}
template<size_t N1, size_t N2>
int_constant<N1 - N2> operator - (const int_constant<N1>&, const int_constant<N2>&) {
    return int_constant<N1 - N2>();
}
template<size_t N1>
int_constant<N1> operator - (const int_constant<N1>&, const int_constant<0>&) {
    return int_constant<N1>();
}
template<size_t N1>
int_constant<N1> operator - (const int_constant<0>&, const int_constant<N1>&) {
    return int_constant<N1>();
}

///mult, divide expressions and operators
///
template<typename T1, char OP, typename T2> struct mul_div_expression;
template<typename T1, char OP, typename T2>
using operator_type_mul_ = typename boost::mpl::if_c<(is_int_constant<T1>::value || is_int_constant<T2>::value),
                                                int_constant<is_constant_value<T1>::value * is_constant_value<T2>::value>,
                            typename boost::mpl::if_c<(is_scalar<T1>::value || is_scalar<T2>::value),
                                                boost::is_same<T1, T2>, mul_div_expression<T1, OP, T2>>::type>::type;
template<typename T1, char OP, typename T2>
using operator_type_div_ = typename boost::mpl::if_c<(is_int_constant<T1>::value || is_int_constant<T2>::value),
                                                int_constant<boost::mpl::divides<typename is_constant_value<T1>::type, typename is_constant_value<T2>::type>::value>,
                            typename boost::mpl::if_c<(is_scalar<T1>::value || is_scalar<T2>::value),
                                                boost::is_same<T1, T2>, mul_div_expression<T1, OP, T2>>::type>::type;

template<typename T1, char OP, typename T2>
struct mul_div_expression_type {
    typedef typename boost::mpl::if_c<OP == '*', operator_type_mul_<T1, OP, T2>,
                     boost::mpl::if_c<OP == '/', operator_type_div_<T1, OP, T2>, boost::mpl::na>>::type type;
};

template<char OP, typename T1, typename T2>
typename boost::enable_if_t<OP == '*', typename mul_div_expression_type<typename T1::diff_type, OP, typename T2::diff_type>::type>
mul_div_expression_diff(const T1& x1, const T2& x2) {
    return x1.diff() * x2.diff();
}
template<char OP, typename T1, typename T2>
typename boost::enable_if_t<OP == '/', typename mul_div_expression_type<typename T1::diff_type, OP, typename T2::diff_type>::type>
mul_div_expression_diff(const T1& x1, const T2& x2) {
    return /*(x1.diff() * x2 - x1 * x2.diff())/(x2*x2)*/x1.diff()/x2.diff();
}

template<typename T1, char OP, typename T2>
struct mul_div_expression : expression<mul_div_expression<T1, OP, T2>> {
    typedef typename mul_div_expression_type<typename T1::diff_type, OP, typename T2::diff_type>::type diff_type;

    const T1 x1;
    const T2 x2;

    mul_div_expression(const expression<T1>& x1, const expression<T2>& x2) : x1(x1.self()), x2(x2.self()) {}

    diff_type diff() const {
        return mul_div_expression_diff<OP>(x1, x2);
    }
    template<typename T>
    typename boost::enable_if_t<OP == '*', T> operator() (const T&x) const {
        return x1(x) * x2(x);
    }
    template<typename T>
    typename boost::enable_if_t<OP == '/', T> operator() (const T&x) const {
        return x1(x) / x2(x);
    }
};

template<typename T1, typename T2>
mul_div_expression<T1, '*', T2> operator * (const expression<T1>& x1, const expression<T2>& x2) {
    return mul_div_expression<T1, '*', T2>(x1, x2);
}
template<typename T>
const T& operator * (const expression<T>& x, const int_constant<0>&) {
    return T(0);
}
template<typename T>
const T& operator * (const int_constant<0>&, const expression<T>& x) {
    return T(0);
}
template<typename T1, typename T2>
mul_div_expression<T1, '/', T2> operator / (const expression<T1>& x1, const expression<T2>& x2) {
    return mul_div_expression<T1, '/', T2>(x1, x2);
}


template<size_t N1, size_t N2>
int_constant<N1 * N2> operator * (const int_constant<N1>&, const int_constant<N2>&) {
    return int_constant<N1 * N2>();
}

template<size_t N1, size_t N2>
int_constant<N1 / N2> operator / (const int_constant<N1>&, const int_constant<N2>&) {
    return int_constant<N1 / N2>();
}

///math_function
///
/*!
  * Sinus
  */

template<typename T> struct cos_expression;

template<typename T>
struct sin_expression : expression<sin_expression<T>>
{
    typedef typename mul_div_expression_type<cos_expression<T>, '*', typename T::diff_type>::type diff_type;
    const T f;

    sin_expression(const expression<T>& x) : f(x.self()) {}
    diff_type diff() const {
        return cos_(f) * f.diff();
    }
    template<typename T1>
    T1 operator() (const T1&x) const {
        return std::sin(f(x));
    }
};
template<typename T>
sin_expression<T> sin_(const expression<T>& x) {
    return sin_expression<T>(x);
}
/*!
  * Cosinus
  */
template<typename T>
struct cos_expression : expression<cos_expression<T>>
{
    typedef typename mul_div_expression_type<sin_expression<T>, '*', typename T::diff_type>::type diff_type;
    const T f;

    cos_expression(const expression<T>& x) : f(x.self()) {}
    diff_type diff() const {
        return sin(f) * f.diff();
    }
    template<typename T1>
    T1 operator() (const T1&x) const {
        return std::cos(f(x));
    }
};
template<typename T>
cos_expression<T> cos_(const expression<T>& x) {
    return cos_expression<T>(x);
}

/*!
  * Log naturaln
  */

template<typename T>
struct log_expression : expression<log_expression<T>>
{
    typedef typename mul_div_expression_type<typename T::diff_type, '/', T>::type diff_type;
    const T f;

    log_expression(const expression<T>& x) : f(x.self()) {}
    diff_type diff() const {
        return f.diff() / f;
    }
    template<typename T1>
    T1 operator() (const T1&x) const {
        return std::log(f(x));
    }
};
template<typename T>
log_expression<T> log_(const expression<T>& x) {
    return log_expression<T>(x);
}

/*!
 *
 * Variable
 */
template <size_t N>
struct variable : expression<variable<N>> {
    template <size_t M>
    struct diff_type
    {
        using type = int_constant<M == N>;
    };

    template <size_t M>
    typename diff_type<M>::type diff() const {
        return diff_type<M>::type();
    }
    template<typename T, size_t SIZE>
    T operator() (const boost::array<T, SIZE> &arr) const {
        return arr[N];
    }
};
template <>
struct variable<0> : expression<variable<0>> {
    typedef int_constant<1> diff_type;
    diff_type diff() const {
        return diff_type();
    }
    template<typename T>
    T operator() (const T& x) const {
        return x;
    }
};

template<typename T, typename F = double>
F metod_newton(const expression<T> & expr, F x, size_t MAXITER, F eps) {
    const T &func = expr.self();
    const typename T::diff_type fd = func.diff();
    size_t N_iter = 0;
    F func_x = x*x;
    while (std::abs(func_x - func(x)) < eps) {
        x -= func_x / fd(x);
        if(++N_iter > MAXITER) throw std::out_of_range("OUT N_iter > MAXITER");
    }
    return x;
}
//Осталось опеределить типы
