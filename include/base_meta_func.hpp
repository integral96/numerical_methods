#pragma once

#include <iostream>
#include <boost/mpl/if.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/foreach.hpp>


/*!
 * struct Factorial
 */
template <size_t N>
struct factorial {
    static constexpr size_t value = N*(factorial<N - 1>::value);
};
template <>
struct factorial<0> : boost::mpl::int_<1> {};

/*!
 * struct Vector
 */

template<size_t... N> struct int_vector {};
template<typename Int_vector, size_t i> struct push_back;

template<size_t... N, size_t i> struct push_back<int_vector<N...>, i> {
    typedef int_vector<N..., i> type;
};
template<typename Int_vector, size_t j> struct at;

template<size_t i, size_t... N, size_t j> struct at<int_vector<i, N...>, j> : at<int_vector<N...>, j - 1> {};
template<size_t i, size_t... N> struct at<int_vector<i, N...>, 0> : boost::mpl::integral_c<size_t, i> {};


/*!
 * struct Binom C
 */
template <size_t N, size_t K>
struct binom_C
{
    static constexpr size_t value = factorial<N>::value / factorial<K>::value /factorial<N - K>::value;
};

/*!
 * struct Pow
 */
    template <size_t N, class T>
    constexpr T pow_(const T& x)
    {
        return N > 1 ? x*pow_<(N-1)*(N > 1)>(x)
                     : N < 0 ? T(1)/pow_<(-N)*(N < 0)>(x)
                             : N == 1 ? x
                                      : T(1);
    }

template<typename O>
struct math_obj_base {
    O& self() {
        return static_cast<O&>(*this);
    }
    const O& self() const {
        return static_cast<const O&>(*this);
    }
};

template<typename T> struct expression : math_obj_base<T> {};


