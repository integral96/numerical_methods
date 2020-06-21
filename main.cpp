#include <iostream>
#include <iomanip>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/spirit/include/karma.hpp>


#include "include/simply_meta.hpp"
#include "include/integrate_meta.hpp"
#include "include/diff_equation.hpp"
#include "include/Runge_Kutta.hpp"

namespace ubls = boost::numeric::ublas;
namespace krm  = boost::spirit::karma;

static Simpson::Func<double> func;
static size_t i{3};
template<size_t N>
using solver = Simpson::Simpson2<N, decltype(func)>;

typedef boost::mpl::vector<solver<3>, solver<6>, solver<12>, solver<24>, solver<48>, solver<96>, solver<192>,
                           solver<384>, solver<768>, solver<1536>, solver<3072>, solver<6144>> vector_Sympson;

struct print_partial_sum
{
    print_partial_sum(std::ostream& s) : out(s) {}
    template<typename T>
    void operator()(boost::mpl::identity<T>) const {
        static T a;
        out <<"Partial sum<" << i << "> = " << a.value(func, 0, 6) << "\tx0 = .0, xn = 6."<< std::endl;
        i *= 2;
    }
private:
    std::ostream& out;
};
////Классический подход
template<size_t N, typename Func>
struct classic_Sympson
{
    typedef typename Func::value_type value_type;
    value_type value(Func& func, value_type x0, value_type xn) const {
        value_type h = (xn - x0)/N;
        value_type s = func(x0) + func(xn) + 4*func(x0 + h);
        for(size_t i = 3; i <= N - 1; i += 2) s += 4*func(x0 + i*h) + 2*func(x0 + (i - 1)*h);
        return s * h / 3;
    }
};
template<size_t N>
using solver_class = classic_Sympson<N, decltype(func)>;
typedef boost::mpl::vector<solver_class<3>, solver_class<6>, solver_class<12>, solver_class<24>, solver_class<48>, solver_class<96>, solver_class<192>,
                           solver_class<384>, solver_class<768>, solver_class<1536>, solver_class<3072>, solver_class<6144>> vector_Sympson_class;


void binom_partial_sum(double a, double b) {
     std::cout << "Integral_2D from Binom. " << std::endl;
     abstr_sums<0>::type ass;
     abstr_sums<1>::type ass1;
     std::cout << std::setw(6) << "partial sum<1>: " << ass.value(a, b) << "\t\tpartial sum<2>: " << ass1.value(a, b) << std::endl;
     abstr_sums<2>::type ass2;
     abstr_sums<3>::type ass3;
     std::cout << std::setw(6)<< "partial sum<3>: " << ass2.value(a, b) << "\t\tpartial sum<4>: " << ass3.value(a, b) << std::endl;
     abstr_sums<4>::type ass4;
     abstr_sums<5>::type ass5;
     std::cout << std::setw(6)<< "partial sum<5>: " << ass4.value(a, b) << "\t\tpartial sum<6>: " << ass5.value(a, b) << std::endl;
     abstr_sums<6>::type ass6;
     abstr_sums<7>::type ass7;
     std::cout << std::setw(6)<< "partial sum<7>: " << ass6.value(a, b) << "\t\tpartial sum<8>: " << ass7.value(a, b) << std::endl;
     abstr_sums<8>::type ass8;
     abstr_sums<9>::type ass9;
     std::cout << std::setw(6) << "partial sum<9>: " << ass8.value(a, b) << "\t\tpartial sum<20>: " << ass9.value(a, b) << std::endl;
     std::cout << " ============================================= " << std::endl;
}

template <typename T = double>
struct Vec_simplix
{
private:
    T x, y, z;
public:
    Vec_simplix(T x, T y, T z) : x(x), y(y), z(z) {}
    T value_x() const {
        return x;
    }
    T value_y() const {
        return y;
    }
    T value_z() const {
        return z;
    }
};

int main(int argc, char **argv)
{
     std::cout << "Двумерный Интеграл от Бинома(Мета) f(x, y) = (ax + by)^N: " << std::endl;

     double a = -5.5;
     double b = -5.6;
     std::cout << "if (a < -1 && b < - 1) : " << std::endl;
     std::cout << "Последовательность частичных сумм в этих точках:" << std::endl;
     binom_partial_sum(a, b);
     double a1 = -0.5;
     double b1 = -0.6;
     std::cout << "if (a < 0 && b < 0) && (a > -1 && b > - 1) : " << std::endl;
     std::cout << "Последовательность частичных сумм в этих точках:" << std::endl;
     binom_partial_sum(a1, b1);
     double a2 = 0.5;
     double b2 = 0.6;
     std::cout << "if (a > 0 && b > 0) && (a < 1 && b < 1) : " << std::endl;
     std::cout << "Последовательность частичных сумм в этих точках:" << std::endl;
     binom_partial_sum(a2, b2);
     double a3 = 2.5;
     double b3 = 2.6;
     std::cout << "if (a > 1 && b > 1) : " << std::endl;
     std::cout << "Последовательность частичных сумм в этих точках:" << std::endl;
     binom_partial_sum(a3, b3);
     double a4 = -5.5;
     double b4 = 5.6;
     std::cout << "if (a < -1 && b > 1) : " << std::endl;
     std::cout << "Последовательность частичных сумм в этих точках:" << std::endl;
     binom_partial_sum(a4, b4);

     std::cout << "трёхмерные тетраэдральныесет-ки, разрывный метод Галёркина" << std::endl;
     std::cout << "Размер тетраэдра dx(1.2, 2.1, 1.3), Якобиан = 1,25 Получим матрицу масс: " << std::endl;
     ubls::matrix<double> A(10, 10);
     double aa[10], ba[10], ca[10];
     for(size_t i = 0; i < 10; ++i) {
         aa[i] = ba[i] = ca[i] = (double)i/2.;
     }
     Vec_simplix dx(1.2, 2.1, 1.3);
     met_Galerkin::calc_matrix<ubls::matrix<double>, Vec_simplix<double>>(A, 1.25, aa, ba, ca, dx);
     std::cout << krm::format_delimited(krm::columns(A.size2()) [krm::auto_], '\t', A.data()) << std::endl;
     std::cout << "==========================================================" << std::endl;
    std::cout << "Интеграл Симпсона(Мета) f(x) = 1/(1 + x*x): " << std::endl;
    std::cout << "Последовательность частичных сумм:" << std::endl;
    boost::mpl::for_each<vector_Sympson, boost::mpl::make_identity<boost::mpl::_1> >(print_partial_sum(std::cout));
    i = 3;
    std::cout << "Интеграл Симпсона(Классический) f(x) = 1/(1 + x*x): " << std::endl;
    std::cout << "Последовательность частичных сумм:" << std::endl;
    boost::mpl::for_each<vector_Sympson_class, boost::mpl::make_identity<boost::mpl::_1> >(print_partial_sum(std::cout));
    
    std::cout << std::boolalpha << int_constant<12>::value << " " <<  is_int_constant<int_constant<12>>::value <<
        " " << is_constant_value<int_constant<12>>::value << std::endl;
    variable<0> x;
    double result = metod_newton(x*x*x+x, .0058, 200, 1e-12);
    // Нужно определить функциональные типы
    std::cout << result << std::endl;
    std::cout << "Ну и напоследок самая известная траектория системы\n "
                 << "дифференциальныx уравнений Лоренца(Метод Рунге Кутта 4-го порядка): " << std::endl;
    const char* win {"GRAPH 2D"};
        glutInit(&argc, argv);
        glutInitWindowSize(winw, winh);
        glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
        glutCreateWindow(win);
        glutDisplayFunc(Draw);
        init();
        glutMainLoop();

    return 0;
}
