#pragma once

#include "base_meta_func.hpp"

#include <GL/glut.h>

#include <boost/array.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/bind.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/function.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/lambda.hpp>

#include <vector>
#include <memory>
#include <utility>

static int winw = 1480;
static int winh = 800;

static const int steps {15000};
static const double dt {.01};
static const float Y {45.};
static const float X {25.};

template <typename Elem, typename Alloc>
struct is_alloc : boost::mpl::false_ {};
template <typename Elem>
struct is_alloc<Elem, std::allocator<Elem> > : boost::mpl::true_ {};

template <typename Elem, typename Alloc = std::allocator<Elem> >
using check_mem = typename boost::mpl::if_c<is_alloc<Elem, Alloc>::value, std::vector<Elem, std::allocator<Elem> >,
                    boost::array<Elem, boost::is_integral<Alloc>::value > >::type;

template <typename Cont>
void resize(const Cont& in, Cont& out) {
    out.resize(std::size(in));
}
template <typename T, size_t N>
void resize(const boost::array<T, N>& in, const boost::array<T, N>& out) {}

template <typename value_type, typename Alloc>
void resize2(const check_mem<value_type, Alloc>& in, const check_mem<value_type, Alloc>& out) {
    if constexpr (is_alloc<value_type, Alloc>::value)
        out.resize(std::size(in));
    else {

    }
}
typedef check_mem<double> Cont_arr;
typedef std::vector<double> Cont_vec;

template<typename System_equ>
void Euler(System_equ sys_eq, Cont_vec& func, const double t, const double dt) {
    Cont_vec vec(func.size());
    sys_eq(func, vec, t);
    for(size_t i = 0; i < func.size(); ++i) func[i] += dt*vec[i];
}

struct container_algebra
{
    template<typename S1, typename S2, typename S3, typename Op>
    void for_each3(S1& s1, S2& s2, S3& s3, Op op) const
    {
        auto first1 = std::begin(s1);
        auto last1 = std::end(s1);
        auto first2 = std::begin(s2);
        auto first3 = std::begin(s3);
        for(; first1 != last1;)
            op(*first1++, *first2++, *first3++);
    }
    template<typename S1, typename S2, typename S3, typename S4, typename S5, typename S6>
    using vector_each = boost::mpl::vector<S1, S2, S3, S4, S5, S6>;
    struct for_eachs
    {
        for_eachs(std::ostream& s) : out(s) {}
        template<typename T>
        void operator()(boost::mpl::identity<T>) const {
            static T a;
            auto first1 = std::begin(a);
            auto last1 = std::end(a);
            out <<"Partial sum = " << /**first1 <<  */std::endl;
        }
    private:
        std::ostream& out;
    };

    template<typename S1, typename S2, typename S3, typename S4, typename S5, typename S6, typename Op>
    void for_each6(S1& s1, S2& s2, S3& s3, S4& s4, S5& s5, S6& s6, Op op) const
    {
        auto first1 = std::begin(s1);
        auto last1  = std::end(s1);
        auto first2 = std::begin(s2);
        auto first3 = std::begin(s3);
        auto first4 = std::begin(s4);
        auto first5 = std::begin(s5);
        auto first6 = std::begin(s6);
        for(; first1 != last1;)
            op(*first1++, *first2++, *first3++, *first4++, *first5++, *first6++);
    }
};

struct default_operation
{
    template<typename F1 = double, typename F2 = F1>
    struct scale_sum2
    {
        typedef void result_type;
        const F1 alpha1;
        const F2 alpha2;
        scale_sum2(F1 a1, F2 a2):alpha1(a1), alpha2(a2) {}
        template<typename T0, typename T1, typename T2>
        void operator() (T0& t0, const T1& t1, const T2& t2) const
        {
            t0 = alpha1*t1 + alpha2 *t2;
        }
    };
    template<typename F1 = double, typename F2 = F1, typename F3 = F1, typename F4 = F1, typename F5 = F1>
    struct scale_sum5
    {
        typedef void result_type;
        const F1 alpha1;
        const F2 alpha2;
        const F3 alpha3;
        const F4 alpha4;
        const F5 alpha5;
        scale_sum5(F1 a1, F2 a2, F3 a3, F4 a4, F5 a5):alpha1(a1), alpha2(a2), alpha3(a3), alpha4(a4), alpha5(a5) {}
        template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
        void operator() (T0& t0, const T1& t1, const T2& t2, const T3& t3, const T4& t4, const T5& t5) const
        {
            t0 = alpha1*t1 + alpha2 *t2 + alpha3*t3 + alpha4 *t4 + alpha5*t5;
        }
    };
    template<typename S1 = double, typename S2 = S1, typename S3 = S1, typename S4 = S1, typename S5 = S1>
    using vector_each = boost::mpl::vector<S1, S2, S3, S4, S5>;
    struct scale_sums
    {
        scale_sums(std::ostream& s) : out(s) {}
        template<typename T>
        void operator()(boost::mpl::identity<T>) const {
            static T a;
            static const T alpha;
            out <<"Partial sum = " << /**first1 <<  */std::endl;
        }
    private:
        std::ostream& out;
    };
};

template<typename value_type, template<typename Elem, typename Alloc = std::allocator<Elem> > class Cont>
class Runge_Kutta4
{
public:
    typedef value_type time_type;


    template<typename System_eq>
    void do_step(System_eq& sys_eq, Cont<value_type>& func, time_type t, time_type dt)
    {
        adjust_size(func);
        const value_type one = 1;
        const time_type dt2 = dt/2, dt3 = dt/3, dt6 = dt/6;
        typedef typename default_operation::template scale_sum2<value_type, time_type> scale_sum2;
        typedef typename default_operation::template scale_sum5<value_type, time_type, time_type, time_type, time_type> scale_sum5;
        using vector_algebra5 = boost::mpl::vector5<Cont<value_type>, Cont<value_type>, Cont<value_type>, Cont<value_type>, Cont<value_type>>;
        sys_eq(func, K1, t);
        m_algebra.for_each3(func_tmp, func, K1, scale_sum2(one, dt2));
        sys_eq(func_tmp, K2, t+dt2);
        m_algebra.for_each3(func_tmp, func, K2, scale_sum2(one, dt2));
        sys_eq(func_tmp, K3, t+dt2);
        m_algebra.for_each3(func_tmp, func, K3, scale_sum2(one, dt));
        sys_eq(func_tmp, K4, t+dt);
        m_algebra.for_each6(func, func, K1, K2, K3, K4, scale_sum5(one, dt6, dt3, dt3, dt6));
    }
private:
    Cont<value_type> func_tmp, K1, K2, K3, K4;
    container_algebra m_algebra;
    void adjust_size(Cont<value_type>& vec) {
        resize(vec, func_tmp);
        resize(vec, K1);
        resize(vec, K2);
        resize(vec, K3);
        resize(vec, K4);
    }
};


using runge_solver = Runge_Kutta4<double, std::vector>;


struct lorenz
{
    const double sigma, R, b;
    lorenz(const double sigma, const double R, const double b):sigma(sigma), R(R), b(b) {}
    void operator() (const Cont_vec& x, Cont_vec& dxdt, double t)
    {
        dxdt[0] = sigma*(x[1] - x[0]);
        dxdt[1] = R*x[0] - x[1] - x[0]*x[2];
        dxdt[2] = -b*x[2] + x[0]*x[1];
    }
};
void init() {
    glClearColor(0.25, 0.0, 0.2, 1.0);
    glEnable(GL_LIGHT0);
}


///----------------------------------------------------------------------------------------
void glWrite(float x, float y, int *font, const char* text, size_t kls) {
    glRasterPos2f(x, y);
    for (size_t i  =0; i < kls; i++)
    glutBitmapCharacter(font, text[i]);
}

void Draw_Graph()
{
    glColor3f(1,1,0);
        glWrite(-.003, -.001, (int*)GLUT_BITMAP_HELVETICA_18, (char*)"(0, 0)", 6);
        glWrite(-.003, Y-.8, (int*)GLUT_BITMAP_HELVETICA_18, (char*)"X(t)", 5);
        glWrite(X-.5, -.01, (int*)GLUT_BITMAP_HELVETICA_18, (char*)"t", 2);
        runge_solver stepper;
        lorenz system(10., 28., 8./3.);
        Cont_vec x(3, 1.);
        x[0] = 10.;
        for (size_t i = 0; i < steps; i++)
        {
            stepper.do_step(system, x, i*dt, dt);
            // std::cout << i*dt << " ";
            // std::cout << x[0] << " = " << x[1] << " = " << x[2] << std::endl;
        }

        glBegin(GL_LINES);
            glColor3f(1, 0, 0); glVertex2f(-200, 0); glVertex2f(400, 0);
            glColor3f(0, 0, 1); glVertex2f(0, -100); glVertex2f(0, 100);
        glEnd();
        glBegin(GL_LINE_STRIP);
            for (size_t i = 0; i < steps; i++)
            {
                stepper.do_step(system, x, i*dt, dt);
                glColor3f(1, 1, 1); glVertex2f(x[0], x[2]);
            }
        glEnd();
    glFlush();
}

//----------------------------------------------------------------------------------------
void Draw()
{
    glClear( GL_COLOR_BUFFER_BIT |
            GL_DEPTH_BUFFER_BIT );
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-25.,X,-25.,Y,.0,10.);
    Draw_Graph();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glutSwapBuffers();
}
