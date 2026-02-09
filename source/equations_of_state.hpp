#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <utility>
#include <functional>

#include "./constants.hpp"

using namespace std;
using namespace constants;

namespace EoS
{
    // Polytropic equation state
    constexpr long double polytropic_equation_of_state(long double P, long double n, long double K)
    {
        return powl(P/K, n / (n+1.0L)); // Basic polytropic equation of state
    }
    // Non-relativistic white dwarf polytropic equation state
    const long double mu_e_wd = 2.0L;
    const long double n_wd_nonR = 1.5L;
    const long double K_wd_nonR = 1.0L / 20.0L * powl((3.0L / __pi), (2.0L / 3.0L)) * ((__planck * __planck) / __mass_electron) * powl((1.0L / (mu_e_wd * __mass_proton)), (5.0L / 3.0L));
    constexpr long double equation_of_state_wd_non_r(long double P)
    {
        return polytropic_equation_of_state(P, n_wd_nonR, K_wd_nonR);
    }
    // Relativistic white dwarf polytropic equation state
    const long double n_wd_R = 3.0L;
    const long double K_wd_R = 1.0L / 8.0L * powl((3.0L / __pi), (2.0L / 3.0L)) * (__planck * __c) * powl((1.0L / (mu_e_wd * __mass_proton)), (4.0L / 3.0L));
    constexpr long double equation_of_state_wd_r(long double P)
    {
        return polytropic_equation_of_state(P, n_wd_R, K_wd_R);
    }
    // Interpolated numerical equations of state
    constexpr long double SLy4(long double P)
    {
        // csv data format:
        // a,b,d,c,x_min,discard
        // linear sorted x_min
        // goto minimum x_min < x
        // return a+b*x+c*x^2+d*x^3
        return 0;
    }
    constexpr long double NPS(long double P)
    {
        return 0;
    }
}