#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <utility>
#include <functional>

#include "./constants.hpp"
#include "../gen/Runge_Kutta.hpp"

using namespace std;
using namespace project;

namespace project
{
    //// Wrapper to use the generic runge kutta scheme in the fourth order to solve TOV

    // Equations
    constexpr long double equation_of_state(long double P, long double n, long double K)
    {
        return powl(P/K, n / (n+1.0L)); // Basic polytropic equation of state
    }
    constexpr long double mass_continuity(long double r, long double rho)
    {
        return __4pi * r * r * rho;
    }
    constexpr long double tolman_oppenheimer_volkoff(long double r, long double P, long double M, long double rho)
    {
        long double term1 = (rho + P / __c2);
        long double term2 = (M + 4.0L * __pi * r * r * r * P / __c2);
        long double term3 = r * (r - 2.0L * __gravitation * M / __c2);
        return -(__gravitation * term1 * term2) / term3;
    }

    // Apply the equation of state to create generic first order couple ODEs for the solver
    // Solver must take ODEs in the form: F_n = dy_n/dt(r, {M, P})

    // Needed for equation of state
    const long double K_white_dwarf_non_rel = 1.0L / 20.0L * powl((3.0L / __pi), (2.0L / 3.0L)) * ((__planck * __planck) / __mass_electron) * powl((1.0L / (2.0L * __mass_proton)), (5.0L / 3.0L));
    const long double n_white_dwarf_non_rel = 1.5L;

    // Re-wrap functions to be in the desired generic form
    long double generic_first_order_ODE_1(long double r, array<long double, 2> y)
    {
        return mass_continuity(r, equation_of_state(y[1], n_white_dwarf_non_rel, K_white_dwarf_non_rel));
    }
    long double generic_first_order_ODE_2(long double r, array<long double, 2> y)
    {
        return tolman_oppenheimer_volkoff(r, y[1], y[0], equation_of_state(y[1], n_white_dwarf_non_rel, K_white_dwarf_non_rel));
    }

    // Solved criterion must also take same input form as ODEs but must be boolean valued
    bool TOV_criterion( long double r, array<long double, 2> y)
    {
        if (y[0] <= 0 || isnan(y[0]) || r > 1000000000) { return true; }
        else { return false; }
    }

    // Use generic solver to solve TOV
    pair< vector<long double>, array<vector<long double>, 2> >  TOV_RK_wrapper(long double h, long double P_c)
    {
        // Initial conditions to turn ODE into IVP
        long double M_c = 0.0L;
        long double init[2] = { M_c, P_c };

        // List of ODEs
        array<function<long double( long double r, array<long double, 2> y )>, 2> ODES = {&generic_first_order_ODE_1, &generic_first_order_ODE_2};

        // Solved criterion
        function<bool( long double r, array<long double, 2> y )> solved_criterion = TOV_criterion;

        // Apply RK scheme to problem
        auto t = Runge_Kutta::RK4_2( h, init, ODES, solved_criterion);

        return t;
    }

    void TOV_mass_radius(long double P_c_min, long double P_c_max, int runs)
    {
        vector<long double> P_c_list;

        for (int i = 0; i < runs; i++)
        {
            P_c_list.push_back(P_c_min * powl(P_c_max / P_c_min, static_cast<long double>(i) / (runs - 1)));
        }

        
        for (long double P_c : P_c_list)
        {
            vector<long double> r;
            array<vector<long double>, 2> y;

            // Unwrap answer to ODE
            tie( r, y ) = TOV_RK_wrapper(1000000.0L, P_c);

            // Grab desired values
            long double R = r.back();
            long double M = y[0].back();

            cout << M << "," << R << endl;
        }
    }

    void TOV_Pr(long double h, long double P_c)
    {
        vector<long double> r;
        array<vector<long double>, 2> y;

        // Unwrap answer to ODE
        tie(r, y) = TOV_RK_wrapper(1000000.0L, P_c);

        // Grab desired values
        for (int j = 0; j < size(y[0]); j++)
        {
            cout << to_string(r[j]) + "," + to_string(y[0][j]) + "," + to_string(y[1][j]) << endl;
        }
    }
}