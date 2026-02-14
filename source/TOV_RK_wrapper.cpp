#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <utility>
#include <functional>
#include <string>
#include <fstream>

#include "./equations_of_state.hpp"
#include "../gen/Runge_Kutta.hpp"

using namespace std;
using namespace constants;
using namespace EoS;

namespace project
{
    //// Wrapper to use the generic runge kutta scheme in the fourth order to solve TOV



    // Equations
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

    // Re-wrap functions to be in the desired generic form
    std::function<long double(long double, std::array<long double, 2>)> generic_first_order_ODE_1(function<long double(long double)> equation_of_state)
    {
        return [equation_of_state](long double r, std::array<long double, 2> y) -> long double
        {
            return mass_continuity(r, equation_of_state(y[1]));
        };
    }
    std::function<long double(long double, std::array<long double, 2>)> generic_first_order_ODE_2(function<long double(long double)> equation_of_state)
    {
        return [equation_of_state](long double r, std::array<long double, 2> y) -> long double
        {
            return tolman_oppenheimer_volkoff(r, y[1], y[0], equation_of_state(y[1]));
        };
    }



    // Solved criterion must also take same input form as ODEs but must be boolean valued
    bool TOV_criterion(long double r, array<long double, 2> y)
    {
        if (y[1] <= 0 || isnan(y[1]) || r > 1000000000000)
        {
            return true;
        }
        else
        {
            return false;
        }
    }



    // Use generic solver to solve TOV
    pair<vector<long double>, array<vector<long double>, 2>> TOV_RK_wrapper(long double h, long double P_c, function<long double(long double)> equation_of_state)
    {
        // Initial conditions to turn ODE into IVP
        long double M_c = 0.0L;
        long double init[2] = {M_c, P_c};

        // List of ODEs
        array<function<long double(long double r, array<long double, 2> y)>, 2> ODES = {generic_first_order_ODE_1(equation_of_state), generic_first_order_ODE_2(equation_of_state)};

        // Solved criterion
        function<bool(long double r, array<long double, 2> y)> solved_criterion = TOV_criterion;

        // Apply RK scheme to problem
        auto t = Runge_Kutta::RK4_2(h, init, ODES, solved_criterion);

        return t;
    }



    // Find mass-radius trend for many stars
    void TOV_mass_radius(long double h, long double P_c_min, long double P_c_max, int runs, function<long double(long double)> equation_of_state, string f_name)
    {
	ofstream file(f_name);

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
            tie(r, y) = TOV_RK_wrapper(h, P_c, equation_of_state);

            // Grab desired values
            long double R = r.back();
            long double M = y[0].back();

            file << M << "," << R << endl;
        }
    }



    // Find pressure radius plot for a single star
    void TOV_Pr(long double h, long double P_c, function<long double(long double)> equation_of_state, string f_name)
    {
	ofstream file(f_name);

        vector<long double> r;
        array<vector<long double>, 2> y;

        // Unwrap answer to ODE
        tie(r, y) = TOV_RK_wrapper(h, P_c, equation_of_state);

        // Grab desired values
        for (int j = 0; j < size(y[0]); j++)
        {
            file << to_string(r[j]) + "," + to_string(y[0][j]) + "," + to_string(y[1][j]) << endl;
        }
    }
}
