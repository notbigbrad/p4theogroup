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
    //// Wrapper to use the generic runge kutta scheme in the fourth order to solve LE

    // Polytropic equation state
    constexpr long double equation_of_state(long double P, long double n, long double K)
    {
        return powl(P/K, n / (n+1.0L)); // Basic polytropic equation of state
    }

    // Re-wrap functions to be in the desired generic form
    long double theta(long double xi, array<long double, 2> y)
    {
        return y[0];
    }
    std::function<long double(long double, std::array<long double, 2>)> get_dtheta(long double n)
    {
        return [n](long double xi, std::array<long double, 2> y) -> long double {
            return -2.0L / xi * y[0] - (powl(y[1], n));
            // return -2.0L / xi * (powl(y[1], n) + xi * y[0]);
        };
    }

    // Solved criterion must also take same input form as ODEs but must be boolean valued
    bool LE_criterion( long double xi, array<long double, 2> y)
    {
        if (xi > 20.0L) { return true; }
        else { return false; }
    }
    bool LE_criterion2( long double xi, array<long double, 2> y)
    {
        if (y[1] <= 0 || isnan(y[0]) || xi > 20.0L) { return true; }
        else { return false; }
    }

    // Use generic solver to solve LE
    pair< vector<long double>, array<vector<long double>, 2> >  LE_RK_wrapper(int n)
    {
        // Initial conditions to turn ODE into IVP
        long double h = 0.01L;
        long double theta_0 = {1.0L - (h * h) / 6.0L}; // Avoid singularity by using xi0=h
        long double dtheta_0 = -h / 3.0L;
        long double init[2] = { dtheta_0, theta_0 };

        // List of ODEs
        array<function<long double( long double r, array<long double, 2> y )>, 2> ODES = {get_dtheta(n), &theta};

        // Solved criterion
        function<bool( long double r, array<long double, 2> y )> solved_criterion = LE_criterion;

        // Apply RK scheme to problem
        auto t = Runge_Kutta::runge_kutta_scheme( h, init, ODES, solved_criterion);

        return t;
    }
    pair< vector<long double>, array<vector<long double>, 2> >  LE_RK_wrapper2(int n)
    {
        // Initial conditions to turn ODE into IVP
        long double h = 0.01L;
        long double theta_0 = {1.0L - (h * h) / 6.0L}; // Avoid singularity by using xi0=h
        long double dtheta_0 = -h / 3.0L;
        long double init[2] = { dtheta_0, theta_0 };

        // List of ODEs
        array<function<long double( long double r, array<long double, 2> y )>, 2> ODES = {get_dtheta(n), &theta};

        // Solved criterion
        function<bool( long double r, array<long double, 2> y )> solved_criterion = LE_criterion2;

        // Apply RK scheme to problem
        auto t = Runge_Kutta::runge_kutta_scheme( h, init, ODES, solved_criterion);

        return t;
    }

    void LE_mass_radius(long double n, long double K, long double P_c_min, long double P_c_max, int runs)
    {
        vector<long double> _xi;
        array<vector<long double>, 2> _thetas;

        tie(_xi, _thetas) = LE_RK_wrapper2(n);

        vector<long double> _theta = _thetas[0];


        vector<long double> P_c_list;

        for (int i = 0; i < runs; i++)
        {
            P_c_list.push_back(P_c_min * powl(P_c_max / P_c_min, static_cast<long double>(i) / (runs - 1)));
        }

        
        for (long double P_c : P_c_list)
        {
            long double alpha = sqrtl((n + 1.0L) * K * powl(equation_of_state(P_c, n, K), (1.0L / n - 1)) / (4.0L * __pi * __gravitation));
            long double R = alpha * _xi.back();
            long double M = 0.0L;
            for (int j = 1; j < _xi.size(); j++)
            {
                long double r = alpha * _xi[j];
                long double dr = alpha * (_xi[j] - _xi[j - 1]);
                long double rho_r = equation_of_state(P_c, n, K) * powl(_theta[j], n);

                M += 4.0L * __pi * (r * r) * rho_r * dr;
            }

            cout << M << "," << R << endl;
        }
    }

    void LE_n(long double h, long double n_min, long double n_max, long double n_h, long double xi_max)
    {
        vector<string> csv;

        vector<long double> r;
        array<vector<long double>, 2> y;

        // Unwrap answer to ODE
        tie(r, y) = LE_RK_wrapper(n_min);

        // Resize csv string array
        csv.resize(y[0].size());

        // Grab desired values
        for (int j = 0; j < size(y[0]); j++)
        {
            csv[j] += to_string(r[j]) + "," + to_string(y[1][j]);
        }
        
        for (long double i = n_min + n_h; i <= n_max; i += n_h)
        {
            vector<long double> r;
            array<vector<long double>, 2> y;

            // Unwrap answer to ODE
            tie( r, y ) = LE_RK_wrapper(i);

            // Grab desired values
            for (int j = 0; j < size(y[1]); j++)
            {
                csv[j] += "," + to_string(y[1][j]);
            }
        }
        for (auto& entry : csv)
        {
            cout << entry << endl;
        }
    }
}