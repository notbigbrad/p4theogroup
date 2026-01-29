#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <utility>

#include "./constants.hpp"
#include "./butcher_tableau.hpp"

// this must be removed
#include "./project.hpp"

using namespace std;
using namespace project;

namespace Runge_Kutta
{
    pair<vector<long double>, vector<long double>> runge_kutta_scheme(vector<long double> inital_conditions, long double h, long double n, long double P_c, long double K)
    {
        //// GENERIC RUNGE-KUTTA SCHEME FOR COUPLED OR NON COUPLED ODE

        // Inputs
        //
        // initial_conditions - vector of initial conditions
        //



        // <TEMPORARY ZONE FOR THINGS THAT MUST BE MADE MODULAR>

        // template params
        int RK_num = 4;
        int ODE_num = 2;
        auto tableau = butcher_tableau::tableau4;
        auto stepping_factors = butcher_tableau::stepping4;
        auto weighting_factors = butcher_tableau::weighting4;

        // needs to be broken out, see lab p1 functionals or create custom wrapper that fits requirements
        auto ODES = {generic_first_order_ODE_1, generic_first_order_ODE_2};

        // </TEMPORARY ZONE FOR THINGS THAT MUST BE MADE MODULAR>



        // Boolean turns true when pressure drops below 0 and this determines the radial extent of star
        bool solved = false;

        // Initial conditions
        // Radius vector starting from initial step size to avoid singularity CHECK IF REQUIRED IN TOV
        vector<long double> _t = {h};
        vector<long double> _values[ODE_num];

        // Step forward in time using RK scheme
        while (!solved)
        {
            // Grab radius and pressure for THIS step
            long double t = _t.back();

            // Append radius of THIS step
            _t.push_back(t + h);

            // Initialise Runge-Kutta variables
            long double RK_values[RK_num][ODE_num];

            // Perform RK
            for (int step = 0; step < RK_num; step++)
            {
                // Run for each ODE for each step keeping coupled
                for (int ODE = 0; ODE < ODE_num; ODE++)
                {
                    // Get sum from butcher tableau and prior steps
                    long double sums[ODE_num] = { 0 };
                    for (int sub_ODE = 0; sub_ODE < ODE_num; sub_ODE++)
                    {
                        for (int sub_step = 0; sub_step < step - 1; sub_step++)
                        {
                            sums[sub_ODE] += tableau[step][sub_step] * RK_values[sub_step][sub_ODE];
                        };
                    };
                    // Apply applicable ODE, tableau, and stepping to get next value
                    RK_values[step][ODE] = generic_first_order_ODE_1(
                        t + h * stepping_factors[step],
                        _values[0].back() + h * sums[0],
                        _values[1].back() + h * sums[1]);
                };
            };

            // Obtain results for this time step
            for (int ODE = 0; ODE < ODE_num; ODE++)
            {
                long double y_prior = _values[ODE].back();
                long double sum = 0;
                for (int step = 0; step < RK_num; step++)
                {
                    sum += weighting_factors[step] * RK_values[step][ODE];
                };
                _values[ODE].push_back(y_prior + h * sum);
            };

            // THIS IS NOT FULLY MODULAR SUCCESS CRITERION, MAYBE MOVE TO CALLER INSTEAD OF SOLVER
            // If pressure is negative then set to solved and remove non-physical step
            if (_values[1].back() <= 0 || isnan(_values[1].back()))
            {
                solved = true;
                // Pop non-physical step
                for (int ODE = 0; ODE < ODE_num; ODE++)
                {
                    _values[ODE].pop_back();
                };
            }
        }

        // return {_r, _P};
    };   
}