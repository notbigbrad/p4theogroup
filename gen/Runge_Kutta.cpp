#include <vector>
#include <array>
#include <functional>
#include <iostream>

// #include "./butcher_tableau.hpp"

using namespace std;

namespace Runge_Kutta
{
    //// GENERIC RUNGE-KUTTA SCHEME
    // template<int RK_num, int ODE_num>                                                               // Template for differing RK schemes and ODE counts
    constexpr int ODE_num = 2;
    constexpr int RK_num = 4;
    pair< vector<long double>, array<vector<long double>, ODE_num> >                                // Returns value of each solved ODE for all time steps
    runge_kutta_scheme
    (
        long double initial_conditions[ODE_num],                                                    // Initial conditions
        array<function<long double( long double r, array<long double, ODE_num> y )>, ODE_num> ODES, // Array of ODEs
        long double tableau[RK_num][RK_num],                                                        // Butcher tableau
        long double stepping_factors[RK_num],                                                       // Time stepping factors
        long double weighting_factors[RK_num],                                                      // Weighted sum factors
        long double h,                                                                              // Stepping size
        function<bool( long double r, array<long double, ODE_num> y )> solved                       // Condition for ODE being solved
    )
    {
        // Boolean state determined by solved() function will stop solver
        bool _solved = false;

        // Initial conditions
        // Radius vector starting from initial step size to avoid singularity CHECK IF REQUIRED IN TOV
        vector<long double> _t = {h};
        array<vector<long double>, ODE_num> _values;

        // Apply initial conditions to the values array
        for (int ODE = 0; ODE < ODE_num; ODE++)
        {
            _values[ODE].push_back(initial_conditions[ODE]);
        };

        // Step forward in time using RK scheme
        while (!_solved)
        {
            // Grab radius and pressure for THIS step
            long double t = _t.back();

            // Append radius of THIS step
            _t.push_back(t + h);

            // Initialise Runge-Kutta variables
            long double RK_values[RK_num][ODE_num] = {};
            array<long double, ODE_num> yS = {};

            // Perform RK
            for (int step = 0; step < RK_num; step++)
            {
                // Run for each ODE for each step keeping coupled
                for (int ODE = 0; ODE < ODE_num; ODE++)
                {
                    // Get sum from butcher tableau and prior steps
                    for (int sub_ODE = 0; sub_ODE < ODE_num; sub_ODE++)
                    {
                        array<long double, ODE_num> sums = {0};
                        // Sum for this ODE to encapsulate all prior steps applied via butcher tableau
                        for (int sub_step = 0; sub_step < step; sub_step++)
                        {
                            sums[sub_ODE] += tableau[step][sub_step] * RK_values[sub_step][sub_ODE];
                        };
                        // Add new value for each y to list
                        yS[sub_ODE] = _values[sub_ODE].back() + h * sums[sub_ODE];
                    };
                    // Apply applicable ODE, tableau, and stepping to get next value
                    RK_values[step][ODE] = ODES[ODE]( t + h * stepping_factors[step], yS );
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

            // Make list of latest values
            array<long double, ODE_num> yV;
            for (int ODE = 0; ODE < ODE_num; ODE++)
            {
                // Grab last point in ODE
                yV[ODE] = _values[ODE].back();
            };

            // Solved criterion
            if ( solved( t, yV ) )
            {
                // Stop solver
                _solved = true;
                // Remove step beyond solved point CHECK THIS
                for (int ODE = 0; ODE < ODE_num; ODE++)
                {
                    _values[ODE].pop_back();
                };
                _t.pop_back();
            };
        };

        return { _t, _values };
    };
}
