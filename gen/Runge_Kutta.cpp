#include <vector>
#include <array>
#include <functional>

using namespace std;

namespace Runge_Kutta
{
    //// GENERIC RUNGE-KUTTA SCHEME CLASS
    template <const int RK_num, const int ODE_num>
    class Runge_Kutta_Scheme
    {
    public:
        using Tableau = array<array<long double, RK_num>, RK_num>;
        using Chart = array<long double, RK_num>;
        Runge_Kutta_Scheme(
            const Tableau& given_tableau,                                                               // Butcher tableau
            const Chart& given_stepping_factors,                                                        // Time stepping factors
            const Chart& given_weighting_factors                                                        // Weighted sum factors
        ) : _tableau(given_tableau),
            _stepping_factors(given_stepping_factors),
            _weighting_factors(given_weighting_factors)
        {}
        pair<vector<long double>, array<vector<long double>, ODE_num>>                                   // Returns value of each solved ODE for all time steps
        solve(
            long double h,                                                                               // Stepping size
            long double initial_conditions[ODE_num],                                                     // Initial conditions
            array<function<long double(long double r, array<long double, ODE_num> y)>, ODE_num> ODES,    // Array of ODEs
            function<bool(long double r, array<long double, ODE_num> y)> solved                          // Condition for ODE being solved
        )
        {
            // Boolean state determined by solved() function will stop solver
            bool _solved = false;

            // Initial conditions
            // Radius vector starting from initial step size to avoid singularity CHECK IF REQUIRED IN TOV
            vector<long double> _t = {h};
            array<vector<long double>, ODE_num> _y_n;

            // Apply initial conditions to the values array
            for (int ODE = 0; ODE < ODE_num; ODE++)
            {
                _y_n[ODE].push_back(initial_conditions[ODE]);
            };

            // Step forward in time using RK scheme
            while (!_solved)
            {
                // Grab radius and pressure for THIS step
                long double t = _t.back();

                // Append radius of THIS step
                _t.push_back(t + h);

                // Initialise Runge-Kutta variables
                long double K_values[RK_num][ODE_num] = {0};

                // Make list of latest values, transposer
                array<long double, ODE_num> yT;
                for (int ODE = 0; ODE < ODE_num; ODE++)
                {
                    // Grab last point in ODE
                    yT[ODE] = _y_n[ODE].back();
                };

                // For k_1 each y value input into function must be the current y value
                array<long double, ODE_num> y_for_each_ODE = yT;

                // Perform RK
                for (int step = 0; step < RK_num; step++)
                {
                    // Obtain all K values
                    for (int ODE = 0; ODE < ODE_num; ODE++)
                    {
                        // Get sum from butcher tableau and prior steps
                        for (int sub_ODE = 0; sub_ODE < ODE; sub_ODE++)
                        {
                            // Sum for this ODE to encapsulate all prior steps applied via butcher tableau
                            for (int sub_step = 0; sub_step < step; sub_step++)
                            {
                                // Add new value for each y to list
                                y_for_each_ODE[sub_ODE] += h * _tableau[step][sub_step] * K_values[sub_step][sub_ODE];
                            };
                        };
                        // Apply applicable ODE, tableau, and stepping to get next value
                        K_values[step][ODE] = ODES[ODE](t + h * _stepping_factors[step], y_for_each_ODE);
                    };
                };

                // Obtain results for this time step for each ODE
                for (int ODE = 0; ODE < ODE_num; ODE++)
                {
                    long double y_prior = _y_n[ODE].back();
                    long double sum = 0;
                    for (int step = 0; step < RK_num; step++)
                    {
                        sum += _weighting_factors[step] * K_values[step][ODE];
                    };
                    _y_n[ODE].push_back(y_prior + h * sum); // Find y_n+1
                };

                // Solved criterion
                if (solved(t, yT))
                {
                    // Stop solver
                    _solved = true;
                    // Remove step beyond solved point CHECK THIS
                    for (int ODE = 0; ODE < ODE_num; ODE++)
                    {
                        _y_n[ODE].pop_back();
                    };
                    _t.pop_back();
                };
            };

            return {_t, _y_n};
        };
        private:
            Tableau _tableau;           // Butcher tableau
            Chart _stepping_factors;    // Time stepping factors
            Chart _weighting_factors;   // Weighted sum factors
    };

    // Default Tableau for RK4 Scheme
    constexpr Runge_Kutta_Scheme<4, 1>::Tableau default_tableau4 =
    {{
        {{ 0.0L , 0.0L , 0.0L , 0.0L }},
        {{ 0.5L , 0.0L , 0.0L , 0.0L }},
        {{ 0.0L , 0.5L , 0.0L , 0.0L }},
        {{ 0.0L , 0.0L , 1.0L , 0.0L }}
    }};
    constexpr Runge_Kutta_Scheme<4, 1>::Chart default_stepping4 =
    {{
        0.0L , 0.5L , 0.5L , 0.0L
    }};
    constexpr Runge_Kutta_Scheme<4, 1>::Chart default_weighting4 =
    {{
        1.0L/6.0L , 1.0L/3.0L , 1.0L/3.0L , 1.0L/6.0L
    }};

    // Default Solver for RK4 for Single and Pair of ODE Cases
    static Runge_Kutta_Scheme<4, 1> _RK4(default_tableau4, default_stepping4, default_weighting4);
    static Runge_Kutta_Scheme<4, 2> _RK4_2(default_tableau4, default_stepping4, default_weighting4);

    pair<vector<long double>, array<vector<long double>, 1>> RK4(
        long double h,                                                            // Stepping size
        long double initial_conditions[1],                                        // Initial conditions
        array<function<long double(long double, array<long double, 1>)>, 1> ODES, // Array of ODEs
        function<bool(long double, array<long double, 1>)> solved                 // Condition for ODE being solved
    )
    {
        return _RK4.solve( h, initial_conditions, ODES, solved);                  // Return value from predefined solver
    };
    pair<vector<long double>, array<vector<long double>, 2>> RK4_2(
        long double h,                                                            // Stepping size
        long double initial_conditions[2],                                        // Initial conditions
        array<function<long double(long double, array<long double, 2>)>, 2> ODES, // Array of ODEs
        function<bool(long double, array<long double, 2>)> solved                 // Condition for ODE being solved
    )
    {
        return _RK4_2.solve( h, initial_conditions, ODES, solved);                // Return value from predefined solver
    };
}
