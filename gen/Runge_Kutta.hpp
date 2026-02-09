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
        );
        pair<vector<long double>, array<vector<long double>, ODE_num>>                                   // Returns value of each solved ODE for all time steps
        solve(
            long double h,                                                                               // Stepping size
            long double initial_conditions[ODE_num],                                                     // Initial conditions
            array<function<long double(long double r, array<long double, ODE_num> y)>, ODE_num> ODES,    // Array of ODEs
            function<bool(long double r, array<long double, ODE_num> y)> solved                          // Condition for ODE being solved
        );
        private:
            Tableau _tableau;           // Butcher tableau
            Chart _stepping_factors;    // Time stepping factors
            Chart _weighting_factors;   // Weighted sum factors
    };
    // Default Solver for RK4 for Single and Pair of ODE Cases
    pair<vector<long double>, array<vector<long double>, 1>>
    RK4(
        long double h,                                                            // Stepping size
        long double initial_conditions[1],                                        // Initial conditions
        array<function<long double(long double, array<long double, 1>)>, 1> ODES, // Array of ODEs
        function<bool(long double, array<long double, 1>)> solved                 // Condition for ODE being solved
    );
    pair<vector<long double>, array<vector<long double>, 2>>
    RK4_2(
        long double h,                                                            // Stepping size
        long double initial_conditions[2],                                        // Initial conditions
        array<function<long double(long double, array<long double, 2>)>, 2> ODES, // Array of ODEs
        function<bool(long double, array<long double, 2>)> solved                 // Condition for ODE being solved
    );
}