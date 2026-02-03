#include <vector>
#include <functional>

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
        long double inital_conditions[ODE_num],                                                     // Initial conditions
        array<function<long double( long double r, array<long double, ODE_num> y)>, ODE_num> ODES, // Array of ODEs
        long double tableau[RK_num][RK_num],                                                        // Butcher tableau
        long double stepping_factors[RK_num],                                                       // Time stepping factors
        long double weighting_factors[RK_num],                                                      // Weighted sum factors
        long double h,                                                                              // Stepping size
        function<bool( long double r, array<long double, 2> y )> solved                              // Condition for ODE being solved
    );
}