#include <vector>
#include <functional>

// #include "./butcher_tableau.hpp"

using namespace std;

namespace Runge_Kutta
{
    //// GENERIC RUNGE-KUTTA SCHEME
    template<int RK_num, int ODE_num>                                                 // Template for differing RK schemes and ODE counts
    pair< vector<long double>, array<vector<long double>, ODE_num> >                  // Returns value of each solved ODE for all time steps
    runge_kutta_scheme
    (
        long double inital_conditions[ODE_num],                                       // Initial conditions
        function<long double( long double r, long double y[ODE_num] )> ODES[ODE_num], // Array of ODEs
        long double tableau[RK_num][RK_num],                                          // Butcher tableau
        long double stepping_factors[RK_num],                                         // Time stepping factors
        long double weighting_factors[RK_num],                                        // Weighted sum factors
        long double h,                                                                // Stepping size
        function<bool( long double r, long double y[ODE_num] )> solved                // Condition for ODE being solved
    );
}