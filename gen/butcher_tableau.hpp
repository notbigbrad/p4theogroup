#include <cmath>
using namespace std;

namespace butcher_tableau
{
    // RK4 Butcher Tableau, Stepping Factors, and Summing Weights
    // [step][weight]
    long double tableau4[4][4] =
    {
        { 0   , 0   , 0   , 0   },
        { 1/2 , 0   , 0   , 0   },
        { 0   , 1/2 , 0   , 0   },
        { 0   , 0   , 1   , 0   }
    };
    long double stepping4[4] =
    {
        0   , 1/2 , 1/2 , 0
    };
    long double weighting4[4] =
    {
        1/6 , 1/3 , 1/3 , 1/6
    };
}