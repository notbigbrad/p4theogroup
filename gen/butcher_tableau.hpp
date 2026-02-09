#include <cmath>
using namespace std;

namespace butcher_tableau
{
    // RK4 Butcher Tableau, Stepping Factors, and Summing Weights
    // [step][weight]
    const long double tableau4[4][4] =
    {
        { 0.0L , 0.0L , 0.0L , 0.0L },
        { 0.5L , 0.0L , 0.0L , 0.0L },
        { 0.0L , 0.5L , 0.0L , 0.0L },
        { 0.0L , 0.0L , 1.0L , 0.0L }
    };
    const long double stepping4[4] =
    {
        0.0L , 0.5L , 0.5L , 0.0L
    };
    const long double weighting4[4] =
    {
        1.0L/6.0L , 1.0L/3.0L , 1.0L/3.0L , 1.0L/6.0L
    };
}