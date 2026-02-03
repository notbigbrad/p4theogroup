# clear
# g++ ./main.cpp -O3 -ffast-math -march=native -pthread --std=c++26 -o ./o.exe
g++ ./main.cpp ./source/constants.hpp ./source/project.hpp ./source/TOV_RK_wrapper.cpp ./source/LE_RK_wrapper.cpp ./gen/butcher_tableau.hpp ./gen/Runge_Kutta.hpp ./gen/Runge_Kutta.cpp --std=c++26 -o ./.o.exe
./.o.exe
