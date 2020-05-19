#pragma once 
#include "libs.h"

std::array<int, 5> get_energy_neighbourhood(int p) ;
std::array<int, 3> get_neighbourhood(int p);
float compute_energy(unsigned char* lat, int x, int y, const int mx, const int my) ;

