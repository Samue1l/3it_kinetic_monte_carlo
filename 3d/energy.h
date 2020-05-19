#pragma once 
#include "libs.h"
/*
float e_rest = 0.5f;
float e_100 = 1.2f;
float e_111 = 1.2f;
float anisotropy[25] = {
                                        e_100, e_rest, e_111, e_rest, e_100, 
                                        e_rest, e_100, e_111, e_100, e_rest, 
                                        e_111, e_111, 0.0f, e_111, e_111,
                                        e_rest, e_100, e_111, e_100, e_rest,
                                        e_100, e_rest, e_111, e_rest, e_100};

*/
std::array<int, 5> get_energy_neighbourhood(int p) ;

std::array<int, 3> get_neighbourhood(int p);
float compute_energy(unsigned char* lat, int x, int y, int z, const int mx, const int my, const int mz);
