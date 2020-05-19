#pragma once
#include "libs.h"
#include "energy.h"


double update_lattice(unsigned char* lat, float* arr,std::vector<double> &rates, std::vector<double> &c_rates,std::unordered_map<float, std::vector<int>> &dict,
                             const int mx, const int my, const int n, const float beta, int i);
