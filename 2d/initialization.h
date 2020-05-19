#pragma once
#include "libs.h"
#include "energy.h"

void init_lat(unsigned char* lat, const int mx, const int my);
void init_dict(std::unordered_map<float, std::vector<int>> &dict, float* rates, const int size);
void calculate_prob(float* arr, std::vector<double> &rates, std::vector<double> &count, int size) ;
float* init_rates(unsigned char* lat,float*arr,  const int mx, const int my, const float beta, const int n) ;