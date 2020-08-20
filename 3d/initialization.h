#pragma once
#include "libs.h"
#include "energy.h"
void init_lat(unsigned char* lat,unsigned char* copied_lat,  const int mx, const int my, const int mz);

void init_dict(std::unordered_map<float, std::vector<int>> &dict, float* rates, const int size);
void compute_rates(unsigned char* lat, float* arr, int x, int y, int z, const int mx, const int my, const int mz, const float beta, const int n);
void calculate_prob(float* arr, std::vector<double> &rates, std::vector<double> &count, int size) ;
float* init_rates(unsigned char* lat, const int mx, const int my, const int mz, const float beta, const int n);

void recompute_epitaxy(unsigned char* lat, int* epi, int xo, int yo, const int mx, const int my, const int mz) ;
void init_epitaxy(unsigned char* lat, float* arr, int* epi, float epitaxial_rate, const int mx, const int my, const int mz, const int n) ;
